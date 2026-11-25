#include "measure.hpp"
#include "InitialData.hpp"
#include <omp.h>
#include <algorithm>
#include <filesystem>
#include <limits>
#include <set>
#include <sstream>

namespace fs = std::filesystem;
namespace {
struct Options {
    fs::path config, input, output;
    double periods = 2, dt = 0;
    std::vector<std::string> cases{"C0","D0"};
    std::vector<double> radii;
};

double positive(const std::string &text)
{
    size_t end;
    const double value = std::stod(text,&end);
    if (end != text.size() || !std::isfinite(value) || value <= 0)
        throw std::invalid_argument("Expected a positive finite number: "+text);
    return value;
}

std::vector<std::string> split(const std::string &value)
{
    std::vector<std::string> result;
    std::istringstream input(value);
    std::string item;
    while (std::getline(input,item,',')) result.push_back(item);
    if (result.empty() || value.back()==',') throw std::invalid_argument("Empty list item");
    return result;
}

Options options(int argc, char **argv)
{
    if (argc < 3) throw std::invalid_argument("Expected config.toml and a MODESBIN input");
    Options o;
    o.config=fs::canonical(argv[1]); o.input=fs::canonical(argv[2]);
    for (int i=3; i<argc; i+=2) {
        if (i+1>=argc) throw std::invalid_argument("Missing option value");
        const std::string key=argv[i], value=argv[i+1];
        if (key=="--cases") o.cases=split(value);
        else if (key=="--periods") o.periods=positive(value);
        else if (key=="--dt") o.dt=positive(value);
        else if (key=="--radii") for (const auto &r:split(value)) o.radii.push_back(positive(r));
        else if (key=="--output") o.output=fs::weakly_canonical(fs::absolute(value));
        else throw std::invalid_argument("Unknown option: "+key);
    }
    std::set<std::string> seen;
    for (const auto &name:o.cases)
        if ((name!="C0" && name!="D0" && name!="D1") || !seen.insert(name).second)
            throw std::invalid_argument("Expected distinct cases from C0,D0,D1");
    if (seen.count("D1") && seen.size()!=1)
        throw std::invalid_argument("Run D1 separately with its corrected input file");
    if (o.output.empty()) throw std::invalid_argument("--output requires a fresh directory");
    return o;
}

// The production reader deliberately fills absent modes with zero. Inspect the
// header/records before allowing it to allocate, and reject ambiguous labels.
void verify_input(const fs::path &file, int nr)
{
    std::ifstream in(file,std::ios::binary);
    in.exceptions(std::ios::badbit|std::ios::failbit);
    char magic[8]; int32_t version; int64_t nodes,count;
    auto read = [&](auto &v){in.read(reinterpret_cast<char*>(&v),sizeof(v));};
    read(magic); read(version); read(nodes); read(count);
    const auto bytes=fs::file_size(file);
    const auto stride=8+16ULL*nr;
    if (std::string(magic,8)!="MODESBIN" || version!=1 || nodes!=nr || count<=0 ||
        (bytes-28)%stride || (bytes-28)/stride!=static_cast<uint64_t>(count))
        throw std::runtime_error("Invalid MODESBIN format, radial size or length");
    std::set<std::pair<int,int>> labels;
    for (int64_t i=0;i<count;++i) {
        int32_t l,m; read(l); read(m);
        if (l<0 || m < -l || m>l || !labels.emplace(l,m).second)
            throw std::runtime_error("Invalid or duplicate input mode");
        in.seekg(16LL*nr,std::ios::cur);
    }
    if (!labels.count({1,1})) throw std::runtime_error("Input lacks the (1,+1) mode");
}

Eigen::VectorXcd make_centered_211(const CrankNicolsonModes &s, double input_norm)
{
    Eigen::VectorXcd state=Eigen::VectorXcd::Zero(s.Nr*s.n_modes);
    const double bohr=1/(s.mu*s.mu);
    for (int j=0;j<s.n_modes;++j) if (s.modes[j].l==1 && s.modes[j].m==1) {
        state.segment(j*s.Nr,s.Nr)=
            (s.r.array().square()*(-s.r.array()/(2*bohr)).exp()/
             (2*std::sqrt(6.)*std::pow(bohr,2.5))).matrix().cast<dcomplex>();
    }
    const double norm=s.dr*state.squaredNorm();
    if (!(norm>0)) throw std::runtime_error("Unresolved centered 211");
    return state*std::sqrt(input_norm/norm);
}

void run_case(const Options &o, const Config &c, const std::string &name)
{
    const FixedKeplerOrbit orbit(1+c.q,c.R,c.e,c.counter);
    const double period=orbit.period();
    const double dt=o.dt ? o.dt : (c.geometry.normalized ? period/c.geometry.startup_steps : c.dt);
    const double duration=o.periods*period, count=duration/dt;
    if (!std::isfinite(count) || count<=0 || count>=std::numeric_limits<int>::max()-1)
        throw std::invalid_argument("Invalid number of diagnostic steps");
    const int steps=std::max(1,static_cast<int>(std::ceil(count-8*std::numeric_limits<double>::epsilon()*count)));
    CrankNicolsonModes s(c.N,c.L,dt,c.epsilon,c.l_max,c.mu,c.M_c,c.q,c.a,c.R,c.counter,
        false,c.save_every,0,c.epsilon_smear,(o.output/"source.toml").string(),c.e);
    const double source_x=s.orbit.R1(), x=name=="C0" ? 0 : source_x;
    std::cout << name << " stationary primary at x=" << x << ", source period=" << period
              << ", fixed dt=" << dt << ", final time=" << duration << '\n';
    s.configure_isolated_black_hole(x);
    s.use_block_preconditioner=c.preconditioner=="radial";
    auto state=load_initial_modes_binary((o.output/"source.bin").string(),s.modes);
    const double raw_norm=s.dr*state.squaredNorm();
    if (!(raw_norm>0) || !std::isfinite(raw_norm)) throw std::runtime_error("Invalid input norm");
    if (name=="C0") state=make_centered_211(s,raw_norm);
    s.set_initial_state_full(state);
    s.initialize();
    if (s.V2.norm()!=0 || s.dVdt.norm()!=0 || s.absorption_sink_active() ||
        s.orbit.M1!=1 || s.orbit.M2!=0 || s.orbit.R!=std::abs(x) || s.orbit.theta_dot!=0)
        throw std::runtime_error("Failed stationary single-primary check");
    const auto case_dir=o.output/name;
    fs::create_directory(case_dir);
    fs::current_path(case_dir);
    cloud211::Measurements measurements(s,o.radii);
    s.run(1,o.radii.front());
    s.t.conservativeResize(steps+1);
    toml::array modes, radii;
    for (const auto &mode:s.modes) modes.push_back(toml::array{mode.l,mode.m});
    for (double radius:o.radii) radii.push_back(radius);
    const toml::table manifest{{"case",name},{"source_config",o.config.string()},
        {"source_input",o.input.string()},{"source_x",source_x},{"x",x},{"period",period},
        {"dt",dt},{"steps",steps},{"duration",duration},{"N",c.N},{"dr",s.dr},
        {"epsilon",c.epsilon},{"mu",c.mu},{"M_c",c.M_c},{"input_norm",raw_norm},
        {"control_integrator","static_cn"},{"M1",s.orbit.M1},{"M2",s.orbit.M2},
        {"ramp_periods",0.0},{"stationary",true},{"absorption",false},{"cap",false},{"mask",false},
        {"preconditioner",c.preconditioner},{"thread_limit",omp_get_thread_limit()},
        {"modes",std::move(modes)},{"radii",std::move(radii)}};
    std::ofstream metadata("manifest.json");
    metadata.exceptions(std::ios::badbit|std::ios::failbit);
    metadata << toml::json_formatter{manifest}; metadata.close();
    const Eigen::VectorXcd hpsi=s.H*s.psi;
    const double mean=s.psi.dot(hpsi).real()/s.psi.squaredNorm();
    const double variance=(hpsi-mean*s.psi).squaredNorm()/s.psi.squaredNorm();
    const double e211=-std::pow(c.mu,3)/8;
    const toml::table stats{{"Q0",s.dr*s.psi.squaredNorm()},{"mean_energy",mean},
        {"variance",variance},{"rms_residual",std::sqrt(variance)},{"E211",e211},
        {"analytical_residual",(hpsi-e211*s.psi).norm()/s.psi.norm()}};
    std::ofstream initial("initial.json");
    initial.exceptions(std::ios::badbit|std::ios::failbit);
    initial << toml::json_formatter{stats}; initial.close();
    cloud211::save_operator(s);
    cloud211::save_state("state_0.c128",s.psi);
    measurements.sample(s,s.psi,0,0);
    const int snapshot_every=std::max(1,steps/12);
    for (int i=1;i<=steps;++i) {
        const Eigen::VectorXcd before=s.psi;
        s.dt=i==steps ? duration-s.t[i-1] : dt;
        const double elapsed=s.dt;
        s.step(i);
        measurements.sample(s,before,i,elapsed);
        if (i%snapshot_every==0 || i==steps)
            cloud211::save_state("state_"+std::to_string(i)+".c128",s.psi);
    }
    measurements.finish();
    std::ofstream done("completed.txt");
    done.exceptions(std::ios::badbit|std::ios::failbit); done << steps << '\n'; done.close();
}
}

int main(int argc,char **argv) try {
    if (argc==2 && std::string(argv[1])=="--help") {
        std::cout << "Cloud211Diagnostic config.toml initial_modes.bin --output NEW_DIRECTORY "
                     "[--cases C0,D0|D1] [--periods 2] [--dt DT] [--radii R1,R2,...]\n";
        return 0;
    }
    auto o=options(argc,argv);
    const Config c=load_config(o.config.string());
    if (c.n!=2 || c.l!=1 || c.m!=1 || c.l_max<1 || c.mu<=0 || c.M_c<=0 ||
        c.N<3 || c.L<=0 || c.q<0)
        throw std::invalid_argument("Expected a positive-mass primary 211 config and valid grid");
    verify_input(o.input,c.N);
    if (c.threads>0) omp_set_num_threads(c.threads);
    Eigen::setNbThreads(omp_get_max_threads());
    if (o.radii.empty()) o.radii={c.r_ex};
    for (auto parent=o.output.parent_path(); !parent.empty(); parent=parent.parent_path()) {
        if (fs::exists(parent/"checkpoint.bin")) throw std::invalid_argument("Output inside a simulation archive");
        if (parent==parent.root_path()) break;
    }
    if (!fs::create_directories(o.output)) throw std::runtime_error("Output directory already exists");
    fs::copy_file(o.config,o.output/"source.toml");
    fs::copy_file(o.input,o.output/"source.bin");
    for (const auto &name:o.cases) run_case(o,c,name);
    return 0;
} catch (const std::exception &e) {
    std::cerr << "Cloud211Diagnostic: " << e.what() << '\n';
    return 1;
}
