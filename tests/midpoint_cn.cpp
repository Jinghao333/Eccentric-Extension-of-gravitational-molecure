// Compare the propagation and measurement operators independently, then
// refine only time on the same radial grid and compare exact restart outputs.
#include "CrankNicolsonModes.hpp"
#include "TestUtils.hpp"
#include <omp.h>
#include <memory>

namespace fs = std::filesystem;
using test::close;
using test::require;
using Solver = CrankNicolsonModes;
namespace {
constexpr double extraction = 12.0;
struct Case {
    double e, ramp_duration;
    bool counter = false, absorption = false, cap = false;
};

std::unique_ptr<Solver> fixture(Case c, double h,
                                FieldIntegrator method = FieldIntegrator::MidpointCN)
{
    const double period = 2 * constants::pi * std::sqrt(std::pow(3.25/(1-c.e), 3)/1.4);
    std::ofstream("fixture.toml") << "[simulation]\nintegrator = \""
        << (method == FieldIntegrator::MidpointCN ? "midpoint_cn" : "legacy_endpoint") << "\"\n";
    auto cn = std::make_unique<Solver>(25, 24., h, 0.17, 1, 0.8, 0.8, 0.4, 0.,
        3.25, c.counter, c.absorption, 4, c.ramp_duration/period, 3., "fixture.toml", c.e, method);
    cn->evolve_orbit = false;
    Eigen::VectorXcd initial = Eigen::VectorXcd::Zero(cn->Nr * cn->n_modes);
    for (int j = 0; j < cn->n_modes; ++j)
        for (int i = 1; i < cn->Nr-1; ++i)
            initial[j*cn->Nr+i] = std::polar(std::sin(constants::pi*i/(cn->Nr-1))/(1+j), 0.27*i+0.31*j);
    cn->set_initial_state_full(initial.normalized());
    cn->initialize();
    if (c.cap) cn->build_cap(0.7, 0.12, 2);
    return cn;
}

// Seed the desired anomaly in the legacy integrator. Its first step has
// zero clock rate for T>0, so it builds exactly the operator at that anomaly
// and the requested endpoint time, independently of midpoint step().
std::unique_ptr<Solver> reference(Case c, double time)
{
    test::ScratchDir dir("legacy_operator");
    auto ref = fixture(c, time, FieldIntegrator::LegacyEndpoint);
    const double x = std::min(time/c.ramp_duration, 1.);
    ref->orbit.phase = ref->mean_motion() *
        (c.ramp_duration*std::pow(x,3)*(1.-0.5*x) + std::max(0., time-c.ramp_duration));
    ref->run(2, extraction);
    return ref;
}

void operators_and_balance(Case c)
{
    test::ScratchDir dir("operators");
    const double h = 1.;
    auto cn = fixture(c, h);
    const auto before = cn->psi;
    cn->run(2, extraction);
    auto mid = reference(c, h/2);
    auto end = reference(c, h);
    const auto expected_A = (cn->Id + I*(h/2)*mid->H).eval();
    const auto expected_B = (cn->Id - I*(h/2)*mid->H).eval();
    require((cn->A-expected_A).norm()<1e-13 && (cn->B-expected_B).norm()<1e-13,
            "A/B did not use the midpoint Hamiltonian");
    require((end->H-mid->H).norm()>1e-4, "Fixture did not resolve the two operator times");
    require((cn->H-end->H).norm()<1e-13 && (cn->dVdt-end->dVdt).norm()<1e-13 &&
            (cn->dVdt_switch-end->dVdt_switch).norm()<1e-13, "Endpoint operators differ");
    close(cn->orbit.phase,end->orbit.phase,1e-15,"Committed mean anomaly");
    close(cn->orbit.theta,end->orbit.theta,1e-15,"Committed azimuth");
    require((cn->apply_sink(before)-end->apply_sink(before)).norm()<1e-13,"Endpoint sink differs");
    Eigen::VectorXcd avg = (cn->psi+before)/2.;
    const double residual = (cn->psi-before-h*mid->schrodinger_rhs(avg)).norm();
    require(residual < 1e-9, "Field did not solve the full midpoint CN equation");

    double loss = h*avg.dot(mid->apply_sink(avg)).real();
    if (c.cap) loss -= 2*h*avg.dot(mid->CAP*avg).imag();
    const double balance = before.squaredNorm()-cn->psi.squaredNorm()-loss;
    require(std::abs(balance)<2e-10, "Midpoint norm balance failed");
    if (c.absorption || c.cap) require(loss>1e-8, "Dissipative balance fixture has negligible loss");

    const auto folder = test::simulation_folder(dir.path());
    const int cut = static_cast<int>(extraction/cn->dr);
    auto expectation = [&](const SparseComplex &op) {
        const Eigen::VectorXcd applied = op*cn->psi;
        dcomplex sum = 0.;
        for (int j=0; j<cn->n_modes; ++j)
            sum += cn->psi.segment(j*cn->Nr,cut+1).dot(applied.segment(j*cn->Nr,cut+1));
        return cn->dr*sum.real();
    };
    auto scalar = [&](const char *name) { return test::read_rows(folder/name).at(0).at(0); };
    close(scalar("t.txt"), h, 0., "Diagnostic timestamp");
    close(scalar("energy.txt"),expectation(end->H),1e-13,"Endpoint energy");
    close(scalar("p_bin_inst.txt"),expectation(end->dVdt),1e-13,"Endpoint orbital power");
    close(scalar("p_ramp_inst.txt"),expectation(end->dVdt_switch),1e-13,"Endpoint switching power");
    const auto flux = test::read_rows(folder/"energy_flux.txt").at(0);
    const auto rhs = end->schrodinger_rhs(cn->psi);
    for (int j=0; j<cn->n_modes; ++j) {
        const int k=j*cn->Nr+cut;
        const dcomplex radial=(cn->psi[k+1]-cn->psi[k-1])/(2.*cn->dr)-cn->psi[k]/cn->r[cut];
        close(flux.at(j),(std::conj(rhs[k])*radial).real()/cn->mu,1e-13,"Endpoint energy flux");
    }
    std::cout << "e=" << c.e << " absorption=" << c.absorption << " CAP=" << c.cap
              << " CN residual=" << residual << " norm balance error=" << balance << '\n';
}

void convergence()
{
    for (double e : {0., 0.4}) {
        auto field = [&](int steps) -> Eigen::VectorXcd {
            test::ScratchDir dir("refinement");
            auto cn = fixture({e,1.23,e>0.}, 4./steps);
            cn->run(steps+1,extraction);
            close(cn->t[steps],4.,0.,"Common final time");
            return cn->psi;
        };
        const auto coarse=field(32), medium=field(64), fine=field(128);
        const double d1=(coarse-medium).norm(), d2=(medium-fine).norm();
        std::cout << "Refinement e=" << e << " differences=" << d1 << ", " << d2
                  << " ratio=" << d1/d2 << '\n';
        require(d1/d2>3.8 && d1/d2<4.2,"Expected second-order complex field convergence");
    }
}

std::string bytes(const fs::path &p)
{
    std::ifstream in(p,std::ios::binary);
    require(bool(in),"Missing "+p.string());
    return {std::istreambuf_iterator<char>(in),{}};
}

void restart()
{
    auto run = [](Case c, int split) {
        test::ScratchDir dir("restart");
        auto cn=fixture(c,0.5);
        cn->run(split ? split+1 : 13,extraction);
        auto folder=test::simulation_folder(dir.path());
        if (split) {
            cn.reset();
            cn=fixture(c,0.5);
            cn->resume(folder.string(),13,extraction);
        }
        return folder;
    };
    for (double e : {0.,0.4}) {
        const Case c{e,2.2,e>0.,true,true};
        const auto full=run(c,0);
        for (int split : {2,7}) {
            const auto resumed=run(c,split);
            for (const auto &entry : fs::recursive_directory_iterator(full)) {
                if (!entry.is_regular_file()) continue;
                const auto name=fs::relative(entry.path(),full);
                require(bytes(entry.path())==bytes(resumed/name),"Restart differs: "+name.string());
            }
        }
    }
}

void reject_integrator_change()
{
    for (auto from : {FieldIntegrator::MidpointCN, FieldIntegrator::LegacyEndpoint}) {
        test::ScratchDir dir("mismatch");
        Case c{0.,0.};
        fixture(c,0.5,from)->run(3,extraction);
        const auto folder=test::simulation_folder(dir.path());
        const auto saved=bytes(folder/"checkpoint.bin");
        auto other=fixture(c,0.5,from==FieldIntegrator::MidpointCN ?
            FieldIntegrator::LegacyEndpoint : FieldIntegrator::MidpointCN);
        try { other->resume(folder.string(),4,extraction); }
        catch (const std::runtime_error &e) {
            require(std::string(e.what()).find("integrator")!=std::string::npos,"Wrong resume error");
            require(bytes(folder/"checkpoint.bin")==saved,"Mismatch changed the checkpoint");
            continue;
        }
        throw std::runtime_error("Integrator change accepted on resume");
    }
}
} // namespace

int main()
{
    omp_set_num_threads(1);
    Eigen::setNbThreads(1);
    using Test = std::pair<const char *, void(*)()>;
    return test::run_all(
        Test{"circular ramp crossing", [] { operators_and_balance({0.,0.75}); }},
        Test{"counter-rotating eccentric startup", [] { operators_and_balance({0.4,1.6,true}); }},
        Test{"midpoint absorption", [] { operators_and_balance({0.4,0.1,true,true}); }},
        Test{"midpoint CAP", [] { operators_and_balance({0.,0.75,false,false,true}); }},
        Test{"temporal convergence",convergence}, Test{"exact restart",restart},
        Test{"integrator mismatch",reject_integrator_change});
}
