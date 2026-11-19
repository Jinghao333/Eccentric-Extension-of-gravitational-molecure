#include "CrankNicolsonModes.hpp"
#include "TestUtils.hpp"
#include <omp.h>
#include <memory>

namespace fs = std::filesystem;
using test::require;
using test::close;
using Solver = CrankNicolsonModes;
namespace {

std::unique_ptr<Solver> fixture(double e, double eta, double ramp, bool counter,
                                double rp=35.71, double maximum=50.)
{
    GeometryTimestep geometry{true,eta,2,3,0.02};
    auto cn=std::make_unique<Solver>(25,24.,maximum,0.1,1,0.8,0.8,1e-5,0.,rp,
        counter,false,4,ramp,3.,"fixture.toml",e,FieldIntegrator::MidpointCN,geometry);
    cn->evolve_orbit=false;
    cn->initialize();
    return cn;
}

void orbit_schedules()
{
    for (double e : {0.,0.4,0.9,0.999}) {
        for (double eta : {0.5,0.25,0.125}) {
            const double maximum=e>0.8 ? 1e12 : 50.;
            auto cn=fixture(e,eta,2.,false,35.71,maximum);
            auto reverse=fixture(e,eta,2.,true,35.71,maximum);
            const double P=cn->initial_orbital_period();
            double t=2*P;
            OrbitState state=cn->evaluate_prescribed_orbit_at(cn->orbit,0.,t);
            OrbitState retro=reverse->evaluate_prescribed_orbit_at(reverse->orbit,0.,t);
            const double first=cn->propose_geometry_timestep(t,state,3*P);
            if (e>0) {
                auto apo=cn->evaluate_prescribed_orbit_at(state,t,t+P/2);
                require(first<=cn->propose_geometry_timestep(t+P/2,apo,3*P),"Pericenter underresolved");
            } else {
                close(first,std::min(maximum,eta/(2*cn->mean_motion())),1e-12,"Circular step");
            }
            int count=0;
            double smallest=maximum,largest=0.,ratio=0.;
            while (t<3*P) {
                const double h=cn->propose_geometry_timestep(t,state,3*P);
                close(h,reverse->propose_geometry_timestep(t,retro,3*P),1e-13,"Counter rotation");
                const auto end=cn->evaluate_prescribed_orbit_at(state,t,t+h);
                for (const auto &s : {state,cn->evaluate_prescribed_orbit_at(state,t,t+h/2),end})
                    ratio=std::max(ratio,h*(2*std::abs(s.theta_dot)+3*std::abs(s.R_dot/s.R))/eta);
                require(2*std::abs(end.theta-state.theta)<=eta*(1+1e-12),"Skipped angular sweep");
                require(t+h>t && t+h<=3*P && h<=maximum,"Invalid physical time");
                retro=reverse->evaluate_prescribed_orbit_at(retro,t,t+h);
                t+=h; state=end;
                smallest=std::min(smallest,h); largest=std::max(largest,h);
                require(++count<10000,"Excessive orbital samples");
            }
            require(ratio<=1+1e-14,"Sampled geometry budget violated");
            if (e==0.4)
                std::cout<<"Orbit eta="<<eta<<" steps="<<count<<" min="<<smallest
                         <<" max="<<largest<<" max(h*rate/eta)="<<ratio<<'\n';
        }
    }
}

void startup()
{
    for (bool switching : {false,true}) {
        auto cn=fixture(0.4,0.25,0.001,false);
        cn->ramp_secondary_potential=switching;
        double time=0, ramp=0.001*cn->initial_orbital_period();
        auto state=cn->orbit;
        int count=0;
        while (time<ramp) {
            double h=cn->propose_geometry_timestep(time,state,2*ramp);
            if (switching) require(h<=0.02*ramp/1.5*(1+1e-14),"Companion weight budget");
            state=cn->evaluate_prescribed_orbit_at(state,time,time+h);
            time+=h; ++count;
            require(time<=ramp,"Straddled ramp boundary");
        }
        require(time==ramp && (switching ? count>=75 : count==1),"Ramp landing or unwanted weight bound");
        close(state.phase,cn->mean_motion()*ramp/2,1e-13,"Integrated startup clock");
    }
}

void radial_only_alias()
{
    // A two-period trial starting at apoapsis samples apoapsis three times.
    // Endpoint log-radius differences alone would miss both close passages.
    GeometryTimestep g{true,0.25,0,3,0.02};
    Solver cn(9,20.,1e12,0.1,0,0.8,0.8,1e-5,0.,4.,false,false,1,
              0.,2.,"unused",0.999,FieldIntegrator::MidpointCN,g);
    cn.evolve_orbit=false;
    cn.initialize();
    double P=cn.initial_orbital_period(), time=P/2;
    auto apo=cn.evaluate_prescribed_orbit_at(cn.orbit,0.,time);
    double h=cn.propose_geometry_timestep(time,apo,time+2*P);
    require(h>0 && h<P/2,"Radial-only step skipped periapsis");
}

void config_validation()
{
    test::ScratchDir dir("config");
    auto load=[](std::string sim,std::string extra=""){
        std::ofstream("input.toml")<<"[simulation]\ntimestep_rule='orbital_geometry'\n"<<sim<<"\n"<<extra;
        return load_config("input.toml");
    };
    std::string valid="integrator='midpoint_cn'\nnum_periods=1\n";
    require(load(valid).dt==50., "Default geometric timestep ceiling");
    require(load(valid+"l_max=0").geometry.azimuthal_order==0,"Monopole default");
    for (auto bad : {"num_steps=4","dt='wrong'","geometry_tolerance=nan","ramp_weight_tolerance=0.0",
                     "effective_radial_order=-1","effective_azimuthal_order=20","dt=inf",
                     "evolve_orbit=true","[mask]\nL_factor=0.7"}) {
        try { load(valid+bad); } catch (const std::invalid_argument &) { continue; }
        throw std::runtime_error("Invalid config accepted: "+std::string(bad));
    }
    for (auto bad : {"integrator='midpoint_cn'","num_periods=2"}) {
        try { load(bad); } catch (const std::invalid_argument &) { continue; }
        throw std::runtime_error("Missing required settings accepted");
    }
}

std::string bytes(const fs::path &p) {
    std::ifstream in(p,std::ios::binary); require(bool(in),"Missing "+p.string());
    return {std::istreambuf_iterator<char>(in),{}};
}

struct Result { fs::path folder; Eigen::VectorXcd field; Eigen::VectorXd times; };

Result field_run(double eta,int split)
{
    test::ScratchDir dir("field");
    std::ofstream("fixture.toml")<<"# geometry workflow\n";
    auto make=[&]() {
        auto cn=fixture(0.4,eta,0.03,true,3.25,50.);
        cn->build_cap(0.8,0.03,2);
        Eigen::VectorXcd initial=Eigen::VectorXcd::Zero(cn->Nr*cn->n_modes);
        for (int j=0;j<cn->n_modes;++j)
            for (int i=1;i<cn->Nr-1;++i)
                initial[j*cn->Nr+i]=std::polar(std::sin(constants::pi*i/(cn->Nr-1))/(1+j),0.2*i);
        cn->set_initial_state_full(initial.normalized());
        return cn;
    };
    auto cn=make();
    // Pause before the first step, then advance the ordinary step routine
    // to an interior checkpoint without shortening the original time target.
    volatile std::sig_atomic_t stop=1;
    if (split) cn->set_stop_flag(&stop);
    constexpr double final_time=25.13;
    cn->run_until(final_time,12.);
    auto folder=test::simulation_folder(dir.path());
    if (split) {
        cn->t.conservativeResize(split+1);
        for (int i=1;i<=split;++i) cn->step(i);
        cn.reset();
        // I/O after the durable checkpoint must not survive resume.
        std::ofstream(folder/"dt_step.txt",std::ios::app)<<"99\n1.2";
        cn=make();
        cn->resume_until(folder.string(),final_time,12.);
    }
    require(cn->t[cn->t.size()-1]==final_time,"Did not hit physical target");
    auto rows=test::read_rows(folder/"dt_step.txt");
    require(rows.size()+1==size_t(cn->t.size()),"Accepted duration count");
    for (size_t k=0;k<rows.size();++k)
        require(cn->t[k]+rows[k][0]==cn->t[k+1],"Accepted time history inconsistent");
    std::cout<<"Field eta="<<eta<<" steps="<<rows.size()
             <<" final solve residual="<<cn->last_solver_residual<<'\n';
    return {folder,cn->psi,cn->t};
}

void fields()
{
    auto coarse=field_run(0.5,0), med=field_run(0.25,0), fine=field_run(0.125,0);
    double d1=(coarse.field-med.field).norm(),d2=(med.field-fine.field).norm();
    std::cout<<"Complex-field differences "<<d1<<" "<<d2<<" ratio="<<d1/d2<<'\n';
    require(d2>0 && d2<d1,"Field failed to converge");
    for (int split : {4,80}) {
        auto resumed=field_run(0.25,split);
        require((med.times.array()==resumed.times.array()).all(),"Resumed time sequence");
        for (const auto &entry:fs::recursive_directory_iterator(med.folder)) {
            if (!entry.is_regular_file()) continue;
            const auto name=fs::relative(entry.path(),med.folder);
            require(bytes(entry.path())==bytes(resumed.folder/name),"Restart differs: "+name.string());
        }
    }
    auto wrong=fixture(0.4,0.125,0.03,true,3.25,50.);
    wrong->build_cap(0.8,0.03,2);
    const auto checkpoint=bytes(med.folder/"checkpoint.bin");
    try { wrong->resume_until(med.folder.string(),50.,12.); }
    catch (const std::runtime_error &err) {
        require(std::string(err.what()).find("physics")!=std::string::npos,"Wrong compatibility error");
        require(bytes(med.folder/"checkpoint.bin")==checkpoint,"Rejected restart modified checkpoint");
        auto extend=fixture(0.4,0.25,0.03,true,3.25,50.);
        extend->build_cap(0.8,0.03,2);
        extend->resume_until(med.folder.string(),27.,12.);
        require(extend->t.tail(1)[0]==27. && extend->t[med.times.size()-1]==25.13,
                "Extension lost original shortened endpoint");
        return;
    }
    throw std::runtime_error("Changed scheduler accepted on resume");
}
}

int main() {
    omp_set_num_threads(1);
    using Test=std::pair<const char *,void(*)()>;
    return test::run_all(Test{"orbit schedules",orbit_schedules},Test{"startup",startup},
        Test{"radial alias",radial_only_alias},
        Test{"configuration",config_validation},Test{"field and restart",fields});
}
