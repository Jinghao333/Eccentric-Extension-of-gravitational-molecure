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
                                double rp=35.71, double maximum=50., int count=0, int startup_count=120)
{
    GeometryTimestep geometry{true,eta,2,3,0.02};
    geometry.normalized = count>0;
    geometry.orbit_steps = count;
    geometry.startup_steps = startup_count;
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

void matched_step_counts()
{
    for (double fixed_dt : {50.,25.}) {
        const double eta = fixed_dt == 50. ? 0.326 : 0.1594;
        auto cn = fixture(0.4,eta,2.,false,35.71,2*fixed_dt);
        const double period = cn->initial_orbital_period();
        double t = 2*period;
        auto state = cn->evaluate_prescribed_orbit_at(cn->orbit,0.,t);
        int count = 0;
        while (t < 3*period) {
            const double h = cn->propose_geometry_timestep(t,state,3*period);
            state = cn->evaluate_prescribed_orbit_at(state,t,t+h);
            t += h;
            ++count;
        }
        require(count == static_cast<int>(std::ceil(period/fixed_dt)), "Matched orbit step count");
        std::cout << "Fixed dt=" << fixed_dt << " and geometry eta=" << eta
                  << " steps per calibration orbit=" << count << '\n';
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
    auto load=[](std::string sim,std::string extra="", std::string rule="orbital_geometry"){
        std::ofstream("input.toml")<<"[simulation]\ntimestep_rule='"<<rule<<"'\n"<<sim<<"\n"<<extra;
        return load_config("input.toml");
    };
    std::string valid="integrator='midpoint_cn'\nnum_periods=1\n";
    auto norm=load(valid+"startup_steps_per_period=7\norbit_steps_per_period=33","","normalized_geometry");
    require(norm.geometry.normalized && norm.geometry.startup_steps==7 && norm.geometry.orbit_steps==33,
            "Normalized config");
    for (auto bad : {"startup_steps_per_period=0","orbit_steps_per_period=1.5","num_steps=1",
                     "evolve_orbit=true","[mask]\nL_factor=0.7"}) {
        try { load(valid+bad,"","normalized_geometry"); } catch (const std::invalid_argument &) { continue; }
        throw std::runtime_error("Invalid normalized setting accepted");
    }
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

Result field_run(double eta,int split, int count=0, bool stop_at_ramp=false)
{
    test::ScratchDir dir("field");
    std::ofstream("fixture.toml")<<"# geometry workflow\n";
    const double ramp_periods=count ? 0.06 : 0.03;
    auto make=[&]() {
        auto cn=fixture(0.4,eta,ramp_periods,true,3.25,50.,count);
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
    const double final_time=stop_at_ramp ? ramp_periods*cn->initial_orbital_period() : 25.13;
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
        require(count ? cn->t[k+1]-cn->t[k]==rows[k][0] : cn->t[k]+rows[k][0]==cn->t[k+1],
                "Accepted time history inconsistent");
    std::cout<<"Field eta="<<eta<<" steps="<<rows.size()
             <<" final solve residual="<<cn->last_solver_residual<<'\n';
    return {folder,cn->psi,cn->t};
}

void normalized_orbits()
{
    for (double e : {0.,0.4,0.9}) for (bool retro : {false,true})
        for (int count : {1,60,120,240}) {
            const double ramps = count==1 ? 0.375 : 2.;
            auto cn=fixture(e,0.25,ramps,retro,35.71,1.,count);
            const double P=cn->initial_orbital_period(), ramp=ramps*P, pi=constants::pi;
            double time=ramp;
            auto start=cn->evaluate_prescribed_orbit_at(cn->orbit,0.,time);
            double minimum=P, maximum=0.;
            for (int cycle=1; cycle<=3; ++cycle) {
                for (int k=0; k<count; ++k) {
                    const double h=cn->propose_normalized_timestep(time,(ramps+3)*P);
                    const auto end=cn->evaluate_prescribed_orbit_at(start,time,time+h);
                    // Independent radial primitive on successive half revolutions.
                    const auto radial=[&](const OrbitState &s) {
                        const int half=int(std::floor(s.phase/pi));
                        const double span=std::log((1+e)/(1-e));
                        const double offset=std::log(s.R/35.71);
                        return half*span+(half%2 ? span-offset : offset);
                    };
                    const double dV=2*std::abs(end.theta-start.theta)+3*(radial(end)-radial(start));
                    close((h/P+dV/(4*pi+6*std::log((1+e)/(1-e))))/2,
                          1./count,3e-12,"normalized coordinate increment");
                    if (!e) close(h/P,1./count,2e-15,"uniform circular limit");
                    require(h>0 && std::isfinite(h),"invalid normalized h");
                    minimum=std::min(minimum,h); maximum=std::max(maximum,h);
                    time+=h; start=end;
                }
                require(time==(ramps+cycle)*P,"wrong number of normalized steps in period");
            }
            if(e==0.4 && !retro) std::cout<<"Normalized orbit N="<<count
                <<" min="<<minimum<<" max="<<maximum<<'\n';
        }
    auto a=fixture(0.4,0.25,0.,false,35.71,1.,120,17);
    auto b=fixture(0.4,0.25,0.,false,35.71,50.,120,300);
    double t=0., P=a->initial_orbital_period();
    for (int i=0; i<120; ++i) {
        double h=a->propose_normalized_timestep(t,P);
        require(h==b->propose_normalized_timestep(t,P),"startup resolution changed orbital grid");
        t+=h;
    }
    require(t==P,"zero-ramp count");
    auto sparse=fixture(0.4,.25,.375,false,35.71,1.,60,9);
    auto dense=fixture(0.4,.25,.375,false,35.71,1.,60,120);
    const double seam=.375*P, target=1.375*P;
    for(auto *solver : {sparse.get(),dense.get()}) {
        t=0.;
        while(t<seam) {
            t+=solver->propose_normalized_timestep(t,target);
            require(t<=seam,"overshot fractional ramp");
        }
    }
    for(int k=0;k<60;++k) {
        double h=sparse->propose_normalized_timestep(t,target);
        require(h==dense->propose_normalized_timestep(t,target),"post-startup grid changed");
        t+=h;
    }
    require(t==target,"fractional ramp period count");
}

void normalized_fields()
{
    auto at_ramp=field_run(.25,0,60,true);
    auto same=field_run(.25,0,240,true);
    require(at_ramp.times==same.times && at_ramp.field==same.field,"startup field or mesh differs");
    for (auto name : {"phase_orb.txt","theta_orb.txt","R_orb.txt","R_dot_orb.txt","Omega_orb.txt"})
        require(bytes(at_ramp.folder/name)==bytes(same.folder/name),"startup orbit differs");
    auto a=field_run(.25,0,60), b=field_run(.25,0,120), c=field_run(.25,0,240);
    const double coarse=(a.field-b.field).norm(), fine=(b.field-c.field).norm();
    std::cout<<"Normalized field differences "<<coarse<<" "<<fine<<" ratio="<<coarse/fine<<'\n';
    require(fine>0 && fine<coarse,"normalized field refinement");
    for (int split : {4,12}) {
        auto resumed=field_run(.25,split,120);
        require(b.field==resumed.field && b.times==resumed.times,"normalized restart state");
        for (const auto &entry:fs::recursive_directory_iterator(b.folder))
            if(entry.is_regular_file()) {
                const auto name=fs::relative(entry.path(),b.folder);
                require(bytes(entry.path())==bytes(resumed.folder/name),"normalized restart: "+name.string());
            }
    }
    auto changed=fixture(0.4,.25,.06,true,3.25,50.,121);
    changed->build_cap(.8,.03,2);
    try { changed->resume_until(b.folder.string(),30.,12.); }
    catch (const std::runtime_error &error) {
        require(std::string(error.what()).find("physics")!=std::string::npos,"wrong fingerprint failure");
        return;
    }
    throw std::runtime_error("Changed normalized resolution accepted");
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
    return test::run_all(Test{"orbit schedules",orbit_schedules},Test{"matched step counts",matched_step_counts},Test{"startup",startup},
        Test{"radial alias",radial_only_alias},Test{"normalized orbits",normalized_orbits},
        Test{"normalized fields",normalized_fields},
        Test{"configuration",config_validation},Test{"field and restart",fields});
}
