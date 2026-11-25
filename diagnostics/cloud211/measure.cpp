#include "measure.hpp"
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <set>

namespace cloud211 {
namespace {
std::ofstream output(const std::filesystem::path &name, bool binary = false)
{
    std::ofstream out;
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out.open(name, binary ? std::ios::binary : std::ios::out);
    out << std::setprecision(17);
    return out;
}
}

void save_state(const std::filesystem::path &filename, const Eigen::VectorXcd &state)
{
    auto out = output(filename, true);
    out.write(reinterpret_cast<const char *>(state.data()), state.size() * sizeof(dcomplex));
    out.close();
}

// Binary COO records on the dynamical (interior) subspace, in mode-major order.
void save_operator(const CrankNicolsonModes &s)
{
    auto out = output("interior.coo", true);
    struct Entry { std::int32_t row, col; double re, im; };
    static_assert(sizeof(Entry) == 24);
    for (int row = 0; row < s.H.rows(); ++row) {
        const int i = row % s.Nr;
        if (i == 0 || i == s.Nr-1) continue;
        for (SparseComplex::InnerIterator it(s.H,row); it; ++it) {
            const int col = static_cast<int>(it.col()), j = col % s.Nr;
            if (j == 0 || j == s.Nr-1 || it.value() == dcomplex{}) continue;
            const Entry entry{(row/s.Nr)*(s.Nr-2)+i-1, (col/s.Nr)*(s.Nr-2)+j-1,
                              it.value().real(), it.value().imag()};
            out.write(reinterpret_cast<const char *>(&entry), sizeof(entry));
        }
    }
    out.close();
}

Measurements::Measurements(const CrankNicolsonModes &s, const std::vector<double> &radii)
    : flux(output("flux.csv")), balance(output("balance.csv"))
{
    std::set<int> nodes;
    for (double radius : radii) {
        if (!(radius >= s.r[1] && radius < s.r[s.Nr-1]))
            throw std::invalid_argument("Extraction radius must select an interior node");
        const int k = static_cast<int>(std::upper_bound(s.r.data(),s.r.data()+s.Nr,radius)-s.r.data())-1;
        if (!nodes.insert(k).second) throw std::invalid_argument("Duplicate extraction node");
        const auto zeros = Eigen::VectorXd::Zero(s.n_modes);
        surfaces.push_back({radius,k,0.,zeros,zeros,zeros});
    }
    flux << "step,t,r_requested,r_node,r_face,mode,l,m,JQ,JE,Jexcess,deltaQ_node,t_mid,JQ_face,deltaQ_face\n";
    balance << "step,t,r_node,Q_full,Q_inside,deltaQ_face,continuity,solve_residual\n";
}

void Measurements::sample(const CrankNicolsonModes &s, const Eigen::VectorXcd &before,
                          int step, double dt)
{
    const auto rhs = s.schrodinger_rhs(s.psi);
    const double time = s.t[step], full = s.dr*s.psi.squaredNorm();
    const double binding = -std::pow(s.mu,3)/8.;
    for (auto &surface : surfaces) {
        const int k = surface.node;
        double inside = 0, total_q = 0, total_e = 0, total_face = 0;
        const auto write = [&](const char *label, int l, int m, double q, double e,
                               double cumulative, double face, double cumulative_face) {
            flux << step << ',' << time << ',' << surface.requested << ',' << s.r[k] << ','
                 << (s.r[k]+s.r[k+1])/2 << ',' << label << ',' << l << ',' << m << ','
                 << q << ',' << e << ',' << e-binding*q << ',' << cumulative << ',';
            if (step) flux << time-dt/2 << ',' << face;
            else flux << ',';
            flux << ',' << cumulative_face << '\n';
        };
        for (int j=0; j<s.n_modes; ++j) {
            const int start = j*s.Nr, p = start+k;
            const auto u = s.psi[p];
            const auto derivative = (s.psi[p+1]-s.psi[p-1])/(2*s.dr)-u/s.r[k];
            const double q = std::imag(std::conj(u)*derivative)/s.mu;
            const double e = -std::real(std::conj(rhs[p])*derivative)/s.mu;
            double face = 0;
            if (step) {
                const auto left = (before[p]+u)*0.5, right = (before[p+1]+s.psi[p+1])*0.5;
                face = std::imag(std::conj(left)*right)/(s.mu*s.dr);
                surface.integrated_face[j] += dt*face;
                surface.integrated_node[j] += dt*(q+surface.last_current[j])*0.5;
            }
            surface.last_current[j] = q;
            inside += s.dr*s.psi.segment(start,k+1).squaredNorm();
            total_q += q;
            total_e += e;
            total_face += face;
            write("mode",s.modes[j].l,s.modes[j].m,q,e,surface.integrated_node[j],
                  face,surface.integrated_face[j]);
        }
        const double transported = surface.integrated_face.sum();
        write("total",-1,0,total_q,total_e,surface.integrated_node.sum(),total_face,transported);
        if (!step) surface.initial_inside = inside;
        balance << step << ',' << time << ',' << s.r[k] << ',' << full << ',' << inside << ','
                << transported << ',' << (inside-surface.initial_inside)+transported << ','
                << s.last_solver_residual << '\n';
    }
}
}
