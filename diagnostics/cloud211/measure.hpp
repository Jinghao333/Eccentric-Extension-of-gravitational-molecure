#pragma once

#include "CrankNicolsonModes.hpp"
#include <filesystem>

namespace cloud211 {

struct Surface {
    double requested;
    int node;
    double initial_inside;
    Eigen::VectorXd last_current, integrated_node, integrated_face;
};

class Measurements {
public:
    Measurements(const CrankNicolsonModes &solver, const std::vector<double> &radii);
    void sample(const CrankNicolsonModes &solver, const Eigen::VectorXcd &before,
                int step, double dt);
    void finish() { flux.close(); balance.close(); }
private:
    std::vector<Surface> surfaces;
    std::ofstream flux, balance;
};

void save_operator(const CrankNicolsonModes &solver);
void save_state(const std::filesystem::path &filename, const Eigen::VectorXcd &state);
}
