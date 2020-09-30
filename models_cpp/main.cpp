// Author: Francesco Regazzoni - MOX, Politecnico di Milano
// Email:  francesco.regazzoni@polimi.it
// Date:   2020

#include <cmath>
#include <map>
#include <memory>
#include <vector>

#include "model_RDQ18.hpp"
#include "model_RDQ20_MF.hpp"
#include "model_RDQ20_SE.hpp"

void print_help() {
  std::cout << "Usage:" << std::endl;
  std::cout << "   run_model [-h] <MODEL_NAME>" << std::endl;
  std::cout << "Positional arguments:" << std::endl;
  std::cout
      << "   <MODEL_NAME>   choice of the model (RDQ18, RDQ20-MF, RDQ20-SE)"
      << std::endl;
  std::cout << "Optional arguments:" << std::endl;
  std::cout << "   -h, --help     print this help and exit" << std::endl;
}

int main(int argc, char **argv) {

  // Model choice
  std::string model_name;
  if (argc == 2)
    model_name = argv[1];
  else {
    std::cout << "Wrong number of input arguments." << std::endl << std::endl;
    print_help();
    return 0;
  }

  std::unique_ptr<sarcomere> model;
  if (model_name == "RDQ20-MF")
    model = std::make_unique<model_RDQ20_MF>();
  else if (model_name == "RDQ20-SE")
    model = std::make_unique<model_RDQ20_SE>();
  else if (model_name == "RDQ18")
    model = std::make_unique<model_RDQ18>();
  else if (model_name == "-h" || model_name == "--help") {
    print_help();
    return 0;
  } else {
    std::cout << "Unknown model " << model_name << "." << std::endl
              << std::endl;
    print_help();
    return 0;
  }

  // Time interval
  double Tmax = .6; // [s]

  // Calcium transient
  double c0 = .1;    // [micro M]
  double cmax = 0.9; // [micro M]
  double tau1 = .02; // [s]
  double tau2 = .05; // [s]
  double t0 = 0.01;  // [s]
  double beta = std::pow(tau1 / tau2, -1 / (tau1 / tau2 - 1)) -
                std::pow(tau1 / tau2, -1 / (1 - tau2 / tau1));
  auto Ca = [&](double t) {
    if (t <= t0)
      return c0;
    else
      return c0 + ((cmax - c0) / beta *
                   (std::exp(-(t - t0) / tau1) - std::exp(-(t - t0) / tau2)));
  };

  // SL transient
  double SL0 = 2.2;       // [micro m]
  double SL1 = SL0 * .97; // [micro m]
  double SLt0 = .05;      // [s]
  double SLt1 = .35;      // [s]
  double SLtau0 = .05;    // [s]
  double SLtau1 = .02;    // [s]

  auto SL = [&](double t) {
    return SL0 +
           (SL1 - SL0) * (std::max(0.0, 1.0 - std::exp((SLt0 - t) / SLtau0)) -
                          std::max(0.0, 1.0 - std::exp((SLt1 - t) / SLtau1)));
  };

  // Simulation
  std::map<std::string, std::vector<double>> results =
      model->solve(Ca, SL, Tmax);

  // Writing results to CSV file
  model->write_csv(results, "output_" + model_name + ".csv");

  return 0;
}