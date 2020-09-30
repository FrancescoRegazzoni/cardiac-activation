// Author: Francesco Regazzoni - MOX, Politecnico di Milano
// Email:  francesco.regazzoni@polimi.it
// Date:   2020

#include <chrono>
#include <fstream>
#include <iomanip>

#include "sarcomere.hpp"

std::map<std::string, std::vector<double>>
sarcomere::solve(const std::function<double(const double &time)> &Ca,
                 const std::function<double(const double &time)> &SL,
                 const double &Tmax, const double &time_step) {
  double calcium;
  double sarcomere_length;
  double dSL_dt;

  std::map<std::string, std::vector<double>> results;

  double time = 0.0;
  std::vector<double> state(initial_state);

  std::cout << model_name << " model. Computing... " << std::flush;
  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();
  // Time loop
  while (time <= Tmax + 1e-10) {
    // Evaluate inputs
    calcium = Ca(time);
    sarcomere_length = SL(time);
    dSL_dt = (SL(time) - SL(time - time_step)) / time_step;

    // Update state
    solve_time_step(state, calcium, sarcomere_length, dSL_dt, time_step);

    // Store the results
    results["time"].push_back(time);
    results["Ca"].push_back(calcium);
    results["SL"].push_back(sarcomere_length);
    results["dSL_dt"].push_back(dSL_dt);
    results["Ta"].push_back(get_active_tension(state, sarcomere_length));

    // Increase time
    time += time_step;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  double duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
          .count() *
      1e-6;
  std::cout << "done. Time elapsed: " << duration << " s" << std::endl;

  return results;
}

void sarcomere::write_csv(std::map<std::string, std::vector<double>> results,
                          std::string file_name) {
  std::ofstream csvfile;
  csvfile.open(file_name);
  csvfile << "t,Ca,SL,dSL_dt,Ta" << std::endl;
  csvfile << std::scientific << std::setprecision(16);
  for (unsigned int i = 0; i < results["time"].size(); ++i)
    csvfile << results["time"][i] << "," << results["Ca"][i] << ","
            << results["SL"][i] << "," << results["dSL_dt"][i] << ","
            << results["Ta"][i] << std::endl;

  csvfile.close();
}