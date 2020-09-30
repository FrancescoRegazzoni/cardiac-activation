// Author: Francesco Regazzoni - MOX, Politecnico di Milano
// Email:  francesco.regazzoni@polimi.it
// Date:   2020

#ifndef SARCOMERE_HPP
#define SARCOMERE_HPP

#include <functional>
#include <iostream>
#include <map>
#include <vector>

class sarcomere {
  // Abstract class representing a generic sarcomere model.

public:
  // Default constructor.
  sarcomere(std::string model_name_) : model_name(model_name_){};

  // Perform a simulation with the model.
  // Inputs:
  //    Ca         time evolution of calcium concentration [micro M]
  //    SL         time evolution of sarcomere length [micro m]
  //    Tmax       total simulation time [s]
  //    time_step  time step length [s]
  // Outputs:
  //    results    map containing the results
  std::map<std::string, std::vector<double>>
  solve(const std::function<double(const double &time)> &Ca,
        const std::function<double(const double &time)> &SL, const double &Tmax,
        const double &time_step = 1e-3);

  // Write the results of a simulation to a csv file.
  void write_csv(std::map<std::string, std::vector<double>> results,
                 std::string file_name);

protected:
  size_t n_variables;                // Total number of variables.
  std::vector<double> initial_state; // Initial state.
  std::string model_name;            // Model name.

  // Solve one time step, updating the state passed as reference.
  // Inputs:
  //    state             state variables
  //    calcium           intracellular calcium concentration [micro M]
  //    sarcomere_length  sarcomere length [micro m]
  //    dSL_dt            time derivative of sarcomere length [micro m / s]
  //    time_step         time step length [s]
  virtual void solve_time_step(std::vector<double> &state,
                               const double &calcium,
                               const double &sarcomere_length,
                               const double &dSL_dt,
                               const double &time_step) = 0;

  // Compute the active tension.
  // Inputs:
  //    state             state variables
  //    sarcomere_length  sarcomere length [micro m]
  // Outputs:
  //    active_tension    active tension [kPa]
  virtual double get_active_tension(const std::vector<double> &state,
                                    const double &sarcomere_length) = 0;
};

#endif /* SARCOMERE_HPP */