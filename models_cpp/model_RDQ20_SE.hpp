// Author: Francesco Regazzoni - MOX, Politecnico di Milano
// Email:  francesco.regazzoni@polimi.it
// Date:   2020

#ifndef MODEL_RDQ20_SE_HPP
#define MODEL_RDQ20_SE_HPP

#include <vector>

#include <Eigen/Dense>

#include "sarcomere.hpp"

class model_RDQ20_SE : public sarcomere {

  //    Class implementing the spatially-explicit ODE model (SE-ODE) for
  //    cardiomyocytes force generation presented in [1, 2].
  //
  //    References:
  //
  //    [1] F. Regazzoni "Mathematical modeling and Machine Learning for the
  //        numerical simulation of cardiac electromechanics", PhD Thesis -
  //        Politecnico di Milano (2020)
  //        http://hdl.handle.net/10589/152617
  //    [2] F. Regazzoni, L. Dede', A. Quarteroni "Biophysically detailed
  //        mathematical models of multiscale cardiac active mechanics",
  //        PLOS Computational Biology (2020)
  //        https://doi.org/10.1371/journal.pcbi.1008294

public:
  // Constructor.
  model_RDQ20_SE(
      std::string parameters_file =
          "../../params/params_RDQ20-SE_human_body-temperature.json");

  // Solve one time step, updating the state passed as reference.
  void solve_time_step(std::vector<double> &state, const double &calcium,
                       const double &sarcomere_length, const double &dSL_dt,
                       const double &dt);

  // Compute the active tension.
  double get_active_tension(const std::vector<double> &state,
                            const double &sarcomere_length);

  // Compute the permissivity.
  double get_permissivity(const std::vector<double> &state);

  // Compute the active stiffness.
  double get_active_stiffness(const std::vector<double> &state,
                              const double &sarcomere_length);

private:
  // Allocate memory to store variables
  void allocate_variables();

  // Initialize transition rates.
  void initialize_rates();

  // Deserialize the model state into state_RU and state_XB.
  void deserialize_state(const std::vector<double> &state);

  // Serialize the model state from state_RU and state_XB.
  void serialize_state(std::vector<double> &state);

  // Update the RU transition rates.
  void RU_update_rates(const double &calcium, const double &sarcomere_length);

  // Update the RU state variables.
  void RU_update_state(const double &dt);

  // Update the XB state variables.
  void XB_update_state(const double &sarcomere_length, const double &dSL_dt,
                       const double &dt);

  // Coordinate of the i-th RU.
  inline double y_j(const unsigned int &i_RU) {
    return prm_LA * (i_RU + 0.5) / prm_n_RU;
  }

  // Coordinate of the left end of the AF.
  inline double yLA(const double &sarcomere_length) {
    return 2 * prm_LA - sarcomere_length;
  }

  // First coordinate of the MF.
  inline double yM0(const double &sarcomere_length) {
    return (2 * prm_LA - sarcomere_length + prm_LB) * 0.5;
  }

  // Second coordinate of the MF.
  inline double yM1(const double &sarcomere_length) {
    return (2 * prm_LA - sarcomere_length + prm_LM) * 0.5;
  }

  // Smoothed indicator function of the MF.
  inline double ChiMF(const double &sarcomere_length,
                      const unsigned int &i_RU) {
    return 0.5 * std::tanh((y_j(i_RU) - yM0(sarcomere_length)) / prm_Lsmooth) +
           0.5 * std::tanh(-(y_j(i_RU) - yM1(sarcomere_length)) / prm_Lsmooth);
  }

  // Smoothed indicator function of the single-overlap region.
  inline double ChiSF(const double &sarcomere_length,
                      const unsigned int &i_RU) {
    return 0.5 +
           0.5 * std::tanh((y_j(i_RU) - yLA(sarcomere_length)) / prm_Lsmooth);
  }

  // Model Parameters.
  unsigned int prm_n_RU; // [-]
  double prm_LA;         // [micro m]
  double prm_LM;         // [micro m]
  double prm_LB;         // [micro m]
  double prm_SL0;        // [micro m]
  double prm_Lsmooth;    // [micro m]
  double prm_Q;          // [-]
  double prm_Kd0;        // [micro M]
  double prm_alphaKd;    // [micro M / micro m]
  double prm_mu;         // [-]
  double prm_gamma;      // [-]
  double prm_Koff;       // [s^-1]
  double prm_Kbasic;     // [s^-1]
  double prm_r0;         // [s^-1]
  double prm_alpha;      // [-]
  double prm_mu0_fP;     // [s^-1]
  double prm_mu1_fP;     // [s^-1]
  double prm_a_XB;       // [kPa]

  // Numerical Parameters.
  double prm_time_step_update_RU_state =
      2.5e-5; // Time step used to update the RU state. [s]

  // Additional variables.
  std::vector<std::array<std::array<std::array<double, 4>, 4>, 4>>
      state_RU;                                // State variables of RU.
  std::vector<std::array<double, 4>> state_XB; // State variables of XB.

  std::vector<
      std::array<std::array<std::array<std::array<double, 4>, 4>, 4>, 4>>
      rates_RU; // RU rates.

  int permissivity_of_state[4] = {0, 0, 1,
                                  1}; // Permissivity state (0 = N, 1 = P)
                                      // associated with the basic RU states.

  std::vector<
      std::array<std::array<std::array<std::array<double, 4>, 4>, 4>, 4>>
      flux_RU_L; // Probability flux of left-ward RUs of the triplets.
  std::vector<
      std::array<std::array<std::array<std::array<double, 4>, 4>, 4>, 4>>
      flux_RU_C; // Probability flux of central RUs of the triplets.
  std::vector<
      std::array<std::array<std::array<std::array<double, 4>, 4>, 4>, 4>>
      flux_RU_R; // Probability flux of right-ward RUs of the triplets.

  Eigen::Matrix<double, 4, 4> XB_A;   // Matrix defining the local XB system.
  Eigen::Matrix<double, 4, 1> XB_rhs; // Right-hand side of the local XB system.
  Eigen::Matrix<double, 4, 1>
      XB_sol; // Vector used to store the solution of the local XB system.

  unsigned int n_states_RU; // Number of RU states.
  unsigned int n_states_XB; // Number of XB states.
  unsigned int i_RU;        // Index of regulatory_unit.
  unsigned int RU_L;        // Index of left regulatory unit.
  unsigned int RU_C;        // Index of central regulatory unit.
  unsigned int RU_R;        // Index of right regulatory unit.
  unsigned int RU_new;      // Index of regulatory unit (new state).
  unsigned int RU_dummy;    // Index of regulatory unit external to the
                            // considered triplet.
  unsigned int i_XB;        // Index of XB state.
};

#endif /* MODEL_RDQ20_SE_HPP */