// Author: Francesco Regazzoni - MOX, Politecnico di Milano
// Email:  francesco.regazzoni@polimi.it
// Date:   2020

#ifndef MODEL_RDQ18_HPP
#define MODEL_RDQ18_HPP

#include <vector>

#include "sarcomere.hpp"

class model_RDQ18 : public sarcomere {

  //    Class implementing the ODE model for sarcomere activation proposed in
  //    [1].
  //
  //    References:
  //
  //    [1] F. Regazzoni, L. Dedè, A. Quarteroni "Active contraction of cardiac
  //        cells: a reduced model for sarcomere dynamics with cooperative
  //        interactions", Biomechanics and Modeling in Mechanobiology (2018)
  //        https://doi.org/10.1007/s10237-018-1049-0

public:
  // Constructor.
  model_RDQ18(std::string parameters_file = "../../params/params_RDQ18.json");

  // Solve one time step, updating the state passed as reference.
  void solve_time_step(std::vector<double> &state, const double &calcium,
                       const double &sarcomere_length, const double &dSL_dt,
                       const double &dt);

  // Compute the active tension.
  double get_active_tension(const std::vector<double> &state,
                            const double &sarcomere_length);

  // Compute the permissivity.
  double get_permissivity(const std::vector<double> &state);

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

  /// Coordinate of the i-th RU.
  inline double x_i(const unsigned int &i_RU) {
    return (prm_LM - prm_LB) * .5 * (i_RU + 1) / prm_n_RU;
  }

  /// Coordinate of the end of the MF.
  inline double xAZ(const double &sarcomere_length) {
    return (sarcomere_length - prm_LB) * .5;
  }

  /// Coordinate of the end of the double overlap region.
  inline double xLA(const double &sarcomere_length) {
    return prm_LA - xAZ(sarcomere_length) - prm_LB;
  }

  /// Coordinate of the beginning of the double overlap region.
  inline double xRA(const double &sarcomere_length) {
    return xAZ(sarcomere_length) - prm_LA;
  }

  /// Smoothed indicator function of the overlap region.
  inline double ChiRA(const double &sarcomere_length,
                      const unsigned int &i_RU) {
    if (x_i(i_RU) <= xRA(sarcomere_length))
      return std::exp(-std::pow(xRA(sarcomere_length) - x_i(i_RU), 2) /
                      std::pow(prm_Lsmooth, 2));
    else if (x_i(i_RU) > xRA(sarcomere_length) &&
             x_i(i_RU) < xAZ(sarcomere_length))
      return 1.0;
    else
      return std::exp(-std::pow(x_i(i_RU) - xAZ(sarcomere_length), 2) /
                      std::pow(prm_Lsmooth, 2));
  }

  /// Smoothed indicator function of the single overlap region.
  inline double ChiLA(const double &sarcomere_length,
                      const unsigned int &i_RU) {
    if (x_i(i_RU) <= xLA(sarcomere_length))
      return std::exp(-std::pow(xLA(sarcomere_length) - x_i(i_RU), 2) /
                      std::pow(prm_Lsmooth, 2));
    else
      return 1.0;
  }

  /// Function Q.
  inline double Q(const double &sarcomere_length) {
    if (sarcomere_length >= prm_SLQ)
      return prm_Q0;
    else
      return prm_Q0 - prm_alphaQ * (prm_SLQ - sarcomere_length);
  }

  // Model Parameters.
  unsigned int prm_n_RU; // [-]
  double prm_LA;         // [micro m]
  double prm_LM;         // [micro m]
  double prm_LB;         // [micro m]
  double prm_Lsmooth;    // [micro m]
  double prm_Q0;         // [-]
  double prm_SLQ;        // [micro m]
  double prm_alphaQ;     // [1 / micro m]
  double prm_mu;         // [-]
  double prm_gamma;      // [-]
  double prm_Kon;        // [micro M^-1 * s^-1]
  double prm_Koff;       // [s^-1]
  double prm_Kbasic;     // [s^-1]
  double prm_TaMax;      // [Pa]

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

  unsigned int i_RU;     // Index of regulatory_unit.
  unsigned int RU_L;     // Index of left regulatory unit.
  unsigned int RU_C;     // Index of central regulatory unit.
  unsigned int RU_R;     // Index of right regulatory unit.
  unsigned int RU_new;   // Index of regulatory unit (new state).
  unsigned int RU_dummy; // Index of regulatory unit external to the
                         // considered triplet.
};

#endif /* MODEL_RDQ18 */