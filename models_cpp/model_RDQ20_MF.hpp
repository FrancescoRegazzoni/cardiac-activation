// Author: Francesco Regazzoni - MOX, Politecnico di Milano
// Email:  francesco.regazzoni@polimi.it
// Date:   2020

#ifndef MODEL_RDQ20_MF_HPP
#define MODEL_RDQ20_MF_HPP

#include <vector>

#include "sarcomere.hpp"

class model_RDQ20_MF : public sarcomere {

  //    Class implementing the mean-field ODE model (MF-ODE) for cardiomyocytes
  //    force generation presented in [1, 2].
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
  model_RDQ20_MF(
      std::string parameters_file =
          "../../params/params_RDQ20-MF_human_body-temperature.json");

  // Solve one time step, updating the state passed as reference.
  void solve_time_step(std::vector<double> &state, const double &calcium,
                       const double &sarcomere_length, const double &dSL_dt,
                       const double &dt);

  // Compute the active tension.
  double get_active_tension(const std::vector<double> &state,
                            const double &sarcomere_length);

  // Compute the active stiffness.
  double get_active_stiffness(const std::vector<double> &state,
                              const double &sarcomere_length);

private:
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
  void XB_update_state(const double &dSL_dt, const double &dt);

  // Single-overlap ratio function at a given sarcomere length.
  double fraction_single_overlap(const double &sarcomere_length) const;

  // Model Parameters.
  double prm_LA;      // [micro m]
  double prm_LM;      // [micro m]
  double prm_LB;      // [micro m]
  double prm_SL0;     // [micro m]
  double prm_Q;       // [-]
  double prm_Kd0;     // [micro M]
  double prm_alphaKd; // [micro M / micro m]
  double prm_mu;      // [-]
  double prm_gamma;   // [-]
  double prm_Koff;    // [s^-1]
  double prm_Kbasic;  // [s^-1]
  double prm_r0;      // [s^-1]
  double prm_alpha;   // [-]
  double prm_mu0_fP;  // [s^-1]
  double prm_mu1_fP;  // [s^-1]
  double prm_a_XB;    // [kPa]

  // Numerical Parameters.
  double prm_time_step_update_RU_state =
      2.5e-5; // Time step used to update the RU state. [s]

  // Additional variables.
  double state_RU[2][2][2][2]; // State variables of RU.
  double state_XB[4];          // State variables of XB.

  double PhiT_L[2][2][2][2]; // Probability fluxes of left tropomyosin.
  double PhiT_C[2][2][2][2]; // Probability fluxes of central tropomyosin.
  double PhiT_R[2][2][2][2]; // Probability fluxes of right tropomyosin.
  double PhiC_C[2][2][2][2]; // Probability fluxes of central troponin.

  double rates_T_L[2][2]; // Transition rates of left tropomyosin.
  double rates_T_R[2][2]; // Transition rates of right tropomyosin.

  double rates_C[2][2];       // Transition rates of troponin.
  double rates_T[2][2][2][2]; // Transition rates of tropomyosin.

  unsigned int TL;   // Index of left tropomyosin.
  unsigned int TC;   // Index of central tropomyosin.
  unsigned int TR;   // Index of right tropomyosin.
  unsigned int CC;   // Index of central troponin.
  unsigned int i_XB; // Index of XB state.
};

#endif /* MODEL_RDQ20_MF_HPP */