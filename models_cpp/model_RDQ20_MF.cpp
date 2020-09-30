// Author: Francesco Regazzoni - MOX, Politecnico di Milano
// Email:  francesco.regazzoni@polimi.it
// Date:   2020

#include <cmath>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <Eigen/Dense>

#include "model_RDQ20_MF.hpp"

model_RDQ20_MF::model_RDQ20_MF(std::string parameters_file)
    : sarcomere("RDQ20-MF") {
  // Read JSON options file
  boost::property_tree::ptree root;
  boost::property_tree::read_json(parameters_file, root);

  prm_LA = root.get_child("geometry")
               .get_child("LA")
               .get_value<double>(); // [micro m]
  prm_LM = root.get_child("geometry")
               .get_child("LM")
               .get_value<double>(); // [micro m]
  prm_LB = root.get_child("geometry")
               .get_child("LB")
               .get_value<double>(); // [micro m]
  prm_SL0 = root.get_child("geometry")
                .get_child("SL0")
                .get_value<double>(); // [micro m]
  prm_Q = root.get_child("RU_steady_state")
              .get_child("Q")
              .get_value<double>(); // [-]
  prm_Kd0 = root.get_child("RU_steady_state")
                .get_child("Kd0")
                .get_value<double>(); // [micro M]
  prm_alphaKd = root.get_child("RU_steady_state")
                    .get_child("alphaKd")
                    .get_value<double>(); // [micro M / micro m]
  prm_mu = root.get_child("RU_steady_state")
               .get_child("mu")
               .get_value<double>(); // [-]
  prm_gamma = root.get_child("RU_steady_state")
                  .get_child("gamma")
                  .get_value<double>(); // [-]
  prm_Koff = root.get_child("RU_kinetics")
                 .get_child("Koff")
                 .get_value<double>(); // [s^-1]
  prm_Kbasic = root.get_child("RU_kinetics")
                   .get_child("Kbasic")
                   .get_value<double>(); // [s^-1]
  prm_r0 = root.get_child("XB_cycling")
               .get_child("r0")
               .get_value<double>(); // [s^-1]
  prm_alpha = root.get_child("XB_cycling")
                  .get_child("alpha")
                  .get_value<double>(); // [-]
  prm_mu0_fP = root.get_child("XB_cycling")
                   .get_child("mu0_fP")
                   .get_value<double>(); // [s^-1]
  prm_mu1_fP = root.get_child("XB_cycling")
                   .get_child("mu1_fP")
                   .get_value<double>(); // [s^-1]
  prm_a_XB = root.get_child("upscaling")
                 .get_child("a_XB")
                 .get_value<double>(); // [kPa]

  initialize_rates();

  n_variables = 20;
  initial_state.push_back(1.0);
  for (unsigned int i = 0; i < n_variables - 1; ++i)
    initial_state.push_back(0.0);
}

void model_RDQ20_MF::solve_time_step(std::vector<double> &state,
                                     const double &calcium,
                                     const double &sarcomere_length,
                                     const double &dSL_dt, const double &dt) {

  // Deserialize state
  deserialize_state(state);

  // Update RU transition rates
  RU_update_rates(calcium, sarcomere_length);

  // Advance RU state
  double RU_dt = 0.0;
  double time_advanced = 0.0;
  while (time_advanced <=
         dt - 1e-10) // Cover the time-step up to a given tolerance
  {
    RU_dt = std::min(prm_time_step_update_RU_state, dt - time_advanced);
    RU_update_state(RU_dt);
    time_advanced += RU_dt;
  }

  // Advance XB state
  XB_update_state(dSL_dt, dt);

  // Re-serialize state
  serialize_state(state);
}

void model_RDQ20_MF::deserialize_state(const std::vector<double> &state) {
  unsigned int i_current = 0;
  for (TL = 0; TL < 2; ++TL)
    for (TC = 0; TC < 2; ++TC)
      for (TR = 0; TR < 2; ++TR)
        for (CC = 0; CC < 2; ++CC) {
          state_RU[TL][TC][TR][CC] = state[i_current];
          i_current++;
        }

  for (i_XB = 0; i_XB < 4; ++i_XB) {
    state_XB[i_XB] = state[i_current];
    i_current++;
  }
}

void model_RDQ20_MF::serialize_state(std::vector<double> &state) {
  unsigned int i_current = 0;
  for (TL = 0; TL < 2; ++TL)
    for (TC = 0; TC < 2; ++TC)
      for (TR = 0; TR < 2; ++TR)
        for (CC = 0; CC < 2; ++CC) {
          state[i_current] = state_RU[TL][TC][TR][CC];
          i_current++;
        }

  for (i_XB = 0; i_XB < 4; ++i_XB) {
    state[i_current] = state_XB[i_XB];
    i_current++;
  }
}

double model_RDQ20_MF::get_active_tension(const std::vector<double> &state,
                                          const double &sarcomere_length) {
  return prm_a_XB * (state[17] + state[19]) *
         fraction_single_overlap(sarcomere_length);
}

double model_RDQ20_MF::get_active_stiffness(const std::vector<double> &state,
                                            const double &sarcomere_length) {
  return prm_a_XB * (state[16] + state[18]) *
         fraction_single_overlap(sarcomere_length);
}

void model_RDQ20_MF::initialize_rates() {
  rates_C[1][0] = prm_Koff;
  rates_C[1][1] = prm_Koff / prm_mu;

  int permissive_neighbors;
  for (TL = 0; TL < 2; ++TL)
    for (TR = 0; TR < 2; ++TR) {
      permissive_neighbors = TL + TR;
      rates_T[TL][1][TR][0] =
          prm_Kbasic * std::pow(prm_gamma, 2 - permissive_neighbors);
      rates_T[TL][1][TR][1] =
          prm_Kbasic * std::pow(prm_gamma, 2 - permissive_neighbors);
      rates_T[TL][0][TR][0] = prm_Q * prm_Kbasic *
                              std::pow(prm_gamma, permissive_neighbors) /
                              prm_mu;
      rates_T[TL][0][TR][1] =
          prm_Q * prm_Kbasic * std::pow(prm_gamma, permissive_neighbors);
    }
}

void model_RDQ20_MF::RU_update_rates(const double &calcium,
                                     const double &sarcomere_length) {
  for (unsigned int i = 0; i < 2; ++i)
    rates_C[0][i] = prm_Koff /
                    (prm_Kd0 - prm_alphaKd * (2.15 - sarcomere_length)) *
                    calcium;
}

void model_RDQ20_MF::RU_update_state(const double &dt) {
  // Compute fluxes associated with center unit
  for (TL = 0; TL < 2; ++TL)
    for (TC = 0; TC < 2; ++TC)
      for (TR = 0; TR < 2; ++TR)
        for (CC = 0; CC < 2; ++CC) {
          // PhiT_C(a,b,c,B) =  P( T_i^t+dt =~ b,
          // (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) ) / dt
          PhiT_C[TL][TC][TR][CC] =
              state_RU[TL][TC][TR][CC] * rates_T[TL][TC][TR][CC];
          // PhiC_C(a,b,c,B) =  P( C_i^t+dt =~ B,
          // (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) ) / dt
          PhiC_C[TL][TC][TR][CC] = state_RU[TL][TC][TR][CC] * rates_C[CC][TC];
        }

  double flux_tot;
  double prob_tot;

  // Compute rates associated with left unit
  // rates_T_L(a,b) ~= P( T_{i-1}^t+dt =~ a | (T_{i-1},T_i)^t = (a,b) ) / dt
  for (TL = 0; TL < 2; ++TL)
    for (TC = 0; TC < 2; ++TC) {
      flux_tot = 0.0;
      prob_tot = 0.0;
      for (TR = 0; TR < 2; ++TR)
        for (CC = 0; CC < 2; ++CC) {
          flux_tot += PhiT_C[TL][TC][TR][CC];
          prob_tot += state_RU[TL][TC][TR][CC];
        }
      if (prob_tot > 1e-12) {
        rates_T_L[TL][TC] = flux_tot / prob_tot;
      } else {
        rates_T_L[TL][TC] = 0.0;
      }
    }

  // Compute rates associated with right unit
  // rates_T_L(c,b) ~= P( T_{i+1}^t+dt =~ c | (T_i,T_{i+1})^t = (b,c) ) / dt
  for (TR = 0; TR < 2; ++TR)
    for (TC = 0; TC < 2; ++TC) {
      flux_tot = 0.0;
      prob_tot = 0.0;
      for (TL = 0; TL < 2; ++TL)
        for (CC = 0; CC < 2; ++CC) {
          flux_tot += PhiT_C[TL][TC][TR][CC];
          prob_tot += state_RU[TL][TC][TR][CC];
        }
      if (prob_tot > 1e-12) {
        rates_T_R[TR][TC] = flux_tot / prob_tot;
      } else {
        rates_T_R[TR][TC] = 0.0;
      }
    }

  // Compute fluxes associated with external units
  for (TL = 0; TL < 2; ++TL)
    for (TC = 0; TC < 2; ++TC)
      for (TR = 0; TR < 2; ++TR)
        for (CC = 0; CC < 2; ++CC) {
          // PhiT_L(a,b,c,B) ~= P( T_{i-1}^t+dt =~ a,
          // (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) ) / dt
          PhiT_L[TL][TC][TR][CC] = state_RU[TL][TC][TR][CC] * rates_T_L[TC][TL];
          // PhiT_R(a,b,c,B) ~= P( T_{i+1}^t+dt =~ c,
          // (T_{i-1},T_i,T_{i+1},C_i)^t = (a,b,c,B) ) / dt
          PhiT_R[TL][TC][TR][CC] = state_RU[TL][TC][TR][CC] * rates_T_R[TC][TR];
        }

  // Forward Euler advance
  for (TL = 0; TL < 2; ++TL)
    for (TC = 0; TC < 2; ++TC)
      for (TR = 0; TR < 2; ++TR)
        for (CC = 0; CC < 2; ++CC) {
          state_RU[TL][TC][TR][CC] +=
              dt * (-PhiT_L[TL][TC][TR][CC] + PhiT_L[1 - TL][TC][TR][CC] -
                    PhiT_C[TL][TC][TR][CC] + PhiT_C[TL][1 - TC][TR][CC] -
                    PhiT_R[TL][TC][TR][CC] + PhiT_R[TL][TC][1 - TR][CC] -
                    PhiC_C[TL][TC][TR][CC] + PhiC_C[TL][TC][TR][1 - CC]);
        }
}

void model_RDQ20_MF::XB_update_state(const double &dSL_dt, const double &dt) {
  double v = -dSL_dt / prm_SL0;

  double permissivity = 0.0;
  double flux_PN = 0.0;
  double flux_NP = 0.0;

  for (TL = 0; TL < 2; ++TL)
    for (TR = 0; TR < 2; ++TR)
      for (CC = 0; CC < 2; ++CC) {
        permissivity += state_RU[TL][1][TR][CC];
        flux_PN += state_RU[TL][1][TR][CC] * rates_T[TL][1][TR][CC];
        flux_NP += state_RU[TL][0][TR][CC] * rates_T[TL][0][TR][CC];
      }

  double k_PN = 0.0;
  double k_NP = 0.0;
  if (permissivity >= 1e-12)
    k_PN = flux_PN / permissivity;
  if (1 - permissivity >= 1e-12)
    k_NP = flux_NP / (1 - permissivity);

  double r = prm_r0 + prm_alpha * std::abs(v);
  double diag_P = r + k_PN;
  double diag_N = r + k_NP;

  Eigen::Matrix<double, 4, 4> XB_A;
  Eigen::Matrix<double, 4, 1> XB_rhs, XB_sol;

  // Fill matrix
  XB_A(0, 0) = -diag_P;
  XB_A(1, 1) = -diag_P;
  XB_A(2, 2) = -diag_N;
  XB_A(3, 3) = -diag_N;
  XB_A(0, 2) = k_NP;
  XB_A(1, 3) = k_NP;
  XB_A(2, 0) = k_PN;
  XB_A(3, 1) = k_PN;
  XB_A(1, 0) = -v;
  XB_A(3, 2) = -v;

  XB_A *= -dt;
  for (i_XB = 0; i_XB < 4; ++i_XB)
    XB_A(i_XB, i_XB) += 1.0;

  // Fill rhs
  XB_rhs(0) = state_XB[0] + dt * permissivity * prm_mu0_fP;
  XB_rhs(1) = state_XB[1] + dt * permissivity * prm_mu1_fP;
  XB_rhs(2) = state_XB[2];
  XB_rhs(3) = state_XB[3];

  // Implicit Euler advance
  XB_sol = XB_A.colPivHouseholderQr().solve(XB_rhs);

  for (i_XB = 0; i_XB < 4; ++i_XB)
    state_XB[i_XB] = XB_sol(i_XB);
}

double model_RDQ20_MF::fraction_single_overlap(const double &SL) const {
  double LMh = (prm_LM - prm_LB) * 0.5;
  if (SL > prm_LA && SL <= prm_LM)
    return (SL - prm_LA) / LMh;
  if (SL > prm_LM && SL <= 2 * prm_LA - prm_LB)
    return (SL + prm_LM - 2 * prm_LA) * 0.5 / LMh;
  if (SL > 2 * prm_LA - prm_LB && SL <= 2 * prm_LA + prm_LB)
    return 1.0;
  if (SL > 2 * prm_LA + prm_LB && SL <= 2 * prm_LA + prm_LM)
    return (prm_LM + 2 * prm_LA - SL) * 0.5 / LMh;
  return 0.0;
}