// Author: Francesco Regazzoni - MOX, Politecnico di Milano
// Email:  francesco.regazzoni@polimi.it
// Date:   2020

#include <cmath>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "model_RDQ18.hpp"

model_RDQ18::model_RDQ18(std::string parameters_file) : sarcomere("RDQ18") {
  // Read JSON options file
  boost::property_tree::ptree root;
  boost::property_tree::read_json(parameters_file, root);

  prm_n_RU = root.get_child("geometry")
                 .get_child("n_RU")
                 .get_value<unsigned int>(); // [-]
  prm_LA = root.get_child("geometry")
               .get_child("LA")
               .get_value<double>(); // [micro m]
  prm_LM = root.get_child("geometry")
               .get_child("LM")
               .get_value<double>(); // [micro m]
  prm_LB = root.get_child("geometry")
               .get_child("LB")
               .get_value<double>(); // [micro m]
  prm_Lsmooth = root.get_child("geometry")
                    .get_child("Lsmooth")
                    .get_value<double>(); // [micro m]
  prm_Q0 = root.get_child("RU_steady_state")
               .get_child("Q0")
               .get_value<double>(); // [-]
  prm_SLQ = root.get_child("RU_steady_state")
                .get_child("SLQ")
                .get_value<double>(); // [micro m]
  prm_alphaQ = root.get_child("RU_steady_state")
                   .get_child("alphaQ")
                   .get_value<double>(); // [1 / micro m]
  prm_mu = root.get_child("RU_steady_state")
               .get_child("mu")
               .get_value<double>(); // [-]
  prm_gamma = root.get_child("RU_steady_state")
                  .get_child("gamma")
                  .get_value<double>(); // [-]
  prm_Kon = root.get_child("RU_kinetics")
                .get_child("Kon")
                .get_value<double>(); // [micro M^-1 * s^-1]
  prm_Koff = root.get_child("RU_kinetics")
                 .get_child("Koff")
                 .get_value<double>(); // [s^-1]
  prm_Kbasic = root.get_child("RU_kinetics")
                   .get_child("Kbasic")
                   .get_value<double>(); // [s^-1]
  prm_TaMax = root.get_child("upscaling")
                  .get_child("TaMax")
                  .get_value<double>(); // [kPa]

  allocate_variables();

  initialize_rates();
}

void model_RDQ18::allocate_variables() {
  // Variable numbers
  n_variables = (prm_n_RU - 2) * 4 * 4 * 4;

  // Allocation of state_RU
  std::array<std::array<std::array<double, 4>, 4>, 4> base_RU_state;
  for (RU_L = 0; RU_L < 4; ++RU_L)
    for (RU_C = 0; RU_C < 4; ++RU_C)
      for (RU_R = 0; RU_R < 4; ++RU_R)
        base_RU_state[RU_L][RU_C][RU_R] = 0.0;
  base_RU_state[0][0][0] = 1.0;

  for (i_RU = 0; i_RU < prm_n_RU - 2; ++i_RU)
    state_RU.push_back(base_RU_state);

  // Allocation of initial_state
  for (unsigned int i = 0; i < n_variables; ++i)
    initial_state.push_back(0.0);
  serialize_state(initial_state);

  // Allocation of rates_RU
  std::array<std::array<std::array<std::array<double, 4>, 4>, 4>, 4>
      base_RU_rates_or_flux;
  for (RU_L = 0; RU_L < 4; ++RU_L)
    for (RU_C = 0; RU_C < 4; ++RU_C)
      for (RU_R = 0; RU_R < 4; ++RU_R)
        for (RU_new = 0; RU_new < 4; ++RU_new)
          base_RU_rates_or_flux[RU_L][RU_C][RU_R][RU_new] = 0.0;

  for (i_RU = 0; i_RU < prm_n_RU; ++i_RU) {
    rates_RU.push_back(base_RU_rates_or_flux);
  }

  // Allocation of flux_RU_*
  for (i_RU = 0; i_RU < prm_n_RU - 2; ++i_RU) {
    flux_RU_L.push_back(base_RU_rates_or_flux);
    flux_RU_C.push_back(base_RU_rates_or_flux);
    flux_RU_R.push_back(base_RU_rates_or_flux);
  }
}

void model_RDQ18::solve_time_step(std::vector<double> &state,
                                  const double &calcium,
                                  const double &sarcomere_length,
                                  const double & /*dSL_dt*/, const double &dt) {
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

  // Re-serialize state
  serialize_state(state);
}

void model_RDQ18::deserialize_state(const std::vector<double> &state) {
  unsigned int i_current = 0;
  for (i_RU = 0; i_RU < prm_n_RU - 2; ++i_RU)
    for (RU_L = 0; RU_L < 4; ++RU_L)
      for (RU_C = 0; RU_C < 4; ++RU_C)
        for (RU_R = 0; RU_R < 4; ++RU_R) {
          state_RU[i_RU][RU_L][RU_C][RU_R] = state[i_current];
          i_current++;
        }
}

void model_RDQ18::serialize_state(std::vector<double> &state) {
  unsigned int i_current = 0;
  for (i_RU = 0; i_RU < prm_n_RU - 2; ++i_RU)
    for (RU_L = 0; RU_L < 4; ++RU_L)
      for (RU_C = 0; RU_C < 4; ++RU_C)
        for (RU_R = 0; RU_R < 4; ++RU_R) {
          state[i_current] = state_RU[i_RU][RU_L][RU_C][RU_R];
          i_current++;
        }
}

double model_RDQ18::get_active_tension(const std::vector<double> &state,
                                       const double & /*sarcomere_length*/) {
  return prm_TaMax * get_permissivity(state);
}

double model_RDQ18::get_permissivity(const std::vector<double> &state) {
  deserialize_state(state);

  double permissivity = 0;
  for (i_RU = 0; i_RU < prm_n_RU; ++i_RU) {
    if (i_RU == 0) // First RU
    {
      for (RU_C = 0; RU_C < 4; ++RU_C)
        for (RU_R = 0; RU_R < 4; ++RU_R)
          permissivity +=
              state_RU[0][2][RU_C][RU_R] + state_RU[0][3][RU_C][RU_R];
    } else if (i_RU == prm_n_RU - 1) // Last RU
    {
      for (RU_L = 0; RU_L < 4; ++RU_L)
        for (RU_C = 0; RU_C < 4; ++RU_C)
          permissivity += state_RU[prm_n_RU - 3][RU_L][RU_C][2] +
                          state_RU[prm_n_RU - 3][RU_L][RU_C][3];
    } else // Intermediate RUs
    {
      for (RU_L = 0; RU_L < 4; ++RU_L)
        for (RU_R = 0; RU_R < 4; ++RU_R)
          permissivity += state_RU[i_RU - 1][RU_L][2][RU_R] +
                          state_RU[i_RU - 1][RU_L][3][RU_R];
    }
  }
  return permissivity / prm_n_RU;
}

void model_RDQ18::initialize_rates() {
  int permissive_neighbors;
  for (i_RU = 0; i_RU < prm_n_RU; ++i_RU)
    for (RU_L = 0; RU_L < 4; ++RU_L)
      for (RU_R = 0; RU_R < 4; ++RU_R) {
        permissive_neighbors =
            permissivity_of_state[RU_L] + permissivity_of_state[RU_R];
        rates_RU[i_RU][RU_L][1][RU_R][0] = prm_Koff;
        rates_RU[i_RU][RU_L][2][RU_R][3] = prm_Koff / prm_mu;
        rates_RU[i_RU][RU_L][3][RU_R][0] =
            prm_Kbasic * std::pow(prm_gamma, 2 - permissive_neighbors);
        rates_RU[i_RU][RU_L][2][RU_R][1] =
            prm_Kbasic * std::pow(prm_gamma, 2 - permissive_neighbors);
      }
}

void model_RDQ18::RU_update_rates(const double &calcium,
                                  const double &sarcomere_length) {
  int permissive_neighbors;
  double ChiRA_i_RU, ChiLA_i_RU;
  double Q_SL = Q(sarcomere_length);
  for (i_RU = 0; i_RU < prm_n_RU; ++i_RU) {
    ChiRA_i_RU = ChiRA(sarcomere_length, i_RU);
    ChiLA_i_RU = ChiLA(sarcomere_length, i_RU);
    for (RU_L = 0; RU_L < 4; ++RU_L)
      for (RU_R = 0; RU_R < 4; ++RU_R) {
        permissive_neighbors =
            permissivity_of_state[RU_L] + permissivity_of_state[RU_R];
        rates_RU[i_RU][RU_L][0][RU_R][1] = prm_Kon * calcium * ChiRA_i_RU;
        rates_RU[i_RU][RU_L][1][RU_R][2] =
            ChiRA_i_RU * ChiLA_i_RU *
            std::pow(prm_gamma, permissive_neighbors) * Q_SL * prm_Kbasic;
        rates_RU[i_RU][RU_L][3][RU_R][2] = prm_Kon * calcium * ChiRA_i_RU;
        rates_RU[i_RU][RU_L][0][RU_R][3] =
            ChiRA_i_RU * ChiLA_i_RU *
            std::pow(prm_gamma, permissive_neighbors) * Q_SL * prm_Kbasic /
            prm_mu;
      }
  }
}

void model_RDQ18::RU_update_state(const double &dt) {
  // Compute fluxes associated with center units
  for (i_RU = 0; i_RU < prm_n_RU - 2; ++i_RU)
    for (RU_L = 0; RU_L < 4; ++RU_L)
      for (RU_C = 0; RU_C < 4; ++RU_C)
        for (RU_R = 0; RU_R < 4; ++RU_R)
          for (RU_new = 0; RU_new < 4; ++RU_new) {
            flux_RU_C[i_RU][RU_L][RU_C][RU_R][RU_new] =
                rates_RU[i_RU + 1][RU_L][RU_C][RU_R][RU_new] *
                state_RU[i_RU][RU_L][RU_C][RU_R];
          }

  double prob_tot;
  double rate_tot;

  // Compute fluxes associated with left units
  // --- Most left-ward triplet
  for (RU_L = 0; RU_L < 4; ++RU_L)
    for (RU_C = 0; RU_C < 4; ++RU_C)
      for (RU_R = 0; RU_R < 4; ++RU_R)
        for (RU_new = 0; RU_new < 4; ++RU_new) {
          flux_RU_L[0][RU_L][RU_C][RU_R][RU_new] =
              rates_RU[0][0][RU_L][RU_C][RU_new] *
              state_RU[0][RU_L][RU_C][RU_R];
        }

  // --- Other triplets
  for (i_RU = 1; i_RU < prm_n_RU - 2; ++i_RU)
    for (RU_L = 0; RU_L < 4; ++RU_L)
      for (RU_C = 0; RU_C < 4; ++RU_C) {
        prob_tot = 0.0;
        for (RU_dummy = 0; RU_dummy < 4; ++RU_dummy)
          prob_tot += state_RU[i_RU - 1][RU_dummy][RU_L][RU_C];

        for (RU_new = 0; RU_new < 4; ++RU_new) {
          rate_tot = 0.0;
          if (prob_tot > 1e-12) {
            for (RU_dummy = 0; RU_dummy < 4; ++RU_dummy)
              rate_tot += flux_RU_C[i_RU - 1][RU_dummy][RU_L][RU_C][RU_new];

            rate_tot /= prob_tot;
          }

          for (RU_R = 0; RU_R < 4; ++RU_R) {
            flux_RU_L[i_RU][RU_L][RU_C][RU_R][RU_new] =
                rate_tot * state_RU[i_RU][RU_L][RU_C][RU_R];
          }
        }
      }

  // Compute fluxes associated with right units
  // --- Most right-ward triplet
  for (RU_L = 0; RU_L < 4; ++RU_L)
    for (RU_C = 0; RU_C < 4; ++RU_C)
      for (RU_R = 0; RU_R < 4; ++RU_R)
        for (RU_new = 0; RU_new < 4; ++RU_new) {
          flux_RU_R[prm_n_RU - 3][RU_L][RU_C][RU_R][RU_new] =
              rates_RU[prm_n_RU - 1][RU_C][RU_R][0][RU_new] *
              state_RU[prm_n_RU - 3][RU_L][RU_C][RU_R];
        }

  // --- Other triplets
  for (i_RU = 0; i_RU < prm_n_RU - 3; ++i_RU)
    for (RU_R = 0; RU_R < 4; ++RU_R)
      for (RU_C = 0; RU_C < 4; ++RU_C) {
        prob_tot = 0.0;
        for (RU_dummy = 0; RU_dummy < 4; ++RU_dummy)
          prob_tot += state_RU[i_RU + 1][RU_C][RU_R][RU_dummy];

        for (RU_new = 0; RU_new < 4; ++RU_new) {
          rate_tot = 0.0;
          if (prob_tot > 1e-12) {
            for (RU_dummy = 0; RU_dummy < 4; ++RU_dummy)
              rate_tot += flux_RU_C[i_RU + 1][RU_C][RU_R][RU_dummy][RU_new];

            rate_tot /= prob_tot;
          }

          for (RU_L = 0; RU_L < 4; ++RU_L) {
            flux_RU_R[i_RU][RU_L][RU_C][RU_R][RU_new] =
                rate_tot * state_RU[i_RU][RU_L][RU_C][RU_R];
          }
        }
      }

  // Forward Euler advance
  for (i_RU = 0; i_RU < prm_n_RU - 2; ++i_RU)
    for (RU_L = 0; RU_L < 4; ++RU_L)
      for (RU_C = 0; RU_C < 4; ++RU_C)
        for (RU_R = 0; RU_R < 4; ++RU_R)
          for (RU_new = 0; RU_new < 4; ++RU_new) {
            state_RU[i_RU][RU_L][RU_C][RU_R] +=
                dt * (flux_RU_L[i_RU][RU_new][RU_C][RU_R][RU_L] +
                      flux_RU_C[i_RU][RU_L][RU_new][RU_R][RU_C] +
                      flux_RU_R[i_RU][RU_L][RU_C][RU_new][RU_R] -
                      flux_RU_L[i_RU][RU_L][RU_C][RU_R][RU_new] -
                      flux_RU_C[i_RU][RU_L][RU_C][RU_R][RU_new] -
                      flux_RU_R[i_RU][RU_L][RU_C][RU_R][RU_new]);
          }
}