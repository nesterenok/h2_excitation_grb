#pragma once

#include <cstring>
#include <vector>

#include "dynamic_array.h"
#include "spectroscopy.h"
#include "h2_elspectrum_eqs.h"


// Saving parameters of the model
void save_model_parameters(const std::string& output_path,  double conc_h_tot, double op_ratio_h2, double ioniz_fract, double dust_gas_mass_ratio, 
    double grain_radius, double grain_nb_density, double grb_cloud_distance = 0., double grb_distance = 0., double hcolumn_dens = 0., int layer_nb = -1);


// Evolution of electron spectrum with time,  
// Note, the data unit is the number of electrons in the energy interval per cm3 (is not divided by energy bin),
void save_electron_spectrum_evolution(const std::string& output_path, const std::vector<double>& electron_energies_grid, const std::vector<double>& electron_energy_bin_size, 
    const std::vector<double>& time_moments, const std::vector<dynamic_array>& spectrum_data, double conc_h_tot);

// Specimen concentration [cm-3] as a function of time,
void save_specimen_conc_evolution(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<dynamic_array>& conc_data, double conc_h_tot);

// H2 population densities [cm-3] as a function of time,
void save_h2_populations_evolution(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<dynamic_array>& h2_popdens_data, const std::vector<dynamic_array>& h2_popdens_v_data, double conc_h_tot, int nb_lev_h2);

// He level population densities [cm-3] as a function of time,
void save_hei_populations_evolution(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<dynamic_array>& hei_popdens_data, double conc_h_tot, int nb_lev_hei);

void save_electron_energy_loss_rates(const std::string& output_path, double conc_h_tot, const std::vector<double>& time_moments,
    const std::vector<double>& enloss_rates_mt, 
    const std::vector<double>& enloss_rates_h2_rot,
    const std::vector<double>& enloss_rates_h2_rot_pos,
    const std::vector<double>& enloss_rates_h2_vibr,
    const std::vector<double>& enloss_rates_h2_vibr_pos,
    const std::vector<double>& enloss_rates_h2_singlet,
    const std::vector<double>& enloss_rates_h2_triplet,
    const std::vector<double>& enloss_rates_ioniz,
    const std::vector<double>& enloss_rates_hei,
    const std::vector<double>& enloss_rates_coulomb_el,
    const std::vector<double>& diss_decay_heating_rate,
    const std::vector<double>& neut_heat_coll_rates);


void save_electron_energy_losses(const std::string& output_path, double conc_h_tot, const std::vector<double>& time_moments,
    const std::vector<double>& enloss_mt,
    const std::vector<double>& enloss_h2_rot,
    const std::vector<double>& enloss_h2_vibr,
    const std::vector<double>& enloss_h2_singlet,
    const std::vector<double>& enloss_h2_triplet,
    const std::vector<double>& enloss_ioniz,
    const std::vector<double>& enloss_hei,
    const std::vector<double>& enloss_coulomb_el,
    const std::vector<double>& diss_decay_heating,
    const std::vector<double>& neut_heat_coll);


// Dissociation or excition rate of species, [cm-3 s-1], as a function of time,
void save_excit_rates(const std::string& output_path, double conc_h_tot, const std::vector<double>& time_moments,
    const std::vector<double>& h2_excit_electr_rate_arr, 
    const std::vector<double>& h2_excit_vibr_rate_arr, 
    const std::vector<double>& h2_excit_rot_rate_arr);


void save_diss(const std::string& output_path, double conc_h_tot, const std::vector<double>& time_moments,
    const std::vector<double>& h2_solomon_diss_arr,
    const std::vector<double>& h2_diss_exc_singlet_arr,
    const std::vector<double>& h2_diss_excit_electr_triplet_arr,
    const std::vector<double>& h2_excit_electr_triplet_arr);


// Dissociated or excited number of species up to a given time, [cm-3], as a function of time,
void save_excit(const std::string& output_path, double conc_h_tot, const std::vector<double>& time_moments,
    const std::vector<double>& h2_excit_electr_arr,
    const std::vector<double>& h2_excit_vibr_arr,
    const std::vector<double>& h2_excit_rot_arr,
    const std::vector<double>& h2_excit_electr_bs_arr,
    const std::vector<double>& h2_excit_electr_bps_arr,
    const std::vector<double>& h2_excit_electr_cp_arr,
    const std::vector<double>& h2_excit_electr_dp_arr,
    const std::vector<double>& h2_excit_vibr_v1_arr,
    const std::vector<double>& h2_excit_vibr_v2_arr, 
    const std::vector<double>& h2_excit_vibr_v3_arr, 
    const std::vector<double>& h2_excit_el_vstate_arr,
    const std::vector<double>& hei_exc_arr);

//
void save_phys_parameters(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<double>& neut_temp_arr, const std::vector<double>& ion_temp_arr, const std::vector<double>& dust_charge_arr,
    const std::vector<dynamic_array>& conc_data, double conc_h_tot);
