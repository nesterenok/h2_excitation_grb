#pragma once

#include <cstring>
#include <vector>

#include "dynamic_array.h"
#include "spectroscopy.h"
#include "h2_elspectrum_eqs.h"

// Evolution of electron spectrum with time,  
// Note, the data unit is the number of electrons in the energy interval per cm3 (is not divided by energy bin),
void save_electron_spectrum_evolution(const std::string& output_path, const std::vector<double>& electron_energies_grid, const std::vector<double>& electron_energy_bin_size, 
    const std::vector<double>& time_moments, const std::vector<dynamic_array>& spectrum_data, double grb_distance, double hcolumn_dens, double conc_h_tot);

// Specimen concentration [cm-3] as a function of time,
void save_specimen_conc_evolution(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<dynamic_array>& conc_data, double grb_distance, double hcolumn_dens, double conc_h_tot);

// H2 population densities [cm-3] as a function of time,
void save_h2_populations_evolution(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<dynamic_array>& h2_popdens_data, const std::vector<dynamic_array>& h2_popdens_v_data, double grb_distance, double hcolumn_dens, 
    double conc_h_tot, int nb_lev_h2);

// He level population densities [cm-3] as a function of time,
void save_hei_populations_evolution(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<dynamic_array>& hei_popdens_data, double grb_distance, double hcolumn_dens, double conc_h_tot, int nb_lev_hei);

void save_electron_energy_loss_rates(const std::string& output_path, double grb_distance, double hcolumn_dens, double conc_h_tot, const std::vector<double>& time_moments,
    const std::vector<double>& enloss_rates_mt, 
    const std::vector<double>& enloss_rates_h2_rot,
    const std::vector<double>& enloss_rates_h2_vibr,
    const std::vector<double>& enloss_rates_h2_singlet,
    const std::vector<double>& enloss_rates_h2_triplet,
    const std::vector<double>& enloss_rates_ioniz,
    const std::vector<double>& enloss_rates_coloumb_ions,
    const std::vector<double>& enloss_rates_coloumb_el,
    const std::vector<double>& enloss_rates_hei);

void save_electron_energy_losses(const std::string& output_path, double grb_distance, double hcolumn_dens, double conc_h_tot, const std::vector<double>& time_moments,
    const std::vector<double>& enloss_mt,
    const std::vector<double>& enloss_h2_rot,
    const std::vector<double>& enloss_h2_vibr,
    const std::vector<double>& enloss_h2_singlet,
    const std::vector<double>& enloss_h2_triplet,
    const std::vector<double>& enloss_ioniz,
    const std::vector<double>& enloss_coloumb_ions,
    const std::vector<double>& enloss_coloumb_el,
    const std::vector<double>& enloss_hei);

// Dissociated or excited number of species per cm3 up to a given time, as a function of time,
void save_diss_excit_data(const std::string& output_path, const std::vector<double>& time_moments, 
    const std::vector<double>& h2_solomon_diss_arr, const std::vector<double>& h2_diss_exc_singlet_arr,
    const std::vector<double>& h2_diss_exc_triplet_arr, const std::vector<double>& hei_exc_arr, double grb_distance, double hcolumn_dens, double conc_h_tot);

//
void save_phys_parameters(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<double>& neut_temp_arr, const std::vector<double>& ion_temp_arr, const std::vector<double>& dust_charge_arr,
    const std::vector<dynamic_array>& conc_data, double grb_distance, double hcolumn_dens, double conc_h_tot);
