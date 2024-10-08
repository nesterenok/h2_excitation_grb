#pragma once

#include <cstring>
#include <vector>

#include "dynamic_array.h"
#include "spectroscopy.h"
#include "h2_elspectrum_eqs.h"


//
// Evolution of electron energy distribution
//

// Be careful, the data point is the number of electrons in the energy interval per cm3 (is not divided by energy bin),
void save_electron_spectrum_evolution(const std::string& output_path, const elspectra_evolution_data* user_data,
    const std::vector<double>& time_moments, const std::vector<dynamic_array>& spectrum_data);

void save_specimen_conc_evolution(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<dynamic_array>& conc_data);

void save_h2_populations_evolution(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<dynamic_array>& h2_pop_data, const energy_diagram* h2_di);

void save_hei_populations_evolution(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<dynamic_array>& hei_pop_data, const energy_diagram* hei_di);

void save_electron_energy_losses(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<double>& deriv_mt_arr, const std::vector<double>& mt_arr,
    const std::vector<double>& deriv_h2_rot_arr, const std::vector<double>& h2_rot_arr,
    const std::vector<double>& deriv_h2_vibr_arr, const std::vector<double>& h2_vibr_arr,
    const std::vector<double>& deriv_h2_electr_arr, const std::vector<double>& h2_electr_arr,
    const std::vector<double>& deriv_h2_electr_tr_arr, const std::vector<double>& h2_electr_tr_arr,
    const std::vector<double>& deriv_ioniz_arr, const std::vector<double>& ioniz_arr,
    const std::vector<double>& deriv_coloumb_ions_arr, const std::vector<double>& coloumb_ions_arr,
    const std::vector<double>& deriv_coloumb_el_arr, const std::vector<double>& coloumb_el_arr,
    const std::vector<double>& deriv_hei_arr, const std::vector<double>& hei_arr);

void save_h2_data(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<double>& h2_solomon_diss_arr, const std::vector<double>& h2_diss_exc_arr,
    const std::vector<double>& h2_diss_exc_triplet_arr, const std::vector<double>& hei_exc_arr);

void save_phys_parameters(const std::string& output_path, const std::vector<double>& time_moments,
    const std::vector<double>& neut_temp_arr, const std::vector<double>& ion_temp_arr, const std::vector<double>& dust_charge_arr,
    const std::vector<dynamic_array>& conc_data);


