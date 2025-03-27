
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>
#include <ctime>
#include <limits>

#include "constants.h"
#include "h2_saving_data.h"
#include "h2_parameters.h"
using namespace std;


void save_model_parameters(const std::string& output_path, double grb_cloud_distance, double grb_distance, double hcolumn_dens, double conc_h_tot, 
	double op_ratio_h2, double dust_gas_mass_ratio, double grain_radius, double grain_nb_density, int layer_nb)
{
	time_t timer;
	char* ctime_str;

	string fname;
	ofstream output;

	fname = output_path + "h2grb_parameters_input.txt";
	output.open(fname.c_str());
	
	output << scientific;
	output.precision(5);

	timer = time(NULL);
	ctime_str = ctime(&timer);

	output << "!H2-GRB model parameters, date " << ctime_str;
	output << left << "!Molecular cloud parameters:" << endl
		<< "!Cloud layer number:" << endl << layer_nb << endl
		<< "!Distance from the GRB source to the cloud boundary [pc]:" << endl << grb_cloud_distance / PARSEC << endl
		<< "!Distance from the GRB source to the centre of the cloud layer [pc]:" << endl << grb_distance / PARSEC << endl
		<< "!H column density passed by the GRB, from the cloud boundary to the cloud layer centre [cm-2]:" << endl << hcolumn_dens << endl
		<< "!H nuclei concentration [cm-3]:" << endl << conc_h_tot << endl
		<< "!Initial H2 concentration [cm-3]:" << endl << 0.5 * conc_h_tot << endl
		<< "!Initial H2 OPR:" << endl << op_ratio_h2 << endl;
	
	output << left << "!Coulomb losses (0/1):" << endl << CALC_EL_LOSSES_THERMAL_EL << endl
		<< "!Ionization fraction of thermal electrons:" << endl << IONIZATION_FRACTION_THERMAL_EL << endl
		<< "!Thermal electron temperature [K]:" << endl << THERMAL_EL_TEMPERATURE << endl;
	
	output << left << "!Collisions of H2 with H2, He included (0/1):" << endl << H2_COLLISIONS_WITH_H2_HE << endl
		<< "!Temperature of the neutral gas in the cloud [K]:" << endl << THERMAL_NEUTRAL_TEMPERATURE << endl;

	output << left << "!Nb of H2 vibrational states taking into account in electron-impact excitation:" << endl << NB_OF_H2_VSTATES_X1SU << endl;
	
	output << left << "!Dust parameters:" << endl
		<< "!Dust-gas mass ratio:" << endl << dust_gas_mass_ratio << endl
		<< "!Initial grain radius [cm]:" << endl << grain_radius << endl
		<< "!Grain number density [cm-3]:" << endl << grain_nb_density << endl;		
	output.close();
}

void save_electron_spectrum_evolution(const string& output_path, const vector<double>& electron_energies_grid, const vector<double>& electron_energy_bin_size,
	const vector<double>& time_moments, const vector<dynamic_array>& spectrum_data, double grb_distance, double hcolumn_dens, double conc_h_tot)
{
	int i, l, nb_of_times, nb_of_el_energies;
	double en, d_en, tot_en, tot_n;

	string fname;
	ofstream output;

	nb_of_el_energies = (int)electron_energy_bin_size.size();
	nb_of_times = (int)spectrum_data.size();

	fname = output_path + "h2grb_electron_spectrum.txt";
	output.open(fname.c_str());

	output << "!Electron spectrum evolution," << endl
		<< "!first row: distance from the GRB source to cloud layer centre [pc], H column density from the cloud boundary to layer centre [cm-2], H nuclei concentration [cm-3], nb of rows (= nb of el. energies) and nb of times," << endl
		<< "!second row: time (s)," << endl
		<< "!data: col.1 - centre of the energy bin [eV], col.2 - bin semi-length [eV], col.3 - spectrum N(E) [cm-3 eV-1]," << endl;

	output << scientific;
	output.precision(6);

	output << left << "! " << setw(15) << grb_distance / PARSEC << setw(15) << hcolumn_dens << setw(15) << conc_h_tot 
		<< setw(7) << nb_of_el_energies << setw(7) << nb_of_times << endl;

	output.precision(4);
	output << left << setw(28) << "!" << setw(13) << time_moments[0];
	for (l = 1; l < nb_of_times; l += NB_OF_TIME_STEPS) {
		output << left << setw(13) << time_moments[l];
	}
	output << endl;

	for (i = 0; i < nb_of_el_energies; i++) {
		en = electron_energies_grid[i];
		d_en = electron_energy_bin_size[i];

		output.precision(6);
		output << left << setw(14) << en + 0.5 * d_en << setw(14) << 0.5 * d_en;

		output.precision(4);
		output << left << setw(13) << spectrum_data[0].arr[i] / d_en;
		
		for (l = 1; l < nb_of_times; l += NB_OF_TIME_STEPS) {
			// spectrum in the data array - nb of electrons in the energy interval per cm3 [cm-3]
			output << left << setw(13) << spectrum_data[l].arr[i] / d_en; 
		}
		output << endl;
	}
	output.close();

	fname = output_path + "h2grb_electron_energy.txt";
	output.open(fname.c_str());

	output << "!Kinetic energy of electrons," << endl
		<< "!first row: distance from the GRB source to cloud layer centre [pc], H column density from the cloud boundary to layer centre [cm-2], H nuclei concentration [cm-3], nb of rows (= nb of times)," << endl
		<< "!data: time (s), conc of electrons (from spectrum) [cm-3], total energy [eV cm-3], median energy [eV]," << endl;

	output << scientific;
	output.precision(6);

	output << left << "! " << setw(15) << grb_distance / PARSEC << setw(15) << hcolumn_dens << setw(15) << conc_h_tot 
		<< setw(7) << nb_of_times << endl;

	output.precision(4);
	for (l = 0; l < nb_of_times; l++) {
		tot_n = tot_en = 0.;

		for (i = 0; i < nb_of_el_energies; i++) {
			en = electron_energies_grid[i];
			d_en = electron_energy_bin_size[i];
			en += 0.5 * d_en;

			tot_n += spectrum_data[l].arr[i];  // [cm-3]
			tot_en += spectrum_data[l].arr[i] * en;
		}
		output << left << setw(13) << time_moments[l] << setw(13) << tot_n << setw(13) << tot_en
			<< setw(13) << tot_en / tot_n << endl;
	}
	output.close();
}


void save_specimen_conc_evolution(const string& output_path, const vector<double>& time_moments,
	const vector<dynamic_array>& conc_data, double grb_distance, double hcolumn_dens, double conc_h_tot)
{
	int i, l;
	double w;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_specimen_conc.txt";
	output.open(fname.c_str());

	output << left << "!Evolution of concentrations of chemical species," << endl
		<< "!first row - distance from the GRB source to cloud layer centre [pc], H column density from the cloud boundary to layer centre [cm-2], H nuclei concentration [cm-3], nb of times, number of species," << endl
		<< "!second row - specimen list," << endl;

	output << scientific;
	output.precision(6);

	output << left << "! " << setw(15) << grb_distance / PARSEC << setw(15) << hcolumn_dens << setw(15) << conc_h_tot 
		<< setw(7) << (int)conc_data.size() << setw(7) << NB_OF_CHEM_SPECIES << endl;
	
	output << left << setw(13) << "!time(s)";
	for (i = 0; i < NB_OF_CHEM_SPECIES; i++) {
		output << left << setw(12) << chemical_species[i];
	}
	output << endl;

	for (l = 0; l < (int)conc_data.size(); l++) {
		output.precision(4);
		output << left << setw(13) << time_moments[l];

		output.precision(3);
		for (i = 0; i < NB_OF_CHEM_SPECIES; i++) 
		{
			w = conc_data[l].arr[i];
			if (w < MINIMAL_ABUNDANCE)
				w = MINIMAL_ABUNDANCE;
			output << left << setw(12) << w;
		}
		output << endl;
	}
	output.close();
}


void save_h2_populations_evolution(const string& output_path, const vector<double>& time_moments,
	const vector<dynamic_array>& h2_popdens_data, const vector<dynamic_array>& h2_popdens_v_data, double grb_distance, double hcolumn_dens, 
	double conc_h_tot, int nb_lev_h2)
{
	int i, l, v;
	double w;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_popdens.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);

	output << left << "!Evolution of population densities [cm-3] of H2 energy levels," << endl
		<< "!first row - distance from the GRB source to cloud layer centre [pc], H column density from the cloud boundary to layer centre [cm-2], H nuclei concentration [cm-3], nb of times, number of levels," << endl
		<< "!second row - level nb," << endl;

	output << left << "! " << setw(15) << grb_distance / PARSEC << setw(15) << hcolumn_dens << setw(15) << conc_h_tot 
		<< setw(7) << (int)h2_popdens_data.size() << setw(7) << nb_lev_h2 << endl;
	
	output << left << setw(13) << "!time(s)";
	for (i = 1; i <= nb_lev_h2; i++) {
		output << left << setw(12) << i;
	}
	output << endl;

	for (l = 0; l < (int)h2_popdens_data.size(); l++) {
		output.precision(4);
		output << left << setw(13) << time_moments[l];

		output.precision(3);
		for (i = 0; i < nb_lev_h2; i++) 
		{
			w = h2_popdens_data[l].arr[i];
			if (w < MINIMAL_ABUNDANCE)
				w = MINIMAL_ABUNDANCE;
			output << left << setw(12) << w;
		}
		output << endl;
	}
	output.close();

	fname = output_path + "h2grb_popdens_v.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);

	output << left << "!Evolution of population densities [cm-3] of vibrational states of H2," << endl
		<< "!first row - distance from the GRB source to cloud layer centre [pc], H column density from the cloud boundary to layer centre [cm-2], H nuclei concentration [cm-3], nb of times, number of vibr states," << endl
		<< "!second row - summed value, v state nb," << endl;

	output << left << "! " << setw(15) << grb_distance / PARSEC << setw(15) << hcolumn_dens << setw(15) << conc_h_tot 
		<< setw(7) << (int)h2_popdens_v_data.size() << setw(7) << MAX_H2_VSTATES_X1SU << endl;

	output << left << setw(13) << "!time(s)" << setw(13) << "sum";
	for (v = 0; v < MAX_H2_VSTATES_X1SU; v++) {
		output << left << setw(12) << v;
	}
	output << endl;

	for (l = 0; l < (int)h2_popdens_v_data.size(); l++) {
		output.precision(4);
		output << left << setw(13) << time_moments[l];

		// summed over the all vibrational states, is equal to H2 concentration,
		w = 0.;
		for (v = 0; v < MAX_H2_VSTATES_X1SU; v++) {
			w += h2_popdens_v_data[l].arr[v];
		}
		output.precision(3);
		output << left << setw(13) << w;

		for (v = 0; v < MAX_H2_VSTATES_X1SU; v++)
		{
			w = h2_popdens_v_data[l].arr[v];
			if (w < MINIMAL_ABUNDANCE)
				w = MINIMAL_ABUNDANCE;
			output << left << setw(12) << w;
		}
		output << endl;
	}
	output.close();
}


void save_hei_populations_evolution(const string& output_path, const vector<double>& time_moments,
	const vector<dynamic_array>& hei_popdens_data, double grb_distance, double hcolumn_dens, double conc_h_tot, int nb_lev_hei)
{
	int i, l;
	double w;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_hei_populations.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);

	output << left << "!Evolution of population densities [cm-3] of HeI energy levels," << endl
		<< "!first row - distance from the GRB source to cloud layer centre [pc], H column density from the cloud boundary to layer centre [cm-2], H nuclei concentration [cm-3], nb of times, number of HeI levels," << endl
		<< "!second row - level nb," << endl;

	output << left << "! " << setw(15) << grb_distance / PARSEC << setw(15) << hcolumn_dens << setw(15) << conc_h_tot 
		<< setw(7) << (int)hei_popdens_data.size() << setw(7) << nb_lev_hei << endl;

	output << left << setw(13) << "!time(s)";
	for (i = 1; i <= nb_lev_hei; i++) {
		output << left << setw(12) << i;
	}
	output << endl;

	for (l = 0; l < (int)hei_popdens_data.size(); l++) {
		output.precision(4);
		output << left << setw(13) << time_moments[l];

		output.precision(3);
		for (i = 0; i < nb_lev_hei; i++)
		{
			w = hei_popdens_data[l].arr[i];
			if (w < MINIMAL_ABUNDANCE)
				w = MINIMAL_ABUNDANCE;
			output << left << setw(12) << w;
		}
		output << endl;
	}
	output.close();
}


void save_electron_energy_loss_rates(const string& output_path, double grb_distance, double hcolumn_dens, double conc_h_tot, const vector<double>& time_moments,
	const vector<double>& enloss_rates_mt,
	const vector<double>& enloss_rates_h2_rot, 
	const vector<double>& enloss_rates_h2_rot_pos,
	const vector<double>& enloss_rates_h2_vibr, 
	const vector<double>& enloss_rates_h2_vibr_pos,
	const vector<double>& enloss_rates_h2_singlet, 
	const vector<double>& enloss_rates_h2_triplet, 
	const vector<double>& enloss_rates_ioniz, 
	const vector<double>& enloss_rates_hei,
	const vector<double>& enloss_rates_coulomb_el,
	const vector<double>& diss_decay_heating_rate,
	const vector<double>& neut_heat_coll_rates)
{
	int i;
	double tot;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_electron_enloss_rates.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);

	output << left << "!Energy loss rates of electrons, parameter < 0 if electrons lose energy, rate in [eV cm-3 s-1]," << endl
		<< "!tot - total rate, " << endl
		<< "!mt - momentum transfer, " << endl
		<< "!rot      - H2 rotational excitation," << endl
		<< "!rot_pos  - H2 rotational excitation, electrons gain energy," << endl
		<< "!vibr     - H2 ro-vibrational excitation (within the ground electronic state, including pure rotational transitions)," << endl
		<< "!vibr_pos - H2 ro-vibrational excitation, electrons gain energy" << endl
		<< "!elec_s - H2 electronic states excitation (singlet)," << endl
		<< "!elec_t - H2 electronic states excitation (triplet)," << endl
		<< "!ion - molecule/atom ionization," << endl
		<< "!hei - HeI excitation by electron impact" << endl
		<< "!ecoul - Coulomb loses on electrons," << endl
		<< "!ndiss - neutral heating rate due to H2 dissociation, > 0" << endl
		<< "!ncoll - neutral heating rate due to collisions H2-H2, H2-He," << endl;

	output << left << "!distance from the GRB source to cloud layer centre [pc], H column density from the cloud boundary to layer centre [cm-2], H nuclei concentration [cm-3], nb of times," << endl
		<< "! " << setw(15) << grb_distance / PARSEC << setw(15) << hcolumn_dens << setw(15) << conc_h_tot 
		<< setw(7) << (int)enloss_rates_mt.size() << endl;
	
	output << left << setw(13) << "!time(s)"
		<< setw(11) << "tot" << setw(11) << "mt" << setw(11) << "rot" << setw(11) << "rot_pos" << setw(11) << "vibr" << setw(11) << "vibr_pos" 
		<< setw(11) << "elec_s" << setw(11) << "elec_t" << setw(11) << "ion" << setw(11) << "hei" << setw(11) << "ecoul" << setw(11) << "ndiss" 
		<< setw(11) << "ncoll" << endl;

	// losses on vibrational excitation of H2 include pure rotational excitation, 
	for (i = 0; i < (int)enloss_rates_mt.size(); i++) {
		tot = enloss_rates_mt[i] + enloss_rates_h2_vibr[i] + enloss_rates_h2_singlet[i] + enloss_rates_h2_triplet[i]
			+ enloss_rates_ioniz[i] + enloss_rates_coulomb_el[i] + enloss_rates_hei[i];
	
		output.precision(3);
		output << left << setw(13) << time_moments[i];

		output.precision(2);
		output << left << setw(11) << tot 
			<< setw(11) << enloss_rates_mt[i]
			<< setw(11) << enloss_rates_h2_rot[i]
			<< setw(11) << enloss_rates_h2_rot_pos[i]
			<< setw(11) << enloss_rates_h2_vibr[i]
			<< setw(11) << enloss_rates_h2_vibr_pos[i]
			<< setw(11) << enloss_rates_h2_singlet[i]
			<< setw(11) << enloss_rates_h2_triplet[i]
			<< setw(11) << enloss_rates_ioniz[i]
			<< setw(11) << enloss_rates_hei[i] 
			<< setw(11) << enloss_rates_coulomb_el[i]
			<< setw(11) << diss_decay_heating_rate[i]
			<< setw(11) << neut_heat_coll_rates[i] << endl;
	}
	output.close();
}


void save_electron_energy_losses(const string& output_path, double grb_distance, double hcolumn_dens, double conc_h_tot, const vector<double>& time_moments,
	const vector<double>& enloss_mt,
	const vector<double>& enloss_h2_rot,
	const vector<double>& enloss_h2_vibr,
	const vector<double>& enloss_h2_singlet,
	const vector<double>& enloss_h2_triplet,
	const vector<double>& enloss_ioniz,
	const vector<double>& enloss_hei,
	const vector<double>& enloss_coulomb_el,
	const vector<double>& diss_decay_heating,
	const vector<double>& neut_heat_coll)
{
	int i;
	double tot_int;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_electron_enloss_int.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);

	output << left << "!Energy losses of electrons, parameter < 0 if electrons lose energy, integrated up to a given time [eV cm-3]," << endl
		<< "!tot_int - total integrated losses, " << endl
		<< "!mt_int - momentum transfer, " << endl
		<< "!rot_int - H2 rotational excitation," << endl
		<< "!vibr_int - H2 ro-vibrational excitation (within the ground electronic state, including pure rotational transitions)," << endl
		<< "!elec_s_int - H2 electronic states excitation (singlet)," << endl
		<< "!elec_t_int - H2 electronic states excitation (triplet)," << endl
		<< "!ion_int - molecule/atom ionization," << endl
		<< "!hei_int - HeI excitation by electron impact" << endl
		<< "!ecoul_int - Coulomb loses on electrons," << endl
		<< "!ndiss_int - neutral heating rate due to H2 dissociation, > 0" << endl
		<< "!ncoll_int - neutral heating rate due to collisions H2-H2, H2-He," << endl;

	output << left << "!distance from the GRB source to cloud layer centre [pc], H column density from the cloud boundary to layer centre [cm-2], H nuclei concentration [cm-3], nb of times," << endl
		<< "! " << setw(15) << grb_distance / PARSEC << setw(15) << hcolumn_dens << setw(15) << conc_h_tot 
		<< setw(7) << (int)enloss_mt.size() << endl;

	output << left << setw(13) << "!time(s)" 
		<< setw(11) << "tot_int"
		<< setw(11) << "mt_int"
		<< setw(11) << "rot_int"
		<< setw(11) << "vibr_int"
		<< setw(11) << "elec_s_int"
		<< setw(11) << "elec_t_int"
		<< setw(11) << "ion_int"
		<< setw(11) << "hei_int" 
		<< setw(11) << "ecoul_int" 
		<< setw(11) << "ndiss_int"
		<< setw(11) << "ncoll_int" << endl;

	// losses on vibrational excitation of H2 include pure rotational excitation, 
	for (i = 0; i < (int)enloss_mt.size(); i++) {
		tot_int = enloss_mt[i] + enloss_h2_vibr[i] + enloss_h2_singlet[i] + enloss_h2_triplet[i] + enloss_ioniz[i] 
			+ enloss_coulomb_el[i] + enloss_hei[i];

		output.precision(3);
		output << left << setw(13) << time_moments[i];

		output.precision(2);
		output << left << setw(11) << tot_int
			<< setw(11) << enloss_mt[i]
			<< setw(11) << enloss_h2_rot[i]
			<< setw(11) << enloss_h2_vibr[i]
			<< setw(11) << enloss_h2_singlet[i]
			<< setw(11) << enloss_h2_triplet[i]
			<< setw(11) << enloss_ioniz[i]
			<< setw(11) << enloss_hei[i] 
			<< setw(11) << enloss_coulomb_el[i]
			<< setw(11) << diss_decay_heating[i]
			<< setw(11) << neut_heat_coll[i] << endl;
	}
	output.close();
}


void save_diss_excit_rates(const string& output_path, double grb_distance, double hcolumn_dens, double conc_h_tot, const vector<double>& time_moments,
	const vector<double>& h2_excit_electr_rate_arr, 
	const vector<double>& h2_excit_vibr_rate_arr, 
	const vector<double>& h2_excit_rot_rate_arr)
{
	int i;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_diss_excit_rates.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);

	output << left << "!Dissociation / excitation rates," << endl
		<< "!h2_exc_el- H2 electronic excitation rate, in [cm-3 s-1]," << endl
		<< "!h2_exc_v - H2 vibrational excitation rate, vi = 0 -> vf > 0 (within the ground electronic state), in [cm-3 s-1]," << endl
		<< "!h2_exc_r - H2 pure rotational excitation within v = 0, ji = 0,1 -> jf > 1, in [cm-3 s-1]," << endl;

	output << left << "!distance from the GRB source to cloud layer centre [pc], H column density from the cloud boundary to layer centre [cm-2], H nuclei concentration [cm-3], nb of times," << endl
		<< "! " << setw(15) << grb_distance / PARSEC << setw(15) << hcolumn_dens << setw(15) << conc_h_tot
		<< setw(7) << (int)h2_excit_electr_rate_arr.size() << endl;

	output << left << setw(13) << "!time(s)" << setw(12) << "h2_exc_el" << setw(12) << "h2_exc_v" << setw(12) << "h2_exc_r" << endl;

	for (i = 0; i < (int)h2_excit_electr_rate_arr.size(); i++) {
		output.precision(3);
		output << left << setw(13) << time_moments[i];

		output.precision(2);
		output << left << setw(12) << h2_excit_electr_rate_arr[i]
			<< setw(12) << h2_excit_vibr_rate_arr[i]
			<< setw(12) << h2_excit_rot_rate_arr[i] << endl;
	}
	output.close();
}


void save_diss_excit(const string& output_path, double grb_distance, double hcolumn_dens, double conc_h_tot, 
	const vector<double>& time_moments,
	const vector<double>& h2_solomon_diss_arr, 
	const vector<double>& h2_diss_exc_singlet_arr, 
	const vector<double>& h2_diss_exc_triplet_arr,
	const vector<double>& hei_exc_arr,
	const vector<double>& h2_excit_electr_arr,
	const vector<double>& h2_excit_electr_bs_arr,
	const vector<double>& h2_excit_electr_cp_arr,
	const vector<double>& h2_excit_vibr_arr,
	const vector<double>& h2_excit_vibr_1_arr,
	const vector<double>& h2_excit_vibr_2_arr)
{
	int i;
	double tot_int;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_diss_excit_int.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);

	output << left << "!Number of dissociations/excitations up to a given time [cm-3]," << endl
		<< "!tot_int(H2) - integrated amount of H2 dissociated [cm-3]," << endl
		<< "!sol_int  - dissociation following excitation to the singlet states (Solomon) [cm-3], " << endl
		<< "!de_s_int - dissociative excitation of H2 (via singlet states)," << endl
		<< "!de_t_int - dissociative excitation of H2 (via triplet states)," << endl
		<< "!hei_int  - HeI excitation [cm-3];" << endl
		<< "!el_int   - excitation to all H2 electronic states (that are taken into account in the model);" << endl
		<< "!el_bs_int- excitation to B1S+u H2 electronic state;" << endl
		<< "!el_cp_int- excitation to C1Pu H2 electronic state;" << endl
		<< "!vibr_int - excitation to vibrational states within the ground electronic state, vi = 0 -> vf > 0;" << endl
		<< "!v1_int   - excitation vi = 0 -> vf = 1 within the ground electronic state;" << endl
		<< "!v2_int   - excitation vi = 0,1 -> vf = 2 within the ground electronic state;" << endl;

	output << left << "!distance from the GRB source to cloud layer centre [pc], H column density from the cloud boundary to layer centre [cm-2], H nuclei concentration [cm-3], nb of times," << endl
		<< "! " << setw(15) << grb_distance / PARSEC << setw(15) << hcolumn_dens << setw(15) << conc_h_tot 
		<< setw(7) << (int)h2_solomon_diss_arr.size() << endl;

	output << left << setw(13) << "!time(s)" << setw(12) << "tot_int(H2)"
		<< setw(12) << "sol_int" << setw(12) << "de_s_int" << setw(12) << "de_t_int" << setw(12) << "hei_int"
		<< setw(12) << "el_int" << setw(12) << "el_bs_int" << setw(12) << "el_cp_int" 
		<< setw(12) << "vibr_int" << setw(12) << "v1_int" << setw(12) << "v2_int" << endl;

	for (i = 0; i < (int)h2_solomon_diss_arr.size(); i++) {
		tot_int = h2_solomon_diss_arr[i] + h2_diss_exc_singlet_arr[i] + h2_diss_exc_triplet_arr[i];

		output.precision(3);
		output << left << setw(13) << time_moments[i];

		output.precision(2);
		output << left << setw(12) << tot_int
			<< setw(12) << h2_solomon_diss_arr[i]
			<< setw(12) << h2_diss_exc_singlet_arr[i]
			<< setw(12) << h2_diss_exc_triplet_arr[i]
			<< setw(12) << hei_exc_arr[i] 
			<< setw(12) << h2_excit_electr_arr[i]
			<< setw(12) << h2_excit_electr_bs_arr[i]
			<< setw(12) << h2_excit_electr_cp_arr[i] 
			<< setw(12) << h2_excit_vibr_arr[i] 
			<< setw(12) << h2_excit_vibr_1_arr[i] 
			<< setw(12) << h2_excit_vibr_2_arr[i] << endl;
	}
	output.close();
}


void save_phys_parameters(const string& output_path, const vector<double>& time_moments,
	const vector<double>& neut_temp_arr, const vector<double>& ion_temp_arr, const vector<double>& dust_charge_arr,
	const vector<dynamic_array>& conc_data, double grb_distance, double hcolumn_dens, double conc_h_tot)
{
	int i;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_physparam_data.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(6);

	output << left << "!Physical parameters," << endl
		<< "!temp_n - temperature of neutral component [K]," << endl
		<< "!temp_i - ion temperature [K]," << endl
		<< "!dcharge - charge of dust grains [cm-3]," << endl
		<< "!conc_e - concentration of electrons [cm-3]" << endl;

	output << left << "!distance from the GRB source cloud layer centre [pc], H column density from the cloud boundary to layer centre [cm-2], H nuclei concentration [cm-3], nb of times," << endl
		<< "! " << setw(15) << grb_distance / PARSEC << setw(15) << hcolumn_dens << setw(15) << conc_h_tot 
		<< setw(7) << (int)neut_temp_arr.size() << endl;

	output << left << setw(13) << "!time(s)" << setw(12) << "temp_n(K)"
		<< setw(12) << "temp_i(K)" << setw(12) << "dcharge" << setw(12) << "conc_e" << endl;

	for (i = 0; i < (int)neut_temp_arr.size(); i++) {
		output.precision(3);
		output << left << setw(13) << time_moments[i];

		output.precision(2);
		output << left << setw(12) << neut_temp_arr[i]
			<< setw(12) << ion_temp_arr[i]
			<< setw(12) << dust_charge_arr[i]
			<< setw(12) << conc_data[i].arr[EL_NB] << endl;
	}
	output.close();
}
