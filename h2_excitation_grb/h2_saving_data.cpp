
#define _CRT_SECURE_NO_WARNINGS

#include <omp.h>
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


void save_model_parameters(const std::string& output_path,  double conc_h_tot, double op_ratio_h2, double ioniz_fract, double dust_gas_mass_ratio, 
	double grain_radius, double grain_nb_density, double grb_cloud_distance, double grb_distance, double hcolumn_dens, int layer_nb)
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

	output << "!H2-electron propagation parameters, date " << ctime_str;
	if (layer_nb >= 0) {
		output << left << "!Molecular cloud parameters:" << endl
			<< "!Cloud layer number:" << endl << layer_nb << endl
			<< "!Distance from the GRB source to the cloud boundary [pc]:" << endl << grb_cloud_distance / PARSEC << endl
			<< "!Distance from the GRB source to the centre of the cloud layer [pc]:" << endl << grb_distance / PARSEC << endl
			<< "!H column density passed by the GRB, from the cloud boundary to the cloud layer centre [cm-2]:" << endl << hcolumn_dens << endl;
	}
	output << left << "!Gas parameters:" << endl
		<< "!H nuclei concentration [cm-3]:" << endl << conc_h_tot << endl
		<< "!Initial H2 concentration [cm-3]:" << endl << 0.5 * conc_h_tot << endl
		<< "!Initial H2 OPR:" << endl << op_ratio_h2 << endl;
	
	output << left << "!Coulomb losses (0/1):" << endl << CALC_EL_LOSSES_THERMAL_EL << endl
		<< "!Ionization fraction of thermal electrons:" << endl << ioniz_fract << endl
		<< "!Thermal electron temperature [K]:" << endl << THERMAL_EL_TEMPERATURE << endl;
	
	output << left << "!Collisions of H2 with H2, He included (0/1):" << endl << H2_COLLISIONS_WITH_H2_HE << endl
		<< "!Temperature of the neutral gas in the cloud [K]:" << endl << THERMAL_NEUTRAL_TEMPERATURE << endl;

	output << left << "!Nb of H2 vibrational states taking into account in electron-impact excitation:" << endl << USED_NB_OF_H2_VSTATES_X1SG << endl;
	
	output << left << "!Dust parameters:" << endl
		<< "!Dust-gas mass ratio:" << endl << dust_gas_mass_ratio << endl
		<< "!Initial grain radius [cm]:" << endl << grain_radius << endl
		<< "!Grain number density [cm-3]:" << endl << grain_nb_density << endl;	

	output << left << "!Nb of bins per order in the electron energy grid (default 100):" << endl
		<< NB_OF_BINS_PER_ORDER_EL << endl;

#ifdef _OPENMP
	output << left << "!Nb of processors in OpenMP:" << endl 
		<< omp_get_num_threads() << endl;
#else
	output << left << "!No OpenMP" << endl << "1" << endl;
#endif

	output.close();
}

void save_electron_spectrum_evolution(const string& output_path, const vector<double>& electron_energies_grid, const vector<double>& electron_energy_bin_size,
	const vector<double>& time_moments, const vector<dynamic_array>& spectrum_data, double conc_h_tot)
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
		// << "!first row: distance from the GRB source to cloud layer centre [pc], H column density from the cloud boundary to layer centre [cm-2] <<
		<< "!first row: H nuclei concentration [cm-3], nb of rows (= nb of el. energies) and nb of times," << endl
		<< "!second row: time (s)," << endl
		<< "!data: col.1 - centre of the energy bin [eV], col.2 - bin semi-length [eV], col.3 - spectrum N(E) [cm-3 eV-1]," << endl;

	output << scientific;
	output.precision(6);
	//output << left << "! " << setw(15) << grb_distance / PARSEC << setw(15) << hcolumn_dens;
	output << left << "! " << setw(15) << conc_h_tot << setw(7) << nb_of_el_energies << setw(7) << nb_of_times << endl;

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
		<< "!first row: H nuclei concentration [cm-3], nb of rows (= nb of times)," << endl
		<< "!data: time (s), conc of electrons (from spectrum) [cm-3], total energy [eV cm-3], median energy [eV]," << endl;

	output << scientific;
	output.precision(6);
	output << left << "! " << setw(15) << conc_h_tot << setw(7) << nb_of_times << endl;

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


void save_specimen_conc_evolution(const string& output_path, const vector<double>& time_moments, const vector<dynamic_array>& conc_data, 
	double conc_h_tot)
{
	int i, l;
	double w;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_specimen_conc.txt";
	output.open(fname.c_str());

	output << left << "!Evolution of concentrations of chemical species," << endl
		<< "!first row: H nuclei concentration [cm-3], nb of times, number of species," << endl
		<< "!second row: specimen list," << endl;

	output << scientific;
	output.precision(6);
	output << left << "! " << setw(15) << conc_h_tot << setw(7) << (int)conc_data.size() << setw(7) << NB_OF_CHEM_SPECIES << endl;
	
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
	const vector<dynamic_array>& h2_popdens_data, const vector<dynamic_array>& h2_popdens_v_data, double conc_h_tot, int nb_lev_h2)
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
		<< "!first row: H nuclei concentration [cm-3], nb of times, number of levels," << endl
		<< "!second row: level nb," << endl;

	output << left << "! " << setw(15) << conc_h_tot << setw(7) << (int)h2_popdens_data.size() << setw(7) << nb_lev_h2 << endl;
	
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
		<< "!first row: H nuclei concentration [cm-3], nb of times, number of vibr states," << endl
		<< "!second row: summed value, v state nb," << endl;

	output << left << "! " << setw(15) << conc_h_tot << setw(7) << (int)h2_popdens_v_data.size() << setw(7) << MAX_NB_H2_VSTATES_X1SG << endl;

	output << left << setw(13) << "!time(s)" << setw(13) << "sum";
	for (v = 0; v < MAX_NB_H2_VSTATES_X1SG; v++) {
		output << left << setw(12) << v;
	}
	output << endl;

	for (l = 0; l < (int)h2_popdens_v_data.size(); l++) {
		output.precision(4);
		output << left << setw(13) << time_moments[l];

		// summed over the all vibrational states, is equal to H2 concentration,
		w = 0.;
		for (v = 0; v < MAX_NB_H2_VSTATES_X1SG; v++) {
			w += h2_popdens_v_data[l].arr[v];
		}
		output.precision(3);
		output << left << setw(13) << w;

		for (v = 0; v < MAX_NB_H2_VSTATES_X1SG; v++)
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
	const vector<dynamic_array>& hei_popdens_data, double conc_h_tot, int nb_lev_hei)
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
		<< "!first row: H nuclei concentration [cm-3], nb of times, number of HeI levels," << endl
		<< "!second row: level nb," << endl;

	output << left << "! " << setw(15) << conc_h_tot << setw(7) << (int)hei_popdens_data.size() << setw(7) << nb_lev_hei << endl;

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


void save_electron_energy_loss_rates(const string& output_path, double conc_h_tot, const vector<double>& time_moments,
	const vector<double>& enloss_rates_mt,
	const vector<double>& enloss_rates_h2_rot, 
	const vector<double>& enloss_rates_h2_rot_pos,
	const vector<double>& enloss_rates_h2_vibr, 
	const vector<double>& enloss_rates_h2_vibr_pos,
	const vector<double>& enloss_rates_ioniz, 
	const vector<double>& enloss_rates_hei,
	const vector<double>& enloss_rates_coulomb_el,
	const vector<double>& neut_heat_coll_rates, 
	vector< array<electronic_excitation_data_unit, NB_EXC_ELECTRONIC_STATES>>& h2_state_data_rate_arr)
{
	int i, j;
	double tot, enloss_sin, enloss_tri, diss_heat_input;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_electron_enloss_rates.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(3);

	output << left << "!Energy loss rates of electrons, parameter < 0 if electrons lose energy, rate in [eV cm-3 s-1]," << endl
		<< "!tot - total rate, " << endl
		<< "!rot      - H2 rotational excitation, only transitions with the same vibrational quantum nb," << endl
		<< "!rot_pos  - H2 rotational excitation, electrons gain energy," << endl
		<< "!vibr     - H2 ro-vibrational excitation (within the ground electronic state, including pure rotational transitions)," << endl
		<< "!vibr_pos - H2 ro-vibrational excitation, electrons gain energy" << endl
		<< "!elec_s - H2 electronic states excitation (singlet)," << endl
		<< "!elec_t - H2 electronic states excitation (triplet)," << endl
		<< "!ion - molecule/atom ionization," << endl
		<< "!hei - HeI excitation by electron impact" << endl
		
		<< "!mt - momentum transfer, " << endl
		<< "!ecoul - Coulomb loses on electrons," << endl
		<< "!ncoll - neutral heating rate due to collisions H2-H2, H2-He," << endl 
		<< "!ndiss - neutral heating rate due to H2 dissociation (two-step dissociation, triplet excitation), > 0" << endl;

	output << left << "!H nuclei concentration [cm-3], nb of times," << endl
		<< "! " << setw(15) << conc_h_tot << setw(7) << (int)enloss_rates_mt.size() << endl;
	
	output << left << setw(12) << "!time(s)"
		<< setw(12) << "tot" << setw(12) << "rot" << setw(12) << "rot_pos" << setw(12) << "vibr" << setw(12) << "vibr_pos" 
		<< setw(12) << "elec_s" << setw(12) << "elec_t" << setw(12) << "ion" << setw(12) << "hei" << setw(12) << "mt" << setw(12) << "ecoul" 
		<< setw(12) << "ncoll" << setw(12) << "ndiss" << endl;

	// losses on vibrational excitation of H2 include pure rotational excitation, 
	for (i = 0; i < (int)enloss_rates_mt.size(); i++) 
	{
		enloss_sin = enloss_tri = diss_heat_input = 0.;
		for (j = 0; j < NB_EXC_ELECTRONIC_SINGLET_STATES; j++) {
			enloss_sin += (h2_state_data_rate_arr[i])[j].enloss_excit + (h2_state_data_rate_arr[i])[j].enloss_direct_diss;
		}

		for (j = NB_EXC_ELECTRONIC_SINGLET_STATES; j < NB_EXC_ELECTRONIC_STATES; j++) {
			enloss_tri += (h2_state_data_rate_arr[i])[j].enloss_excit + (h2_state_data_rate_arr[i])[j].enloss_direct_diss;
		}
		
		for (j = 0; j < NB_EXC_ELECTRONIC_STATES; j++) {
			diss_heat_input += (h2_state_data_rate_arr[i])[j].diss_heat_input;
		}

		tot = enloss_rates_mt[i] + enloss_rates_h2_vibr[i] + enloss_rates_ioniz[i] + enloss_rates_coulomb_el[i] + enloss_rates_hei[i]
			+ enloss_sin + enloss_tri;

		output << left << setw(12) << time_moments[i];

		output << left << setw(12) << tot 
			<< setw(12) << enloss_rates_h2_rot[i]
			<< setw(12) << enloss_rates_h2_rot_pos[i]
			<< setw(12) << enloss_rates_h2_vibr[i]
			<< setw(12) << enloss_rates_h2_vibr_pos[i]
			<< setw(12) << enloss_sin
			<< setw(12) << enloss_tri
			<< setw(12) << enloss_rates_ioniz[i]
			<< setw(12) << enloss_rates_hei[i] 
			<< setw(12) << enloss_rates_mt[i]
			<< setw(12) << enloss_rates_coulomb_el[i]
			<< setw(12) << neut_heat_coll_rates[i]
			<< setw(12) << diss_heat_input << endl;
	}
	output.close();
}


void save_electron_energy_losses(const string& output_path, double conc_h_tot, const vector<double>& time_moments,
	const vector<double>& enloss_mt,
	const vector<double>& enloss_h2_rot,
	const vector<double>& enloss_h2_vibr,
	const vector<double>& enloss_ioniz,
	const vector<double>& enloss_hei,
	const vector<double>& enloss_coulomb_el,
	const vector<double>& neut_heat_coll, 
	vector< array<electronic_excitation_data_unit, NB_EXC_ELECTRONIC_STATES> >& h2_state_data_arr)
{
	int i, j;
	double tot_int, enloss_sin, enloss_tri, diss_heat_input;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_electron_enloss_int.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(3);

	output << left << "!Energy losses of electrons, parameter < 0 if electrons lose energy, integrated up to a given time [eV cm-3]," << endl
		<< "!tot_int - total integrated losses," << endl
		<< "!rot_int - pure H2 rotational excitation, only transitions with the same vibrational quantum nb," << endl
		<< "!vibr_int - H2 ro-vibrational excitation (within the ground electronic state, including pure rotational transitions)," << endl
		<< "!elec_s_int - due to H2 electronic states excitation (singlet), including dissociative excitation," << endl
		<< "!elec_t_int - due to H2 electronic states excitation (triplet), including dissociative excitation," << endl
		<< "!ion_int - molecule/atom ionization," << endl
		<< "!hei_int - HeI excitation by electron impact," << endl
		
		<< "!mt_int - momentum transfer, " << endl
		<< "!ecoul_int - Coulomb loses on electrons," << endl
		<< "!ncoll_int - neutral heating rate due to collisions H2-H2, H2-He," << endl
		<< "!ndiss_int - neutral heating rate due to H2 dissociation (two-step dissociation, triplet excitation), > 0" << endl;

	output << left << "! H nuclei concentration [cm-3], nb of times," << endl
		<< "! " << setw(15) << conc_h_tot << setw(7) << (int)enloss_mt.size() << endl;

	output << left << setw(12) << "!time(s)" 
		<< setw(12) << "tot_int"
		<< setw(12) << "rot_int"
		<< setw(12) << "vibr_int"
		<< setw(12) << "elec_s_int"
		<< setw(12) << "elec_t_int"
		<< setw(12) << "ion_int"
		<< setw(12) << "hei_int"

		<< setw(12) << "mt_int"
		<< setw(12) << "ecoul_int"
		<< setw(12) << "ncoll_int" 
		<< setw(12) << "ndiss_int"<< endl;

	// losses on vibrational excitation of H2 include pure rotational excitation, 
	for (i = 0; i < (int)enloss_mt.size(); i++) 
	{
		enloss_sin = enloss_tri = diss_heat_input = 0.;
		for (j = 0; j < NB_EXC_ELECTRONIC_SINGLET_STATES; j++) {
			enloss_sin += (h2_state_data_arr[i])[j].enloss_excit + (h2_state_data_arr[i])[j].enloss_direct_diss;
		}

		for (j = NB_EXC_ELECTRONIC_SINGLET_STATES; j < NB_EXC_ELECTRONIC_STATES; j++) {
			enloss_tri += (h2_state_data_arr[i])[j].enloss_excit + (h2_state_data_arr[i])[j].enloss_direct_diss;
		}

		for (j = 0; j < NB_EXC_ELECTRONIC_STATES; j++) {
			diss_heat_input += (h2_state_data_arr[i])[j].diss_heat_input;
		}

		tot_int = enloss_mt[i] + enloss_h2_vibr[i] + enloss_ioniz[i] + enloss_coulomb_el[i] + enloss_hei[i] 
			+ enloss_sin + enloss_tri;

		output << left << setw(12) << time_moments[i];

		output << left << setw(12) << tot_int
			<< setw(12) << enloss_h2_rot[i]
			<< setw(12) << enloss_h2_vibr[i]
			<< setw(12) << enloss_sin
			<< setw(12) << enloss_tri
			<< setw(12) << enloss_ioniz[i]
			<< setw(12) << enloss_hei[i] 
			<< setw(12) << enloss_mt[i]
			<< setw(12) << enloss_coulomb_el[i]
			<< setw(12) << neut_heat_coll[i] 
			<< setw(12) << diss_heat_input << endl;
	}
	output.close();
}


void save_electronic_states_excit_rates(const string& output_path, double conc_h_tot, const vector<double>& time_moments,
	vector< array<electronic_excitation_data_unit, NB_EXC_ELECTRONIC_STATES>>& h2_state_data_rate_arr)
{
	int i, j;
	double tot_excit;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_electr_excit_rates.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(3);

	output << left << "!Rates of dissociations/excitations for each of the electronic states [cm-3 s-1]," << endl
		<< "!total - rate of excitation to all H2 electronic states (WITH excitations to dissociation continuum);" << endl
		<< "!Electronic states: B, C, Bp, D, EF, Bpp, Dp, b, a, c, d" << endl
		<< "!Parameters: 1. excitation (no direct dissociation); 2. direct dissociation; 3. twostep dissociation" << endl;

	output << left << "!H nuclei concentration [cm-3], nb of times," << endl
		<< "! " << setw(15) << conc_h_tot << setw(7) << (int)h2_state_data_rate_arr.size() << endl;

	output << left << setw(12) << "!time(s)" << setw(12) << "total";

	for (j = 0; j < NB_EXC_ELECTRONIC_STATES; j++) {
		output << left << setw(12) << H2_ELECTRONIC_STATE_NAMES[j] + "_1" 
			<< setw(12) << H2_ELECTRONIC_STATE_NAMES[j] + "_2"
			<< setw(14) << H2_ELECTRONIC_STATE_NAMES[j] + "_3";
	}
	output << endl;

	for (i = 0; i < (int)h2_state_data_rate_arr.size(); i++) 
	{
		tot_excit = 0.;
		for (j = 0; j < NB_EXC_ELECTRONIC_STATES; j++) {
			tot_excit += (h2_state_data_rate_arr[i])[j].excit + (h2_state_data_rate_arr[i])[j].direct_diss;
		}

		output << left << setw(12) << time_moments[i] << setw(12) << tot_excit;

		for (j = 0; j < NB_EXC_ELECTRONIC_STATES; j++) {
			output << left << setw(12) << (h2_state_data_rate_arr[i])[j].excit
				<< setw(12) << (h2_state_data_rate_arr[i])[j].direct_diss
				<< setw(14) << (h2_state_data_rate_arr[i])[j].twostep_diss;
		}
		output << endl;
	}
	output.close();
}


void save_vibrational_states_excit_rates(const string& output_path, double conc_h_tot, const vector<double>& time_moments,
	const vector<dynamic_array>& h2_electr_vstates_rate_arr,
	const vector<dynamic_array>& h2_vibr_vstates_rate_arr, 
	const vector<double>& h2_excit_rot_rate_arr)
{
	int i, j;
	string fname;
	ofstream output;

	fname = output_path + "h2grb_vibr_excit_rates.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(3);

	output << left << "!Excitation of vibrational states of the ground electronic state," << endl
		<< "!First parameter - pure rotational excitation [cm-3 s-1], from v= 0, j = 0 or j = 1, to level with j > 1" << endl
		<< "!For each vibration state (0,1,...14) two parameters:" << endl
		<< "!1. Vibrational excitation rate, vi -> vf, vi < vf (within the ground electronic state) [cm-3 s-1]," << endl
		<< "!2. Excitation rate of vibration state through electronic excitation [cm-3 s-1]," << endl;

	output << left << "!H nuclei concentration [cm-3], nb of times," << endl
		<< "! " << setw(15) << conc_h_tot << setw(7) << (int)h2_vibr_vstates_rate_arr.size() << endl;

	output << left << setw(12) << "!time(s)" << setw(14) << "rot";
	for (j = 0; j < MAX_NB_H2_VSTATES_X1SG; j++) {
		output << left << setw(12) << j << setw(14) << j;
	}
	output << endl;

	for (i = 0; i < (int)h2_vibr_vstates_rate_arr.size(); i++) {
		output << left << setw(12) << time_moments[i] << setw(14) << h2_excit_rot_rate_arr[i];

		for (j = 0; j < MAX_NB_H2_VSTATES_X1SG; j++) {
			output << left << setw(12) << h2_vibr_vstates_rate_arr[i].arr[j] << setw(14) << h2_electr_vstates_rate_arr[i].arr[j];
		}
		output  << endl;
	}
	output.close();
}


void save_output_parameters(const string& output_path, double conc_h_tot, const vector<double>& electron_energies_grid, const vector<double>& electron_energy_bin_size,
	const vector<dynamic_array>& spectrum_data,
	const vector<dynamic_array>& conc_data,
	const vector<double>& enloss_mt_arr,
	const vector<double>& enloss_coulomb_el_arr,
	const vector<double>& neutral_coll_heating_arr,
	const vector<double>& enloss_h2_vibr_arr,
	vector< array<electronic_excitation_data_unit, NB_EXC_ELECTRONIC_STATES>>& h2_state_data_arr,
	const vector<dynamic_array> & h2_electr_vstates_arr,
	const vector<dynamic_array> & h2_vibr_vstates_arr, 
	const std::vector<double>& h2_excit_rot_arr, 
	const std::vector<double>& hei_excit_arr)
{
	bool do_files_exist = false;
	int j, nb_of_el_energies;
	double energy, nb_density, nb_ions, heating_effic, diss_heat_singlet, diss_heat_triplet, rovibr_losses_effic, total_electr_excit, total_diss;

	string fname;
	ofstream output;

	nb_of_el_energies = (int)electron_energy_bin_size.size();

	// Energy in electrons, number density of electrons at the start of simulations,
	energy = nb_density = 0.;
	for (j = 0; j < nb_of_el_energies; j++) {
		energy += spectrum_data.front().arr[j] * 0.5 * (electron_energies_grid[j] + electron_energies_grid[j + 1]);  // [eV cm-3]
		nb_density += spectrum_data.front().arr[j];  // [cm-3]
	}

	// simple species are: "e-", "H", "H+", "H2", "H2+", "He", "He+", "He++"
	// number density of ions at the start of the simulations,
	nb_ions = -(conc_data.front().arr[2] + conc_data.front().arr[4] + conc_data.front().arr[6]
		+ 2. * conc_data.front().arr[7]);  // [cm-3]

	// number density of ions produced, at the end of the simulations,
	nb_ions += conc_data.back().arr[2] + conc_data.back().arr[4] + conc_data.back().arr[6]
		+ 2. * conc_data.back().arr[7];  // [cm-3]

	heating_effic = (fabs(enloss_mt_arr.back()) + fabs(enloss_coulomb_el_arr.back()) + neutral_coll_heating_arr.back()) / energy;  // dimensionless
	rovibr_losses_effic = fabs(enloss_h2_vibr_arr.back()) / energy;  // dimensionless

	diss_heat_singlet = diss_heat_triplet = 0.;
	for (j = 0; j < NB_EXC_ELECTRONIC_SINGLET_STATES; j++) {
		diss_heat_singlet += (h2_state_data_arr.back())[j].diss_heat_input;  // [eV cm-3]
	}
	for (j = NB_EXC_ELECTRONIC_SINGLET_STATES; j < NB_EXC_ELECTRONIC_STATES; j++) {
		diss_heat_triplet += (h2_state_data_arr.back())[j].diss_heat_input;  // [eV cm-3]
	}

	diss_heat_singlet /= energy;  // dimensionless
	diss_heat_triplet /= energy;  // dimensionless
	
	// Saving the general parameters of the electron energy degradation,
	fname = output_path + "h2grb_output_parameters.txt";
	
	// Checking if files exist
	ifstream fcheck(fname.c_str());
	if (fcheck.good()) {
		do_files_exist = true;
	}
	else {
		do_files_exist = false;
	}
	
	output.open(fname.c_str(), ios_base::app);
	output << scientific;
	output.precision(4);

	if (!do_files_exist) {
		output << left << "!Parameters of the electron energy degradation," << endl
			<< "!H nuclei concentration [cm-3]: " << conc_h_tot << "; initial number of electrons [cm-3]: " << nb_density << endl
			<< "!1. Initial electron energy E [eV]" << endl
			<< "!2. Number density of ions N [per electron]" << endl
			<< "!3. Mean energy per ion pair W [eV]" << endl
			<< "!4. 1/W [eV-1]" << endl
			<< "!5. Gas heating efficiency (momentum transfer, Coulomb, collisions with neutrals), dimensionless" << endl
			<< "!6. Gas heating through H2 dissociation (singlet), dimensionless" << endl
			<< "!7. Gas heating through H2 dissociation (triplet), dimensionless" << endl
			<< "!8. Energy fraction lost in ro-vibr., dimensionless" << endl
			<< "!9. Helium excitations [per electron]" << endl;

		output << "!";
		for (j = 1; j <= 9; j++) {
			output << setw(13) << j;
		}
		output << endl;
	}
	output << left << setw(13) << energy / nb_density
		<< setw(13) << nb_ions / nb_density
		<< setw(13) << energy / nb_ions
		<< setw(13) << nb_ions / energy
		<< setw(13) << heating_effic
		<< setw(13) << diss_heat_singlet
		<< setw(13) << diss_heat_triplet
		<< setw(13) << rovibr_losses_effic 
		<< setw(13) << hei_excit_arr.back() / nb_density << endl;
	output.close();

	// rate of excitation to all H2 electronic states(WITH excitations to dissociation continuum) [cm-3]
	total_electr_excit = total_diss = 0.;

	for (j = 0; j < NB_EXC_ELECTRONIC_STATES; j++) {
		total_electr_excit += (h2_state_data_arr.back())[j].excit + (h2_state_data_arr.back())[j].direct_diss;
		total_diss += (h2_state_data_arr.back())[j].direct_diss + (h2_state_data_arr.back())[j].twostep_diss;
	}

	// Excitation of H2 electronic states
	fname = output_path + "h2grb_output_electr_excit.txt";
	output.open(fname.c_str(), ios_base::app);

	output << scientific;
	output.precision(4);

	if (!do_files_exist) {
		output << left << "!Dissociations/excitations for each of the electronic states [per electron]," << endl
			<< "!H nuclei concentration [cm-3]: " << conc_h_tot << "; initial number of electrons [cm-3]: " << nb_density << endl
			<< "!1. Initial electron energy E [eV]" << endl
			<< "!2. total_exc  - Excitations to all H2 electronic states (WITH excitations to dissociation continuum), per primary electron;" << endl
			<< "!3. total_diss - Dissociations (excitations to dissociation continuum, all triplet states, Solomon process), per primary electron;" << endl
			<< "!Electronic states: B, C, Bp, D, EF, Bpp, Dp, b, a, c, d" << endl
			<< "!Parameters: 1. excitation (no direct dissociation); 2. direct dissociation; 3. two-step dissociation" << endl;

		output << left << setw(13) << "!1" << setw(13) << "2" << setw(13) << "3";
		for (j = 0; j < NB_EXC_ELECTRONIC_STATES; j++) {
			output << left << setw(13) << j + 3 << setw(13) << j + 4 << setw(15) << j + 5;
		}
		output << endl;

		output << left << setw(13) << "!Energy" << setw(13) << "total_exc" << setw(13) << "total_diss";
		for (j = 0; j < NB_EXC_ELECTRONIC_STATES; j++) {
			output << left << setw(13) << H2_ELECTRONIC_STATE_NAMES[j] + "_1"
				<< setw(13) << H2_ELECTRONIC_STATE_NAMES[j] + "_2"
				<< setw(15) << H2_ELECTRONIC_STATE_NAMES[j] + "_3";
		}
		output << endl;
	}

	output << left << setw(13) << energy / nb_density << setw(13) << total_electr_excit / nb_density << setw(13) << total_diss / nb_density;
	for (j = 0; j < NB_EXC_ELECTRONIC_STATES; j++) {
		output << left << setw(13) << (h2_state_data_arr.back())[j].excit / nb_density
			<< setw(13) << (h2_state_data_arr.back())[j].direct_diss / nb_density
			<< setw(15) << (h2_state_data_arr.back())[j].twostep_diss / nb_density;
	}
	output << endl;
	output.close();

	// Excitation of vibration states,
	fname = output_path + "h2grb_output_vibr_excit.txt";
	output.open(fname.c_str(), ios_base::app);

	output << scientific;
	output.precision(4);

	if (!do_files_exist) {
		output << left << "!Excitation of vibrational states of the ground electronic state," << endl
			<< "!H nuclei concentration [cm-3]: " << conc_h_tot << "; initial number of electrons [cm-3]: " << nb_density << endl
			<< "!1. Initial electron energy E [eV]" << endl
			<< "!There are two parameters for each vibration state v = 0,1,...14:" << endl
			<< "!1. v1 vibrational excitation, vi -> vf, vi < vf (within the ground electronic state) [per electron]," << endl
			<< "!for v = 0, H2 rotational excitation, from v= 0, j = 0 or j = 1, to level with j > 1 [per electron]," << endl
			<< "!2. v2 Excitation of vibration state through electronic excitation [per electron]," << endl;

		output << left << setw(13) << "!Energy" 
			<< setw(13) << "rotational" << setw(15) << "0";
		
		for (j = 1; j < MAX_NB_H2_VSTATES_X1SG; j++) {
			output << left << setw(13) << j << setw(15) << j;
		}
		output << endl;
	}
	
	output << left << setw(13) << energy / nb_density 
		<< setw(13) << h2_excit_rot_arr.back() / nb_density << setw(15) << h2_electr_vstates_arr.back().arr[0] / nb_density;
	
	for (j = 1; j < MAX_NB_H2_VSTATES_X1SG; j++) {
		output << left << setw(13) << h2_vibr_vstates_arr.back().arr[j] / nb_density 
			<< setw(15) << h2_electr_vstates_arr.back().arr[j] / nb_density;
	}
	output << endl;
	output.close();
}


// Not used
void save_phys_param_evolution(const string& output_path, double conc_h_tot, const vector<double>& time_moments,
	const vector<double>& neut_temp_arr, const vector<double>& ion_temp_arr, const vector<double>& dust_charge_arr,
	const vector<dynamic_array>& conc_data)
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

	output << left << "!H nuclei concentration [cm-3], nb of times," << endl
		<< "! " << setw(15) << conc_h_tot << setw(7) << (int)neut_temp_arr.size() << endl;

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
