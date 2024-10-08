
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>
#include <ctime>
#include <limits>

#include "h2_saving_data.h"
#include "h2_parameters.h"

using namespace std;

#define NB_OF_TIME_STEPS 4

void save_electron_spectrum_evolution(const string& output_path, const elspectra_evolution_data* user_data,
	const vector<double>& time_moments, const vector<dynamic_array>& spectrum_data)
{
	int i, l, n;
	double en, d_en, tot_en, tot_n;

	string fname;
	ofstream output;

	n = (int)spectrum_data.size();

	fname = output_path + "h2exc_electron_spectrum.txt";
	output.open(fname.c_str());

	output << "! Electron spectrum evolution," << endl
		<< "! first row: nb of rows (= nb of el. energies) and nb of times," << endl
		<< "! second row: time (s)," << endl
		<< "! data: col.1 - centre of the energy bin (eV), col.2 - bin semi-length (eV), col.3 - N(E), nb of electrons per eV per cm3," << endl;

	output << scientific;
	output.precision(4);
	output << left << "! " << setw(7) << user_data->get_nb_of_el_en() << setw(7) << n << endl;

	output << left << setw(26) << "!";
	for (l = 0; l < n; l++) {
		output << left << setw(13) << time_moments[l];
	}
	output << endl;

	for (i = 0; i < user_data->get_nb_of_el_en(); i++) {
		en = user_data->get_electron_energy(i);
		d_en = user_data->get_electron_energy_bin(i);
		output << left << setw(13) << en << setw(13) << 0.5 * d_en;

		for (l = 0; l < n; l++) {
			output << left << setw(13) << spectrum_data[l].arr[i] / d_en;
		}
		output << endl;
	}
	output.close();

	fname = output_path + "h2exc_electron_energy.txt";
	output.open(fname.c_str());

	output << "! Kinetic energy of electrons," << endl
		<< "! first row: nb of rows (= nb of times)," << endl
		<< "! data: time (s) - conc of electrons (cm-3) - total energy (eV cm-3) - median energy (eV)" << endl;

	output << scientific;
	output.precision(4);
	output << left << "! " << setw(7) << n << endl;

	for (l = 0; l < n; l++) {
		tot_n = tot_en = 0.;

		for (i = 0; i < user_data->get_nb_of_el_en(); i++) {
			en = user_data->get_electron_energy(i);
			tot_n += spectrum_data[l].arr[i];
			tot_en += spectrum_data[l].arr[i] * en;
		}
		output << left << setw(13) << time_moments[l] << setw(13) << tot_n << setw(13) << tot_en
			<< setw(13) << tot_en / tot_n << endl;
	}
	output.close();
}

void save_specimen_conc_evolution(const string& output_path, const vector<double>& time_moments,
	const vector<dynamic_array>& conc_data)
{
	int i, l;
	string fname;
	ofstream output;

	fname = output_path + "h2exc_specimen_conc.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(3);

	output << left << "! Evolution of abundances of chemical species," << endl
		<< "! first row - nb of times, number of species," << endl
		<< "! second row - specimen list," << endl;

	output << left << "! " << setw(7) << (int)conc_data.size() << setw(7) << NB_OF_CHEM_SPECIES << endl;
	output << left << setw(13) << "!";

	for (i = 0; i < NB_OF_CHEM_SPECIES; i++) {
		output << left << setw(12) << chemical_species[i];
	}
	output << endl;

	for (l = 0; l < (int)conc_data.size(); l++) {
		output.precision(4);
		output << left << setw(13) << time_moments[l];

		output.precision(3);
		for (i = 0; i < NB_OF_CHEM_SPECIES; i++) {
			output << left << setw(12) << conc_data[l].arr[i];
		}
		output << endl;
	}
	output.close();
}

void save_h2_populations_evolution(const string& output_path, const vector<double>& time_moments,
	const vector<dynamic_array>& h2_pop_data, const energy_diagram* h2_di)
{
	const int max_v = 15;  // the number of vibrational states (v = 0,1,2,.., max_v-1 )
	const int max_j = 32;

	int i, l, v, n;
	int index_arr[max_v][max_j];
	double h2_pop_v[max_v];

	string fname;
	ofstream output;

	for (v = 0; v < max_v; v++) {
		for (i = 0; i < max_j; i++) {
			index_arr[v][i] = -1;
		}
	}
	for (i = 0; i < h2_di->nb_lev; i++) {
		index_arr[h2_di->lev_array[i].v][rounding(h2_di->lev_array[i].j)] = i;
	}

	// Saving level population densities,
	fname = output_path + "h2exc_h2_populations.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(3);

	output << left << "! Evolution of population densities of H2 energy levels," << endl
		<< "! first row - nb of times, number of levels," << endl
		<< "! second row - level nb," << endl;

	output << left << "! " << setw(7) << (int)h2_pop_data.size() << setw(7) << h2_di->nb_lev << endl;
	output << left << setw(13) << "!";

	for (v = 0; v < max_v; v++) {
		for (i = 0; i < max_j; i++) {
			n = index_arr[v][i];
			// it is assumed that data size > 0
			if (n != -1 && n < h2_pop_data[0].dim)
				output << left << setw(4) << v << setw(8) << i;
		}
		output << setw(3) << "";
	}
	output << endl;

	for (l = 0; l < (int)h2_pop_data.size(); l++) {
		output.precision(4);
		output << left << setw(13) << time_moments[l];

		output.precision(3);
		for (v = 0; v < max_v; v++) {
			for (i = 0; i < max_j; i++) {
				n = index_arr[v][i];
				if (n != -1 && n < h2_pop_data[l].dim)
					output << left << setw(12) << h2_pop_data[l].arr[n];
			}
			output << setw(3) << "";
		}
		output << endl;
	}
	output.close();

	// Saving the populations of vibrational states,
	fname = output_path + "h2exc_h2_populations_vstates.txt";
	output.open(fname.c_str());

	output << scientific;
	output << left << "! Evolution of population densities of vibrational states of H2," << endl
		<< "! first row - nb of times, number of vibrational states max_v (0,1,2,.., max_v-1)," << endl
		<< "! second row - v state nb," << endl;

	output << left << "! " << setw(7) << (int)h2_pop_data.size() << setw(7) << max_v << endl;
	output << left << setw(13) << "!";
	for (v = 0; v < max_v; v++) {
		output << left << setw(12) << v;
	}
	output << endl;

	for (l = 0; l < (int)h2_pop_data.size(); l++) {
		output.precision(4);
		output << left << setw(13) << time_moments[l];

		memset(h2_pop_v, 0, max_v * sizeof(double));
		for (i = 0; i < h2_pop_data[l].dim; i++) {
			h2_pop_v[h2_di->lev_array[i].v] += h2_pop_data[l].arr[i];
		}

		output.precision(3);
		for (v = 0; v < max_v; v++) {
			output << left << setw(12) << h2_pop_v[v];
		}
		output << endl;
	}
	output.close();
}

void save_hei_populations_evolution(const string& output_path, const vector<double>& time_moments,
	const vector<dynamic_array>& hei_pop_data, const energy_diagram* hei_di)
{
	int i, l;
	string fname;
	ofstream output;

	fname = output_path + "h2exc_hei_populations.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(3);

	output << left << "! Evolution of population densities of HeI energy levels," << endl
		<< "! first row - nb of times, number of levels," << endl
		<< "! second row - level nb," << endl;

	output << left << "! " << setw(7) << (int)hei_pop_data.size() << setw(7) << hei_di->nb_lev << endl;
	output << left << setw(13) << "!";

	for (i = 1; i <= hei_di->nb_lev; i++) {
		output << left << setw(12) << i;
	}
	output << endl;

	for (l = 0; l < (int)hei_pop_data.size(); l++) {
		output.precision(4);
		output << left << setw(13) << time_moments[l];

		output.precision(3);
		for (i = 0; i < hei_pop_data[l].dim; i++) {
			output << left << setw(12) << hei_pop_data[l].arr[i];
		}
		output << endl;
	}
	output.close();
}


void save_electron_energy_losses(const string& output_path, const vector<double>& time_moments,
	const vector<double>& deriv_mt_arr, const vector<double>& mt_arr,
	const vector<double>& deriv_h2_rot_arr, const vector<double>& h2_rot_arr,
	const vector<double>& deriv_h2_vibr_arr, const vector<double>& h2_vibr_arr,
	const vector<double>& deriv_h2_electr_arr, const vector<double>& h2_electr_arr,
	const vector<double>& deriv_h2_electr_tr_arr, const vector<double>& h2_electr_tr_arr,
	const vector<double>& deriv_ioniz_arr, const vector<double>& ioniz_arr,
	const vector<double>& deriv_coloumb_ions_arr, const vector<double>& coloumb_ions_arr,
	const vector<double>& deriv_coloumb_el_arr, const vector<double>& coloumb_el_arr,
	const vector<double>& deriv_hei_arr, const vector<double>& hei_arr)
{
	int i;
	double tot, tot_int;
	string fname;
	ofstream output;

	fname = output_path + "h2exc_electron_enloss.txt";
	output.open(fname.c_str());

	output << scientific;
	output << left << "! Evolution of energy losses of electrons, parameter < 0 if electrons lose energy," << endl
		<< "! rate [eV cm-3 s-1] and integrated [eV cm-3]," << endl
		<< "! total rate and total integrated losses (without scattering on electrons), " << endl
		<< "! mt, mt_int - momentum transfer, " << endl
		<< "! rot, rot_int - H2 rotational excitation," << endl
		<< "! vibr, vibr_int - H2 vibrational excitation," << endl
		<< "! elecs, elecs_int - H2 electronic states excitation (singlet)," << endl
		<< "! elect, elect_int - H2 electronic states excitation (triplet)," << endl
		<< "! ion, ion_int - molecule/atom ionization," << endl
		<< "! coli, coli_int - Coulomb loses on ions," << endl
		<< "! cole, cole_int - Coulomb loses on electrons, |energy loss| + |energy gain|," << endl
		<< "! hei, hei_int - HeI excitation by electron impact" << endl;

	output << left << "! nb of times:" << endl
		<< "! " << setw(7) << (int)deriv_mt_arr.size() << endl;

	output << left << setw(13) << "!t(s)"
		<< setw(11) << "tot" << setw(13) << "tot_int"
		<< setw(11) << "mt" << setw(13) << "mt_int"
		<< setw(11) << "rot" << setw(13) << "rot_int"
		<< setw(11) << "vibr" << setw(13) << "vibr_int"
		<< setw(11) << "elecs" << setw(13) << "elecs_int"
		<< setw(11) << "elect" << setw(13) << "elect_int"
		<< setw(11) << "ion" << setw(13) << "ion_int"
		<< setw(11) << "coli" << setw(13) << "coli_int"
		<< setw(11) << "cole" << setw(13) << "cole_int"
		<< setw(11) << "hei" << setw(13) << "hei_int" << endl;

	for (i = 0; i < (int)deriv_mt_arr.size(); i++) {
		tot = deriv_mt_arr[i] + deriv_h2_rot_arr[i] + deriv_h2_vibr_arr[i] + deriv_h2_electr_arr[i] + deriv_h2_electr_tr_arr[i]
			+ deriv_ioniz_arr[i] + deriv_coloumb_ions_arr[i] + deriv_hei_arr[i];
		tot_int = mt_arr[i] + h2_rot_arr[i] + h2_vibr_arr[i] + h2_electr_arr[i] + h2_electr_tr_arr[i] + ioniz_arr[i]
			+ coloumb_ions_arr[i] + hei_arr[i];

		output.precision(3);
		output << left << setw(13) << time_moments[i];

		output.precision(2);
		output << left << setw(11) << tot << setw(13) << tot_int
			<< setw(11) << deriv_mt_arr[i] << setw(13) << mt_arr[i]
			<< setw(11) << deriv_h2_rot_arr[i] << setw(13) << h2_rot_arr[i]
			<< setw(11) << deriv_h2_vibr_arr[i] << setw(13) << h2_vibr_arr[i]
			<< setw(11) << deriv_h2_electr_arr[i] << setw(13) << h2_electr_arr[i]
			<< setw(11) << deriv_h2_electr_tr_arr[i] << setw(13) << h2_electr_tr_arr[i]
			<< setw(11) << deriv_ioniz_arr[i] << setw(13) << ioniz_arr[i]
			<< setw(11) << deriv_coloumb_ions_arr[i] << setw(13) << coloumb_ions_arr[i]
			<< setw(11) << deriv_coloumb_el_arr[i] << setw(13) << coloumb_el_arr[i]
			<< setw(11) << deriv_hei_arr[i] << setw(13) << hei_arr[i] << endl;
	}
	output.close();
}

void save_h2_data(const string& output_path, const vector<double>& time_moments,
	const vector<double>& h2_solomon_diss_arr, const vector<double>& h2_diss_exc_arr, const vector<double>& h2_diss_exc_triplet_arr,
	const vector<double>& hei_exc_arr)
{
	int i;
	double tot_int;
	string fname;
	ofstream output;

	fname = output_path + "h2exc_h2_data.txt";
	output.open(fname.c_str());

	output << scientific;
	output << left << "! Reaction contribution," << endl
		<< "! integrated losses [cm-3]," << endl
		<< "! sol_int - Solomon process of H2, " << endl
		<< "! de_int  - dissociative excitation of H2 (via singlet states)," << endl
		<< "! det_int - dissociative excitation of H2 (via triplet states)," << endl
		<< "! hei_int - HeI excitation" << endl;

	output << left << "! nb of times:" << endl
		<< "! " << setw(7) << (int)h2_solomon_diss_arr.size() << endl;

	output << left << setw(13) << "!t(s)" << setw(12) << "tot_int(H2)"
		<< setw(12) << "sol_int" << setw(12) << "de_int" << setw(12) << "det_int" << setw(12) << "hei_int" << endl;

	for (i = 0; i < (int)h2_solomon_diss_arr.size(); i++) {
		tot_int = h2_solomon_diss_arr[i] + h2_diss_exc_arr[i] + h2_diss_exc_triplet_arr[i];

		output.precision(3);
		output << left << setw(13) << time_moments[i];

		output.precision(2);
		output << left << setw(12) << tot_int
			<< setw(12) << h2_solomon_diss_arr[i]
			<< setw(12) << h2_diss_exc_arr[i]
			<< setw(12) << h2_diss_exc_triplet_arr[i]
			<< setw(12) << hei_exc_arr[i] << endl;
	}
	output.close();
}

void save_phys_parameters(const string& output_path, const vector<double>& time_moments,
	const vector<double>& neut_temp_arr, const vector<double>& ion_temp_arr, const vector<double>& dust_charge_arr,
	const vector<dynamic_array>& conc_data)
{
	int i;
	string fname;
	ofstream output;

	fname = output_path + "h2exc_physparam_data.txt";
	output.open(fname.c_str());

	output << scientific;
	output << left << "! Physical parameters," << endl
		<< "! temp_n - temperature of neutral component [K]," << endl
		<< "! temp_i - ion temperature [K]," << endl
		<< "! dcharge - charge of dust grains [cm-3]," << endl
		<< "! conc_e - concentration of electrons [cm-3]" << endl;

	output << left << "! nb of times:" << endl
		<< "! " << setw(7) << (int)neut_temp_arr.size() << endl;

	output << left << setw(13) << "!t(s)" << setw(12) << "temp_n(K)"
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
