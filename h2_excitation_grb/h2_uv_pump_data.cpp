// Checks:
// 23.01.2024 - memcpy instead of cycle.

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>

#include "constants.h"
#include "h2_uv_pump_data.h"

#define MAX_TEXT_LINE_WIDTH 240  // maximal size of the comment lines in the files,
#define SOURCE_NAME "h2_uv_pump_data.cpp"
using namespace std;


void h2_band_transitions(const std::string& path, std::vector<transition>& einstein_coeff_vector, 
	const energy_diagram* h2_excited_di, const energy_diagram* h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int	e, v, j, up, low;
	double coeff;

	string fname, name;
	ifstream input;
	stringstream ss;
	transition* trans;

	if (verbosity)
		cout << "H2 molecule radiative coefficients are being initializing..." << endl;

	if (h2_excited_di->electronic_state == 1) {
		name = "B";
		fname = path + "h2_cloudy_data/transprob_B.dat";
	}
	else if (h2_excited_di->electronic_state == 2) {
		name = "C_plus";
		fname = path + "h2_cloudy_data/transprob_C_plus.dat";
	}
	else if (h2_excited_di->electronic_state == 3) {
		name = "C_minus";
		fname = path + "h2_cloudy_data/transprob_C_minus.dat";
	}
	else if (h2_excited_di->electronic_state == 4) {
		name = "B_primed";
		fname = path + "h2_cloudy_data/transprob_B_primed.dat";
	}
	else if (h2_excited_di->electronic_state == 5) {
		name = "D_plus";
		fname = path + "h2_cloudy_data/transprob_D_plus.dat";
	}
	else if (h2_excited_di->electronic_state == 6) {
		name = "D_minus";
		fname = path + "h2_cloudy_data/transprob_D_minus.dat";
	}
	// M. Glass-Maujean, private communication (2025),
	else if (h2_excited_di->electronic_state == 7) {
		name = "B_primed_primed";
		fname = path + "h2_additional_data/transprob_glass_maujean_4ps.txt";
	}
	else if (h2_excited_di->electronic_state == 8) {
		name = "D_primed_plus";
		fname = path + "h2_additional_data/transprob_glass_maujean_4ppi_plus.txt";
	}
	else if (h2_excited_di->electronic_state == 9) {
		name = "D_primed_minus";
		fname = path + "h2_additional_data/transprob_glass_maujean_4ppi_minus.txt";
	}
	// transitions EF -> B
	// Glass-Maujean, Quadrelli & Dressler, Atomic Data and Nuclear Data Tables 30, 273 (1984) - vibrationally resolved Einstein coefficients,
	// Honl-London factors are used to calculated ro-vibrationally resolved Einstein coefficients,
	// Comment by M. Glass-Maujean: 
	// "To get the right non adiabatic transition probabilities, it is not so simple, the Honl-London factors are true in the adiabatic approximation only.
	// The full equations are in the paper showing the ratio A(R)/A(P) depends strongly of the mixing between parallel and perpendicular transition."
	else if (h2_excited_di->electronic_state == 10) {
		name = "EF-B";
		fname = path + "h2_additional_data/transprob_glass_maujean_ef_b.txt";
	}
	else {
		cout << "Error in " << SOURCE_NAME << ": the unknown excited electronic state of H2 molecule" << endl;
		exit(1);
	}
	
	input.open(fname.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	einstein_coeff_vector.clear();
	while (!input.eof()) {
		// comment lines are read (all lines must be commented, no empty lines at the file end)
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		if (text_line[0] == '\0')
			break;

		ss.clear();
		ss.str(text_line);

		ss >> e >> v >> j;  // upper level belongs to the excited electronic state,
		up = h2_excited_di->get_nb(v, j);

		ss >> e >> v >> j;  // lower level belongs to the ground (or excited) electronic state,
		low = h2_di->get_nb(v, j);

		ss >> coeff;  // Einstein coefficients in s-1,
		if (up != -1 and low != -1) {
			trans = new transition(h2_di->lev_array[low], h2_excited_di->lev_array[up], coeff);
			einstein_coeff_vector.push_back(*trans);
			delete trans;
		}
	}
	input.close();

	if (verbosity) {
		cout << "  data have been read from file " << fname << endl;
	}
}

void print_statistics(string output_path, string name, const std::vector<transition> & einstein_coeff_vector)
{
	const int nb_per_order = 100;
	int i, j, nb;
	
	double x, en_min, en_max, trans_en_low(12.), trans_en_high(12.), ground_trans_en_low(13.), ground_trans_en_high(13.);
	int *arr;

	string fname;
	ofstream output;

	x = pow(10., 1. / nb_per_order);
	en_min = pow(x, (int) (log(6.)/log(x)));
	en_max = pow(x, (int) (log(15.)/log(x)) + 1);

	nb = (int)(nb_per_order * log10(en_max / en_min)) + 1;
	arr = new int[nb];
	memset(arr, 0, nb * sizeof(int));

	for (i = 0; i < (int)einstein_coeff_vector.size(); i++) {
		j = (int) (log(einstein_coeff_vector[i].energy * CM_INVERSE_TO_EV / en_min) / log(x));
		if (j < 0)
			j = 0;
		if (j > nb - 1)
			j = nb - 1;
		arr[j]++;

		if (einstein_coeff_vector[i].energy * CM_INVERSE_TO_EV < trans_en_low)
			trans_en_low = einstein_coeff_vector[i].energy * CM_INVERSE_TO_EV;

		if (einstein_coeff_vector[i].energy * CM_INVERSE_TO_EV > trans_en_high)
			trans_en_high = einstein_coeff_vector[i].energy * CM_INVERSE_TO_EV;

		if (einstein_coeff_vector[i].low_lev.el == 0 && einstein_coeff_vector[i].low_lev.v == 0 
			&& (einstein_coeff_vector[i].low_lev.j == 0 || einstein_coeff_vector[i].low_lev.j == 1)) 
		{
			if (einstein_coeff_vector[i].energy * CM_INVERSE_TO_EV < ground_trans_en_low)
				ground_trans_en_low = einstein_coeff_vector[i].energy * CM_INVERSE_TO_EV;

			if (einstein_coeff_vector[i].energy * CM_INVERSE_TO_EV > ground_trans_en_high)
				ground_trans_en_high = einstein_coeff_vector[i].energy * CM_INVERSE_TO_EV;
		}
	}
	
	fname = output_path + "stat_energy_distr_" + name + ".txt";
	output.open(fname.c_str(), ios_base::out);
	
	output << scientific;
	output.precision(6);
	output << "! The nb of transitions in the interval [en, en + d_en]; " << endl
		<< "! lowest transition energy (eV): " << trans_en_low << endl
		<< "! highest transition energy (eV): " << trans_en_high << endl
		<< "! lowest transition energy from the ground level v=0, j=0,1 (eV): " << ground_trans_en_low << endl
		<< "! highest transition energy from the ground level v=0, j=0,1 (eV): " << ground_trans_en_high << endl
		<< "! i, en (eV), number of transitions" << endl;
	
	output.precision(3);
	for (i = 0; i < nb; i++) {
		output << left << setw(6) << i << setw(12) << en_min * pow(x, i) << setw(12) << arr[i] << endl;
	}
	output.close();
}


//
//
//
dissociation_data_unit::dissociation_data_unit()
	: v(0), j(0), prob(0.), kin_energy(0.)
{;}

dissociation_data_unit::dissociation_data_unit(const dissociation_data_unit& obj) {
	v = obj.v;
	j = obj.j;
	prob = obj.prob;
	kin_energy = obj.kin_energy;
}


h2_energy_level_param::h2_energy_level_param(int i, int nbd, double dp, double td, double ke)
	: nb(i), diss_prob(dp), tot_decay(td), kin_energy(ke), nb_of_decays(nbd)
{
	decay_level_nbs = new int[nb_of_decays];
	memset(decay_level_nbs, 0, nb_of_decays*sizeof(int));

	decay_probs = new double[nb_of_decays];
	memset(decay_probs, 0, nb_of_decays * sizeof(double));
}

h2_energy_level_param::~h2_energy_level_param()
{
	delete[] decay_level_nbs;
	delete[] decay_probs;
}

h2_energy_level_param::h2_energy_level_param(const h2_energy_level_param &obj)
{
	nb = obj.nb;
	nb_of_decays = obj.nb_of_decays;
	diss_prob = obj.diss_prob;
	tot_decay = obj.tot_decay;
	kin_energy = obj.kin_energy;

	decay_level_nbs = new int[nb_of_decays];
	decay_probs = new double[nb_of_decays];

	memcpy(decay_level_nbs, obj.decay_level_nbs, nb_of_decays *sizeof(decay_level_nbs[0]));
	memcpy(decay_probs, obj.decay_probs, nb_of_decays * sizeof(decay_probs[0]));
}

h2_energy_level_param& h2_energy_level_param::operator=(const h2_energy_level_param & obj)
{
	if (this == &obj)
		return *this;

	nb = obj.nb;
	nb_of_decays = obj.nb_of_decays;
	diss_prob = obj.diss_prob;
	tot_decay = obj.tot_decay;
	kin_energy = obj.kin_energy;

	delete[] decay_level_nbs;
	delete[] decay_probs;

	decay_level_nbs = new int[nb_of_decays];
	decay_probs = new double[nb_of_decays];

	memcpy(decay_level_nbs, obj.decay_level_nbs, nb_of_decays * sizeof(decay_level_nbs[0]));
	memcpy(decay_probs, obj.decay_probs, nb_of_decays * sizeof(decay_probs[0]));

	return *this;
}


void h2_excited_state_data(const std::string& data_path, std::vector<h2_energy_level_param> & level_data, 
	const energy_diagram* h2_excited_di, const std::vector<transition> & h2_band, int diss_data_exists, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, k, v, j;
	double prob, tot_decay_rate, kin_energy;

	string fname;
	stringstream ss;
	ifstream input;

	h2_energy_level_param* parameters;
	decay_channel dec_ch;
	vector<decay_channel> decay_list;

	dissociation_data_unit diss_fdata;
	vector<dissociation_data_unit> diss_fdata_arr;

	level_data.clear();
	diss_fdata_arr.clear();

	if (verbosity) {
		cout << "H2 molecule dissociation probabilities are being initializing..." << endl;
	}

	if (diss_data_exists)
	{
		if (h2_excited_di->electronic_state == 1) { // B state,
			fname = "h2_cloudy_data/dissprob_B.dat";
		}
		else if (h2_excited_di->electronic_state == 2) { // C+ state
			fname = "h2_cloudy_data/dissprob_C_plus.dat";
		}
		else if (h2_excited_di->electronic_state == 3) { // C- state
			fname = "h2_cloudy_data/dissprob_C_minus.dat";
		}
		else if (h2_excited_di->electronic_state == 4) { // Bp state,
			fname = "h2_cloudy_data/dissprob_B_primed.dat";
		}
		else if (h2_excited_di->electronic_state == 5) { // D+ state
			fname = "h2_cloudy_data/dissprob_D_plus.dat";
		}
		else if (h2_excited_di->electronic_state == 6) { // D- state
			fname = "h2_cloudy_data/dissprob_D_minus.dat";
		}
		else {
			cout << "Error in " << SOURCE_NAME << ": there is no data for the excited electronic state..." << endl;
			exit(1);
		}

		fname = data_path + fname;
		input.open(fname.c_str(), ios_base::in);

		if (!input.is_open()) {
			cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
			exit(1);
		}
		else {
			if (verbosity) {
				cout << "	reading data from file: " << fname << endl;
			}
		}

		while (!input.eof()) {
			// comment lines are read (all lines must be commented, no empty lines at the file end)
			do
				input.getline(text_line, MAX_TEXT_LINE_WIDTH);
			while (text_line[0] == '#');

			if (text_line[0] == '\0')
				break;

			ss.clear();
			ss.str(text_line);

			ss >> v >> j >> prob >> kin_energy;  // dissociation probability here in s-1, energy in eV

			diss_fdata.v = v;
			diss_fdata.j = j;
			diss_fdata.prob = prob;
			diss_fdata.kin_energy = kin_energy;  // in eV;

			diss_fdata_arr.push_back(diss_fdata);
		}
		input.close();

		if (verbosity) {
			cout << "  data have been read from file " << fname << endl;
		}
	}

	for (i = 0; i < h2_excited_di->nb_lev; i++) {
		prob = tot_decay_rate = kin_energy = 0.;
		
		for (k = 0; k < (int)diss_fdata_arr.size(); k++) {
			if (diss_fdata_arr[k].v == h2_excited_di->lev_array[i].v && diss_fdata_arr[k].j == h2_excited_di->lev_array[i].j) 
			{
				prob = diss_fdata_arr[k].prob;
				tot_decay_rate += prob;
				kin_energy = diss_fdata_arr[k].kin_energy;
				break;
			}
		}
		
		decay_list.clear();
		for (k = 0; k < (int)(h2_band.size()); k++) {
			if (h2_band[k].up_lev.j == h2_excited_di->lev_array[i].j && h2_band[k].up_lev.v == h2_excited_di->lev_array[i].v) 
			{
				tot_decay_rate += h2_band[k].einst_coeff;
				dec_ch.nb = h2_band[k].low_lev.nb;
				dec_ch.rate = h2_band[k].einst_coeff;
				decay_list.push_back(dec_ch);
			}
		}
		prob /= tot_decay_rate;  // dissociation probability, making dimensionless

		parameters = new h2_energy_level_param(i, (int)decay_list.size(), prob, tot_decay_rate, kin_energy);

		for (k = 0; k < parameters->nb_of_decays; k++) {
			parameters->decay_level_nbs[k] = decay_list[k].nb;
			parameters->decay_probs[k] = decay_list[k].rate / tot_decay_rate;
		}

		level_data.push_back(*parameters);
		delete parameters;
	}
	sort(level_data.begin(), level_data.end());
}

// A -> B -> X
// level_data_low must already exist,
void h2_excited_state_data(const string& path, vector<h2_energy_level_param>& level_data_up, vector<h2_energy_level_param> level_data_low,
	const energy_diagram* h2_excited_di_up, const energy_diagram* h2_excited_di_low, const vector<transition>& h2_band_up, const vector<transition>& h2_band_low, 
	int verbosity)
{
	int i, k, l, n, m;
	double x, prob, tot_decay_rate, kin_energy;

	h2_energy_level_param* parameters;
	decay_channel dec_ch;
	vector<decay_channel> decay_list;

	for (i = 0; i < h2_excited_di_up->nb_lev; i++) {
		decay_list.clear();
		prob = tot_decay_rate = kin_energy = 0.;

		// Calculation of the total decay rate of the levels of the upper electronic state A
		for (k = 0; k < (int)(h2_band_up.size()); k++) {
			if (h2_band_up[k].up_lev.j == h2_excited_di_up->lev_array[i].j && h2_band_up[k].up_lev.v == h2_excited_di_up->lev_array[i].v) 
			{
				tot_decay_rate += h2_band_up[k].einst_coeff;  // s-1
			}
		}

		for (k = 0; k < (int)(h2_band_up.size()); k++) {
			if (h2_band_up[k].up_lev.j == h2_excited_di_up->lev_array[i].j && h2_band_up[k].up_lev.v == h2_excited_di_up->lev_array[i].v)
			{
				// finding the level of the electronic state B, and using the data on B -> X transitions
				l = h2_excited_di_low->get_nb(h2_band_up[k].low_lev.v, h2_band_up[k].low_lev.j);
				
				for (n = 0; n < level_data_low[l].nb_of_decays; n++) {
					for (m = 0; m < (int)decay_list.size(); m++) {
						if (decay_list[m].nb == level_data_low[l].decay_level_nbs[n]) {
							break;
						}
					}

					if (m < (int)decay_list.size()) {
						decay_list[m].rate += level_data_low[l].decay_probs[n] * h2_band_up[k].einst_coeff;	
					}
					else 
					{
						dec_ch.nb = level_data_low[l].decay_level_nbs[n];
						dec_ch.rate = level_data_low[l].decay_probs[n] * h2_band_up[k].einst_coeff;
						decay_list.push_back(dec_ch);
					}	
				}
				// dissociation probability in level data is dimensionless,
				prob += level_data_low[l].diss_prob * h2_band_up[k].einst_coeff;
				kin_energy += level_data_low[l].kin_energy * h2_band_up[k].einst_coeff;
			}
		}
		prob /= tot_decay_rate;
		kin_energy /= tot_decay_rate;

		parameters = new h2_energy_level_param(i, (int)decay_list.size(), prob, tot_decay_rate, kin_energy);

		x = 0.;
		for (k = 0; k < parameters->nb_of_decays; k++) {
			parameters->decay_level_nbs[k] = decay_list[k].nb;
			parameters->decay_probs[k] = decay_list[k].rate / tot_decay_rate;
			x += parameters->decay_probs[k];
		}

		x += parameters->diss_prob;
		if (x < 1. - 10. * numeric_limits<double>::epsilon()) {
			if (verbosity) {
				cout << "v = " << h2_excited_di_up->lev_array[i].v << "  j = " << h2_excited_di_up->lev_array[i].j << "  total probability " << x << endl;
			}
		}

		level_data_up.push_back(*parameters);
		delete parameters;
	}
	sort(level_data_up.begin(), level_data_up.end());
}



void print_statistics(std::string output_path, std::string name, const energy_diagram* h2_excited_di, 
	const std::vector<h2_energy_level_param>& level_param_vector)
{
	int i;
	string fname;
	ofstream output;

	fname = output_path + "stat_lifetime_" + name + ".txt";
	output.open(fname.c_str(), ios_base::out);

	output << scientific;
	output.precision(6);
	output << "! The lifetime of the levels, s^{-1} " << endl
		<< "! i, v, j, t(s-1)" << endl;

	output.precision(3);
	for (i = 0; i < (int)level_param_vector.size(); i++) {
		output << left << setw(6) << i << setw(6) << h2_excited_di->lev_array[i].v << setw(6) << rounding(h2_excited_di->lev_array[i].j) 
			<< setw(12) << level_param_vector[i].tot_decay << endl;
	}
	output.close();
}


//
// Honl-London factors, normalized, 
// (j, j - 1) + (j, j + 1) = 1
double honl_london_singlet_dl0::operator()(int j, int dj)
{
	double f(0.);
	if (dj == -1) {// P(J) transition, j -> j - 1 (upward)
		f = ((double)j) / (2. * j + 1.);
	}
	else if (dj == 1) {  // R(J) transition, j -> j + 1 (upward)
		f = (j + 1.) / (2. * j + 1.);
	}
	return f;
}

// (j, j - 1) + (j, j) + (j, j + 1) = 1
double honl_london_singlet_dl1::operator()(int j, int dj)
{
	double f(0.);
	if (dj == -1)  // P(j) transition, j -> j - 1 (upward)
		f = 0.5 * (j - 1.);
	else if (dj == 0)  // Q(j) transition, j -> j
		f = 0.5 * (2. * j + 1.);
	else  // R(j) transition, j -> j + 1
		f = 0.5 * (j + 2.);
	
	f /= (2. * j + 1.);
	return f;
}
