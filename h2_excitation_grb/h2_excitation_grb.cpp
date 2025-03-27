// h2_excitation_grb.cpp 
// Important notes:
// The Gauss unit system is adopted elsewhere, exceptions: arguments of the cross sections, 
// concentration [cm-3]


#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#define _CRT_SECURE_NO_WARNINGS

// SUNDIALS CVODE solver headers
// works only with static libraries (shared must not be installed)
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#include <omp.h>
#include <stdio.h>
#include <cmath>
#include <limits>
#include <vector>
#include <cstring>
#include <ctime>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "constants.h"
#include "h2_parameters.h"
#include "dynamic_array.h"
#include "h2_elspectrum_eqs.h"
#include "h2_saving_data.h"

#define MAX_TEXT_LINE_WIDTH 3350  // maximal size of the comment lines in the files in symbols,
#define SOURCE_NAME "h2_excitation_grb.cpp"

using namespace std;


//
// Simulations of the excitation of H2 molecules in the gas (plasma)
void simulate_h2_excitation(const string& data_path, const string& sim_path, const string& output_path, double hcolumn_density, 
	int is_h2_pop_dens_init, double op_ratio_h2,
	int test_fast_electron_yields, double test_electron_conc, double test_electron_energy, double test_he_abund, double test_conc_h_tot, double test_ioniz_fract);

void init_electron_energy_grid(vector<double>& electron_energies_grid, vector<double>& electron_energy_bin_size);

// sim_path - path to the folder with simulation data (GRB propagation results)
// output_path - path to the folder where data of the current simulations must be saved,
void init_specimen_conc(const string& sim_path, double& conc_h_tot, vector<double>& layer_centre_distances, vector<dynamic_array>& d);
void init_dust_abund(const string& sim_path, vector<double>& d, double& grain_radius_init, double& grain_nb_density, double& grain_material_density);
void init_electron_spectra(const string& sim_path, const string& output_path, vector<dynamic_array>& d, const vector<double>& electron_energies_grid);
void init_electron_spectra_test(const string& output_path, dynamic_array & init_data_test, const vector<double>& en_grid, 
	double test_electron_conc, double test_electron_energy);
void init_h2_population_density(const string& sim_path, vector<dynamic_array>& d, vector<int>& qnb_v_arr, vector<int>& qnb_j_arr);


int f_elsp(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
	return static_cast<elspectra_evolution_data*>(user_data)->f(t, y, ydot);
}

int main()
{
	int nb_processors(4);
	int test_fast_electron_yields, is_h2_pop_dens_init;
	double hcolumn_density, op_ratio_h2, test_electron_conc, test_electron_energy, test_he_abund, test_conc_h_tot, test_ioniz_fract;

	char text_line[MAX_TEXT_LINE_WIDTH];
	
	string data_path;    // path to the input_data folder with spectroscopic, collisional data and etc.
	string output_path;  // path to the folder where data of the current simulation must be saved,
	string sim_path;     // path to the folder with simulation data (GRB propagation results)

	stringstream ss;
	ifstream input;

	// Parameters of the H2 excitation in molecular cloud after GRB propagation,
	input.open("input_parameters.txt");
	while (!input.eof())
	{
		do // comment lines are read:
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> nb_processors;

#ifdef _OPENMP
		omp_set_num_threads(nb_processors);

#pragma omp parallel 
		{
#pragma omp master 
			{
				cout << "OpenMP is supported" << endl;
				cout << "Nb of threads: " << omp_get_num_threads() << endl;
			}
		}
#endif
		// path to the directory with data tables(spectroscopic data, collision rates and etc.):
		// linux - "/disk4/nester/input_data/"
		// windows - "C:/Users/Александр/Documents/input_data/"
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> data_path;

		// path to the directory with the simulation data (GRB propagation)
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		// path may contain spaces
		sim_path = text_line;

		// path to the directory for the output
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		output_path = text_line;

		// H column density from the cloud boundary to the cloud layer centre,
		// the layer width is the same for all layers,
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> hcolumn_density;

		// H2 populations are read from file, 0 - no, 1 - yes,
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> is_h2_pop_dens_init;

		// ortho-to-para-H2 ratio, this parameter is used if H2 population densities are not read from file, 
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> op_ratio_h2;

		//
		// test simulations - electrons have fixed energy
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> test_fast_electron_yields;

		// electron concentration in test simulations in [cm-3]
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> test_electron_conc;

		// Electron energy in test simulations in[eV]
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> test_electron_energy;

		// Helium abundance in test simulations(only H2 and He are taken into account),
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> test_he_abund;

		// H nuclei total concentration for test simulations, [cm-3],
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> test_conc_h_tot;

		// Ionization fraction x_e, number density of electrons n_e = n_{H,tot} * x_e,
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> test_ioniz_fract;
		
		input.close();

		// check test mode in the file with parameters, h2_parameters.h
		simulate_h2_excitation(data_path, sim_path, output_path, hcolumn_density, is_h2_pop_dens_init, op_ratio_h2, 
			test_fast_electron_yields, test_electron_conc, test_electron_energy, test_he_abund, test_conc_h_tot, test_ioniz_fract);
	}

	output_path = "C:/Users/Александр/Documents/Данные и графики/paper GRB in molecular cloud/python_scripts_ism/";
	//save_cross_section_table(output_path, data_path);
	//calc_helium_lifetimes(output_path, data_path);
}


// Initialization of the grid for electron energies,
// there IS NO interval for electrons with high energy (E > max_el_energy)
// there IS interval for electrons with low energy [0, min energy], min_energy must be < any excitation process threshold
void init_electron_energy_grid(vector<double>& electron_energies_grid, vector<double>& electron_energy_bin_size)
{
	int i, k, nb_bins;
	double x, y;

	electron_energies_grid.clear();
	electron_energy_bin_size.clear();

	x = pow(10., 1. / NB_OF_BINS_PER_ORDER_EL);
	nb_bins = (int)(1. / (x - 1.)) + 1;

	// equal energy intervals in energy range [0, E_fixed],
	y = ELECTRON_ENERGY_FIXED / nb_bins;
	for (i = 0; i <= nb_bins; i++) {
		electron_energies_grid.push_back(i * y);
	}

	k = (int)(NB_OF_BINS_PER_ORDER_EL * log10((double)MAX_ELECTRON_ENERGY / ELECTRON_ENERGY_FIXED));
	for (i = 0; i < k; i++) {
		electron_energies_grid.push_back(electron_energies_grid.back() * x);
	}

	for (i = 0; i < (int)electron_energies_grid.size() - 1; i++) {
		y = electron_energies_grid[i + 1] - electron_energies_grid[i];
		electron_energy_bin_size.push_back(y);
	}
}


void simulate_h2_excitation(const string& data_path, const string& sim_path, const string& output_path, double hcolumn_density, 
	int is_h2_pop_dens_init, double op_ratio_h2,
	int test_fast_electron_yields, double test_electron_conc, double test_electron_energy, double test_he_abund, double test_conc_h_tot, double test_ioniz_fract)
{
#ifdef __linux__
	stringstream lin_out;
	lin_out << output_path;
	lin_out << "out";
	lin_out << "_screen";
	lin_out << ".txt";
	// lin_out.str("/dev/null");

	//ofstream outerr(lin_out.str().c_str(), ios::app);
	//streambuf* orig_cerr = cerr.rdbuf(); // original cerr;
	//cerr.rdbuf(outerr.rdbuf());

	ofstream out(lin_out.str().c_str(), ios::app);
	streambuf* orig_cout = cout.rdbuf();
	cout.rdbuf(out.rdbuf());
#endif
	const int verbosity = 1;
	
	bool dust_is_presented(true), must_be_saved;
	int i, k, nb_of_time_moments, lay_nb, time_nb, step_nb, flag, nb_of_equat, nb_of_el_energies, nb_lev_h2, nb_lev_hei, 
		h2eq_nb, heieq_nb, physeq_nb;
	long int nb_steps_tot;
	double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, rel_tol, dt, model_time, model_time_aux, model_time_aux_prev, model_time_step, model_time_out,
		conc_h_tot, ioniz_fract, grain_radius_init, grain_nb_density, grain_material_density, grain_radius, dust_mass, gas_mass, grb_distance, 
		dg_ratio, grb_cloud_distance;

	double enloss_rate_mt, enloss_rate_h2_rot, enloss_rate_h2_rot_pos, enloss_rate_h2_vibr, enloss_rate_h2_vibr_pos,
		enloss_rate_h2_singlet, enloss_rate_h2_triplet, enloss_rate_ioniz, enloss_rate_hei, enloss_rate_coulomb_el, diss_decay_heating_rate,
		neutral_coll_heating_rate, enloss_mt, enloss_h2_rot, enloss_h2_vibr, enloss_h2_singlet, enloss_h2_triplet, enloss_ioniz, enloss_hei,
		enloss_coulomb_el, diss_decay_heating, neutral_coll_heating;
		
	double h2_solomon_diss_rate, h2_diss_exc_singlet_rate, h2_diss_exc_triplet_rate, hei_exc_rate, h2_excit_electr_bs_rate, 
		h2_excit_electr_cp_rate, h2_excit_electr_rate, h2_excit_vibr_1_rate, h2_excit_vibr_2_rate, h2_excit_vibr_rate,  h2_excit_rot_rate,
		h2_solomon_diss, h2_diss_exc_singlet, h2_diss_exc_triplet, hei_exc, h2_excit_electr_bs, h2_excit_electr_cp, h2_excit_electr, 
		h2_excit_vibr_1, h2_excit_vibr_2, h2_excit_vibr;

	time_t timer;
	char* ctime_str;

	string fname;
	ofstream output;

	vector<int>    qnb_v_arr, qnb_j_arr;
	vector<double> time_cloud_arr;      // in s
	vector<double> layer_centre_distances;  // in cm
	vector<double> electron_energies_grid, electron_energy_bin_size;  // in eV,

	// Energy loss rates in [eV cm-3 s-1], as a function of time
	vector<double> enloss_rate_mt_arr, enloss_rate_h2_rot_arr, enloss_rate_h2_rot_arr_pos, enloss_rate_h2_vibr_arr, enloss_rate_h2_vibr_arr_pos, 
		enloss_rate_h2_singlet_arr, enloss_rate_h2_triplet_arr, enloss_rate_ioniz_arr, enloss_rate_hei_arr, 
		enloss_rate_coulomb_el_arr, diss_decay_heating_rate_arr, neutral_coll_heating_rate_arr;

	// Excitation rates in [cm-3 s-1], as a function of time
	vector<double> h2_excit_electr_rate_arr, h2_excit_vibr_rate_arr, h2_excit_rot_rate_arr;

	// Energy loss in [eV cm-3] during the time interval [0, t] as a function of time t,
	vector<double> enloss_mt_arr, enloss_h2_rot_arr, enloss_h2_vibr_arr, enloss_h2_singlet_arr, enloss_h2_triplet_arr, enloss_ioniz_arr, 
		enloss_hei_arr, enloss_coulomb_el_arr, diss_decay_heating_arr, neutral_coll_heating_arr;

	// Number of excitations in [cm-3] up to a given time, as a function of time
	vector<double> h2_excit_electr_arr, h2_excit_electr_bs_arr, h2_excit_electr_cp_arr, 
		h2_excit_vibr_arr, h2_excit_vibr_1_arr, h2_excit_vibr_2_arr;

	// Dissociated concentration of H2 in [cm-3] up to a given time, as a function of time
	vector<double> h2_solomon_diss_arr, h2_diss_exc_singlet_arr, h2_diss_exc_triplet_arr, hei_exc_arr;

	// Physical parameters as a function of time,
	vector<double> neutral_temp_arr, ion_temp_arr, dust_charge_arr;

	// Parameter[layer nb], at the fixed time (start of electron spectra simulations),
	vector<double> dg_ratio_cloud_data, ioniz_cloud_data;

	// Spectrum[layer nb][energy bin] in [cm-3 eV-1], vector index - layer nb, array index - energy bin nb (specimen nb, H2 level nb),
 	vector<dynamic_array> el_spectrum_init_data, specimen_conc_init_data; 
	vector<dynamic_array> h2_pop_dens_init_data;
	
	// Electron spectrum[time nb][energy bin] in [cm-3], vector index - time moment nb, array index - energy bin nb 
	vector<dynamic_array> el_spectrum_evol;
	// Concentration[time nb][specimen nb], or [level nb]
	vector<dynamic_array> specimen_conc_evol, h2_popdens_evol, h2_popdens_v_evol, hei_popul_dens_evol;  // in [cm-3]

	cout << scientific;
	cout.precision(3);

	timer = time(NULL);
	ctime_str = ctime(&timer);

	cout << ctime_str << endl
		<< "Start of the simulations of the evolution of electron spectra" << endl
		<< "Initialization of the data..." << endl;

	// initialization of the electron spectrum grid, that can be finer than the energy grid in the file with simulation data,
	init_electron_energy_grid(electron_energies_grid, electron_energy_bin_size);	
	nb_of_el_energies = (int)electron_energies_grid.size() - 1;
	
	dynamic_array el_spectrum_test(nb_of_el_energies);
	if (test_fast_electron_yields) 
	{
		init_electron_spectra_test(output_path, el_spectrum_test, electron_energies_grid, test_electron_conc, test_electron_energy);
		
		ioniz_fract = test_ioniz_fract;
		conc_h_tot = test_conc_h_tot;
		grb_cloud_distance = 0.;
		hcolumn_density = 0.;
		grb_distance = 0.;
		lay_nb = 0;
	}
	else
	{
		ioniz_fract = IONIZATION_FRACTION_THERMAL_EL;
		init_dust_abund(sim_path, dg_ratio_cloud_data, grain_radius_init, grain_nb_density, grain_material_density);

		// gas density and layer coordinates are initialized here,
		init_specimen_conc(sim_path, conc_h_tot, layer_centre_distances, specimen_conc_init_data);

		// electron spectrum in the file has the dimension [cm-3 eV-1],	
		init_electron_spectra(sim_path, output_path, el_spectrum_init_data, electron_energies_grid);

		// population density in the file has dimension [cm-3],
		init_h2_population_density(sim_path, h2_pop_dens_init_data, qnb_v_arr, qnb_j_arr);

		// distance from the GRB source to the cloud boundary (the layer width is the same for all layers),
		grb_cloud_distance = layer_centre_distances[0] - 0.5 * (layer_centre_distances[1] - layer_centre_distances[0]);

		// Simulations are carried out for fixed cloud layer,
		// Calculation of layer nb:
		for (lay_nb = 0; lay_nb < (int)layer_centre_distances.size(); lay_nb++) {
			if (hcolumn_density <= conc_h_tot * (layer_centre_distances[lay_nb] - grb_cloud_distance))
				break;
		}

		// updating column density for given layer number,
		hcolumn_density = conc_h_tot * (layer_centre_distances[lay_nb] - grb_cloud_distance);

		// distance from the GRB source to the centre of the particular cloud layer,
		grb_distance = layer_centre_distances[lay_nb];
	}
	elspectra_evolution_data user_data(data_path, output_path, conc_h_tot, ioniz_fract, electron_energies_grid, verbosity);

	nb_of_equat = user_data.get_nb_of_equat();
	nb_lev_h2 = user_data.get_nb_of_h2_lev();
	nb_lev_hei = user_data.get_nb_of_hei_lev();

	h2eq_nb = user_data.get_h2eq_nb();
	heieq_nb = user_data.get_heieq_nb();
	physeq_nb = user_data.get_physeq_nb();

	dynamic_array el_spectrum(nb_of_el_energies), specimen_conc(NB_OF_CHEM_SPECIES), h2_popul_dens(nb_lev_h2), hei_popul_dens(nb_lev_hei), 
		h2_pop_v(MAX_H2_VSTATES_X1SU);

	time_cloud_arr.clear();
	time_cloud_arr.push_back(0.);
	time_cloud_arr.push_back(MIN_MODEL_TIME);

	x1 = pow(10., 1. / NB_OF_BINS_PER_ORDER_TIME);
	
	for (i = 0; time_cloud_arr.back() < MAX_MODEL_TIME; i++) {
		time_cloud_arr.push_back(time_cloud_arr.back() * x1);
	}
	nb_of_time_moments = (int)time_cloud_arr.size();

	// initialization of vectors and matrices used by solver SUNDIALS CVODE
	N_Vector y, ydot, abs_tol;
	SUNMatrix A(NULL);
	SUNLinearSolver LS(NULL);

	y = N_VNew_Serial(nb_of_equat);
	ydot = N_VNew_Serial(nb_of_equat);
	abs_tol = N_VNew_Serial(nb_of_equat);

	model_time_aux = model_time = 0.;

	// Initialization for tolerances:
	rel_tol = REL_ERROR_SOLVER;
	user_data.set_tolerances(abs_tol);

	// Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula and the use of a Newton iteration 
	void* cvode_mem = CVodeCreate(CV_BDF);

	// Call CVodeInit to initialize the integrator memory and specify the user's right hand side function in y'=f(t,y), 
	// the initial time t0, and the initial dependent variable vector y;
	// y is not defined yet,
	flag = CVodeInit(cvode_mem, f_elsp, model_time, y);

	// Call CVodeSVtolerances to specify the scalar tolerances:
	flag = CVodeSVtolerances(cvode_mem, rel_tol, abs_tol);

	// The maximal number of steps between simulation stops;
	flag = CVodeSetMaxNumSteps(cvode_mem, 10000);

	// specifies the maximum number of error test failures permitted in attempting one step:
	flag = CVodeSetMaxErrTestFails(cvode_mem, MAX_ERR_TEST_FAILS_SOLVER); // default value is 7; 
	flag = CVodeSetMaxConvFails(cvode_mem, MAX_CONV_FAILS_SOLVER); // default value is 10;

	// Create dense SUNMatrix for use in linear solves 
	A = SUNDenseMatrix(nb_of_equat, nb_of_equat);

	// Create dense SUNLinearSolver object for use by CVode 
	LS = SUNDenseLinearSolver(y, A);

	// Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode 
	flag = CVDlsSetLinearSolver(cvode_mem, LS, A);

	// The function attaches the user data block to the solver;
	flag = CVodeSetUserData(cvode_mem, &user_data);

	// Initialization of initial values of variables,
	for (i = 0; i < nb_of_equat; i++) {
		NV_Ith_S(y, i) = 0.;
	}

	if (test_fast_electron_yields) {
		for (i = 0; i < nb_of_el_energies; i++) {
			// the number of electrons in the interval:
			NV_Ith_S(y, i) = el_spectrum_test.arr[i] * user_data.get_electron_energy_bin(i);
		}

		NV_Ith_S(y, nb_of_el_energies + H2_NB) = 0.5 * conc_h_tot;
		NV_Ith_S(y, nb_of_el_energies + HE_NB) = test_he_abund * conc_h_tot;

		// the initial distribution of energy level populations is postulated, [cm-3]
		NV_Ith_S(y, h2eq_nb) = NV_Ith_S(y, nb_of_el_energies + H2_NB) / (op_ratio_h2 + 1.);
		NV_Ith_S(y, h2eq_nb + 1) = NV_Ith_S(y, nb_of_el_energies + H2_NB) * op_ratio_h2 / (op_ratio_h2 + 1.);

		// the initial distribution of HeI levels,
		NV_Ith_S(y, heieq_nb) = NV_Ith_S(y, nb_of_el_energies + HE_NB);

		// Initial gas/dust temperature, in K
		NV_Ith_S(y, physeq_nb) = NV_Ith_S(y, physeq_nb + 1) = 10.;

		// Dust charge
		NV_Ith_S(y, physeq_nb + 2) = 0.;

		dust_is_presented = false;
		grain_radius = 0.;
		grain_nb_density = 0.;
		dg_ratio = 0.;

		user_data.set_dust_parameters(dust_is_presented, grain_radius, grain_nb_density);
	}
	else {
		for (i = 0; i < nb_of_el_energies; i++) {
			// the number of electrons in the interval:
			NV_Ith_S(y, i) = el_spectrum_init_data[lay_nb].arr[i] * user_data.get_electron_energy_bin(i);
		}

		for (i = 0; i < NB_OF_CHEM_SPECIES; i++) {
			NV_Ith_S(y, nb_of_el_energies + i) = specimen_conc_init_data[lay_nb].arr[i];
		}

		if (is_h2_pop_dens_init) {
			// population densities must be in [cm-3],
			x1 = x2 = 0.;
			for (i = 0; i < h2_pop_dens_init_data[lay_nb].dim; i++) 
			{
				k = user_data.get_level_nb_h2(qnb_v_arr[i], qnb_j_arr[i]);
				if (k >= 0) {
					NV_Ith_S(y, h2eq_nb + k) = h2_pop_dens_init_data[lay_nb].arr[i];
				}

				if (qnb_j_arr[i] % 2 == 0) {
					x2 += h2_pop_dens_init_data[lay_nb].arr[i];
				}
				else {
					x1 += h2_pop_dens_init_data[lay_nb].arr[i];
				}
			}
			// calculating of ortho-para-H2 ratio,
			op_ratio_h2 = x1 / x2;
		}
		else {
			// the initial distribution of energy level populations is postulated, [cm-3]
			NV_Ith_S(y, h2eq_nb) = NV_Ith_S(y, nb_of_el_energies + H2_NB) / (op_ratio_h2 + 1.);
			NV_Ith_S(y, h2eq_nb + 1) = NV_Ith_S(y, nb_of_el_energies + H2_NB) * op_ratio_h2 / (op_ratio_h2 + 1.);
		}

		// the initial distribution of HeI levels,
		NV_Ith_S(y, heieq_nb) = NV_Ith_S(y, nb_of_el_energies + HE_NB);

		// Initial gas/dust temperature, in K
		NV_Ith_S(y, physeq_nb) = NV_Ith_S(y, physeq_nb + 1) = 10.;

		// Dust charge
		NV_Ith_S(y, physeq_nb + 2) = 0.;

		gas_mass = specimen_conc_init_data[lay_nb].arr[H_NB] + specimen_conc_init_data[lay_nb].arr[H_P_NB]
			+ 2. * (specimen_conc_init_data[lay_nb].arr[H2_NB] + specimen_conc_init_data[lay_nb].arr[H2_P_NB])
			+ 4. * (specimen_conc_init_data[lay_nb].arr[HE_NB] + specimen_conc_init_data[lay_nb].arr[HE_P_NB] + specimen_conc_init_data[lay_nb].arr[HE_PP_NB]);
		gas_mass *= ATOMIC_MASS_UNIT;

		// dust-gas mass ratio < 1.
		dust_mass = gas_mass * dg_ratio_cloud_data[lay_nb] / (1. - dg_ratio_cloud_data[lay_nb]);

		// calculation of grain cross section,
		if (dg_ratio_cloud_data[lay_nb] <= numeric_limits<double>::epsilon()) {
			dust_is_presented = false;
			grain_radius = 0.;
			dg_ratio = 0.;
		}
		else {
			dust_is_presented = true;
			x1 = dust_mass * 3. / (grain_nb_density * grain_material_density * 4. * M_PI);
			grain_radius = pow(x1, 1. / 3.);   // grain radius, cm
			dg_ratio = dg_ratio_cloud_data[lay_nb];
		}
		user_data.set_dust_parameters(dust_is_presented, grain_radius, grain_nb_density);
	}
	
	// restart of the solver with new values of initial conditions,
	flag = CVodeReInit(cvode_mem, model_time, y);

	// update the members of user_data class,
	f_elsp(model_time, y, ydot, &user_data);

	enloss_rate_mt = enloss_rate_h2_rot = enloss_rate_h2_rot_pos = enloss_rate_h2_vibr = enloss_rate_h2_vibr_pos = enloss_rate_h2_singlet 
		= enloss_rate_h2_triplet = enloss_rate_ioniz = enloss_rate_hei = enloss_rate_coulomb_el = diss_decay_heating_rate = neutral_coll_heating_rate = 0.;

	enloss_mt = enloss_h2_rot = enloss_h2_vibr = enloss_h2_singlet = enloss_h2_triplet = enloss_ioniz = enloss_hei 
		= enloss_coulomb_el = diss_decay_heating = neutral_coll_heating  = 0.;

	h2_solomon_diss_rate = h2_diss_exc_singlet_rate = h2_diss_exc_triplet_rate = hei_exc_rate = 0.;
	h2_solomon_diss = h2_diss_exc_singlet = h2_diss_exc_triplet = hei_exc  = 0.;

	h2_excit_electr_bs_rate = h2_excit_electr_cp_rate = h2_excit_electr_rate 
		= h2_excit_vibr_rate = h2_excit_vibr_1_rate = h2_excit_vibr_2_rate = h2_excit_rot_rate = 0.;
	
	h2_excit_electr_bs = h2_excit_electr_cp = h2_excit_electr = h2_excit_vibr = h2_excit_vibr_1 = h2_excit_vibr_2 = 0.;

	// for initial time moment:
	for (i = 0; i < nb_of_el_energies; i++) {
		el_spectrum.arr[i] = NV_Ith_S(y, i);  // number of electrons in the energy interval (is not divided by energy bin),
	}
	el_spectrum_evol.push_back(el_spectrum);

	for (i = 0; i < NB_OF_CHEM_SPECIES; i++) {
		specimen_conc.arr[i] = NV_Ith_S(y, nb_of_el_energies + i);
	}
	specimen_conc_evol.push_back(specimen_conc);

	for (i = 0; i < nb_lev_h2; i++) {
		h2_popul_dens.arr[i] = NV_Ith_S(y, h2eq_nb + i);
	}
	h2_popdens_evol.push_back(h2_popul_dens);

	memset(h2_pop_v.arr, 0, MAX_H2_VSTATES_X1SU * sizeof(double));
	for (i = 0; i < nb_lev_h2; i++) {
		k = user_data.get_vibr_nb_h2(i);
		if (k >= 0)
			h2_pop_v.arr[k] += NV_Ith_S(y, h2eq_nb + i);
	}
	h2_popdens_v_evol.push_back(h2_pop_v);

	for (i = 0; i < nb_lev_hei; i++) {
		hei_popul_dens.arr[i] = NV_Ith_S(y, heieq_nb + i);
	}
	hei_popul_dens_evol.push_back(hei_popul_dens);

	// Energy loss rates
	enloss_rate_mt_arr.push_back(enloss_rate_mt);
	enloss_rate_h2_rot_arr.push_back(enloss_rate_h2_rot);
	enloss_rate_h2_vibr_arr.push_back(enloss_rate_h2_vibr);
	enloss_rate_h2_singlet_arr.push_back(enloss_rate_h2_singlet);
	enloss_rate_h2_triplet_arr.push_back(enloss_rate_h2_triplet);
	enloss_rate_ioniz_arr.push_back(enloss_rate_ioniz);
	enloss_rate_hei_arr.push_back(enloss_rate_hei);

	enloss_rate_coulomb_el_arr.push_back(enloss_rate_coulomb_el);
	diss_decay_heating_rate_arr.push_back(diss_decay_heating_rate);
	neutral_coll_heating_rate_arr.push_back(neutral_coll_heating_rate);
	
	enloss_rate_h2_rot_arr_pos.push_back(enloss_rate_h2_rot_pos);
	enloss_rate_h2_vibr_arr_pos.push_back(enloss_rate_h2_vibr_pos);
	
	// energy losses, time-integrated
	enloss_mt_arr.push_back(enloss_mt);
	enloss_h2_rot_arr.push_back(enloss_h2_rot);
	enloss_h2_vibr_arr.push_back(enloss_h2_vibr);
	enloss_h2_singlet_arr.push_back(enloss_h2_singlet);
	enloss_h2_triplet_arr.push_back(enloss_h2_triplet);
	enloss_ioniz_arr.push_back(enloss_ioniz);
	enloss_hei_arr.push_back(enloss_hei);

	enloss_coulomb_el_arr.push_back(enloss_coulomb_el);
	diss_decay_heating_arr.push_back(diss_decay_heating);
	neutral_coll_heating_arr.push_back(neutral_coll_heating);
	
	// Dissociation/excitation rates,
	h2_excit_electr_rate_arr.push_back(h2_excit_electr_rate);
	h2_excit_vibr_rate_arr.push_back(h2_excit_vibr_rate);
	h2_excit_rot_rate_arr.push_back(h2_excit_rot_rate);

	// Dissociation/excitation, time-integrated
	h2_solomon_diss_arr.push_back(h2_solomon_diss);
	h2_diss_exc_singlet_arr.push_back(h2_diss_exc_singlet);
	h2_diss_exc_triplet_arr.push_back(h2_diss_exc_triplet);
	hei_exc_arr.push_back(hei_exc);

	h2_excit_electr_bs_arr.push_back(h2_excit_electr_bs);
	h2_excit_electr_cp_arr.push_back(h2_excit_electr_cp);
	h2_excit_electr_arr.push_back(h2_excit_electr);

	h2_excit_vibr_arr.push_back(h2_excit_vibr);
	h2_excit_vibr_1_arr.push_back(h2_excit_vibr_1);
	h2_excit_vibr_2_arr.push_back(h2_excit_vibr_2);

	neutral_temp_arr.push_back(NV_Ith_S(y, physeq_nb));
	ion_temp_arr.push_back(NV_Ith_S(y, physeq_nb + 1));
	dust_charge_arr.push_back(NV_Ith_S(y, physeq_nb + 2) * grain_nb_density);  // dust charge per cm3

	// Start of the simulations: 
	cout << left << "Initialization step time (s): " << setw(12) << (int)(time(NULL) - timer) << endl
		<< "Starting simulations..." << endl;
	timer = time(NULL);

	if (verbosity) {
		cout << left << "lay_nb: " << lay_nb << endl
			<< left << setw(12) << "model_time " << setw(9) << "steps_nb "
			<< setw(12) << "gas_temp(K)" << setw(12) << "ion_temp(K)" << setw(12) << "dust_charge" << endl;
	}

	save_model_parameters(output_path, grb_cloud_distance, grb_distance, hcolumn_density, conc_h_tot,
		op_ratio_h2, dg_ratio, grain_radius, grain_nb_density, lay_nb);

	time_nb = 1;
	// there is integration of the parameters over the time, this time step must be sufficiently small, 
	model_time_step = 1.;

	while (model_time < time_cloud_arr.back()) 
	{	
		flag = CV_SUCCESS;
		must_be_saved = false;
		model_time_aux_prev = model_time_aux;

		if (model_time < time_cloud_arr[time_nb]) 
		{
			step_nb = 0;
			model_time_out = model_time + model_time_step;  // ?
			
			while (step_nb < MAX_NB_STEPS && flag == CV_SUCCESS && model_time < model_time_out) 
			{
				// actual model time is stored in cvode_mem, 
				flag = CVode(cvode_mem, model_time_out, y, &model_time, CV_ONE_STEP);  // CV_NORMAL or CV_ONE_STEP	
				step_nb++;
			}
		}
		
		if (flag == CV_SUCCESS) {	
			if (model_time > time_cloud_arr[time_nb]) {
				// step back to the time grid value, updating model time,
				flag = CVodeGetDky(cvode_mem, time_cloud_arr[time_nb], 0, y);
				model_time_aux = time_cloud_arr[time_nb];

				must_be_saved = true;
				time_nb++;
				
				if (flag != CV_SUCCESS) {
					cout << "Error in CVodeGetDky() " << flag << endl;
				}
			}
			else {
				model_time_aux = model_time;
			}
			
			// update the members of user_data class,
			f_elsp(model_time_aux, y, ydot, &user_data);
			
			dt = 0.5 * (model_time_aux - model_time_aux_prev);

			// Energy losses, 
			x1 = enloss_rate_mt;
			x2 = enloss_rate_h2_rot;
			x3 = enloss_rate_h2_vibr;
			x4 = enloss_rate_h2_singlet;
			x5 = enloss_rate_ioniz;
			x6 = enloss_rate_hei;
			x7 = enloss_rate_h2_triplet;

			x8 = enloss_rate_coulomb_el;
			x9 = diss_decay_heating_rate;
			x10 = neutral_coll_heating_rate;
			
			user_data.get_el_energy_losses(enloss_rate_mt, enloss_rate_h2_rot, enloss_rate_h2_rot_pos, enloss_rate_h2_vibr, enloss_rate_h2_vibr_pos, 
				enloss_rate_h2_singlet, enloss_rate_h2_triplet, enloss_rate_ioniz, enloss_rate_hei, enloss_rate_coulomb_el, 
				diss_decay_heating_rate, neutral_coll_heating_rate);

			enloss_mt += (x1 + enloss_rate_mt) * dt;
			enloss_h2_rot += (x2 + enloss_rate_h2_rot) * dt;
			enloss_h2_vibr += (x3 + enloss_rate_h2_vibr) * dt;
			enloss_h2_singlet += (x4 + enloss_rate_h2_singlet) * dt;
			enloss_ioniz += (x5 + enloss_rate_ioniz) * dt;
			enloss_hei += (x6 + enloss_rate_hei) * dt;
			enloss_h2_triplet += (x7 + enloss_rate_h2_triplet) * dt;

			enloss_coulomb_el += (x8 + enloss_rate_coulomb_el) * dt;
			diss_decay_heating += (x9 + diss_decay_heating_rate) * dt;
			neutral_coll_heating += (x10 + neutral_coll_heating_rate) * dt;
			
			// Excitation and other rates
			x1 = h2_solomon_diss_rate;
			x2 = h2_diss_exc_singlet_rate;
			x3 = h2_diss_exc_triplet_rate;
			x4 = hei_exc_rate;

			x5 = h2_excit_electr_rate;
			x6 = h2_excit_electr_bs_rate;
			x7 = h2_excit_electr_cp_rate;

			x8 = h2_excit_vibr_rate;
			x9 = h2_excit_vibr_1_rate;
			x10 = h2_excit_vibr_2_rate;

			user_data.get_h2_process_rates(h2_excit_electr_rate, h2_excit_electr_bs_rate, h2_excit_electr_cp_rate, 
				h2_excit_vibr_rate, h2_excit_vibr_1_rate, h2_excit_vibr_2_rate, h2_excit_rot_rate, 
				h2_solomon_diss_rate, h2_diss_exc_singlet_rate, h2_diss_exc_triplet_rate, hei_exc_rate);
			
			h2_solomon_diss += (x1 + h2_solomon_diss_rate) * dt;
			h2_diss_exc_singlet += (x2 + h2_diss_exc_singlet_rate) * dt;
			h2_diss_exc_triplet += (x3 + h2_diss_exc_triplet_rate) * dt;
			hei_exc += (x4 + hei_exc_rate) * dt;

			h2_excit_electr += (x5 + h2_excit_electr_rate) * dt;
			h2_excit_electr_bs += (x6 + h2_excit_electr_bs_rate) * dt;
			h2_excit_electr_cp += (x7 + h2_excit_electr_cp_rate) * dt;
			
			h2_excit_vibr += (x8 + h2_excit_vibr_rate) * dt;
			h2_excit_vibr_1 += (x9 + h2_excit_vibr_1_rate) * dt;
			h2_excit_vibr_2 += (x10 + h2_excit_vibr_2_rate) * dt;
			
			flag = CVodeGetNumSteps(cvode_mem, &nb_steps_tot);

			if (verbosity) {
				cout.precision(2);
				cout << left << setw(12) << model_time_aux << setw(9) << nb_steps_tot <<
					setw(12) << NV_Ith_S(y, physeq_nb) << setw(12) << NV_Ith_S(y, physeq_nb + 1) << setw(12) << NV_Ith_S(y, physeq_nb + 2) << endl;
			}
		
			// saving data,
			if (must_be_saved) {
				for (i = 0; i < nb_of_el_energies; i++) {
					el_spectrum.arr[i] = NV_Ith_S(y, i);  // number of electrons in the energy interval (is not divided by energy bin),
				}
				el_spectrum_evol.push_back(el_spectrum);

				for (i = 0; i < NB_OF_CHEM_SPECIES; i++) {
					specimen_conc.arr[i] = NV_Ith_S(y, nb_of_el_energies + i);
				}
				specimen_conc_evol.push_back(specimen_conc);

				for (i = 0; i < nb_lev_h2; i++) {
					h2_popul_dens.arr[i] = NV_Ith_S(y, h2eq_nb + i);
				}
				h2_popdens_evol.push_back(h2_popul_dens);

				memset(h2_pop_v.arr, 0, MAX_H2_VSTATES_X1SU * sizeof(double));
				for (i = 0; i < nb_lev_h2; i++) {
					k = user_data.get_vibr_nb_h2(i);
					if (k >= 0)
						h2_pop_v.arr[k] += NV_Ith_S(y, h2eq_nb + i);
				}
				h2_popdens_v_evol.push_back(h2_pop_v);

				for (i = 0; i < nb_lev_hei; i++) {
					hei_popul_dens.arr[i] = NV_Ith_S(y, heieq_nb + i);
				}
				hei_popul_dens_evol.push_back(hei_popul_dens);

				// Energy losses, rates
				enloss_rate_mt_arr.push_back(enloss_rate_mt);
				enloss_rate_h2_rot_arr.push_back(enloss_rate_h2_rot);
				enloss_rate_h2_vibr_arr.push_back(enloss_rate_h2_vibr);
				enloss_rate_h2_singlet_arr.push_back(enloss_rate_h2_singlet);
				enloss_rate_h2_triplet_arr.push_back(enloss_rate_h2_triplet);
				enloss_rate_ioniz_arr.push_back(enloss_rate_ioniz);
				enloss_rate_hei_arr.push_back(enloss_rate_hei);
				
				enloss_rate_coulomb_el_arr.push_back(enloss_rate_coulomb_el);
				diss_decay_heating_rate_arr.push_back(diss_decay_heating_rate);
				neutral_coll_heating_rate_arr.push_back(neutral_coll_heating_rate);
				
				enloss_rate_h2_rot_arr_pos.push_back(enloss_rate_h2_rot_pos);
				enloss_rate_h2_vibr_arr_pos.push_back(enloss_rate_h2_vibr_pos);
				
				// energy losses, time-integrated,
				enloss_mt_arr.push_back(enloss_mt);
				enloss_h2_rot_arr.push_back(enloss_h2_rot);
				enloss_h2_vibr_arr.push_back(enloss_h2_vibr);
				enloss_h2_singlet_arr.push_back(enloss_h2_singlet);
				enloss_h2_triplet_arr.push_back(enloss_h2_triplet);
				enloss_ioniz_arr.push_back(enloss_ioniz);
				enloss_hei_arr.push_back(enloss_hei);

				enloss_coulomb_el_arr.push_back(enloss_coulomb_el);
				diss_decay_heating_arr.push_back(diss_decay_heating);
				neutral_coll_heating_arr.push_back(neutral_coll_heating);
				
				// Excitations / dissociations, rates
				h2_excit_electr_rate_arr.push_back(h2_excit_electr_rate);
				h2_excit_vibr_rate_arr.push_back(h2_excit_vibr_rate);
				h2_excit_rot_rate_arr.push_back(h2_excit_rot_rate);

				// Excitations / dissociations, time-integrated,
				h2_solomon_diss_arr.push_back(h2_solomon_diss);
				h2_diss_exc_singlet_arr.push_back(h2_diss_exc_singlet);
				h2_diss_exc_triplet_arr.push_back(h2_diss_exc_triplet);
				hei_exc_arr.push_back(hei_exc);
				
				h2_excit_electr_arr.push_back(h2_excit_electr);
				h2_excit_electr_bs_arr.push_back(h2_excit_electr_bs);
				h2_excit_electr_cp_arr.push_back(h2_excit_electr_cp);

				h2_excit_vibr_arr.push_back(h2_excit_vibr);
				h2_excit_vibr_1_arr.push_back(h2_excit_vibr_1);
				h2_excit_vibr_2_arr.push_back(h2_excit_vibr_2);

				neutral_temp_arr.push_back(NV_Ith_S(y, physeq_nb));
				ion_temp_arr.push_back(NV_Ith_S(y, physeq_nb + 1));
				dust_charge_arr.push_back(NV_Ith_S(y, physeq_nb + 2) * grain_nb_density);  // dust charge per cm3

				save_electron_spectrum_evolution(output_path, electron_energies_grid, electron_energy_bin_size, time_cloud_arr, el_spectrum_evol, 
					grb_distance, hcolumn_density, conc_h_tot);

				save_specimen_conc_evolution(output_path, time_cloud_arr, specimen_conc_evol, grb_distance, hcolumn_density, conc_h_tot);

				save_h2_populations_evolution(output_path, time_cloud_arr, h2_popdens_evol, h2_popdens_v_evol,
					grb_distance, hcolumn_density, conc_h_tot, nb_lev_h2);
			
				save_hei_populations_evolution(output_path, time_cloud_arr, hei_popul_dens_evol, grb_distance, hcolumn_density, conc_h_tot, nb_lev_hei);

				save_electron_energy_loss_rates(output_path, grb_distance, hcolumn_density, conc_h_tot, time_cloud_arr,
					enloss_rate_mt_arr, 
					enloss_rate_h2_rot_arr,
					enloss_rate_h2_rot_arr_pos,
					enloss_rate_h2_vibr_arr,
					enloss_rate_h2_vibr_arr_pos,
					enloss_rate_h2_singlet_arr,
					enloss_rate_h2_triplet_arr,
					enloss_rate_ioniz_arr, 
					enloss_rate_hei_arr, 
					enloss_rate_coulomb_el_arr,
					diss_decay_heating_rate_arr,
					neutral_coll_heating_rate_arr);

				save_electron_energy_losses(output_path, grb_distance, hcolumn_density, conc_h_tot, time_cloud_arr,
					enloss_mt_arr,
					enloss_h2_rot_arr,
					enloss_h2_vibr_arr,
					enloss_h2_singlet_arr,
					enloss_h2_triplet_arr,
					enloss_ioniz_arr,
					enloss_hei_arr, 
					enloss_coulomb_el_arr, 
					diss_decay_heating_arr, 
					neutral_coll_heating_arr);

				save_diss_excit_rates(output_path, grb_distance, hcolumn_density, conc_h_tot, time_cloud_arr, 
					h2_excit_electr_rate_arr, h2_excit_vibr_rate_arr, h2_excit_rot_rate_arr);

				save_diss_excit(output_path, grb_distance, hcolumn_density, conc_h_tot, time_cloud_arr,
					h2_solomon_diss_arr, 
					h2_diss_exc_singlet_arr, 
					h2_diss_exc_triplet_arr, 
					hei_exc_arr, 
					h2_excit_electr_arr, 
					h2_excit_electr_bs_arr,
					h2_excit_electr_cp_arr, 
					h2_excit_vibr_arr, 
					h2_excit_vibr_1_arr, 
					h2_excit_vibr_2_arr);

				save_phys_parameters(output_path, time_cloud_arr, neutral_temp_arr, ion_temp_arr, dust_charge_arr, specimen_conc_evol, 
					grb_distance, hcolumn_density, conc_h_tot);
			}
		}
		else {
			cout << "Some unexpected error has been occurred: " << flag << endl;
			break;
		}

	}
	cout << left << "Total calculation time (s): " << setw(12) << (int)(time(NULL) - timer) << endl;
	
	// Saving parameters of the electron energy degradation,
	fname = output_path + "h2grb_parameters_output.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(4);

	// Energy in electrons, number density of electrons at the start of simulations,
	x3 = x4 = 0.;
	for (i = 0; i < nb_of_el_energies; i++) {
		x3 += el_spectrum_evol.front().arr[i] * 0.5 * (electron_energies_grid[i] + electron_energies_grid[i + 1]);  // [eV cm-3]
		x4 += el_spectrum_evol.front().arr[i];  // [cm-3]
	}

	// simple species are: "e-", "H", "H+", "H2", "H2+", "He", "He+", "He++"
	// number density of ions produced, at the start of the simulations,
	x1 = specimen_conc_evol.front().arr[2] + specimen_conc_evol.front().arr[4] + specimen_conc_evol.front().arr[6]
		+ 2. * specimen_conc_evol.front().arr[7];  // [cm-3]

	// number density of ions produced, at the end of the simulations,
	x2 = specimen_conc_evol.back().arr[2] + specimen_conc_evol.back().arr[4] + specimen_conc_evol.back().arr[6]
		+ 2. * specimen_conc_evol.back().arr[7];  // [cm-3]
		
	output << left << "Parameters of the electron energy degradation:" << endl
		<< setw(38) << "Initial number of electrons, [cm-3]: "<< x4 << endl
		<< setw(38) << "Number of ionizations N[cm-3]: "      << x2 - x1 << endl
		<< setw(38) << "Initial electron energy E[eV cm-3]: " << x3 << endl
		<< setw(38) << "Mean energy per ion pair W: "         << x3 / (x2 - x1) << endl
		<< setw(38) << "Heating efficiency of neutral gas: "  << (fabs(enloss_mt_arr.back()) + neutral_coll_heating_arr.back()) / x3 << endl
		<< setw(38) << "Heating efficiency of electrons: "    << fabs(enloss_coulomb_el_arr.back()) / x3 << endl
		<< setw(38) << "Energy fraction lost in ro-vibr.: "   << fabs(enloss_h2_vibr_arr.back()) / x3 << endl;

	x4 = h2_diss_exc_singlet_arr.back() + h2_diss_exc_triplet_arr.back() + h2_solomon_diss_arr.back();

	output << left
		<< setw(38) << "H2 dissociations [cm-3]: "           << x4 << endl
		<< setw(38) << "H2 dissociations, per ion pair: "    << x4 / (x2 - x1) << endl
		<< setw(38) << "Electronic excitation, per ion pair" << h2_excit_electr_arr.back() / (x2 - x1) << endl
		<< setw(38) << "Excitation to B1S+u, per ion pair: " << h2_excit_electr_bs_arr.back() / (x2 - x1) << endl
		<< setw(38) << "Excitation to C1Pu, per ion pair: "  << h2_excit_electr_cp_arr.back() / (x2 - x1) << endl
		<< setw(38) << "Vibrational excitation [cm-3]: "     << h2_excit_vibr_arr.back() << endl
		<< setw(38) << "Excitation to v=1 [cm-3]: "          << h2_excit_vibr_1_arr.back() << endl
		<< setw(38) << "Ratio of excitations v=2/v=1: "      << h2_excit_vibr_2_arr.back() / h2_excit_vibr_1_arr.back() << endl;
	
	output.close();

	if (verbosity) {
		cout << "The memory is freeing up" << endl;
	}
	N_VDestroy(y);
	N_VDestroy(ydot);
	N_VDestroy(abs_tol);

	// Free integrator memory
	CVodeFree(&cvode_mem);

	// Free the linear solver memory
	SUNLinSolFree(LS);

	// Free the matrix memory
	SUNMatDestroy(A);
}


//
// please, check the value of MAX_TEXT_LINE_WIDTH
void init_dust_abund(const string& sim_path, vector<double>& d, double & grain_radius_init, double& grain_nb_density, double& grain_material_density)
{
	int i, j, j0, nb_lay, nb_t;
	char ch, text_line[MAX_TEXT_LINE_WIDTH];
	double x;

	string fname;
	ifstream input;

	// path to the folder with simulation data (GRB propagation results),
	fname = sim_path + "grb_dust_evol.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> ch >> grain_radius_init >> grain_nb_density >> grain_material_density;
	input >> ch >> nb_lay >> nb_t;
	input >> ch; 
	
	// dust abundances at the end of the GRB signal propagation, last column of data,
	j0 = nb_t - 1;

	for (j = 0; j < nb_t; j++) {
		input >> x;
	}
	for (i = 0; i < nb_lay; i++) {
		input >> x;
		for (j = 0; j < nb_t; j++) {
			input >> x;
			if (j == j0) {
				d.push_back(x);
			}
		}
	}
	input.close();
}

// for test simulations,
void init_electron_spectra_test(const string& output_path, dynamic_array & init_data_test, const vector<double>& en_grid, 
	double test_electron_conc, double test_electron_energy)
{
	int i;
	double de;

	for (i = 0; i < (int)en_grid.size(); i++) {
		if (en_grid[i] > test_electron_energy) {
			i--;
			break;
		}
	}

	// in the case when test energy is very high,
	if (i == (int)en_grid.size()) {
		i = (int)en_grid.size() - 2;
	}
	de = en_grid[i + 1] - en_grid[i];

	memset(init_data_test.arr, 0, init_data_test.dim * sizeof(double));
	init_data_test.arr[i] = test_electron_conc / de;
}


void init_electron_spectra(const string& sim_path, const string& output_path, vector<dynamic_array>& init_data, const vector<double>& en_grid)
{ 
	// layer nb for which the electron spectra is saved (initial and after rescaling of energy grid),
	const int lay_nb_test = 0; 
	
	int i, j, l, k, nb_l, en_nb_f, en_nb, nb_bins;
	char ch, text_line[MAX_TEXT_LINE_WIDTH];
	double e, de, x, min_en, max_en, conc_h_tot, distance;

	string fname, str;
	ifstream input;
	ofstream output;
	
	vector<dynamic_array> file_data;
	vector<double> en_grid_file;

	en_nb = (int)en_grid.size() - 1;  // nb of energy bins,

	fname = sim_path + "grb_electron_spectrum.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> ch >> conc_h_tot >> distance >> min_en >> max_en >> nb_bins >> en_nb_f >> nb_l;
	dynamic_array el_spectrum_file(en_nb_f), el_spectrum(en_nb);

	input >> str >> str;
	for (i = 0; i < nb_l; i++) {
		input >> x;
	}

	for (j = 0; j < nb_l; j++) {
		file_data.push_back(el_spectrum_file);  // empty arrays
		init_data.push_back(el_spectrum);
	}

	// in the file the energy spectrum of electrons is given in [cm-3 eV-1], half length of the energy interval is provided,
	en_grid_file.push_back(0.);
	for (i = 0; i < en_nb_f; i++) {
		input >> e >> de;
		en_grid_file.push_back(e + de);

		for (l = 0; l < nb_l; l++) {
			input >> x;
			// vector index - layer nb, array index - electron energy bin nb,	
			file_data[l].arr[i] = x;
		}
	}
	input.close();

	// reformatting the input data,
	for (l = 0; l < nb_l; l++) {
		for (i = 0; i < en_nb; i++) {
			// after cycle en_grid_file[j] >= en_grid[i]	
			for (j = 0; j < (int)en_grid_file.size() - 1 && en_grid_file[j] < en_grid[i]; j++) { ; }

			// after cycle en_grid_file[k] >= en_grid[i+1]
			for (k = j; k < (int)en_grid_file.size() - 1 && en_grid_file[k] < en_grid[i + 1]; k++) { ; }

			init_data[l].arr[i] = 0.;
			if (j > 0 && en_grid_file[j] > en_grid[i]) {
				init_data[l].arr[i] += (en_grid_file[j] - en_grid[i]) * file_data[l].arr[j - 1];
			}

			for (; j < k; j++) {
				init_data[l].arr[i] += file_data[l].arr[j] * (en_grid_file[j + 1] - en_grid_file[j]);
			}

			if (k > 0 && en_grid_file[k] > en_grid[i + 1]) {
				init_data[l].arr[i] -= (en_grid_file[k] - en_grid[i + 1]) * file_data[l].arr[k - 1];
			}
			init_data[l].arr[i] /= en_grid[i + 1] - en_grid[i];
		}
	}

	// Saving electron spectra for test,
	fname = output_path + "electron_spectra_file.txt";
	output.open(fname.c_str());
	output << scientific;
	
	output << left << "!Electron spectra from file with simulation data," << endl
		<< "!Histogram, energy (eV), spectra [cm-3 eV-1] (as in the file with simulation data), layer nb: " << lay_nb_test << endl;

	for (i = 0; i < en_nb_f; i++) {
		output.precision(5);
		output << left << setw(14) << en_grid_file[i] * 1.0001; 
		
		output.precision(3);
		output << left << setw(12) << file_data[lay_nb_test].arr[i] << endl;

		output.precision(5);
		output << left << setw(14) << en_grid_file[i + 1];

		output.precision(3);
		output << left << setw(12) << file_data[lay_nb_test].arr[i] << endl;
	}
	output.close();

	fname = output_path + "electron_spectra_init.txt";
	output.open(fname.c_str());
	output << scientific;

	output << left << "!Initial electron spectra in simulations," << endl
		<< "!Histogram, energy (eV), spectra [cm-3 eV-1], layer nb: " << lay_nb_test << endl;

	for (i = 0; i < en_nb; i++) {
		output.precision(5);
		output << left << setw(14) << en_grid[i]* 1.0001;

		output.precision(3);
		output << left << setw(12) << init_data[lay_nb_test].arr[i] << endl;

		output.precision(5);
		output << left << setw(14) << en_grid[i + 1];

		output.precision(3);
		output << left << setw(12) << init_data[lay_nb_test].arr[i] << endl;
	}
	output.close();

	fname = output_path + "electron_nb_density_test.txt";
	output.open(fname.c_str());
	output << scientific;
	output.precision(3);

	output << left << "!The initial number density of electrons, from initial and rebinned spectra," << endl
		<< "!Layer nb, number density (file), number density (re-binned) [cm-3]:" << endl;

	for (l = 0; l < nb_l; l++) 
	{
		x = 0.;
		for (i = 0; i < en_nb_f; i++) {
			x += file_data[l].arr[i] * (en_grid_file[i + 1] - en_grid_file[i]);
		}
		output << left << setw(5) << l << setw(12) << x;

		x = 0.;
		for (i = 0; i < en_nb; i++) {
			x += init_data[l].arr[i] * (en_grid[i + 1] - en_grid[i]);
		}
		output << left << setw(12) << x << endl;
	}
	output.close();
}

//
void init_specimen_conc(const string& path, double& conc_h_tot, vector<double>& layer_centre_distances, vector<dynamic_array>& d)
{
	int i, j, k, nb_l, nb_sp;
	double x;
	char ch, text_line[MAX_TEXT_LINE_WIDTH];

	string fname, sn;
	stringstream ss;
	ifstream input;

	vector<string> chemical_species_file;
	dynamic_array specimen_conc(NB_OF_CHEM_SPECIES);

	fname = path + "grb_specimen_conc.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> ch >> conc_h_tot >> nb_l >> nb_sp;
	if (nb_sp != NB_OF_CHEM_SPECIES) {
		cout << "Note: nb of species in input file does not coincide with that in the code." << endl;
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);  // reading the end of the previous line
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	ss.clear();
	ss.str(text_line);

	ss >> sn >> sn;  // nb of the layer, distance
	for (i = 0; i < nb_sp; i++) {
		ss >> sn;
		chemical_species_file.push_back(sn);
	}

	for (i = 0; i < nb_l; i++) {
		d.push_back(specimen_conc);
		
		input >> j >> x;
		layer_centre_distances.push_back(x);

		for (j = 0; j < nb_sp; j++) {
			input >> x;
			
			for (k = 0; k < NB_OF_CHEM_SPECIES; k++) {
				if (chemical_species_file[j] == chemical_species[k])
					d[i].arr[k] = x;  // vector index - layer nb, array index - specimen nb,
			}
		}
	}
	input.close();
}

// In the file, the levels are denoted by the number,
// Note, there is no check yet about the consistency of energy levels used in this code and those used in the simulations before, 
void init_h2_population_density(const string& path, vector<dynamic_array>& d, vector<int>& qnb_v_arr, vector<int>& qnb_j_arr)
{
	int i, j, v, nb, nb_l, nb_lev_h2_f;
	double x, conc_h_tot, distance;
	char ch, text_line[MAX_TEXT_LINE_WIDTH];

	string fname, sn;
	stringstream ss;
	ifstream input;

	// Reading the data file with the list of H2 levels of the ground electronic state,
	// these levels were used in GRB propagation simulations,
	qnb_j_arr.clear();
	qnb_v_arr.clear();

	fname = path + "grb_h2_levels.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> nb_lev_h2_f;
	for (i = 0; i < nb_lev_h2_f; i++) {
		input >> nb >> v >> j >> x;
		qnb_v_arr.push_back(v);
		qnb_j_arr.push_back(j);
	}
	input.close();

	// reading the data file with the population densitities of H2 molecule,
	fname = path + "grb_h2_popdens.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> ch >> conc_h_tot >> distance >> nb_l >> nb_lev_h2_f;

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);  // reading the end of the previous line
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	dynamic_array population_densities(nb_lev_h2_f);

	for (i = 0; i < nb_l; i++) {
		d.push_back(population_densities);
		input >> j >> x;

		for (j = 0; j < nb_lev_h2_f; j++) {
			input >> x;
			d[i].arr[j] = x;  // vector index - layer nb, array index - molecule level nb,
		}
	}
	input.close();
}


// calculation of data tables for electron Coulomb scattering
/*
*
void calc_electron_scattering_tables(const string& output_path);

void calc_electron_scattering_tables(const string& output_path)
{
	const int nb_th = 1000;
	const int nb_phi = 100;
	const double theta0 = 0.001;  // radians

	int i, j, k, l, m, p, p1, p2, nb_of_el_energies;
	double x, y, z, w, cs1, cs2, cs3, dcs, en1, en2, en0, dphi, dth, phi(0.), theta(0.);

	double* energies, *energies_grid, *s;
	double** coloumb_indexes;

	string fname;
	ofstream output;

	time_t timer;

	vector<double> energies_grid_vector;
	electron_coulomb_scattering el_scatt_cs;

	init_electron_energy_grid(energies_grid_vector);
	nb_of_el_energies = (int) energies_grid_vector.size() - 1;

	energies_grid = new double[nb_of_el_energies + 1];
	energies = new double[nb_of_el_energies];

	energies_grid[0] = energies_grid_vector[0];
	for (i = 0; i < nb_of_el_energies; i++) {
		energies_grid[i + 1] = energies_grid_vector[i + 1];
		energies[i] = 0.5*(energies_grid[i] + energies_grid[i + 1]);
	}

	s = new double[nb_of_el_energies];
	coloumb_indexes = alloc_2d_array<double>(nb_of_el_energies*(nb_of_el_energies - 1)/2, nb_of_el_energies);

	dth = (0.5 * M_PI - theta0) / nb_th;  // sin(th), cos(th) > 0.
	dphi = M_PI / nb_phi;

	cs1 = 0.;

/*	i = 500;
	j = 100;
	memset(*coloumb_indexes, 0, (nb_of_el_energies * nb_of_el_energies * (nb_of_el_energies - 1) / 2) *sizeof(double));

	for (l = 0; l < nb_phi; l++) {
		phi = dphi * (l + 0.5);
		cout << left << setw(12) << cs1;

		for (k = 0; k < nb_th; k++) {
			theta = theta0 + dth * (k + 0.5);

			// [cm2 *cm/s]
			dcs = el_scatt_cs(energies[i], energies[j], phi, theta) * dth * 0.5 * sin(phi) * dphi;

			en1 = el_scatt_cs.get_energy_gain(energies[i], energies[j], phi, theta, M_PI);  // maximal negative value
			en2 = el_scatt_cs.get_energy_gain(energies[i], energies[j], phi, theta, 0.);

			en0 = 0.5 * (en1 + en2);
			x = 0.5 * (en2 - en1);

			hunt_index(energies_grid, nb_of_el_energies + 1, energies[i] + en1, p1);
			hunt_index(energies_grid, nb_of_el_energies + 1, energies[i] + en2, p2);

			y = (energies_grid[p1] - energies[i] - en0) / x;
			if (y < -1.)
				y = -1.;
			y = acos(y);

			for (p = p1; p <= p2; p++) {
				z = (energies_grid[p + 1] - energies[i] - en0) / x;
				if (z > 1.)
					z = 1.;
				z = acos(z);

				coloumb_indexes[i * (i + 1) / 2 + j][p] += dcs * (y - z) / M_PI;
				y = z;
			}

			hunt_index(energies_grid, nb_of_el_energies + 1, energies[j] - en1, p2);
			hunt_index(energies_grid, nb_of_el_energies + 1, energies[j] - en2, p1);

			y = (energies_grid[p1] - energies[j] + en0) / x;
			if (y < -1.)
				y = -1.;
			y = acos(y);

			for (p = p1; p <= p2; p++) {
				z = (energies_grid[p + 1] - energies[j] + en0) / x;
				if (z > 1.)
					z = 1.;
				z = acos(z);

				coloumb_indexes[i * (i + 1) / 2 + j][p] += dcs * (y - z) / M_PI;
				y = z;
			}
			cs1 += dcs;
		}
	}

	fname = output_path + "el_scatt_tables.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(3);

	output << left << "! " << setw(7) << i << setw(7) << j << endl;
	for (p = 0; p < nb_of_el_energies; p++) {
		output << left << setw(13) << energies[p] << setw(13) << coloumb_indexes[i * (i + 1) / 2 + j][p] << endl;
	}
	output.close();


	timer = time(NULL);
	memset(*coloumb_indexes, 0, (nb_of_el_energies * nb_of_el_energies * (nb_of_el_energies - 1) / 2) * sizeof(double));

	i = 500;
	for (j = 0; j < i; j++) {
		coloumb_scatt_func c_func(energies_grid, nb_of_el_energies);
		cs2 = 0.;
		c_func.set_electron_energies(energies[i], energies[j]);

		for (l = 0; l < nb_phi; l++) {
			phi = dphi * (l + 0.5);

			c_func.set_phi(phi);
			dcs = qromb_mod<coloumb_scatt_func>(c_func, theta0, 0.5 * M_PI, s, nb_of_el_energies, 0.001, false)
				* 0.5 * sin(phi) * dphi;

			for (p = 0; p < nb_of_el_energies; p++) {
				coloumb_indexes[i * (i + 1) / 2 + j][p] += s[p] * 0.5 * sin(phi) * dphi;
			}
			cs2 += dcs;
		}
		cs3 = en0 = 0.;
		for (p = 0; p < nb_of_el_energies; p++) {
			cs3 += coloumb_indexes[i * (i + 1) / 2 + j][p];
			en0 += coloumb_indexes[i * (i + 1) / 2 + j][p] * energies[p];
		}
		cs3 *= 0.5;
		cout << left << setw(12) << j << setw(12) << (int)(time(NULL) - timer) << en0/cs3 - energies[i]- energies[j] << endl;

	}

	fname = output_path + "el_scatt_tables_2.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(3);

	output << left << "! " << setw(7) << i << endl;
	for (j = 0; j < i; j++) {
		output << left << setw(7) << j;
	}
	output << endl;

	for (p = 0; p < nb_of_el_energies; p++) {
		for (j = 0; j < i; j++) {
			output << left << setw(13) << energies[p] << setw(13) << coloumb_indexes[i * (i + 1) / 2 + j][p];
		}
		output << endl;
	}
	output.close();

	cout << endl << "Simple integration: " << cs1 << endl
		<< "qromb/qsimp: " << cs2 << endl
		<< "control sum: " << cs3 << endl;
}
*/