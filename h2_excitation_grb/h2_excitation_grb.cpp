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
#include <sundials/sundials_types.h>   /* definition of type realtype          */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <limits>
#include <vector>
#include <array>
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
const int VERBOSITY = 1;

// sim_path - path to the folder with simulation data (GRB propagation results)
// output_path - path to the folder where data of the current simulations must be saved,
void init_specimen_conc(const string& sim_path, double& conc_h_tot, vector<double>& layer_centre_distances, vector<dynamic_array>& d);
void init_dust_abund(const string& sim_path, vector<double>& d, double& grain_radius_init, double& grain_nb_density, double& grain_material_density);
void init_electron_spectra(const string& sim_path, const string& output_path, vector<dynamic_array>& init_data, const vector<double>& electron_energies_grid);
void init_electron_spectra_mono(dynamic_array & init_data, const vector<double>& en_grid, 
	double electron_conc, double electron_energy);
void init_h2_population_density(const string& sim_path, vector<dynamic_array>& d, vector<int>& qnb_v_arr, vector<int>& qnb_j_arr);


static int f_elsp(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
	static_cast<elspectra_evolution_data*>(user_data)->f(t, y, ydot);
	return (0);
}

//
// internal parameters used to save simulation data,
//
bool dust_is_presented;
int nb_of_time_moments, nb_of_equat, nb_of_el_energies, nb_lev_h2, nb_lev_hei, h2eq_nb, heieq_nb, physeq_nb;

double grain_radius_init, grain_nb_density, grain_material_density, grain_radius, dust_mass, gas_mass, dg_ratio;

double enloss_rate_mt, enloss_rate_h2_rot, enloss_rate_h2_rot_pos, enloss_rate_h2_vibr, enloss_rate_h2_vibr_pos, enloss_rate_ioniz, 
	enloss_rate_hei, enloss_rate_coulomb_el, neutral_coll_heating_rate;

double enloss_mt, enloss_h2_rot, enloss_h2_vibr, enloss_ioniz, enloss_hei, enloss_coulomb_el, neutral_coll_heating;

double hei_excit_rate, h2_excit_rot_rate;
double hei_excit, h2_excit_rot;

// the array with initial values of variables,
double* initial_data;

vector<double> time_cloud_arr;      // in s
vector<double> electron_energies_grid, electron_energy_bin_size;  // in eV,

// array for saving the data on excitation of each of the states, 
array< electronic_excitation_data_unit, NB_EXC_ELECTRONIC_STATES> h2_state_data, h2_state_data_rates, h2_state_data_aux;

// Excitation and dissociation rates [cm-3 s-1], energy loss rate [eV cm-3 s-1], as a function of time
vector< array<electronic_excitation_data_unit, NB_EXC_ELECTRONIC_STATES> > h2_state_data_rate_arr;
vector< dynamic_array> h2_electr_vstates_rate_arr, h2_vibr_vstates_rate_arr;

// Number of excitations and dissociations [cm-3], energy loss [eV cm-3] during the time interval [0, t] as a function of time t,
vector< array<electronic_excitation_data_unit, NB_EXC_ELECTRONIC_STATES> > h2_state_data_arr;
vector< dynamic_array> h2_electr_vstates_arr, h2_vibr_vstates_arr;

// Energy loss rates in [eV cm-3 s-1], as a function of time
vector<double> enloss_rate_mt_arr, enloss_rate_h2_rot_arr, enloss_rate_h2_rot_arr_pos, enloss_rate_h2_vibr_arr, enloss_rate_h2_vibr_arr_pos,
	enloss_rate_ioniz_arr, enloss_rate_hei_arr, enloss_rate_coulomb_el_arr, neutral_coll_heating_rate_arr;

// Energy loss in [eV cm-3] during the time interval [0, t] as a function of time t,
vector<double> enloss_mt_arr, enloss_h2_rot_arr, enloss_h2_vibr_arr, enloss_ioniz_arr, enloss_hei_arr, 
	enloss_coulomb_el_arr, neutral_coll_heating_arr;

// Pure rotational excitations as a function of time, rates
vector<double> h2_excit_rot_rate_arr;

// Number of excitations in [cm-3] up to a given time, as a function of time
vector<double> h2_excit_rot_arr, hei_excit_arr;

// Physical parameters as a function of time,
vector<double> neutral_temp_arr, ion_temp_arr, dust_charge_arr;

// Parameter[layer nb], at the fixed time (start of electron spectra simulations),
vector<double> dg_ratio_cloud_data, ioniz_cloud_data;

// Spectrum[layer nb][energy bin] in [cm-3 eV-1], vector index - layer nb, array index - energy bin nb (specimen nb, H2 level nb),
vector<dynamic_array> el_spectrum_init_data, specimen_conc_init_data;
vector<dynamic_array> h2_pop_dens_init_data;

// Electron spectrum[time nb][energy bin] in [cm-3], vector index - time moment nb, array index - energy bin nb 
vector<dynamic_array> el_spectrum_evol;
// Concentration[time nb][specimen nb], or [time nb][level nb]
vector<dynamic_array> specimen_conc_evol, h2_popdens_evol, h2_popdens_v_evol, hei_popul_dens_evol;  // in [cm-3]

// these are called inside the simulations routine:
void memory_freeing_up();
void init_time_intervals();
void init_electron_energy_grid(double max_electron_energy);

void simulations_monoenergetic(const string& data_path, const string& output_path, double conc_h_tot, double he_abund, double ioniz_fract,
	double op_ratio_h2, double electron_conc, double electron_energy);

void simulations_grb(const string& parameter_fname);


// path to the directory with data tables(spectroscopic data, collision rates and etc.):
// linux - "/disk4/nester/input_data/"
// windows - "C:/Users/Александр/Documents/input_data/"
int main()
{	
	char text_line[MAX_TEXT_LINE_WIDTH];

	int i, nb_processors(6), nb;
	double op_ratio_h2, conc_h_tot, ioniz_fract;
	double x, electron_conc, he_abund;
	
	vector<double> electron_energy_arr;

	string fname;
	string data_path;    // path to the input_data folder with spectroscopic, collisional data and etc.
	string output_path;  // path to the folder where data of the current simulation must be saved,
	string sim_path;     // path to the folder with simulation data (GRB propagation results)

	stringstream ss;
	ifstream input;

	//
	// Initialization of initial data
	//
	fname = "input_parameters_mono.txt";

	input.open(fname.c_str());
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with cross section data " << fname << endl;
		exit(1);
	}

	while (!input.eof())
	{
		do // comment lines are read:
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> nb_processors;

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> data_path;

		// path to the directory for the output
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

#ifdef __linux__
		// if path has not spaces:
		ss.clear();
		ss.str(text_line);
		ss >> output_path;
#else
		// path may have spaces (for linux does not work)
		output_path = text_line;
#endif

		// H nuclei total concentration, in [cm-3],
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> conc_h_tot;

		// Helium abundance (only H2 and He are taken into account),
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> he_abund;

		// Ionization fraction x_e, number density of electrons n_e = n_{H,tot} * x_e,
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> ioniz_fract;

		// ortho-to-para-H2 ratio,
		// test the dependence on this parameter,
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> op_ratio_h2;

		// electron concentration, in [cm-3]
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> electron_conc;

		// Nb of electron energies to calculate:
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> nb;

		// Electron energies, in [eV]
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);

		electron_energy_arr.clear();
		for (i = 0; i < nb; i++) {
			ss >> x;
			electron_energy_arr.push_back(x);
		}
		input.close();
	}

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
#else
	cout << "No OpenMP" << endl;
#endif

	
	for (i = 0; i < nb; i++) {
		simulations_monoenergetic(data_path, output_path, conc_h_tot, he_abund, ioniz_fract, op_ratio_h2, electron_conc, 
			electron_energy_arr[i]);
	}

	output_path = "C:/Users/Александр/Documents/Данные и графики/paper GRB in molecular cloud/python_scripts_ism/";
	data_path   = "C:/Users/Александр/Documents/input_data/";

	//save_cross_section_table(output_path, data_path, IONIZATION_FRACTION_THERMAL_EL, THERMAL_EL_TEMPERATURE);
	//calc_helium_lifetimes(output_path, data_path);
}


/*----------------------------------------------------------------------------------------------------------------------*/
//
/*----------------------------------------------------------------------------------------------------------------------*/

void memory_freeing_up()
{
	int i;
	if (VERBOSITY) {
		cout << "The memory is freeing up" << endl;
	}
	delete[] initial_data;

	enloss_rate_mt = enloss_rate_h2_rot = enloss_rate_h2_rot_pos = enloss_rate_h2_vibr = enloss_rate_h2_vibr_pos 
		= enloss_rate_ioniz = enloss_rate_hei = enloss_rate_coulomb_el = neutral_coll_heating_rate = 0.;

	enloss_mt = enloss_h2_rot = enloss_h2_vibr = enloss_ioniz = enloss_hei = enloss_coulomb_el = neutral_coll_heating = 0.;

	h2_excit_rot_rate = hei_excit_rate = 0.;
	h2_excit_rot = hei_excit = 0.;

	for (i = 0; i < NB_EXC_ELECTRONIC_STATES; i++) {
		h2_state_data[i].excit = h2_state_data[i].direct_diss = h2_state_data[i].twostep_diss =
			h2_state_data[i].diss_heat_input = h2_state_data[i].enloss_excit = h2_state_data[i].enloss_direct_diss = 0;

		h2_state_data_rates[i].excit = h2_state_data_rates[i].direct_diss = h2_state_data_rates[i].twostep_diss =
			h2_state_data_rates[i].diss_heat_input = h2_state_data_rates[i].enloss_excit = h2_state_data_rates[i].enloss_direct_diss = 0;

		h2_state_data_aux[i].excit = h2_state_data_aux[i].direct_diss = h2_state_data_aux[i].twostep_diss =
			h2_state_data_aux[i].diss_heat_input = h2_state_data_aux[i].enloss_excit = h2_state_data_aux[i].enloss_direct_diss = 0;
	}

	h2_state_data_rate_arr.clear();
	h2_electr_vstates_rate_arr.clear();
	h2_vibr_vstates_rate_arr.clear();

	h2_state_data_arr.clear();
	h2_electr_vstates_arr.clear();
	h2_vibr_vstates_arr.clear();

	enloss_rate_mt_arr.clear();
	enloss_rate_h2_rot_arr.clear();
	enloss_rate_h2_rot_arr_pos.clear();
	enloss_rate_h2_vibr_arr.clear();
	enloss_rate_h2_vibr_arr_pos.clear();
	enloss_rate_ioniz_arr.clear(); 
	enloss_rate_hei_arr.clear();
	enloss_rate_coulomb_el_arr.clear();
	neutral_coll_heating_rate_arr.clear();

	enloss_mt_arr.clear(); 
	enloss_h2_rot_arr.clear(); 
	enloss_h2_vibr_arr.clear();
	enloss_ioniz_arr.clear();
	enloss_hei_arr.clear(); 
	enloss_coulomb_el_arr.clear(); 
	neutral_coll_heating_arr.clear();

	h2_excit_rot_rate_arr.clear();
	h2_excit_rot_arr.clear();
	hei_excit_arr.clear();

	neutral_temp_arr.clear();
	ion_temp_arr.clear();
	dust_charge_arr.clear();

	dg_ratio_cloud_data.clear(); 
	ioniz_cloud_data.clear();

	el_spectrum_init_data.clear();
	specimen_conc_init_data.clear();
	h2_pop_dens_init_data.clear();

	el_spectrum_evol.clear();
	specimen_conc_evol.clear(); 
	h2_popdens_evol.clear(); 
	h2_popdens_v_evol.clear(); 
	hei_popul_dens_evol.clear();
}


// Initialization of the grid for electron energies,
// there IS NO interval for electrons with high energy (E > max_el_energy)
// there IS interval for electrons with low energy [0, min energy], min_energy must be < any excitation process threshold
void init_electron_energy_grid(double max_electron_energy)
{
	int i, k, nb_bins;
	double x, y;

	electron_energies_grid.clear();
	electron_energy_bin_size.clear();

	x = pow(10., 1. / NB_OF_BINS_PER_ORDER_EL);
	nb_bins = (int)(1. / (x - 1.)) + 1;

	// equal energy intervals in energy range [0, E_fixed],
	y = ELECTRON_ENERGY_FIXED / nb_bins;
	for (i = 0; i < nb_bins; i++) {
		electron_energies_grid.push_back(i * y);
	}
	electron_energies_grid.push_back(ELECTRON_ENERGY_FIXED);

	// two extra intervals are added at the end,
	k = (int)(NB_OF_BINS_PER_ORDER_EL * log10((double)max_electron_energy / ELECTRON_ENERGY_FIXED)) + 2;
	for (i = 0; i < k; i++) {
		electron_energies_grid.push_back(electron_energies_grid.back() * x);
	}

	for (i = 0; i < (int)electron_energies_grid.size() - 1; i++) {
		y = electron_energies_grid[i + 1] - electron_energies_grid[i];
		electron_energy_bin_size.push_back(y);
	}
	nb_of_el_energies = (int)electron_energies_grid.size() - 1;
}


void init_time_intervals()
{
	int i;
	double x;

	// Initialization of time grid,
	time_cloud_arr.clear();
	time_cloud_arr.push_back(0.);
	time_cloud_arr.push_back(MIN_MODEL_TIME);

	x = pow(10., 1. / NB_OF_BINS_PER_ORDER_TIME);

	for (i = 0; time_cloud_arr.back() < MAX_MODEL_TIME; i++) {
		time_cloud_arr.push_back(time_cloud_arr.back() * x);
	}
	nb_of_time_moments = (int)time_cloud_arr.size();
}


void simulations_monoenergetic(const string& data_path, const string& output_path, double conc_h_tot, double he_abund, double ioniz_fract, 
	double op_ratio_h2, double electron_conc, double electron_energy)
{
	bool must_be_saved;
	int i, k, time_nb, step_nb, flag;
	long int nb_steps_tot;

	double x1, x2, x3, x4, x5, x6, x7, rel_tol, dt, model_time, model_time_aux, model_time_aux_prev, model_time_step, model_time_out;

	stringstream ss;
	string output_path_2;

	stringstream linux_out_fname;
	ofstream linux_out;
	streambuf* orig_cout;

	time_t timer;
	char* ctime_str;

	// the class containing data used by solver,
	elspectra_evolution_data* user_data;

	// directory name for saving the data on the specific electron energy,
	i = (int) log10(conc_h_tot);
	k = rounding(conc_h_tot / pow(10., i));
	
	ss.clear();
	ss.str("");
	ss << output_path << "nh_" << k << "e" << i << "_";

	i = rounding(electron_energy);
	ss << i << "ev" << "/";

	output_path_2 = ss.str();

#ifdef __linux__
	linux_out_fname.str("");
	linux_out_fname << output_path_2;
	linux_out_fname << "out";
	linux_out_fname << "_screen";
	linux_out_fname << ".txt";
	// linux_out_fname.str("/dev/null");

	//ofstream outerr(lin_out.str().c_str(), ios::app);
	//streambuf* orig_cerr = cerr.rdbuf(); // original cerr;
	//cerr.rdbuf(outerr.rdbuf());

	linux_out.open(linux_out_fname.str().c_str(), ios::app);
	orig_cout = cout.rdbuf();
	cout.rdbuf(linux_out.rdbuf());
#endif

	timer = time(NULL);
	ctime_str = ctime(&timer);

	cout << scientific;
	cout.precision(3);

	cout << ctime_str << endl
		<< "Start of the simulations of the evolution of electron spectra" << endl
		<< "Initialization of the data..." << endl;

	memory_freeing_up();
	init_time_intervals();
	init_electron_energy_grid(electron_energy);

	dynamic_array el_spectrum(nb_of_el_energies);
	init_electron_spectra_mono(el_spectrum, electron_energies_grid, electron_conc, electron_energy);

	user_data
		= new elspectra_evolution_data(data_path, output_path, conc_h_tot, ioniz_fract, electron_energies_grid, VERBOSITY);

	nb_of_equat = user_data->get_nb_of_equat();
	nb_lev_h2 = user_data->get_nb_of_h2_lev();
	nb_lev_hei = user_data->get_nb_of_hei_lev();

	h2eq_nb = user_data->get_h2eq_nb();
	heieq_nb = user_data->get_heieq_nb();
	physeq_nb = user_data->get_physeq_nb();

	initial_data = new double[nb_of_equat];
	memset(initial_data, 0, nb_of_equat * sizeof(double));

	for (i = 0; i < nb_of_el_energies; i++) {
		// the number of electrons in the interval:
		initial_data[i] = el_spectrum.arr[i] * user_data->get_electron_energy_bin(i);
	}

	initial_data[nb_of_el_energies + H2_NB] = 0.5 * conc_h_tot;
	initial_data[nb_of_el_energies + HE_NB] = he_abund * conc_h_tot;

	// the initial distribution of energy level populations is postulated, [cm-3]
	initial_data[h2eq_nb] = initial_data[nb_of_el_energies + H2_NB] / (op_ratio_h2 + 1.);
	initial_data[h2eq_nb + 1] = initial_data[nb_of_el_energies + H2_NB] * op_ratio_h2 / (op_ratio_h2 + 1.);

	// the initial distribution of HeI levels,
	initial_data[heieq_nb] = initial_data[nb_of_el_energies + HE_NB];

	// Initial gas/dust temperature, in K
	initial_data[physeq_nb] = initial_data[physeq_nb + 1] = 10.;

	// Dust charge
	initial_data[physeq_nb + 2] = 0.;

	dust_is_presented = false;
	grain_radius = 0.;
	grain_nb_density = 0.;
	dg_ratio = 0.;

	user_data->set_dust_parameters(dust_is_presented, grain_radius, grain_nb_density);
	save_model_parameters(output_path_2, conc_h_tot, op_ratio_h2, ioniz_fract, dg_ratio, grain_radius, grain_nb_density);

	//
	// Simulations
	//

	model_time_aux = model_time = 0.;

	// Auxiliary arrays, 
	// the arrays are initialized by zeros in the class constructor,
	dynamic_array specimen_conc(NB_OF_CHEM_SPECIES), h2_popul_dens(nb_lev_h2), hei_popul_dens(nb_lev_hei), h2_pop_v(MAX_NB_H2_VSTATES_X1SG);

	dynamic_array h2_electr_vstates(MAX_NB_H2_VSTATES_X1SG), h2_electr_vstates_rates(MAX_NB_H2_VSTATES_X1SG),
		h2_electr_vstates_aux(MAX_NB_H2_VSTATES_X1SG);

	dynamic_array h2_vibr_vstates(MAX_NB_H2_VSTATES_X1SG), h2_vibr_vstates_rates(MAX_NB_H2_VSTATES_X1SG),
		h2_vibr_vstates_aux(MAX_NB_H2_VSTATES_X1SG);

	// vectors used by solver SUNDIALS CVODE
	N_Vector y, ydot, abs_tol;

	SUNMatrix A;
	SUNLinearSolver LS;
	void* cvode_mem;

	y = N_VNew_Serial(nb_of_equat);
	ydot = N_VNew_Serial(nb_of_equat);
	abs_tol = N_VNew_Serial(nb_of_equat);

	// Initialization of initial values of variables,
	for (i = 0; i < nb_of_equat; i++) {
		NV_Ith_S(y, i) = initial_data[i];
	}

	// Initialization for tolerances:
	rel_tol = REL_ERROR_SOLVER;
	user_data->set_tolerances(abs_tol);

	// Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula and the use of a Newton iteration 
	cvode_mem = CVodeCreate(CV_BDF);

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
	//LS = SUNDenseLinearSolver(y, A);
	LS = SUNLinSol_Dense(y, A);

	// Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode 
	//flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
	flag = CVodeSetLinearSolver(cvode_mem, LS, A);

	// The function attaches the user data block to the solver;
	flag = CVodeSetUserData(cvode_mem, user_data);

	// restart of the solver with new values of initial conditions,
	// flag = CVodeReInit(cvode_mem, model_time, y);

	//
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

	memset(h2_pop_v.arr, 0, MAX_NB_H2_VSTATES_X1SG * sizeof(double));
	for (i = 0; i < nb_lev_h2; i++) {
		k = user_data->get_vibr_nb_h2(i);
		if (k >= 0) {
			h2_pop_v.arr[k] += NV_Ith_S(y, h2eq_nb + i);
		}
	}
	h2_popdens_v_evol.push_back(h2_pop_v);

	for (i = 0; i < nb_lev_hei; i++) {
		hei_popul_dens.arr[i] = NV_Ith_S(y, heieq_nb + i);
	}
	hei_popul_dens_evol.push_back(hei_popul_dens);

	//
	// update the members of user_data class,
	f_elsp(model_time, y, ydot, (void*)user_data);

	user_data->get_el_energy_losses(enloss_rate_mt, enloss_rate_h2_rot, enloss_rate_h2_rot_pos, enloss_rate_h2_vibr, enloss_rate_h2_vibr_pos,
		enloss_rate_ioniz, enloss_rate_hei, enloss_rate_coulomb_el, neutral_coll_heating_rate);

	user_data->get_h2_state_data(h2_state_data_rates);

	user_data->get_h2_excitation_rates(h2_electr_vstates_rates, h2_vibr_vstates_rates, h2_excit_rot_rate, hei_excit_rate);

	// Energy loss rates
	enloss_rate_mt_arr.push_back(enloss_rate_mt);
	enloss_rate_h2_rot_arr.push_back(enloss_rate_h2_rot);
	enloss_rate_h2_vibr_arr.push_back(enloss_rate_h2_vibr);
	enloss_rate_ioniz_arr.push_back(enloss_rate_ioniz);
	enloss_rate_hei_arr.push_back(enloss_rate_hei);

	enloss_rate_coulomb_el_arr.push_back(enloss_rate_coulomb_el);
	neutral_coll_heating_rate_arr.push_back(neutral_coll_heating_rate);

	enloss_rate_h2_rot_arr_pos.push_back(enloss_rate_h2_rot_pos);
	enloss_rate_h2_vibr_arr_pos.push_back(enloss_rate_h2_vibr_pos);

	// Energy losses, time-integrated
	// Please, check that all variables at this step are set to zero,
	enloss_mt_arr.push_back(enloss_mt);
	enloss_h2_rot_arr.push_back(enloss_h2_rot);
	enloss_h2_vibr_arr.push_back(enloss_h2_vibr);
	enloss_ioniz_arr.push_back(enloss_ioniz);
	enloss_hei_arr.push_back(enloss_hei);

	enloss_coulomb_el_arr.push_back(enloss_coulomb_el);
	neutral_coll_heating_arr.push_back(neutral_coll_heating);

	// Excitation,
	h2_excit_rot_rate_arr.push_back(h2_excit_rot_rate);
	h2_excit_rot_arr.push_back(h2_excit_rot);
	hei_excit_arr.push_back(hei_excit);

	// H2 data
	h2_state_data_rate_arr.push_back(h2_state_data_rates);
	h2_electr_vstates_rate_arr.push_back(h2_electr_vstates_rates);
	h2_vibr_vstates_rate_arr.push_back(h2_vibr_vstates_rates);

	h2_state_data_arr.push_back(h2_state_data);
	h2_electr_vstates_arr.push_back(h2_electr_vstates);
	h2_vibr_vstates_arr.push_back(h2_vibr_vstates);
	
	// Physical parameters
	neutral_temp_arr.push_back(NV_Ith_S(y, physeq_nb));
	ion_temp_arr.push_back(NV_Ith_S(y, physeq_nb + 1));
	dust_charge_arr.push_back(NV_Ith_S(y, physeq_nb + 2) * grain_nb_density);  // dust charge per cm3

	// Start of the simulations: 
	cout << left << "Initialization step time (s): " << setw(12) << (int)(time(NULL) - timer) << endl
		<< "Starting simulations..." << endl;
	timer = time(NULL);

	if (VERBOSITY) {
		cout << left << setw(12) << "model_time " << setw(9) << "steps_nb "
			<< setw(12) << "gas_temp(K)" << setw(12) << "ion_temp(K)" << setw(12) << "dust_charge" << endl;
	}

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
			f_elsp(model_time_aux, y, ydot, user_data);

			// Note the factor 0.5,
			dt = 0.5 * (model_time_aux - model_time_aux_prev);

			// Energy losses, 
			x1 = enloss_rate_h2_rot;
			x2 = enloss_rate_h2_vibr;
			x3 = enloss_rate_ioniz;
			x4 = enloss_rate_hei;
			
			x5 = enloss_rate_mt;
			x6 = enloss_rate_coulomb_el;
			x7 = neutral_coll_heating_rate;

			user_data->get_el_energy_losses(enloss_rate_mt, enloss_rate_h2_rot, enloss_rate_h2_rot_pos, enloss_rate_h2_vibr, enloss_rate_h2_vibr_pos,
				enloss_rate_ioniz, enloss_rate_hei, enloss_rate_coulomb_el, neutral_coll_heating_rate);

			enloss_h2_rot += (x1 + enloss_rate_h2_rot) * dt;
			enloss_h2_vibr += (x2 + enloss_rate_h2_vibr) * dt;
			enloss_ioniz += (x3 + enloss_rate_ioniz) * dt;
			enloss_hei += (x4 + enloss_rate_hei) * dt;
		
			// contribute to thermal heating
			enloss_mt += (x5 + enloss_rate_mt) * dt;
			enloss_coulomb_el += (x6 + enloss_rate_coulomb_el) * dt;
			neutral_coll_heating += (x7 + neutral_coll_heating_rate) * dt;

			// H2 state data
			for (i = 0; i < NB_EXC_ELECTRONIC_STATES; i++) {
				h2_state_data_aux[i] = h2_state_data_rates[i];
			}
			
			user_data->get_h2_state_data(h2_state_data_rates);

			for (i = 0; i < NB_EXC_ELECTRONIC_STATES; i++) {
				h2_state_data[i].excit += (h2_state_data_aux[i].excit + h2_state_data_rates[i].excit) * dt;
				h2_state_data[i].direct_diss += (h2_state_data_aux[i].direct_diss + h2_state_data_rates[i].direct_diss) * dt;
				h2_state_data[i].twostep_diss += (h2_state_data_aux[i].twostep_diss + h2_state_data_rates[i].twostep_diss) * dt;
				h2_state_data[i].diss_heat_input += (h2_state_data_aux[i].diss_heat_input + h2_state_data_rates[i].diss_heat_input) * dt;
				h2_state_data[i].enloss_excit += (h2_state_data_aux[i].enloss_excit + h2_state_data_rates[i].enloss_excit) * dt;
				h2_state_data[i].enloss_direct_diss += (h2_state_data_aux[i].enloss_direct_diss + h2_state_data_rates[i].enloss_direct_diss) * dt;
			}

			// Excitation rates,
			for (i = 0; i < MAX_NB_H2_VSTATES_X1SG; i++) {
				h2_electr_vstates_aux.arr[i] = h2_electr_vstates_rates.arr[i];
				h2_vibr_vstates_aux.arr[i] = h2_vibr_vstates_rates.arr[i];

			}
			x1 = h2_excit_rot_rate;
			x2 = hei_excit_rate;

			user_data->get_h2_excitation_rates(h2_electr_vstates_rates, h2_vibr_vstates_rates, h2_excit_rot_rate, hei_excit_rate);
			
			for (i = 0; i < MAX_NB_H2_VSTATES_X1SG; i++) {
				h2_electr_vstates.arr[i] += (h2_electr_vstates_aux.arr[i] + h2_electr_vstates_rates.arr[i]) *dt;
				h2_vibr_vstates.arr[i] += (h2_vibr_vstates_aux.arr[i] + h2_vibr_vstates_rates.arr[i]) *dt;
			}
			h2_excit_rot += (x1 + h2_excit_rot_rate) * dt;
			hei_excit += (x2 + hei_excit_rate) * dt;

			flag = CVodeGetNumSteps(cvode_mem, &nb_steps_tot);

			if (VERBOSITY) {
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

				memset(h2_pop_v.arr, 0, MAX_NB_H2_VSTATES_X1SG * sizeof(double));
				for (i = 0; i < nb_lev_h2; i++) {
					k = user_data->get_vibr_nb_h2(i);
					if (k >= 0) {
						h2_pop_v.arr[k] += NV_Ith_S(y, h2eq_nb + i);
					}
				}
				h2_popdens_v_evol.push_back(h2_pop_v);

				for (i = 0; i < nb_lev_hei; i++) {
					hei_popul_dens.arr[i] = NV_Ith_S(y, heieq_nb + i);
				}
				hei_popul_dens_evol.push_back(hei_popul_dens);

				// Energy losses, rates
				enloss_rate_h2_rot_arr.push_back(enloss_rate_h2_rot);
				enloss_rate_h2_vibr_arr.push_back(enloss_rate_h2_vibr);
				enloss_rate_ioniz_arr.push_back(enloss_rate_ioniz);
				enloss_rate_hei_arr.push_back(enloss_rate_hei);

				enloss_rate_mt_arr.push_back(enloss_rate_mt);
				enloss_rate_coulomb_el_arr.push_back(enloss_rate_coulomb_el);
				neutral_coll_heating_rate_arr.push_back(neutral_coll_heating_rate);

				enloss_rate_h2_rot_arr_pos.push_back(enloss_rate_h2_rot_pos);
				enloss_rate_h2_vibr_arr_pos.push_back(enloss_rate_h2_vibr_pos);

				// Energy losses, time-integrated,
				enloss_h2_rot_arr.push_back(enloss_h2_rot);
				enloss_h2_vibr_arr.push_back(enloss_h2_vibr);
				enloss_ioniz_arr.push_back(enloss_ioniz);
				enloss_hei_arr.push_back(enloss_hei);

				enloss_mt_arr.push_back(enloss_mt);
				enloss_coulomb_el_arr.push_back(enloss_coulomb_el);
				neutral_coll_heating_arr.push_back(neutral_coll_heating);

				// H2 data
				h2_state_data_rate_arr.push_back(h2_state_data_rates);
				h2_electr_vstates_rate_arr.push_back(h2_electr_vstates_rates);
				h2_vibr_vstates_rate_arr.push_back(h2_vibr_vstates_rates);

				h2_state_data_arr.push_back(h2_state_data);
				h2_electr_vstates_arr.push_back(h2_electr_vstates);
				h2_vibr_vstates_arr.push_back(h2_vibr_vstates);

				// Excitations, rates and time-integrated
				h2_excit_rot_rate_arr.push_back(h2_excit_rot_rate);
				h2_excit_rot_arr.push_back(h2_excit_rot);
				hei_excit_arr.push_back(hei_excit);

				neutral_temp_arr.push_back(NV_Ith_S(y, physeq_nb));
				ion_temp_arr.push_back(NV_Ith_S(y, physeq_nb + 1));
				dust_charge_arr.push_back(NV_Ith_S(y, physeq_nb + 2) * grain_nb_density);  // dust charge per cm3

				save_electron_spectrum_evolution(output_path_2, electron_energies_grid, electron_energy_bin_size, time_cloud_arr, el_spectrum_evol, conc_h_tot);

				save_specimen_conc_evolution(output_path_2, time_cloud_arr, specimen_conc_evol, conc_h_tot);

				save_h2_populations_evolution(output_path_2, time_cloud_arr, h2_popdens_evol, h2_popdens_v_evol, conc_h_tot, nb_lev_h2);

				save_hei_populations_evolution(output_path_2, time_cloud_arr, hei_popul_dens_evol, conc_h_tot, nb_lev_hei);

				save_electron_energy_loss_rates(output_path_2, conc_h_tot, time_cloud_arr,
					enloss_rate_mt_arr,
					enloss_rate_h2_rot_arr,
					enloss_rate_h2_rot_arr_pos,
					enloss_rate_h2_vibr_arr,
					enloss_rate_h2_vibr_arr_pos,
					enloss_rate_ioniz_arr,
					enloss_rate_hei_arr,
					enloss_rate_coulomb_el_arr,
					neutral_coll_heating_rate_arr, 
					h2_state_data_rate_arr);

				save_electron_energy_losses(output_path_2, conc_h_tot, time_cloud_arr,
					enloss_mt_arr,
					enloss_h2_rot_arr,
					enloss_h2_vibr_arr,
					enloss_ioniz_arr,
					enloss_hei_arr,
					enloss_coulomb_el_arr,
					neutral_coll_heating_arr, 
					h2_state_data_arr);

				// Saving H2 electronic excitation data
				save_electronic_states_excit_rates(output_path_2, conc_h_tot, time_cloud_arr, h2_state_data_rate_arr);

				// Saving H2 vibrational excitation data
				save_vibrational_states_excit_rates(output_path_2, conc_h_tot, time_cloud_arr, h2_electr_vstates_rate_arr, h2_vibr_vstates_rate_arr,
					h2_excit_rot_rate_arr);
			}
		}
		else {
			cout << "Some unexpected error has been occurred: " << flag << endl;
			break;
		}
	}
	cout << left << "Total calculation time (s): " << setw(12) << (int)(time(NULL) - timer) << endl;

	N_VDestroy(y);
	N_VDestroy(ydot);
	N_VDestroy(abs_tol);

	// Free integrator memory
	CVodeFree(&cvode_mem);

	// Free the linear solver memory
	SUNLinSolFree(LS);

	// Free the matrix memory
	SUNMatDestroy(A);

	save_output_parameters(output_path, conc_h_tot, electron_energies_grid, electron_energy_bin_size,
		el_spectrum_evol,
		specimen_conc_evol,
		enloss_mt_arr,
		enloss_coulomb_el_arr,
		neutral_coll_heating_arr,
		enloss_h2_vibr_arr,
		h2_state_data_arr,
		h2_electr_vstates_arr,
		h2_vibr_vstates_arr,
		h2_excit_rot_arr,
		hei_excit_arr);
}

/*
void simulations_grb(const string& parameter_fname)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	
	bool must_be_saved;
	int i, k, lay_nb, is_h2_pop_dens_init, time_nb, step_nb, flag;
	long int nb_steps_tot;
	
	double hcolumn_density, grb_distance, grb_cloud_distance;
	double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, rel_tol, dt, model_time, model_time_aux, model_time_aux_prev, model_time_step, model_time_out;

	stringstream ss;
	ifstream input;

	vector<int>    qnb_v_arr, qnb_j_arr;
	vector<double> layer_centre_distances;  // in cm

	// the class containing data used by solver,
	elspectra_evolution_data* user_data;

	//
	// Initialization of initial data
	//
	input.open(parameter_fname.c_str());
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with cross section data " << parameter_fname << endl;
		exit(1);
	}

	while (!input.eof())
	{
		do // comment lines are read:
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		ss.clear();
		ss.str(text_line);
		ss >> nb_processors;

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

		// path may have spaces,
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

		input.close();
	}

#ifdef __linux__
	linux_out_fname.str("");
	linux_out_fname << output_path;
	linux_out_fname << "out";
	linux_out_fname << "_screen";
	linux_out_fname << ".txt";
	// linux_out_fname.str("/dev/null");

	//ofstream outerr(lin_out.str().c_str(), ios::app);
	//streambuf* orig_cerr = cerr.rdbuf(); // original cerr;
	//cerr.rdbuf(outerr.rdbuf());

	linux_out.open(linux_out_fname.str().c_str(), ios::app);
	orig_cout = cout.rdbuf();
	cout.rdbuf(linux_out.rdbuf());
#endif

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

	timer = time(NULL);
	ctime_str = ctime(&timer);

	cout << ctime_str << endl
		<< "Start of the simulations of the evolution of electron spectra" << endl
		<< "Initialization of the data..." << endl;

	memory_freeing_up();
	init_time_intervals();

	ioniz_fract = IONIZATION_FRACTION_THERMAL_EL;
	init_electron_energy_grid(MAX_ELECTRON_ENERGY);

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

	user_data 
		= new elspectra_evolution_data (data_path, output_path, conc_h_tot, ioniz_fract, electron_energies_grid, VERBOSITY);

	nb_of_equat = user_data->get_nb_of_equat();
	nb_lev_h2 = user_data->get_nb_of_h2_lev();
	nb_lev_hei = user_data->get_nb_of_hei_lev();

	h2eq_nb = user_data->get_h2eq_nb();
	heieq_nb = user_data->get_heieq_nb();
	physeq_nb = user_data->get_physeq_nb();

	initial_data = new double[nb_of_equat];
	memset(initial_data, 0, nb_of_equat * sizeof(double));

	for (i = 0; i < nb_of_el_energies; i++) {
		// the number of electrons in the interval:
		initial_data[i] = el_spectrum_init_data[lay_nb].arr[i] * user_data->get_electron_energy_bin(i);
	}

	for (i = 0; i < NB_OF_CHEM_SPECIES; i++) {
		initial_data[nb_of_el_energies + i] = specimen_conc_init_data[lay_nb].arr[i];
	}

	if (is_h2_pop_dens_init) {
		// population densities must be in [cm-3],
		x1 = x2 = 0.;
		for (i = 0; i < h2_pop_dens_init_data[lay_nb].dim; i++)
		{
			k = user_data->get_level_nb_h2(qnb_v_arr[i], qnb_j_arr[i]);
			if (k >= 0) {
				initial_data[h2eq_nb + k] = h2_pop_dens_init_data[lay_nb].arr[i];
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
		initial_data[h2eq_nb] = initial_data[nb_of_el_energies + H2_NB] / (op_ratio_h2 + 1.);
		initial_data[h2eq_nb + 1] = initial_data[nb_of_el_energies + H2_NB] * op_ratio_h2 / (op_ratio_h2 + 1.);
	}

	// the initial distribution of HeI levels,
	initial_data[heieq_nb] = initial_data[nb_of_el_energies + HE_NB];

	// Initial gas/dust temperature, in K
	initial_data[physeq_nb] = initial_data[physeq_nb + 1] = 10.;

	// Dust charge
	initial_data[physeq_nb + 2] = 0.;

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
	user_data->set_dust_parameters(dust_is_presented, grain_radius, grain_nb_density);

	save_model_parameters(output_path, conc_h_tot, op_ratio_h2, ioniz_fract, dg_ratio, grain_radius, grain_nb_density,
		grb_cloud_distance, grb_distance, hcolumn_density, lay_nb);

	//
	// Simulations
	//

	dynamic_array el_spectrum(nb_of_el_energies), specimen_conc(NB_OF_CHEM_SPECIES), h2_popul_dens(nb_lev_h2), hei_popul_dens(nb_lev_hei),
		h2_pop_v(MAX_NB_H2_VSTATES_X1SG);

	// vectors used by solver SUNDIALS CVODE
	N_Vector y, ydot, abs_tol;

	cout << scientific;
	cout.precision(3);

	//elspectra_evolution_data user_data_(data_path, output_path, conc_h_tot, ioniz_fract, electron_energies_grid, VERBOSITY);
	//user_data_.set_dust_parameters(dust_is_presented, grain_radius, grain_nb_density);

	model_time_aux = model_time = 0.;

	SUNMatrix A;
	SUNLinearSolver LS;
	void* cvode_mem;

	y = N_VNew_Serial(nb_of_equat);
	ydot = N_VNew_Serial(nb_of_equat);
	abs_tol = N_VNew_Serial(nb_of_equat);

	// Initialization of initial values of variables,
	for (i = 0; i < nb_of_equat; i++) {
		NV_Ith_S(y, i) = initial_data[i];
	}

	// Initialization for tolerances:
	rel_tol = REL_ERROR_SOLVER;
	user_data->set_tolerances(abs_tol);

	// Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula and the use of a Newton iteration 
	cvode_mem = CVodeCreate(CV_BDF);

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
	//LS = SUNDenseLinearSolver(y, A);
	LS = SUNLinSol_Dense(y, A);

	// Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode 
	//flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
	flag = CVodeSetLinearSolver(cvode_mem, LS, A);

	// The function attaches the user data block to the solver;
	flag = CVodeSetUserData(cvode_mem, user_data);

	// restart of the solver with new values of initial conditions,
	// flag = CVodeReInit(cvode_mem, model_time, y);

	// update the members of user_data class,
	//f_elsp(model_time, y, ydot,  (void*) user_data);

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

	memset(h2_pop_v.arr, 0, MAX_NB_H2_VSTATES_X1SG * sizeof(double));
	for (i = 0; i < nb_lev_h2; i++) {
		k = user_data->get_vibr_nb_h2(i);
		if (k >= 0)
			h2_pop_v.arr[k] += NV_Ith_S(y, h2eq_nb + i);
	}
	h2_popdens_v_evol.push_back(h2_pop_v);

	for (i = 0; i < nb_lev_hei; i++) {
		hei_popul_dens.arr[i] = NV_Ith_S(y, heieq_nb + i);
	}
	hei_popul_dens_evol.push_back(hei_popul_dens);

	// Energy loss rates
	// Please, check that all variables are set to zero,
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

	// Excitation rates,
	h2_excit_electr_rate_arr.push_back(h2_excit_electr_rate);
	h2_excit_vibr_rate_arr.push_back(h2_excit_vibr_rate);
	h2_excit_rot_rate_arr.push_back(h2_excit_rot_rate);

	// Excitation, time-integrated
	h2_excit_electr_arr.push_back(h2_excit_electr);
	h2_excit_vibr_arr.push_back(h2_excit_vibr);
	h2_excit_rot_arr.push_back(h2_excit_rot);

	h2_excit_electr_bs_arr.push_back(h2_excit_electr_bs);
	h2_excit_electr_bps_arr.push_back(h2_excit_electr_bps);
	h2_excit_electr_cp_arr.push_back(h2_excit_electr_cp);
	h2_excit_electr_dp_arr.push_back(h2_excit_electr_dp);
	
	h2_excit_vibr_v1_arr.push_back(h2_excit_vibr_v1);
	h2_excit_vibr_v2_arr.push_back(h2_excit_vibr_v2);
	h2_excit_vibr_v3_arr.push_back(h2_excit_vibr_v3);
	h2_excit_el_vstate_arr.push_back(h2_el_excit_vstate_rate);

	hei_excit_arr.push_back(hei_exc);

	// Dissociation, time-integrated
	h2_solomon_diss_arr.push_back(h2_solomon_diss);
	h2_diss_exc_singlet_arr.push_back(h2_diss_exc_singlet);
	h2_diss_excit_electr_triplet_arr.push_back(h2_diss_excit_electr_triplet);
	h2_excit_electr_triplet_arr.push_back(h2_excit_electr_triplet);

	neutral_temp_arr.push_back(NV_Ith_S(y, physeq_nb));
	ion_temp_arr.push_back(NV_Ith_S(y, physeq_nb + 1));
	dust_charge_arr.push_back(NV_Ith_S(y, physeq_nb + 2) * grain_nb_density);  // dust charge per cm3

	// Start of the simulations: 
	cout << left << "Initialization step time (s): " << setw(12) << (int)(time(NULL) - timer) << endl
		<< "Starting simulations..." << endl;
	timer = time(NULL);

	if (VERBOSITY) {
		cout << left << setw(12) << "model_time " << setw(9) << "steps_nb "
			<< setw(12) << "gas_temp(K)" << setw(12) << "ion_temp(K)" << setw(12) << "dust_charge" << endl;
	}

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
			f_elsp(model_time_aux, y, ydot, user_data);

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

			user_data->get_el_energy_losses(enloss_rate_mt, enloss_rate_h2_rot, enloss_rate_h2_rot_pos, enloss_rate_h2_vibr, enloss_rate_h2_vibr_pos,
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

			// Dissociation rates
			x1 = h2_solomon_diss_rate;
			x2 = h2_diss_exc_singlet_rate;
			x3 = h2_diss_excit_triplet_rate;
			x4 = h2_excit_triplet_rate;

			user_data->get_h2_diss_rates(h2_solomon_diss_rate, h2_diss_exc_singlet_rate, h2_diss_excit_triplet_rate, h2_excit_triplet_rate);

			h2_solomon_diss += (x1 + h2_solomon_diss_rate) * dt;
			h2_diss_exc_singlet += (x2 + h2_diss_exc_singlet_rate) * dt;
			h2_diss_excit_electr_triplet += (x3 + h2_diss_excit_triplet_rate) * dt;
			h2_excit_electr_triplet += (x4 + h2_excit_triplet_rate) * dt;

			// Excitation rates,
			x1 = h2_excit_electr_rate;
			x2 = h2_excit_vibr_rate;
			x3 = h2_excit_rot_rate;

			x4 = h2_excit_electr_bs_rate;
			x5 = h2_excit_electr_bps_rate;
			x6 = h2_excit_electr_cp_rate;
			x7 = h2_excit_electr_dp_rate;
			
			x8 = h2_excit_vibr_v1_rate;
			x9 = h2_excit_vibr_v2_rate;
			x10 = h2_excit_vibr_v3_rate;
			x11 = h2_el_excit_vstate_rate;

			x12 = hei_exc_rate;

			user_data->get_h2_process_rates(h2_excit_electr_rate, h2_excit_vibr_rate, h2_excit_rot_rate,
				h2_excit_electr_bs_rate, h2_excit_electr_bps_rate, h2_excit_electr_cp_rate, h2_excit_electr_dp_rate, 
				h2_excit_vibr_v1_rate, h2_excit_vibr_v2_rate, h2_excit_vibr_v3_rate, h2_el_excit_vstate_rate, hei_exc_rate);

			h2_excit_electr += (x1 + h2_excit_electr_rate) * dt;
			h2_excit_vibr += (x2 + h2_excit_vibr_rate) * dt;
			h2_excit_rot += (x3 + h2_excit_rot_rate) * dt;

			h2_excit_electr_bs += (x4 + h2_excit_electr_bs_rate) * dt;
			h2_excit_electr_bps += (x5 + h2_excit_electr_bps_rate) * dt;
			h2_excit_electr_cp += (x6 + h2_excit_electr_cp_rate) * dt;
			h2_excit_electr_dp += (x7 + h2_excit_electr_dp_rate) * dt;
			
			h2_excit_vibr_v1 += (x8 + h2_excit_vibr_v1_rate) * dt;
			h2_excit_vibr_v2 += (x9 + h2_excit_vibr_v2_rate) * dt;
			h2_excit_vibr_v3 += (x10 + h2_excit_vibr_v3_rate) * dt;
			h2_el_excit_vstate += (x11 + h2_el_excit_vstate_rate) * dt;

			hei_exc += (x12 + hei_exc_rate) * dt;

			flag = CVodeGetNumSteps(cvode_mem, &nb_steps_tot);

			if (VERBOSITY) {
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

				memset(h2_pop_v.arr, 0, MAX_NB_H2_VSTATES_X1SG * sizeof(double));
				for (i = 0; i < nb_lev_h2; i++) {
					k = user_data->get_vibr_nb_h2(i);
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

				// Excitations, rates
				h2_excit_electr_rate_arr.push_back(h2_excit_electr_rate);
				h2_excit_vibr_rate_arr.push_back(h2_excit_vibr_rate);
				h2_excit_rot_rate_arr.push_back(h2_excit_rot_rate);

				// Excitations, time-integrated,
				h2_excit_electr_arr.push_back(h2_excit_electr);
				h2_excit_vibr_arr.push_back(h2_excit_vibr);
				h2_excit_rot_arr.push_back(h2_excit_rot);

				h2_excit_electr_bs_arr.push_back(h2_excit_electr_bs);
				h2_excit_electr_bps_arr.push_back(h2_excit_electr_bps);
				h2_excit_electr_cp_arr.push_back(h2_excit_electr_cp);
				h2_excit_electr_dp_arr.push_back(h2_excit_electr_dp);
				
				h2_excit_vibr_v1_arr.push_back(h2_excit_vibr_v1);
				h2_excit_vibr_v2_arr.push_back(h2_excit_vibr_v2);
				h2_excit_vibr_v3_arr.push_back(h2_excit_vibr_v3);
				h2_excit_el_vstate_arr.push_back(h2_el_excit_vstate);

				hei_excit_arr.push_back(hei_exc);

				// Dissociations, time-integrated,
				h2_solomon_diss_arr.push_back(h2_solomon_diss);
				h2_diss_exc_singlet_arr.push_back(h2_diss_exc_singlet);
				h2_diss_excit_electr_triplet_arr.push_back(h2_diss_excit_electr_triplet);
				h2_excit_electr_triplet_arr.push_back(h2_excit_electr_triplet);

				neutral_temp_arr.push_back(NV_Ith_S(y, physeq_nb));
				ion_temp_arr.push_back(NV_Ith_S(y, physeq_nb + 1));
				dust_charge_arr.push_back(NV_Ith_S(y, physeq_nb + 2) * grain_nb_density);  // dust charge per cm3

				save_electron_spectrum_evolution(output_path, electron_energies_grid, electron_energy_bin_size, time_cloud_arr, el_spectrum_evol, conc_h_tot);

				save_specimen_conc_evolution(output_path, time_cloud_arr, specimen_conc_evol, conc_h_tot);

				save_h2_populations_evolution(output_path, time_cloud_arr, h2_popdens_evol, h2_popdens_v_evol, conc_h_tot, nb_lev_h2);

				save_hei_populations_evolution(output_path, time_cloud_arr, hei_popul_dens_evol, conc_h_tot, nb_lev_hei);

				save_electron_energy_loss_rates(output_path, conc_h_tot, time_cloud_arr,
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

				save_electron_energy_losses(output_path, conc_h_tot, time_cloud_arr,
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

				save_excit_rates(output_path, conc_h_tot, time_cloud_arr, h2_excit_electr_rate_arr, h2_excit_vibr_rate_arr, h2_excit_rot_rate_arr);

				save_excit(output_path, conc_h_tot, time_cloud_arr,
					h2_excit_electr_arr,
					h2_excit_vibr_arr,
					h2_excit_rot_arr,
					h2_excit_electr_bs_arr,
					h2_excit_electr_bps_arr,
					h2_excit_electr_cp_arr,
					h2_excit_electr_dp_arr,
					h2_excit_vibr_v1_arr,
					h2_excit_vibr_v2_arr,
					h2_excit_vibr_v3_arr,
					h2_excit_el_vstate_arr,
					hei_excit_arr);

				save_diss(output_path, conc_h_tot, time_cloud_arr,
					h2_solomon_diss_arr,
					h2_diss_exc_singlet_arr,
					h2_diss_excit_electr_triplet_arr,
					h2_excit_electr_triplet_arr);

				save_phys_parameters(output_path, time_cloud_arr, neutral_temp_arr, ion_temp_arr, dust_charge_arr, specimen_conc_evol, conc_h_tot);
			}
		}
		else {
			cout << "Some unexpected error has been occurred: " << flag << endl;
			break;
		}
	}
	cout << left << "Total calculation time (s): " << setw(12) << (int)(time(NULL) - timer) << endl;

	N_VDestroy(y);
	N_VDestroy(ydot);
	N_VDestroy(abs_tol);

	// Free integrator memory
	CVodeFree(&cvode_mem);

	// Free the linear solver memory
	SUNLinSolFree(LS);

	// Free the matrix memory
	SUNMatDestroy(A);

	save_simulation_report();
}*/


//
// please, check the value of MAX_TEXT_LINE_WIDTH
void init_dust_abund(const string& sim_path, vector<double>& d, double& grain_radius_init, double& grain_nb_density, double& grain_material_density)
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
void init_electron_spectra_mono(dynamic_array& init_spectrum, const vector<double>& en_grid,
	double electron_conc, double electron_energy)
{
	int i;
	double e1, e2, w1, w2;

	// the interval that contains the electron energy is [i, i+1],
	for (i = 0; i < (int)en_grid.size(); i++) {
		if (en_grid[i] * (1. - 1.e-9) > electron_energy) {
			i--;
			break;
		}
	}

	// in the case when test energy is very high,
	if (i == (int)en_grid.size()) {
		i = (int)en_grid.size() - 2;
	}
	memset(init_spectrum.arr, 0, init_spectrum.dim * sizeof(double));

	if (i == (int)en_grid.size() - 2) {
		init_spectrum.arr[i] = electron_conc / (en_grid[i + 1] - en_grid[i]);
	}
	else {
		if (electron_energy > 0.5 * (en_grid[i + 1] + en_grid[i])) {
			e1 = 0.5 * (en_grid[i + 1] + en_grid[i]);
			e2 = 0.5 * (en_grid[i + 2] + en_grid[i + 1]);
			
			w1 = (e2 - electron_energy) / (e2 - e1);
			w2 = 1. - w1;

			init_spectrum.arr[i] = w1 * electron_conc / (en_grid[i + 1] - en_grid[i]);
			init_spectrum.arr[i + 1] = w2 * electron_conc / (en_grid[i + 2] - en_grid[i + 1]);
		}
		else {
			e1 = 0.5 * (en_grid[i] + en_grid[i - 1]);
			e2 = 0.5 * (en_grid[i + 1] + en_grid[i]);
			
			w1 = (e2 - electron_energy) / (e2 - e1);
			w2 = 1. - w1;

			init_spectrum.arr[i - 1] = w1 * electron_conc / (en_grid[i] - en_grid[i - 1]);
			init_spectrum.arr[i] = w2 * electron_conc / (en_grid[i + 1] - en_grid[i]);
		}
	}
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
		output << left << setw(14) << en_grid[i] * 1.0001;

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

