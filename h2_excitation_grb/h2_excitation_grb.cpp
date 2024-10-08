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

#define MAX_TEXT_LINE_WIDTH 240  // maximal size of the comment lines in the files,
#define SOURCE_NAME "h2_excitation_grb.cpp"

using namespace std;


//
// Simulations of the evolution of electron spectra in the gas (plasma)
void simulate_h2_excitation(const string& data_path, const string& sim_path, const string& output_path, double op_ratio_h2,
    bool test_mode);

void init_electron_energy_grid(vector<double>& electron_energies_grid);
void init_dust_abund(const string& path, vector<double>& d, double& grain_nb_density, double& grain_material_density);
void init_electron_spectra(const string& path, vector<dynamic_array>& d, const vector<double>& electron_energies_grid, bool test_mode);
void init_specimen_conc(const string& path, double& conc_h_tot, vector<dynamic_array>& d);


int f_elsp(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
	return static_cast<elspectra_evolution_data*>(user_data)->f(t, y, ydot);
}

int main()
{
	bool test_mode;
	string sim_path;     // path to the folder with simulation data (GRB propagation results)

	// be careful about test simulations,
	test_mode = false;
	sim_path = "C:/Users/Александр/Documents/Данные и графики/GRB in molecular cloud/2022_04-v2/";
	//simulate_h2_excitation(data_path, sim_path, output_path, op_ratio_h2, test_mode);

}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.


// Initialization of the grid for electron energies,
// there IS NO interval for electrons with high energy (E > max_el_energy)
// there IS interval for electrons with low energy [0, min energy], min_energy must be < any excitation process threshold
void init_electron_energy_grid(vector<double>& electron_energies_grid)
{
	int i, k, nb_bins_1, nb_bins_per_order;
	double x1, x2, el_energy_1, max_el_energy;

	nb_bins_per_order = 150;
	el_energy_1 = 1.;  // eV
	max_el_energy = 1.e+4;  // eV

	x1 = pow(10., 1. / nb_bins_per_order);
	nb_bins_1 = (int)(1. / (x1 - 1)) + 1;

	x2 = el_energy_1 / nb_bins_1;
	for (i = 0; i <= nb_bins_1; i++) {
		electron_energies_grid.push_back(i * x2);
	}

	k = (int)(nb_bins_per_order * log10((double)max_el_energy / el_energy_1));
	for (i = 0; i < k; i++) {
		electron_energies_grid.push_back(electron_energies_grid[nb_bins_1 + i] * x1);
	}
}

void simulate_h2_excitation(const string& data_path, const string& sim_path, const string& output_path, double op_ratio_h2, bool test_mode)
{
	const int verbosity = 1;
	const double model_time_fin = 3.e+6;  // in s

	bool dust_is_presented(true);
	int i, k, nb_of_time_moments, lay_nb, step_nb, max_nb_steps, flag, nb_of_equat, nb_of_el_energies,
		nb_lev_h2, nb_lev_hei, h2eq_nb, heieq_nb, physeq_nb;
	long int nb_steps_tot;
	double x1, x2, x3, x4, x5, x6, x7, x8, x9, rel_tol, dt, model_time, model_time_step, model_time_out, conc_h_tot, grain_nb_density, grain_material_density, grain_radius,
		enloss_deriv_mt, enloss_deriv_h2_rot, enloss_deriv_h2_vibr, enloss_deriv_h2_electr, enloss_deriv_h2_electr_tr, enloss_deriv_ioniz, enloss_deriv_coloumb_ions,
		enloss_deriv_coloumb_el, enloss_deriv_hei, enloss_mt, enloss_h2_rot, enloss_h2_vibr, enloss_h2_electr, enloss_h2_electr_tr, enloss_ioniz, enloss_coloumb_ions,
		enloss_coloumb_el, enloss_hei, h2_solomon_diss_rate, h2_diss_exc_rate, h2_diss_exc_triplet_rate, hei_exc_rate,
		h2_solomon_diss_tot, h2_diss_exc_tot, h2_diss_exc_triplet_tot, hei_exc_tot;

	time_t timer;
	char* ctime_str;

	vector<double> time_cloud_arr;  // in s
	vector<double> enloss_deriv_mt_arr, enloss_deriv_h2_rot_arr, enloss_deriv_h2_vibr_arr, enloss_deriv_h2_electr_arr, enloss_deriv_h2_electr_tr_arr,
		enloss_deriv_ioniz_arr, enloss_deriv_coloumb_ions_arr, enloss_deriv_coloumb_el_arr, enloss_deriv_hei_arr;  // eV cm-3 s-1
	vector<double> enloss_mt_arr, enloss_h2_rot_arr, enloss_h2_vibr_arr, enloss_h2_electr_arr, enloss_h2_electr_tr_arr, enloss_ioniz_arr,
		enloss_coloumb_ions_arr, enloss_coloumb_el_arr, enloss_hei_arr;  // eV cm-3
	vector<double> electron_energies_grid;  // electron energy in eV,
	vector<double> h2_solomon_diss_tot_arr, h2_diss_exc_tot_arr, h2_diss_exc_triplet_tot_arr, hei_exc_tot_arr;  // [cm-3]
	vector<double> neutral_temp_arr, ion_temp_arr, dust_charge_arr;

	// vector index - layer nb, array index - electron energy bin nb:
	vector<dynamic_array> el_spectrum_init_data, el_spectrum_evol;
	// vector index - layer nb, array index - specimen nb:
	vector<dynamic_array> specimen_conc_init_data, specimen_conc_evol;
	// vector index - layer nb, array index - energy level nb:
	vector<dynamic_array> h2_popul_dens_evol, hei_popul_dens_evol;

	// dust/gas ratio and ionization fraction at the fixed time (start of electron spectra simulations), index - cloud layer,
	vector<double> dg_ratio_cloud_data, ioniz_cloud_data;

	cout << scientific;
	cout.precision(3);

	timer = time(NULL);
	ctime_str = ctime(&timer);

	cout << ctime_str << endl
		<< "Start of the simulations of the evolution of electron spectra" << endl
		<< "Initialization of the data..." << endl;


	init_electron_energy_grid(electron_energies_grid);
	init_dust_abund(sim_path, dg_ratio_cloud_data, grain_nb_density, grain_material_density);
	init_specimen_conc(sim_path, conc_h_tot, specimen_conc_init_data);  // gas density is initialized here,
	// electron spectrum in the file has the dimension [cm-3 eV-1],
	init_electron_spectra(sim_path, el_spectrum_init_data, electron_energies_grid, test_mode);

	elspectra_evolution_data user_data(data_path, conc_h_tot, electron_energies_grid, verbosity);

	nb_of_el_energies = (int)electron_energies_grid.size() - 1;
	nb_of_equat = user_data.get_nb_of_equat();
	nb_lev_h2 = user_data.get_nb_of_h2_lev();
	nb_lev_hei = user_data.get_nb_of_hei_lev();

	h2eq_nb = user_data.get_h2eq_nb();
	heieq_nb = user_data.get_heieq_nb();
	physeq_nb = user_data.get_physeq_nb();

	dynamic_array el_spectrum(nb_of_el_energies), specimen_conc(NB_OF_CHEM_SPECIES), h2_popul_dens(nb_lev_h2), hei_popul_dens(nb_lev_hei);

	time_cloud_arr.clear();
	for (i = 0; i <= 5; i++) {
		time_cloud_arr.push_back(i * 2.e+3);
	}
	for (i = 1; i <= 5; i++) {
		time_cloud_arr.push_back(i * 2.e+4);
	}
	for (i = 2; time_cloud_arr.back() < model_time_fin; i++) {
		time_cloud_arr.push_back(i * 1.e+5);
	}
	nb_of_time_moments = (int)time_cloud_arr.size();

	// initialization of vectors and matrices used by solver SUNDIALS CVODE
	N_Vector y, ydot, abs_tol;
	SUNMatrix A(NULL);
	SUNLinearSolver LS(NULL);

	y = N_VNew_Serial(nb_of_equat);
	ydot = N_VNew_Serial(nb_of_equat);
	abs_tol = N_VNew_Serial(nb_of_equat);

	// Initialization of y values,
	lay_nb = 1200; // 1070;
	if (lay_nb > (int)el_spectrum_init_data.size())
		lay_nb = (int)el_spectrum_init_data.size() - 1;

	for (i = 0; i < nb_of_equat; i++) {
		NV_Ith_S(y, i) = 0.;
	}
	for (i = 0; i < nb_of_el_energies; i++) {
		// the number of electrons in the interval:
		NV_Ith_S(y, i) = el_spectrum_init_data[lay_nb].arr[i] * user_data.get_electron_energy_bin(i);
	}
	for (i = 0; i < NB_OF_CHEM_SPECIES; i++) {
		NV_Ith_S(y, nb_of_el_energies + i) = specimen_conc_init_data[lay_nb].arr[i];
	}

	// the initial distribution of energy level populations is postulated,
	NV_Ith_S(y, h2eq_nb) = NV_Ith_S(y, nb_of_el_energies + H2_NB) / (op_ratio_h2 + 1.);
	NV_Ith_S(y, h2eq_nb + 1) = NV_Ith_S(y, nb_of_el_energies + H2_NB) * op_ratio_h2 / (op_ratio_h2 + 1.);

	// the initial distribution of HeI levels,
	NV_Ith_S(y, heieq_nb) = NV_Ith_S(y, nb_of_el_energies + HE_NB);

	// Initial gas temperature, in K
	NV_Ith_S(y, physeq_nb) = NV_Ith_S(y, physeq_nb + 1) = 10.;

	// Dust charge
	NV_Ith_S(y, physeq_nb + 2) = 0.;

	// calculation of grain cross section,
	if (dg_ratio_cloud_data[lay_nb] <= numeric_limits<double>::epsilon()) {
		dust_is_presented = false;
		grain_radius = 0.;
	}
	else {
		dust_is_presented = true;  // He?
		x1 = dg_ratio_cloud_data[lay_nb] * conc_h_tot * ATOMIC_MASS_UNIT * 3. / (grain_nb_density * grain_material_density * 4. * M_PI);
		grain_radius = pow(x1, 1. / 3.);   // grain radius, cm
	}
	user_data.set_dust_parameters(dust_is_presented, grain_radius, grain_nb_density);

	// Initialization for tolerances:
	rel_tol = REL_ERROR_SOLVER;
	user_data.set_tolerances(abs_tol);

	// Call CVodeCreate to create the solver memory and specify the Backward Differentiation Formula and the use of a Newton iteration 
	void* cvode_mem = CVodeCreate(CV_BDF);

	// Call CVodeInit to initialize the integrator memory and specify the user's right hand side function in y'=f(t,y), 
	// the initial time t0, and the initial dependent variable vector y;
	flag = CVodeInit(cvode_mem, f_elsp, model_time = 0., y);

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

	cout << left << "Initialization step time (s): " << setw(12) << (int)(time(NULL) - timer) << endl
		<< "Starting simulations..." << endl;
	timer = time(NULL);

	if (verbosity) {
		cout << left << "lay_nb: " << lay_nb << endl
			<< left << setw(12) << "model_time " << setw(9) << "steps_nb "
			<< setw(12) << "gas_temp(K)" << setw(12) << "ion_temp(K)" << setw(12) << "dust_charge" << endl;
	}

	max_nb_steps = 30;
	model_time_step = 1.;

	k = 0;
	enloss_deriv_mt = enloss_deriv_h2_rot = enloss_deriv_h2_vibr = enloss_deriv_h2_electr = enloss_deriv_h2_electr_tr
		= enloss_deriv_ioniz = enloss_deriv_coloumb_ions = enloss_deriv_coloumb_el = enloss_deriv_hei = 0.;

	enloss_mt = enloss_h2_rot = enloss_h2_vibr = enloss_h2_electr = enloss_h2_electr_tr = enloss_ioniz = enloss_coloumb_ions
		= enloss_coloumb_el = enloss_hei = 0.;

	h2_solomon_diss_rate = h2_diss_exc_rate = h2_diss_exc_triplet_rate = hei_exc_rate
		= h2_solomon_diss_tot = h2_diss_exc_tot = h2_diss_exc_triplet_tot = hei_exc_tot = 0.;

	while (model_time < time_cloud_arr.back())
	{
		// saving data,
		if (model_time >= time_cloud_arr[k]) {
			// updating time grid,
			time_cloud_arr[k] = model_time;

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
			h2_popul_dens_evol.push_back(h2_popul_dens);

			for (i = 0; i < nb_lev_hei; i++) {
				hei_popul_dens.arr[i] = NV_Ith_S(y, heieq_nb + i);
			}
			hei_popul_dens_evol.push_back(hei_popul_dens);

			enloss_deriv_mt_arr.push_back(enloss_deriv_mt);
			enloss_deriv_h2_rot_arr.push_back(enloss_deriv_h2_rot);
			enloss_deriv_h2_vibr_arr.push_back(enloss_deriv_h2_vibr);
			enloss_deriv_h2_electr_arr.push_back(enloss_deriv_h2_electr);
			enloss_deriv_h2_electr_tr_arr.push_back(enloss_deriv_h2_electr_tr);
			enloss_deriv_ioniz_arr.push_back(enloss_deriv_ioniz);
			enloss_deriv_coloumb_ions_arr.push_back(enloss_deriv_coloumb_ions);
			enloss_deriv_coloumb_el_arr.push_back(enloss_deriv_coloumb_el);
			enloss_deriv_hei_arr.push_back(enloss_deriv_hei);

			enloss_mt_arr.push_back(enloss_mt);
			enloss_h2_rot_arr.push_back(enloss_h2_rot);
			enloss_h2_vibr_arr.push_back(enloss_h2_vibr);
			enloss_h2_electr_arr.push_back(enloss_h2_electr);
			enloss_h2_electr_tr_arr.push_back(enloss_h2_electr_tr);
			enloss_ioniz_arr.push_back(enloss_ioniz);
			enloss_coloumb_ions_arr.push_back(enloss_coloumb_ions);
			enloss_coloumb_el_arr.push_back(enloss_coloumb_el);
			enloss_hei_arr.push_back(enloss_hei);

			h2_solomon_diss_tot_arr.push_back(h2_solomon_diss_tot);
			h2_diss_exc_tot_arr.push_back(h2_diss_exc_tot);
			h2_diss_exc_triplet_tot_arr.push_back(h2_diss_exc_triplet_tot);
			hei_exc_tot_arr.push_back(hei_exc_tot);

			neutral_temp_arr.push_back(NV_Ith_S(y, physeq_nb));
			ion_temp_arr.push_back(NV_Ith_S(y, physeq_nb + 1));
			dust_charge_arr.push_back(NV_Ith_S(y, physeq_nb + 2) * grain_nb_density);  // dust charge per cm3

			save_electron_spectrum_evolution(output_path, &user_data, time_cloud_arr, el_spectrum_evol);
			save_specimen_conc_evolution(output_path, time_cloud_arr, specimen_conc_evol);
			save_h2_populations_evolution(output_path, time_cloud_arr, h2_popul_dens_evol, user_data.get_pointer_h2_levels());
			save_hei_populations_evolution(output_path, time_cloud_arr, hei_popul_dens_evol, user_data.get_pointer_hei_levels());
			save_electron_energy_losses(output_path, time_cloud_arr,
				enloss_deriv_mt_arr, enloss_mt_arr,
				enloss_deriv_h2_rot_arr, enloss_h2_rot_arr,
				enloss_deriv_h2_vibr_arr, enloss_h2_vibr_arr,
				enloss_deriv_h2_electr_arr, enloss_h2_electr_arr,
				enloss_deriv_h2_electr_tr_arr, enloss_h2_electr_tr_arr,
				enloss_deriv_ioniz_arr, enloss_ioniz_arr,
				enloss_deriv_coloumb_ions_arr, enloss_coloumb_ions_arr,
				enloss_deriv_coloumb_el_arr, enloss_coloumb_el_arr,
				enloss_deriv_hei_arr, enloss_hei_arr);

			save_h2_data(output_path, time_cloud_arr, h2_solomon_diss_tot_arr, h2_diss_exc_tot_arr, h2_diss_exc_triplet_tot_arr,
				hei_exc_tot_arr);

			save_phys_parameters(output_path, time_cloud_arr, neutral_temp_arr, ion_temp_arr, dust_charge_arr, specimen_conc_evol);
			k++;
		}

		step_nb = 0;
		flag = CV_SUCCESS;
		model_time_out = model_time + model_time_step;

		dt = model_time;
		while (step_nb < max_nb_steps && flag == CV_SUCCESS && model_time < model_time_out) {
			step_nb++;
			flag = CVode(cvode_mem, model_time_out, y, &model_time, CV_ONE_STEP);  // CV_NORMAL or CV_ONE_STEP	
		}
		dt = 0.5 * (model_time - dt);

		if (flag == CV_SUCCESS) {
			flag = CVodeGetNumSteps(cvode_mem, &nb_steps_tot);

			x1 = enloss_deriv_mt;
			x2 = enloss_deriv_h2_rot;
			x3 = enloss_deriv_h2_vibr;
			x4 = enloss_deriv_h2_electr;
			x5 = enloss_deriv_ioniz;
			x6 = enloss_deriv_coloumb_ions;
			x7 = enloss_deriv_coloumb_el;
			x8 = enloss_deriv_hei;
			x9 = enloss_deriv_h2_electr_tr;

			user_data.get_el_energy_losses(enloss_deriv_mt, enloss_deriv_h2_rot, enloss_deriv_h2_vibr, enloss_deriv_h2_electr, enloss_deriv_h2_electr_tr,
				enloss_deriv_ioniz, enloss_deriv_coloumb_ions, enloss_deriv_coloumb_el, enloss_deriv_hei);

			enloss_mt += (x1 + enloss_deriv_mt) * dt;
			enloss_h2_rot += (x2 + enloss_deriv_h2_rot) * dt;
			enloss_h2_vibr += (x3 + enloss_deriv_h2_vibr) * dt;
			enloss_h2_electr += (x4 + enloss_deriv_h2_electr) * dt;
			enloss_ioniz += (x5 + enloss_deriv_ioniz) * dt;
			enloss_coloumb_ions += (x6 + enloss_deriv_coloumb_ions) * dt;
			enloss_coloumb_el += (x7 + enloss_deriv_coloumb_el) * dt;
			enloss_hei += (x8 + enloss_deriv_hei) * dt;
			enloss_h2_electr_tr += (x9 + enloss_deriv_h2_electr_tr) * dt;

			x1 = h2_solomon_diss_rate;
			x2 = h2_diss_exc_rate;
			x3 = h2_diss_exc_triplet_rate;
			x4 = hei_exc_rate;

			user_data.get_h2_process_rates(h2_solomon_diss_rate, h2_diss_exc_rate, h2_diss_exc_triplet_rate, hei_exc_rate);

			h2_solomon_diss_tot += (x1 + h2_solomon_diss_rate) * dt;
			h2_diss_exc_tot += (x2 + h2_diss_exc_rate) * dt;
			h2_diss_exc_triplet_tot += (x3 + h2_diss_exc_triplet_rate) * dt;
			hei_exc_tot += (x4 + hei_exc_rate) * dt;

			if (verbosity) {
				cout.precision(2);
				cout << left << setw(12) << model_time << setw(9) << nb_steps_tot <<
					setw(12) << NV_Ith_S(y, physeq_nb) << setw(12) << NV_Ith_S(y, physeq_nb + 1) << setw(12) << NV_Ith_S(y, physeq_nb + 2) << endl;
			}
		}
		else {
			cout << "Some unexpected error has been occurred: " << flag << endl;
			break;
		}
	}
	cout << left << "Total calculation time (s): " << setw(12) << (int)(time(NULL) - timer) << endl;

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

// please, check the value of MAX_TEXT_LINE_WIDTH
void init_dust_abund(const string& path, vector<double>& d, double& grain_nb_density, double& grain_material_density)
{
	int i, j, l, j0, nb_lay, nb_t;
	char ch, text_line[MAX_TEXT_LINE_WIDTH];
	double x;

	string fname;
	ifstream input;

	fname = path + "grb_dust_param.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> ch >> grain_nb_density >> grain_material_density;
	input >> ch >> nb_lay >> nb_t;
	input >> ch;

	j0 = nb_t - 1;  // last column of data,
	for (j = 0; j < nb_t; j++) {
		input >> x;
	}
	for (i = 0; i < nb_lay; i++) {
		input >> l >> x;
		for (j = 0; j < nb_t; j++) {
			input >> x;
			if (j == j0) {
				d.push_back(x);
			}
		}
	}
	input.close();
}

// 
void init_electron_spectra(const string& path, vector<dynamic_array>& init_data, const vector<double>& en_grid, bool test_mode)
{
	int i, j, l, k, nb_l, en_nb_f, en_nb, nb_bins;
	char ch, text_line[MAX_TEXT_LINE_WIDTH];
	double e, de, a, b, min_en, max_en;

	string fname;
	ifstream input;
	vector<dynamic_array> file_data;
	vector<double> en_grid_file;

	en_nb = (int)en_grid.size() - 1;  // nb of energy bins,

	fname = path + "grb_electron_spectrum.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> ch >> min_en >> max_en >> nb_bins >> en_nb_f >> nb_l;
	dynamic_array el_spectrum(en_nb_f), el_spectrum_2(en_nb);

	input >> ch;
	for (j = 0; j < nb_l; j++) {
		input >> a >> b;
		file_data.push_back(el_spectrum);  // empty arrays
		init_data.push_back(el_spectrum_2);
	}

	// in the file the energy spectrum of electrons is given in [cm-3 eV-1], half length of the energy interval is provided,
	en_grid_file.push_back(0.);
	for (i = 0; i < en_nb_f; i++) {
		input >> e >> de;
		en_grid_file.push_back(e + de);

		for (l = 0; l < nb_l; l++) {
			input >> a >> b;

			// vector index - layer nb, array index - electron energy bin nb,	
			if (test_mode) {
				file_data[l].arr[i] = 0.;
				if (i == en_nb_f - 1)
					file_data[l].arr[i] = 0.1 / (2. * de);
			}
			else {
				file_data[l].arr[i] = a;
			}
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
}

void init_specimen_conc(const string& path, double& conc_h_tot, vector<dynamic_array>& d)
{
	int i, j, nb_l, nb_sp;
	double x;
	char ch, text_line[MAX_TEXT_LINE_WIDTH];

	string fname;
	ifstream input;
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

	input >> ch >> conc_h_tot >> nb_l >> nb_sp;
	if (nb_sp != NB_OF_CHEM_SPECIES) {
		cout << "Error in " << SOURCE_NAME << ": nb of species in input file do not coincide with that in the code." << endl;
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);  // reading the end of the previous line
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	for (i = 0; i < nb_l; i++) {
		d.push_back(specimen_conc);
		input >> j >> x;

		for (j = 0; j < nb_sp; j++) {
			input >> x;
			// vector index - layer nb, array index - specimen nb,
			d[i].arr[j] = x;
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