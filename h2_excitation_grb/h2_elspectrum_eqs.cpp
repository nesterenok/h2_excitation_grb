#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "interpolation.h"
#include "h2_elspectrum_eqs.h"
#include "h2_microphysics.h"

#include "constants.h"
#include "spectroscopy.h"
#include "h2_parameters.h"

#include <omp.h>
#include<cmath>
#include<limits>
#include<cstring>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<sstream>

using namespace std;
#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "grb_elspectrum_eqs.cpp"


// Rotational level nb of the initial vibrational state of the ground electronic state for which the excitation vi -> vf is taken into account 
// must be <= nb of levels of the initial vibrational state, maximal for vi = 0 is 29, 
#define MAX_J_H2_VIBR_INIT 15

// Nb of vibrational states of excited electronic state taken into account,
// this nb is reconciled with CLOUDY data (Abgrall et al., A&AS, 141, 297-300, 2000);
#define MAX_H2_VSTATES_B1SU 38  // 0,1,..,37 (including) in CLOUDY
#define MAX_H2_VSTATES_C1PU 14  // 0,1,..,13 in CLOUDY
#define MAX_H2_VSTATES_BP1SU 10  // 0,1,..,9 in CLOUDY
#define MAX_H2_VSTATES_D1PU 19  // 0,1,..,18 in CLOUDY - for D-; for D+ v=0,1,2; for D- levels with j=0 were commented in original file

const double log_coulomb_constant = log(1.5 * pow(BOLTZMANN_CONSTANT, 1.5)
	/ (sqrt(M_PI) * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE)); // *T_e^1.5/n_e^0.5

const double el_ion_scatt_const = 8. * SQRT_PI * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE 
	/(ATOMIC_MASS_UNIT * EV_TO_ERGS);   // [eV cm^4 s^-2]

//
// Structure of the data array:
// electron_spectrum [cm-3] (nb of electrons per cm-3 in energy interval), 
// concentration [cm-3] e-, H, H+, H2, H2+, He, He+, He++,
// H2 level population densities, [cm-3] 
// HeI level populations, [cm-3]
// gas temperature (of neutrals and ions) [K], dust grain charge

el_spectra_data::el_spectra_data(const el_spectra_data& obj)
{
	i0 = obj.i0;
	nb = obj.nb;
	weights = new double[nb];
	memcpy(weights, obj.weights, nb * sizeof(double));
}

el_spectra_data& el_spectra_data::operator=(const el_spectra_data& obj)
{
	i0 = obj.i0;
	nb = obj.nb;

	delete[] weights;
	weights = new double[nb];
	memcpy(weights, obj.weights, nb * sizeof(double));

	return *this;
}

elspectra_evolution_data::elspectra_evolution_data(const string& data_path, double cht, const vector<double>& en_grid, int verb)
	: conc_h_tot(cht), dust_is_presented(false), grain_cs(0.), grain_nb_density(0.),
	enloss_deriv_mt(0.), enloss_deriv_h2_rot(0.), enloss_deriv_h2_vibr(0.), enloss_deriv_h2_electr(0.), enloss_deriv_h2_electr_tr(0.),
	enloss_deriv_ioniz(0.), enloss_deriv_coloumb_ions(0.), enloss_deriv_coloumb_el(0.), enloss_deriv_hei(0.), 
	conc_n(0.), conc_i(0.), energy_gain_n(0.), energy_gain_i(0.), nb_gain_n(0.), nb_gain_i(0.), 
	h2_solomon_diss_rate(0.), h2_diss_exc_rate(0.), h2_diss_exc_triplet_rate(0.), hei_exc_rate(0.)
{
	int i, j0, isotope, electronic_state, verbosity = 1;
	double en, den, vel;
	string fname;
	stringstream ss;

	nb_of_el_energies = (int)en_grid.size() - 1;

	// kinetic energy of electron is considered, E_rel = E_kin + mc2, 
	electron_energies_grid = new double[nb_of_el_energies + 1];
	electron_energy_bin_size = new double[nb_of_el_energies];
	electron_energies = new double[nb_of_el_energies];
	electron_velocities = new double[nb_of_el_energies];

	for (i = 0; i <= nb_of_el_energies; i++) {
		electron_energies_grid[i] = en_grid[i];
	}
	for (i = 0; i < nb_of_el_energies; i++) {
		electron_energy_bin_size[i] = electron_energies_grid[i + 1] - electron_energies_grid[i];
		electron_energies[i] = electron_energies_grid[i] + 0.5 * electron_energy_bin_size[i];

		// relativistic velocity,
		electron_velocities[i] = SPEED_OF_LIGHT * sqrt(2. * ELECTRON_MASS_EV * electron_energies[i] + electron_energies[i] * electron_energies[i])
			/ (ELECTRON_MASS_EV + electron_energies[i]);
		// = sqrt(2. * electron_energies[i] * EV_TO_ERGS / ELECTRON_MASS);
	}
	
	// Initialization of cross section data for elastic scattering:
	elastic_h2_el_cs
		= new cross_section_table_vs1(data_path, "elastic_scattering/e-h2_mt_cs.txt");

	elastic_he_el_cs
		= new cross_section_table_vs1(data_path, "elastic_scattering/e-he_mt_cs.txt");

	el_enloss_h2_mt = new double[nb_of_el_energies];
	el_enloss_he_mt = new double[nb_of_el_energies];
	
	elel_scatt_rates = new double[nb_of_el_energies];
	elion_scatt_rates = new double[nb_of_el_energies];

	el_enloss_h2_mt[0] = el_enloss_he_mt[0] = elel_scatt_rates[0] = elion_scatt_rates[0] = 0.;

	for (i = 1; i < nb_of_el_energies; i++) {
		// it is assumed that electron energy is the kinetic energy, non-relativistic, in eV
		en = electron_energies[i];  // the centre of the interval,
		vel = electron_velocities[i];
		// 
		den = 0.5 * (electron_energy_bin_size[i] + electron_energy_bin_size[i - 1]);
		
		// energy losses are normalized on target (H2, H, He) concentration and energy bin size, [cm2 eV cm/s eV-1]
		// m_h2 = 2 a.m.u., m_he = 4 a.m.u.
		el_enloss_h2_mt[i] = (*elastic_h2_el_cs)(en) *en *vel * ELECTRON_MASS/(ATOMIC_MASS_UNIT * den); 
		el_enloss_he_mt[i] = (*elastic_he_el_cs)(en) * en * vel * 0.5* ELECTRON_MASS / (ATOMIC_MASS_UNIT *den); 

		// energy loss rate coefficient in electron-electron collisions, [eV cm^3 s^-1]
		elel_scatt_rates[i] = 4.*M_PI * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE 
			/ (ELECTRON_MASS * vel * EV_TO_ERGS);

		// rates are normalized on electron energy bin size, [eV cm^3 s^-1 eV-1]
		//elion_scatt_rates[i] = 4. * M_PI * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE
		//	/ (ATOMIC_MASS_UNIT * vel * EV_TO_ERGS * den);
	}

	// H2 molecule data:
	nb_lev_h2 = 298; // the maximal number for which Einstein coefficients are provided - 298 levels,	
	molecule h2_mol("H2", isotope = 1, 2. * ATOMIC_MASS_UNIT);

	h2_di = new h2_diagram(data_path, h2_mol, nb_lev_h2, verbosity);
	h2_einst = new h2_einstein_coeff(data_path, h2_di, verbosity);

	// Spectroscopic data for electronically excited H2 states,
	nb_lev_h2_b = 882;
	h2_di_b = new h2_diagram(electronic_state = 1, data_path, h2_mol, nb_lev_h2_b, verbosity);

	nb_lev_h2_cplus = 248;
	h2_di_cplus = new h2_diagram(electronic_state = 2, data_path, h2_mol, nb_lev_h2_cplus, verbosity);

	nb_lev_h2_cminus = 251;
	h2_di_cminus = new h2_diagram(electronic_state = 3, data_path, h2_mol, nb_lev_h2_cminus, verbosity);

	nb_lev_h2_bp = 108;
	h2_di_bp = new h2_diagram(electronic_state = 4, data_path, h2_mol, nb_lev_h2_bp, verbosity);

	nb_lev_h2_dplus = 28;
	h2_di_dplus = new h2_diagram(electronic_state = 5, data_path, h2_mol, nb_lev_h2_dplus, verbosity);

	nb_lev_h2_dminus = 336;
	h2_di_dminus = new h2_diagram(electronic_state = 6, data_path, h2_mol, nb_lev_h2_dminus, verbosity);

	h2_band_transitions(data_path, lyman_band_h2, h2_di_b, h2_di, verbosity);
	h2_band_transitions(data_path, werner_plus_band_h2, h2_di_cplus, h2_di, verbosity);
	h2_band_transitions(data_path, werner_minus_band_h2, h2_di_cminus, h2_di, verbosity);
	
	h2_band_transitions(data_path, bp_band_h2, h2_di_bp, h2_di, verbosity);
	h2_band_transitions(data_path, dplus_band_h2, h2_di_dplus, h2_di, verbosity);
	h2_band_transitions(data_path, dminus_band_h2, h2_di_dminus, h2_di, verbosity);
	
	h2_excited_state_data(data_path, h2_b_state_data, h2_di_b, lyman_band_h2, verbosity);
	h2_excited_state_data(data_path, h2_cplus_state_data, h2_di_cplus, werner_plus_band_h2, verbosity);
	h2_excited_state_data(data_path, h2_cminus_state_data, h2_di_cminus, werner_minus_band_h2, verbosity);

	h2_excited_state_data(data_path, h2_bp_state_data, h2_di_bp, bp_band_h2, verbosity);
	h2_excited_state_data(data_path, h2_dplus_state_data, h2_di_dplus, dplus_band_h2, verbosity);
	h2_excited_state_data(data_path, h2_dminus_state_data, h2_di_dminus, dminus_band_h2, verbosity);

	// HeI spectroscopic data
	nb_lev_hei = 31;
	molecule ion_hei("HeI", isotope = 1, 4.* ATOMIC_MASS_UNIT);

	hei_di = new ion_diagram(data_path, ion_hei, nb_lev_hei, verbosity);
	hei_einst = new ion_einstein_coeff(data_path, hei_di, verbosity);


	// Excitation of H2 rotational levels,
	rot_h2_j02_cs
		= new cross_section_table_vs1(data_path, "coll_h2/cross_sections/h2_e_J_0_2.txt");

	rot_h2_j24_cs
		= new cross_section_table_vs1(data_path, "coll_h2/cross_sections/h2_e_J_2_4.txt");

	rot_h2_j13_cs
		= new cross_section_table_vs1(data_path, "coll_h2/cross_sections/h2_e_J_1_3.txt");

	rot_h2_j35_cs
		= new cross_section_table_vs1(data_path, "coll_h2/cross_sections/h2_e_J_3_5.txt");

	init_tables_h2_rotational_exc(rot_h2_j02_cs, h2_j02_rates, h2_j02_indexes, h2_j20_rates, h2_j20_indexes, j0 = 0);
	init_tables_h2_rotational_exc(rot_h2_j24_cs, h2_j24_rates, h2_j24_indexes, h2_j42_rates, h2_j42_indexes, j0 = 2);
	init_tables_h2_rotational_exc(rot_h2_j13_cs, h2_j13_rates, h2_j13_indexes, h2_j31_rates, h2_j31_indexes, j0 = 1);
	init_tables_h2_rotational_exc(rot_h2_j35_cs, h2_j35_rates, h2_j35_indexes, h2_j53_rates, h2_j53_indexes, j0 = 3);

	// Excitation of H2 vibrational states of the ground electronic state,
	// Check the cross sections with MCCC calculations published by Padovani et al. A&A (2022),
	vibr_h2_v01_cs 
		= new cross_section_table_vs1(data_path, "coll_h2/cross_sections/h2_e_v_0_1.txt");

	vibr_h2_v02_cs
		= new h2_excitation_vibr02_cs();

	init_tables_h2_vibrational_exc(vibr_h2_v01_cs, h2_v01_rates, h2_v01_indexes, 1, verbosity);
	init_tables_h2_vibrational_exc(vibr_h2_v02_cs, h2_v02_rates, h2_v02_indexes, 2, verbosity);

	// Excitation of electronically excited states of H2,
	// S1g(X) -> S1u(B), v0 -> vf
	h2_v0_bstate_cs = new cross_section_table_mccc * [MAX_H2_VSTATES_B1SU];
	for (i = 0; i < MAX_H2_VSTATES_B1SU; i++) {
		ss.clear();
		ss.str("");
		ss << "MCCC-el-H2-B1Su_vf=";
		ss << i;
		ss << ".X1Sg_vi=0.txt";
		fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=0/" + ss.str();

		h2_v0_bstate_cs[i] = new cross_section_table_mccc(data_path, fname);
	}
	init_tables_h2_electronic_exc(h2_di_b, h2_v0_bstate_cs, h2_v0_bstate_rates, h2_v0_bstate_indexes, 
		0, MAX_H2_VSTATES_B1SU, h2_b1su_min_nb, verbosity);
	
	// S1g(X) -> S1u(B), v1 -> vf
	h2_v1_bstate_cs = new cross_section_table_mccc * [MAX_H2_VSTATES_B1SU];
	for (i = 0; i < MAX_H2_VSTATES_B1SU; i++) {
		ss.clear();
		ss.str("");
		ss << "MCCC-el-H2-B1Su_vf=";
		ss << i;
		ss << ".X1Sg_vi=1.txt";
		fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=1/" + ss.str();

		h2_v1_bstate_cs[i] = new cross_section_table_mccc(data_path, fname);
	}
	// Note: look at the minimal nb of electron energy bin, is lower at higher v_i,
	init_tables_h2_electronic_exc(h2_di_b, h2_v1_bstate_cs, h2_v1_bstate_rates, h2_v1_bstate_indexes, 
		1, MAX_H2_VSTATES_B1SU, h2_b1su_min_nb, verbosity);

	// S1g(X) -> P1u(C-/+)
	// the level energy difference of C- and C+ is neglected, C- state have higher number of levels,
	// v0 -> vf
	h2_v0_cstate_cs = new cross_section_table_mccc * [MAX_H2_VSTATES_C1PU];
	for (i = 0; i < MAX_H2_VSTATES_C1PU; i++) {
		ss.clear();
		ss.str("");
		ss << "MCCC-el-H2-C1Pu_vf=";
		ss << i;
		ss << ".X1Sg_vi=0.txt";
		fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=0/" + ss.str();

		h2_v0_cstate_cs[i] = new cross_section_table_mccc(data_path, fname);
	}
	init_tables_h2_electronic_exc(h2_di_cminus, h2_v0_cstate_cs, h2_v0_cstate_rates, h2_v0_cstate_indexes, 
		0, MAX_H2_VSTATES_C1PU, h2_c1pu_min_nb, verbosity);

	// v1 -> vf
	h2_v1_cstate_cs = new cross_section_table_mccc * [MAX_H2_VSTATES_C1PU];
	for (i = 0; i < MAX_H2_VSTATES_C1PU; i++) {
		ss.clear();
		ss.str("");
		ss << "MCCC-el-H2-C1Pu_vf=";
		ss << i;
		ss << ".X1Sg_vi=1.txt";
		fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=1/" + ss.str();

		h2_v1_cstate_cs[i] = new cross_section_table_mccc(data_path, fname);
	}
	init_tables_h2_electronic_exc(h2_di_cminus, h2_v1_cstate_cs, h2_v1_cstate_rates, h2_v1_cstate_indexes, 
		1, MAX_H2_VSTATES_C1PU, h2_c1pu_min_nb, verbosity);

	// S1g(X) -> S1u(Bp), v0 -> vf
	h2_v0_bpstate_cs = new cross_section_table_mccc * [MAX_H2_VSTATES_BP1SU];
	for (i = 0; i < MAX_H2_VSTATES_BP1SU; i++) {
		ss.clear();
		ss.str("");
		ss << "MCCC-el-H2-Bp1Su_vf=";
		ss << i;
		ss << ".X1Sg_vi=0.txt";
		fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=0/" + ss.str();

		h2_v0_bpstate_cs[i] = new cross_section_table_mccc(data_path, fname);
	}
	init_tables_h2_electronic_exc(h2_di_bp, h2_v0_bpstate_cs, h2_v0_bpstate_rates, h2_v0_bpstate_indexes, 
		0, MAX_H2_VSTATES_BP1SU, h2_bp1su_min_nb, verbosity);

	// S1g(X) -> S1u(Bp), v1 -> vf
	h2_v1_bpstate_cs = new cross_section_table_mccc * [MAX_H2_VSTATES_BP1SU];
	for (i = 0; i < MAX_H2_VSTATES_BP1SU; i++) {
		ss.clear();
		ss.str("");
		ss << "MCCC-el-H2-Bp1Su_vf=";
		ss << i;
		ss << ".X1Sg_vi=1.txt";
		fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=1/" + ss.str();

		h2_v1_bpstate_cs[i] = new cross_section_table_mccc(data_path, fname);
	}
	init_tables_h2_electronic_exc(h2_di_bp, h2_v1_bpstate_cs, h2_v1_bpstate_rates, h2_v1_bpstate_indexes, 
		1, MAX_H2_VSTATES_BP1SU, h2_bp1su_min_nb, verbosity);

	// S1g(X) -> P1u(D-/+)
	// the level energy difference of D- and D+ is neglected, C- state have much higher number of levels,
	// v0 -> vf
	h2_v0_dstate_cs = new cross_section_table_mccc * [MAX_H2_VSTATES_D1PU];
	for (i = 0; i < MAX_H2_VSTATES_D1PU; i++) {
		ss.clear();
		ss.str("");
		ss << "MCCC-el-H2-D1Pu_vf=";
		ss << i;
		ss << ".X1Sg_vi=0.txt";
		fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=0/" + ss.str();

		h2_v0_dstate_cs[i] = new cross_section_table_mccc(data_path, fname);
	}
	init_tables_h2_electronic_exc(h2_di_dminus, h2_v0_dstate_cs, h2_v0_dstate_rates, h2_v0_dstate_indexes,
		0, MAX_H2_VSTATES_D1PU, h2_d1pu_min_nb, verbosity);

	// v1 -> vf
	h2_v1_dstate_cs = new cross_section_table_mccc * [MAX_H2_VSTATES_D1PU];
	for (i = 0; i < MAX_H2_VSTATES_D1PU; i++) {
		ss.clear();
		ss.str("");
		ss << "MCCC-el-H2-D1Pu_vf=";
		ss << i;
		ss << ".X1Sg_vi=1.txt";
		fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=1/" + ss.str();

		h2_v1_dstate_cs[i] = new cross_section_table_mccc(data_path, fname);
	}
	init_tables_h2_electronic_exc(h2_di_dminus, h2_v1_dstate_cs, h2_v1_dstate_rates, h2_v1_dstate_indexes,
		1, MAX_H2_VSTATES_D1PU, h2_d1pu_min_nb, verbosity);


	// Electronic dissociative excitation of H2
	// S1g(X) -> S1u(B), v0,v1 -> dissociative continuum
	fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=0/MCCC-el-H2-B1Su_DE.X1Sg_vi=0.txt";
	h2_v0_bstate_diss_cs = new cross_section_table_mccc(data_path, fname);
	
	init_tables_h2_electronic_diss(h2_v0_bstate_diss_cs, h2_v0_bstate_diss_rates, h2_v0_bstate_diss_indexes, 
		0, h2_b1su_diss_min_nb, verbosity);

	fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=1/MCCC-el-H2-B1Su_DE.X1Sg_vi=1.txt";
	h2_v1_bstate_diss_cs = new cross_section_table_mccc(data_path, fname);

	init_tables_h2_electronic_diss(h2_v1_bstate_diss_cs, h2_v1_bstate_diss_rates, h2_v1_bstate_diss_indexes,
		1, h2_b1su_diss_min_nb, verbosity);

	// S1g(X) -> P1u(C), v0,v1 -> dissociative continuum
	fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=0/MCCC-el-H2-C1Pu_DE.X1Sg_vi=0.txt";
	h2_v0_cstate_diss_cs = new cross_section_table_mccc(data_path, fname);

	init_tables_h2_electronic_diss(h2_v0_cstate_diss_cs, h2_v0_cstate_diss_rates, h2_v0_cstate_diss_indexes,
		0, h2_c1pu_diss_min_nb, verbosity);

	fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=1/MCCC-el-H2-C1Pu_DE.X1Sg_vi=1.txt";
	h2_v1_cstate_diss_cs = new cross_section_table_mccc(data_path, fname);

	init_tables_h2_electronic_diss(h2_v1_cstate_diss_cs, h2_v1_cstate_diss_rates, h2_v1_cstate_diss_indexes,
		1, h2_c1pu_diss_min_nb, verbosity);

	// S1g(X) -> S1u(Bp), v0,v1 -> dissociative continuum
	fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=0/MCCC-el-H2-Bp1Su_DE.X1Sg_vi=0.txt";
	h2_v0_bpstate_diss_cs = new cross_section_table_mccc(data_path, fname);

	init_tables_h2_electronic_diss(h2_v0_bpstate_diss_cs, h2_v0_bpstate_diss_rates, h2_v0_bpstate_diss_indexes,
		0, h2_bp1su_diss_min_nb, verbosity);

	fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=1/MCCC-el-H2-Bp1Su_DE.X1Sg_vi=1.txt";
	h2_v1_bpstate_diss_cs = new cross_section_table_mccc(data_path, fname);

	init_tables_h2_electronic_diss(h2_v1_bpstate_diss_cs, h2_v1_bpstate_diss_rates, h2_v1_bpstate_diss_indexes,
		1, h2_bp1su_diss_min_nb, verbosity);

	// S1g(X) -> P1u(D), v0,v1 -> dissociative continuum
	fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=0/MCCC-el-H2-D1Pu_DE.X1Sg_vi=0.txt";
	h2_v0_dstate_diss_cs = new cross_section_table_mccc(data_path, fname);

	init_tables_h2_electronic_diss(h2_v0_dstate_diss_cs, h2_v0_dstate_diss_rates, h2_v0_dstate_diss_indexes,
		0, h2_d1pu_diss_min_nb, verbosity);

	fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=1/MCCC-el-H2-D1Pu_DE.X1Sg_vi=1.txt";
	h2_v1_dstate_diss_cs = new cross_section_table_mccc(data_path, fname);

	init_tables_h2_electronic_diss(h2_v1_dstate_diss_cs, h2_v1_dstate_diss_rates, h2_v1_dstate_diss_indexes,
		1, h2_d1pu_diss_min_nb, verbosity);

	// S1g(X) -> S3u(b), v0,v1 -> dissociative continuum
	fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=0/MCCC-el-H2-b3Su_DE.X1Sg_vi=0.txt";
	h2_v0_3bstate_diss_cs = new cross_section_table_mccc(data_path, fname);

	init_tables_h2_electronic_diss(h2_v0_3bstate_diss_cs, h2_v0_3bstate_diss_rates, h2_v0_3bstate_diss_indexes,
		0, h2_b3su_diss_min_nb, verbosity);

	fname = "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=1/MCCC-el-H2-b3Su_DE.X1Sg_vi=1.txt";
	h2_v1_3bstate_diss_cs = new cross_section_table_mccc(data_path, fname);

	init_tables_h2_electronic_diss(h2_v1_3bstate_diss_cs, h2_v1_3bstate_diss_rates, h2_v1_3bstate_diss_indexes,
		1, h2_b3su_diss_min_nb, verbosity);

	// Initialization of HeI excitation cross sections,
	init_tables_hei_electronic_exc(data_path);


	// Initialization of indexes and rates of e-e Coulomb scattering,



	// Electron impact ionization
	h2_ioniz_cs
		= new h2_electron_impact_ionization(verbosity);

	he_ioniz_cs
		= new he_electron_impact_ionization(verbosity);

	h_ioniz_cs
		= new h_electron_impact_ionization(verbosity);

	h2p_ioniz_cs
		= new h2p_electron_impact_ionization(verbosity);

	hep_ioniz_cs
		= new hep_electron_impact_ionization(verbosity);

	// electron impact ionization of atoms and molecules
	verbosity = 1;
#pragma omp parallel 
	{
#pragma omp sections nowait
		{
#pragma omp section
			{
				init_tables_ionization(h2_ioniz_cs, h2_ioniz_rates_tot, h2_ioniz_rates, h2_ioniz_indexes, h2_ioniz_min_nb, verbosity);
			}
#pragma omp section
			{
				init_tables_ionization(he_ioniz_cs, he_ioniz_rates_tot, he_ioniz_rates, he_ioniz_indexes, he_ioniz_min_nb, verbosity);
			}
#pragma omp section
			{
				init_tables_ionization(h_ioniz_cs, h_ioniz_rates_tot, h_ioniz_rates, h_ioniz_indexes, h_ioniz_min_nb, verbosity);
			}
#pragma omp section
			{
				init_tables_ionization(h2p_ioniz_cs, h2p_ioniz_rates_tot, h2p_ioniz_rates, h2p_ioniz_indexes, h2p_ioniz_min_nb, verbosity);
			}
#pragma omp section
			{
				init_tables_ionization(hep_ioniz_cs, hep_ioniz_rates_tot, hep_ioniz_rates, hep_ioniz_indexes, hep_ioniz_min_nb, verbosity);
			}
		}
	}

	h2eq_nb = nb_of_el_energies + NB_OF_CHEM_SPECIES;
	heieq_nb = h2eq_nb + nb_lev_h2;
	physeq_nb = heieq_nb + nb_lev_hei;
	nb_of_equat = nb_of_el_energies + NB_OF_CHEM_SPECIES + nb_lev_h2 + nb_lev_hei + 3;

	el_nb = EL_NB + nb_of_el_energies;
	h_nb = H_NB + nb_of_el_energies;
	hp_nb = H_P_NB + nb_of_el_energies;
	h2_nb = H2_NB + nb_of_el_energies;
	h2p_nb = H2_P_NB + nb_of_el_energies;
	he_nb = HE_NB + nb_of_el_energies;
	hep_nb = HE_P_NB + nb_of_el_energies;
	hepp_nb = HE_PP_NB + nb_of_el_energies;
}

elspectra_evolution_data::~elspectra_evolution_data()
{
	delete elastic_h2_el_cs;
	delete elastic_he_el_cs;

	delete h2_ioniz_cs;
	delete he_ioniz_cs;
	delete h_ioniz_cs;
	delete h2p_ioniz_cs;
	delete hep_ioniz_cs;

	delete rot_h2_j02_cs;
	delete rot_h2_j24_cs;
	delete rot_h2_j13_cs;
	delete rot_h2_j35_cs;

	delete vibr_h2_v01_cs;
	delete vibr_h2_v02_cs;

	delete h2_v0_bstate_diss_cs;
	delete h2_v1_bstate_diss_cs;
	delete h2_v0_cstate_diss_cs;
	delete h2_v1_cstate_diss_cs;
	delete h2_v0_bpstate_diss_cs;
	delete h2_v1_bpstate_diss_cs;
	delete h2_v0_dstate_diss_cs;
	delete h2_v1_dstate_diss_cs;
	delete h2_v0_3bstate_diss_cs;
	delete h2_v1_3bstate_diss_cs;

	for (int i = 0; i < (int)hei_cs.size(); i++) {
		delete hei_cs[i];
	}

	delete[] h2_v0_bstate_cs;
	delete[] h2_v1_bstate_cs;
	delete[] h2_v0_bpstate_cs;
	delete[] h2_v1_bpstate_cs;
	delete[] h2_v0_cstate_cs;
	delete[] h2_v1_cstate_cs;
	delete[] h2_v0_dstate_cs;
	delete[] h2_v1_dstate_cs;

	delete h2_di;
	delete h2_di_b; 
	delete h2_di_cminus; 
	delete h2_di_cplus;
	delete h2_di_bp; 
	delete h2_di_dminus; 
	delete h2_di_dplus;
	delete h2_einst;

	delete[] electron_energies_grid;
	delete[] electron_energy_bin_size;
	delete[] electron_energies;
	delete[] electron_velocities;
	
	delete[] el_enloss_h2_mt;
	delete[] el_enloss_he_mt;
	delete[] elel_scatt_rates;
	delete[] elion_scatt_rates;

	delete[] h2_ioniz_rates_tot;
	delete[] he_ioniz_rates_tot;
	delete[] h_ioniz_rates_tot;
	delete[] h2p_ioniz_rates_tot;
	delete[] hep_ioniz_rates_tot;
	
	delete[] h2_j02_rates; 
	delete[] h2_j24_rates;
	delete[] h2_j13_rates; 
	delete[] h2_j35_rates;

	delete[] h2_j20_rates;
	delete[] h2_j42_rates;
	delete[] h2_j31_rates;
	delete[] h2_j53_rates;

	delete[] h2_v01_rates;
	delete[] h2_v02_rates;

	delete[] h2_v0_bstate_diss_rates;
	delete[] h2_v1_bstate_diss_rates;
	delete[] h2_v0_cstate_diss_rates;
	delete[] h2_v1_cstate_diss_rates;
	delete[] h2_v0_bpstate_diss_rates;
	delete[] h2_v1_bpstate_diss_rates;
	delete[] h2_v0_dstate_diss_rates;
	delete[] h2_v1_dstate_diss_rates;
	delete[] h2_v0_3bstate_diss_rates;
	delete[] h2_v1_3bstate_diss_rates;

	free_2d_array(h2_ioniz_rates);
	free_2d_array(he_ioniz_rates);
	free_2d_array(h_ioniz_rates);
	free_2d_array(h2p_ioniz_rates);
	free_2d_array(hep_ioniz_rates);

	free_2d_array(h2_v0_bstate_rates);
	free_2d_array(h2_v1_bstate_rates);
	free_2d_array(h2_v0_bpstate_rates);
	free_2d_array(h2_v1_bpstate_rates);
	free_2d_array(h2_v0_cstate_rates);
	free_2d_array(h2_v1_cstate_rates);
	free_2d_array(h2_v0_dstate_rates);
	free_2d_array(h2_v1_dstate_rates);

	free_2d_array(hei_rates);

	delete[] h2_j02_indexes; 
	delete[] h2_j13_indexes; 
	delete[] h2_j24_indexes; 
	delete[] h2_j35_indexes;

	delete[] h2_j20_indexes;
	delete[] h2_j31_indexes;
	delete[] h2_j42_indexes;
	delete[] h2_j53_indexes;

	free_2d_array(h2_ioniz_indexes);
	free_2d_array(he_ioniz_indexes);
	free_2d_array(h_ioniz_indexes);
	free_2d_array(h2p_ioniz_indexes);
	free_2d_array(hep_ioniz_indexes);

	free_2d_array(h2_v01_indexes);
	free_2d_array(h2_v02_indexes);

	free_2d_array(h2_v0_bstate_indexes); 
	free_2d_array(h2_v1_bstate_indexes);
	free_2d_array(h2_v0_bpstate_indexes);
	free_2d_array(h2_v1_bpstate_indexes);
	free_2d_array(h2_v0_cstate_indexes);
	free_2d_array(h2_v1_cstate_indexes);
	free_2d_array(h2_v0_dstate_indexes);
	free_2d_array(h2_v1_dstate_indexes);

	free_2d_array(h2_v0_bstate_diss_indexes);
	free_2d_array(h2_v1_bstate_diss_indexes);
	free_2d_array(h2_v0_bpstate_diss_indexes);
	free_2d_array(h2_v1_bpstate_diss_indexes);
	free_2d_array(h2_v0_cstate_diss_indexes);
	free_2d_array(h2_v1_cstate_diss_indexes);
	free_2d_array(h2_v0_dstate_diss_indexes);
	free_2d_array(h2_v1_dstate_diss_indexes);
	free_2d_array(h2_v0_3bstate_diss_indexes);
	free_2d_array(h2_v1_3bstate_diss_indexes);

	free_2d_array(hei_indexes);
}

void elspectra_evolution_data::init_tables_ionization(electron_impact_ionization * cs,
	double *&rates_tot, double **&rates, spectra_data **& indexes, int &ioniz_min_nb, int verbosity)
{
	const double err = 1.e-2;
	int i, j, j0, k, n;
	double en0, y, z, fluct;

	cross_section_integral_2d ion_cs_int_2d(cs, 1.e-6);
	cross_section_integral_weight ion_cs_int_weight(cs, 1.e-6);

	n = 0;
#ifdef _OPENMP
	n = omp_get_thread_num();
#endif

	if (verbosity) {
		cout << left << "Calculation of e-ionization rates, " << cs->get_name() << ", thread nb " << n << endl;
	}
	ioniz_min_nb = get_nb_electron_energy_array(cs->get_binding_energy(), 0);

	rates_tot = new double[nb_of_el_energies];
	memset(rates_tot, 0, nb_of_el_energies * sizeof(double));

	rates = alloc_2d_array<double>(nb_of_el_energies, nb_of_el_energies);
	memset(*rates, 0, nb_of_el_energies * nb_of_el_energies * sizeof(double));

	// constructors of the el_spectra_data must be called,
	indexes = alloc_2d_array<spectra_data>(nb_of_el_energies, nb_of_el_energies);

	// incident electron energy en_i, ejected electron energy en_j, 
	for (i = ioniz_min_nb; i < nb_of_el_energies; i++) {
		// upper bound of the energy of the ejected electron, (t - 1)/2 (limit in the integration)
		en0 = 0.5 * (electron_energies[i] - cs->get_binding_energy());
		// the number of the energy interval where the limiting energy is located,
		j0 = get_nb_electron_energy_array(en0, 0);

		for (j = 0; j <= j0; j++) {
			z = electron_energies[j];
			y = electron_energies[i] - z - cs->get_binding_energy();
			
			k = get_nb_electron_energy_array(y, 0);

			if (k == 0 && y < electron_energies[k]) {
				indexes[i][j].i0 = 0;
				indexes[i][j].w1 = 1.;
			}
			else {
				if (y < electron_energies[k])
					k--;

				indexes[i][j].i0 = k;

				indexes[i][j].w1 = (electron_energies[k + 1] - y) / (electron_energies[k + 1] - electron_energies[k]);
				if (indexes[i][j].w1 > 1.)
					indexes[i][j].w1 = 1.;
			}
			indexes[i][j].w2 = 1. - indexes[i][j].w1;
			
			if (j < j0)
				rates[i][j] = cs->get_int_cs(electron_energies[i], electron_energies_grid[j], electron_energies_grid[j+1]);
			else
				rates[i][j] = cs->get_int_cs(electron_energies[i], electron_energies_grid[j], en0);

			rates_tot[i] += rates[i][j];  // [cm2]
			rates[i][j] *= electron_velocities[i];           // [cm2 cm/s], 
		}

		// check the integrated cross section,
		// the large departure appears when interval length becomes higher than binding energy of electron in atom (molecule),
		fluct = fabs(cs->get_int_cs(electron_energies[i]) / (rates_tot[i] + 1.e-99));
		if (fabs(fluct - 1.) > err && verbosity) {
			cout << left << n << "    total cs (relative) difference for interval " << i << ": " << fluct << endl;
		}
		rates_tot[i] *= electron_velocities[i];  // [cm2 cm/s],
	}
}

void elspectra_evolution_data::init_tables_coulomb(electron_impact_ionization *elel_cs, 
	double*& el_ion_rates_tot, double**& el_ion_rates, el_spectra_data**& el_ion_indexes, int verbosity)
{
}


void elspectra_evolution_data::init_tables_h2_rotational_exc(const cross_section* cs, double *& rates_up, spectra_data*& indexes_up, 
	double*& rates_down, spectra_data*& indexes_down, int j0, int verbosity)
{
	int i, k;
	double z;

	rates_up = new double[nb_of_el_energies];
	indexes_up = new spectra_data[nb_of_el_energies];

	rates_down = new double[nb_of_el_energies];
	indexes_down = new spectra_data[nb_of_el_energies];

	if (verbosity) {
		cout << "Calculation of e-H2 rotational excitation rates," << endl;
	}
	for (i = 0; i < nb_of_el_energies; i++) {
		// the cross section at the centre of the interval i (not averaged over the interval),
		rates_up[i] = (*cs)(electron_energies[i]) * electron_velocities[i]; // [cm2 *cm/s]
		
		z = electron_energies[i] - cs->get_threshold_energy();
		if (z < 0.)
			z = 0.;

		// the energy z must be E_k < z < E_k+1 (E is the median energy in the bin)
		k = get_nb_electron_energy_array(z, 0);
		if (k == 0 && z < electron_energies[k]) {
			indexes_up[i].i0 = 0;
			indexes_up[i].w1 = 1.;
		}
		else {
			if (z < electron_energies[k])
				k--;

			indexes_up[i].i0 = k;
			indexes_up[i].w1 = (electron_energies[k + 1] - z) / (electron_energies[k + 1] - electron_energies[k]);
			if (indexes_up[i].w1 > 1.)
				indexes_up[i].w1 = 1.;	
		}
		indexes_up[i].w2 = 1. - indexes_up[i].w1;
	}
	if (verbosity) {
		cout << "Calculation of e-H2 rotational de-excitation rates," << endl;
	}
	for (i = 0; i < nb_of_el_energies; i++) {
		// calculation of de-excitation cross section:
		// incident electron gains energy in this case,
		rates_down[i] = (*cs)(electron_energies[i] + cs->get_threshold_energy()) * (2. * j0 + 1) / (2. * j0 + 5)
			* (electron_energies[i] + cs->get_threshold_energy())/ electron_energies[i] * electron_velocities[i]; // [cm2 *cm/s]
		
		z = electron_energies[i] + cs->get_threshold_energy();
		
		// the energy z must be E_k < z < E_k+1 (E is the median energy in the bin)
		k = get_nb_electron_energy_array(z, 0);
		if (k == 0 && z < electron_energies[k]) {
			indexes_down[i].i0 = 0;
			indexes_down[i].w1 = 1.;
		}
		else {
			if (z < electron_energies[k])
				k--;

			indexes_down[i].i0 = k;
			indexes_down[i].w1 = (electron_energies[k + 1] - z) / (electron_energies[k + 1] - electron_energies[k]);
			if (indexes_down[i].w1 > 1.)
				indexes_down[i].w1 = 1.;
		}
		indexes_down[i].w2 = 1. - indexes_down[i].w1;
	}
}

void elspectra_evolution_data::init_tables_h2_vibrational_exc(const cross_section* cs, double*& rates, 
	spectra_data**& indexes, int vibr_qnb_final, int verbosity)
{
	int i, j, dj, k,l, n, v;
	double en, z;

	rates = new double[nb_of_el_energies];
	indexes = alloc_2d_array<spectra_data>(3 * MAX_J_H2_VIBR_INIT, nb_of_el_energies);

	for (i = 0; i < nb_of_el_energies; i++) {
		rates[i] = (*cs)(electron_energies[i]) * electron_velocities[i];  // [cm2 *cm/s]
	}
	
	for (j = 0; j < MAX_J_H2_VIBR_INIT; j++) {
		for (dj = -2; dj <= 2; dj += 2) {
			k = h2_di->get_nb(v = 0, j);
			l = h2_di->get_nb(vibr_qnb_final, j + dj);

			if (l != -1) {
				en = (h2_di->lev_array[l].energy - h2_di->lev_array[k].energy)* CM_INVERSE_TO_EV;

				for (i = 0; i < nb_of_el_energies; i++) {
					z = electron_energies[i] - en;
					if (z < 0.)
						z = 0.;

					k = get_nb_electron_energy_array(z, 0);
					n = 3 * j + dj / 2 + 1;

					if (k == 0 && en < electron_energies[k]) {
						indexes[n][i].i0 = 0;
						indexes[n][i].w1 = 1.;
					}
					else {
						if (en < electron_energies[k])
							k--;

						indexes[n][i].i0 = k;
						indexes[n][i].w1 = (electron_energies[k + 1] - en) / (electron_energies[k + 1] - electron_energies[k]);
						if (indexes[n][i].w1 > 1.)
							indexes[n][i].w1 = 1.;
					}
					indexes[n][i].w2 = 1. - indexes[n][i].w1;
				}
			}
		}
	}
}

// vi fixed, vf = [0,1,.., max_h2_vstates-1]
// S1g(X) -> S1u(B), S1u(Bp) dJ = +/-1, electronic states
// S1g(X) -> P1u(C), P1u(D), dJ = -1, 0, 1
void elspectra_evolution_data::init_tables_h2_electronic_exc(const energy_diagram* h2_di_elexc, cross_section_table_mccc** cs, 
	double**& rates, spectra_data**& indexes, int vibr_qnb_init, int max_h2_vstates, int & el_min_nb, int verbosity)
{
	int i, j, l, dj, k, p, vf, nb_of_branches;
	double en, z, en_min;
 
	// eV, ionization energy of H2
	en_min = 15.43; 
	
	// for S(X)->S(B,Bp) states only transitions with dj = -/+1 are allowed,
	if (h2_di_elexc->electronic_state == 1 || h2_di_elexc->electronic_state == 4) {
		nb_of_branches = 2;
	}
	else {
		nb_of_branches = 3;
	}

	// the cross sections are given for vi - > vf	
	rates = alloc_2d_array<double>(max_h2_vstates, nb_of_el_energies);
	indexes = alloc_2d_array<spectra_data>(nb_of_branches * max_h2_vstates * MAX_J_H2_VIBR_INIT, nb_of_el_energies);

	for (i = 0; i < nb_of_el_energies; i++) {
		for (vf = 0; vf < max_h2_vstates; vf++) {
			rates[vf][i] = (*cs[vf])(electron_energies[i]) * electron_velocities[i];  // [cm2 *cm/s]
		}
	}

	// loop over vibrational nb of excited electronic state,
	for (vf = 0; vf < max_h2_vstates; vf++) {
		for (j = 0; j < MAX_J_H2_VIBR_INIT; j++) {
			for (dj = -1; dj <= 1; dj++) {
				if (nb_of_branches == 3 || (nb_of_branches == 2 && dj != 0))
				{
					k = h2_di->get_nb(vibr_qnb_init, j);
					l = h2_di_elexc->get_nb(vf, j + dj);

					if (l != -1) {
						en = (h2_di_elexc->lev_array[l].energy - h2_di->lev_array[k].energy) * CM_INVERSE_TO_EV;
						if (en < en_min)
							en_min = en;

						for (i = 0; i < nb_of_el_energies; i++) {
							z = electron_energies[i] - en;
							if (z < 0.)
								z = 0.;

							k = get_nb_electron_energy_array(z, 0);
							p = nb_of_branches * (vf * MAX_J_H2_VIBR_INIT + j) + dj + 1;
							
							if (k == 0 && z < electron_energies[k]) {
								indexes[p][i].i0 = 0;
								indexes[p][i].w1 = 1.;
							}
							else {
								if (z < electron_energies[k])
									k--;

								indexes[p][i].i0 = k;
								indexes[p][i].w1 = (electron_energies[k + 1] - z) / (electron_energies[k + 1] - electron_energies[k]);
								if (indexes[p][i].w1 > 1.)
									indexes[p][i].w1 = 1.;
							}
							indexes[p][i].w2 = 1. - indexes[p][i].w1;
						}
					}
				}
			}
		}
	}
	el_min_nb = get_nb_electron_energy_array(en_min, 0);
}

void elspectra_evolution_data::init_tables_h2_electronic_diss(cross_section_table_mccc* cs, double*& rates, 
	spectra_data**& indexes, int vibr_qnb_init, int& el_min_nb, int verbosity)
{
	int i, j, k, l;
	double z, en, en_min;

	rates = new double[nb_of_el_energies];
	indexes = alloc_2d_array<spectra_data>(MAX_J_H2_VIBR_INIT, nb_of_el_energies);

	if (verbosity) {
		cout << "Calculation of H2 dissociation rates through electronic excitation," << endl;
	}
	
	for (i = 0; i < nb_of_el_energies; i++) {
		// the cross section at the centre of the interval i (not averaged over the interval),
		rates[i] = (*cs)(electron_energies[i]) * electron_velocities[i]; // [cm2 *cm/s]
	}

	// = E(excited state, vibration continuum) - E(X, initial vibr qnb, J=0)
	en_min = cs->get_threshold_energy();

	for (j = 0; j < MAX_J_H2_VIBR_INIT; j++) {
		l = h2_di->get_nb(vibr_qnb_init, 0);
		k = h2_di->get_nb(vibr_qnb_init, j);

		// energy loss by the electron at the collision,
		en = cs->get_threshold_energy() 
			- (h2_di->lev_array[k].energy - h2_di->lev_array[l].energy) * CM_INVERSE_TO_EV;

		if (en < en_min)
			en_min = en;

		for (i = 0; i < nb_of_el_energies; i++) {
			z = electron_energies_grid[i] - en;
			if (z < 0.)
				z = 0.;

			k = get_nb_electron_energy_array(z, 0);
			if (k == 0 && z < electron_energies[k]) {
				indexes[j][i].i0 = 0;
				indexes[j][i].w1 = 1.;
			}
			else {
				if (z < electron_energies[k])
					k--;

				indexes[j][i].i0 = k;
				indexes[j][i].w1 = (electron_energies[k + 1] - z) / (electron_energies[k + 1] - electron_energies[k]);
				if (indexes[j][i].w1 > 1.)
					indexes[j][i].w1 = 1.;
			}
			indexes[j][i].w2 = 1. - indexes[j][i].w1;
		}
	}
	el_min_nb = get_nb_electron_energy_array(en_min, 0);
}


void elspectra_evolution_data::init_tables_hei_electronic_exc(const string& data_path)
{
	int i, m, j, l, s, k, i_low, i_fin, g;
	double a1, a2, a3, a4, a5, a6, en_thr, en_min, en, z;
	char text_line[MAX_TEXT_LINE_WIDTH];

	stringstream ss;
	ifstream input;
	string fname, conf, term, name, type;

	fname = data_path + "coll_he/coll_he_e_1s2_1S.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with cross section data " << fname << endl;
		exit(1);
	}

	hei_cs.clear();
	hei_level_nbs.clear();

	while (!input.eof()) {
		// comment lines are read:
		do {
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		} while (text_line[0] == '#');

		if (text_line[0] == '\0') // check for empty line at the file end;
			break;

		ss.clear();
		ss.str(text_line);

		ss >> type;
		ss >> conf >> term >> l >> s >> j;  // angular momentum and spin are integers for HeI
		name = conf + term;
		
		g = 2 * j + 1;   // statistical weight of the initial state,
		i_low = hei_di->get_nb(name, g);

		ss >> conf >> term >> l >> s >> j;
		name = conf + term;
		i_fin = hei_di->get_nb(name, 2 * j + 1);
		
		if (i_low == 0 && i_fin != -1) {   // only transitions from the ground level are considered,
			for (i = 0; i < 6; i++) {
				ss >> a1 >> a2 >> a3 >> a4 >> a5 >> a6;
			}
			en_thr = (hei_di->lev_array[i_fin].energy - hei_di->lev_array[i_low].energy) * CM_INVERSE_TO_EV;
			if (type == "DA") {
				hei_cs.push_back(new hei_electron_excitation_dipole_allowed(a1, a2, a3, a4, a5, a6, en_thr, g));
				hei_level_nbs.push_back(i_fin);
			}
			else if (type == "DF") {
				hei_cs.push_back(new hei_electron_excitation_dipole_forbidden(a1, a2, a3, a4, a5, a6, en_thr, g));
				hei_level_nbs.push_back(i_fin);
			}
			else if (type == "SF") {
				// relative weight of the final level
				// in ion level object total angular momentum is saved in variable spin,
				z = (2. * hei_di->lev_array[i_fin].spin + 1.) / ((2.*l + 1.)*(2.*s + 1.));
				a1 *= z;
				a2 *= z;
				a3 *= z;
				a4 *= z;

				hei_cs.push_back(new hei_electron_excitation_spin_forbidden(a1, a2, a3, a4, a5, a6, en_thr, g));
				hei_level_nbs.push_back(i_fin);
			}
		}
	}
	input.close();

	nb_coll_trans_hei = (int)hei_level_nbs.size();
	hei_rates = alloc_2d_array<double>(nb_coll_trans_hei, nb_of_el_energies);
	hei_indexes = alloc_2d_array<spectra_data>(nb_coll_trans_hei, nb_of_el_energies);
	
	en_min = 24.59;  // ionization potential of HeI in eV;

	for (m = 0; m < nb_coll_trans_hei; m++) {
		for (i = 0; i < nb_of_el_energies; i++) {
			// the cross section at the centre of the interval i (not averaged over the interval),
			hei_rates[m][i] = (*(hei_cs[m]))(electron_energies[i]) * electron_velocities[i];  // [cm2 *cm/s]
		}

		// energy loss by the electron at the collision,
		en = hei_cs[m]->get_threshold_energy();
		if (en < en_min)
			en_min = en;

		for (i = 0; i < nb_of_el_energies; i++) {
			z = electron_energies_grid[i] - en;
			if (z < 0.)
				z = 0.;

			k = get_nb_electron_energy_array(z, 0);
			if (k == 0 && z < electron_energies[k]) {
				hei_indexes[m][i].i0 = 0;
				hei_indexes[m][i].w1 = 1.;
			}
			else {
				if (z < electron_energies[k])
					k--;

				hei_indexes[m][i].i0 = k;
				hei_indexes[m][i].w1 = (electron_energies[k + 1] - z) / (electron_energies[k + 1] - electron_energies[k]);
				if (hei_indexes[m][i].w1 > 1.)
					hei_indexes[m][i].w1 = 1.;
			}
			hei_indexes[m][i].w2 = 1. - hei_indexes[m][i].w1;
		}
	}
	hei_min_nb = get_nb_electron_energy_array(en_min, 0);
}

void elspectra_evolution_data::set_dust_parameters(bool dip, double r, double nbd)
{
	dust_is_presented = dip;
	grain_cs = M_PI * r * r;
	grain_nb_density = nbd;
	
	// According to Weingartner & Draine, ApJSS 134, 263 (2001)
	// for silicate grains, 1 A  = 1.e-8 cm, radius in cm,
	min_grain_charge = (int)(-(2.5 + 0.07 * 1.e+8 * r) * r * 1.e+8 / 14.4) + 1;
}

double elspectra_evolution_data::get_electron_energy(int i) const {
	return electron_energies[i];  // returns the centre of the interval [i, i+1]
}

double elspectra_evolution_data::get_electron_energy_bin(int i) const {
	return electron_energy_bin_size[i];
}

int elspectra_evolution_data::get_nb_electron_energy_array(double energy, int i) const
{
	hunt_index(electron_energies_grid, nb_of_el_energies + 1, energy, i);
	if (i < 0)
		i = 0;
	else if (i > nb_of_el_energies - 1)
		i = nb_of_el_energies - 1;
	return i;
}

void elspectra_evolution_data::get_el_energy_losses(double& mt, double& h2_rot, double &h2_vibr, double &h2_electr, double &h2_electr_tr,
	double& ion, double& col_ions, double& col_el, double& hei) const
{
	mt = enloss_deriv_mt;
	h2_rot = enloss_deriv_h2_rot;
	h2_vibr = enloss_deriv_h2_vibr;
	h2_electr = enloss_deriv_h2_electr;        // via excitation of singlet H2 states (including dissociation)
	h2_electr_tr = enloss_deriv_h2_electr_tr;  // triplet H2 states (only dissociation yet)
	ion = enloss_deriv_ioniz;
	col_ions = enloss_deriv_coloumb_ions;
	col_el = enloss_deriv_coloumb_el;
	hei = enloss_deriv_hei;
}

void elspectra_evolution_data::get_h2_process_rates(double& sol_diss, double& diss_exc, double& diss_exc_tr, double & hei_exc)
{
	sol_diss = h2_solomon_diss_rate;
	diss_exc = h2_diss_exc_rate;
	diss_exc_tr = h2_diss_exc_triplet_rate;
	hei_exc = hei_exc_rate;
}

int elspectra_evolution_data::f(realtype t, N_Vector y, N_Vector ydot)
{
	int i, j;
	double x, z, rate, log_coulomb, thermal_energy, vel_i;

	realtype* y_data, * ydot_data;
	y_data = NV_DATA_S(y);
	ydot_data = NV_DATA_S(ydot);

	for (i = 0; i < nb_of_equat; i++) {
		ydot_data[i] = 0.;
	}

	log_coulomb = 20.;
	energy_gain_n = energy_gain_i = 0.;

	// The derivative of the electron spectrum, 
	// the spectrum is the nb of electrons in the energy bin per cm3, [cm-3]
	// momentum transfer in collisions with neutral species:
	enloss_deriv_mt = 0.;
	for (i = 1; i < nb_of_el_energies; i++) {
		// [cm2 eV cm/s eV-1] * [cm-3] *[cm-3]
		x = (el_enloss_h2_mt[i]* y_data[h2_nb] + el_enloss_he_mt[i] * y_data[he_nb])* y_data[i];

		ydot_data[i] -= x;
		ydot_data[i-1] += x;

		// [cm-3 s-1 eV], electrons lose energy, < 0
		enloss_deriv_mt -= x * (electron_energies[i] - electron_energies[i - 1]);
	}
	energy_gain_n -= EV_TO_ERGS * enloss_deriv_mt;  // [erg cm-3 s-1], neutrals gain energy in this process, must be > 0.

	// Coulomb losses
	// Electron scattering on free electrons,
	enloss_deriv_coloumb_el = 0.;
#if CALC_COLOUMB_EL_LOSSES	
	for (i = 0; i < nb_of_el_energies; i++) {
		x = 0.;
		for (j = 0; j < i; j++) {
			x -= elel_scatt_rates[i] * y_data[j]; // [eV cm^3 s^-1] * [cm-3]
		}

		for (j = i + 1; j < nb_of_el_energies; j++) {
			x += elel_scatt_rates[j] * y_data[j];
		}
		
		if (i > 0 && x < 0.) {
			// [eV s^-1] * [cm^-3 eV^-1]
			x *= 2. * y_data[i] * log_coulomb / (electron_energy_bin_size[i] + electron_energy_bin_size[i - 1]);  
			
			ydot_data[i] += x;  // x < 0
			ydot_data[i - 1] -= x;
			enloss_deriv_coloumb_el -= x * (electron_energies[i] - electron_energies[i - 1]);
		}
		if (i < nb_of_el_energies - 1 && x > 0.) {
			x *= 2. * y_data[i] * log_coulomb / (electron_energy_bin_size[i] + electron_energy_bin_size[i + 1]);
			
			ydot_data[i] -= x;
			ydot_data[i + 1] += x;
			enloss_deriv_coloumb_el += x * (electron_energies[i + 1] - electron_energies[i]);
		}
	}
#endif

	// Electron scattering on ions,
	/*z = log_coulomb * (0.5 * y_data[h2p_nb] + 0.25 * y_data[hep_nb] + y_data[hepp_nb] + y_data[hp_nb]);
	enloss_deriv_coloumb_ions = 0.;
	for (i = 1; i < nb_of_el_energies; i++) {
		// rates are normalized on electron energy bin size,
		x = y_data[i] * elion_scatt_rates[i] * z;
		ydot_data[i] -= x;
		ydot_data[i - 1] += x;

		// [cm-3 s-1 eV], < 0.
		enloss_deriv_coloumb_ions -= x * (electron_energies[i] - electron_energies[i - 1]);
	}*/

	enloss_deriv_coloumb_ions = 0.;
	for (i = 0; i < nb_of_el_energies; i++) {
		// H+
		vel_i = sqrt(2. * BOLTZMANN_CONSTANT * y_data[physeq_nb + 1] / ATOMIC_MASS_UNIT);
		x = electron_velocities[i] / vel_i;
		z = y_data[hp_nb] / vel_i 
			* (0.5 * SQRT_PI * error_function.f(x) / x - (1. + ATOMIC_MASS_UNIT / ELECTRON_MASS) * exp(-x * x));
		
		// H2+
		vel_i = sqrt(BOLTZMANN_CONSTANT * y_data[physeq_nb + 1] / ATOMIC_MASS_UNIT);
		x = electron_velocities[i] / vel_i;
		z += 0.5 * y_data[h2p_nb] / vel_i
			* (0.5 * SQRT_PI * error_function.f(x) / x - (1. + 2. * ATOMIC_MASS_UNIT / ELECTRON_MASS) * exp(-x * x));

		// He+
		vel_i = sqrt(0.5 * BOLTZMANN_CONSTANT * y_data[physeq_nb + 1] / ATOMIC_MASS_UNIT);
		x = electron_velocities[i] / vel_i;
		z += 0.25 * y_data[hep_nb] / vel_i
			* (0.5 * SQRT_PI * error_function.f(x) / x - (1. + 4.*ATOMIC_MASS_UNIT / ELECTRON_MASS) * exp(-x * x));

		// He++
		vel_i = sqrt(0.5 * BOLTZMANN_CONSTANT * y_data[physeq_nb + 1] / ATOMIC_MASS_UNIT);
		x = electron_velocities[i] / vel_i;
		z += y_data[hepp_nb] / vel_i
			* (0.5 * SQRT_PI * error_function.f(x) / x - (1. + 4. * ATOMIC_MASS_UNIT / ELECTRON_MASS) * exp(-x * x));

		z *= -log_coulomb * el_ion_scatt_const * y_data[i];   // [eV cm^4 s^-2] [cm-3 cm-1 s] [cm-3] = [eV cm-3 s-1]
		
		if (i > 0 && z < 0.) {
			z /= 0.5 * (electron_energy_bin_size[i] + electron_energy_bin_size[i - 1]);
			ydot_data[i] += z;
			ydot_data[i - 1] -= z;
			
			// [cm-3 s-1 eV], < 0
			enloss_deriv_coloumb_ions += z*(electron_energies[i] - electron_energies[i-1]);
		}
		if (i < nb_of_el_energies - 1 && z > 0.) {
			z /= 0.5 * (electron_energy_bin_size[i] + electron_energy_bin_size[i + 1]);
			ydot_data[i] -= z;
			ydot_data[i + 1] += z;

			// [cm-3 s-1 eV], > 0
			enloss_deriv_coloumb_ions += z * (electron_energies[i + 1] - electron_energies[i]);
		}
	}
	energy_gain_i -= EV_TO_ERGS * enloss_deriv_coloumb_ions;  // [erg cm-3 s-1], opposite sign (electrons loose energy, ions gain),
	
	// Electron impact ionization of H2, He, H, H2+, He+
	enloss_deriv_ioniz = 0.;
	derivatives_ionization(y_data, ydot_data, h2_ioniz_cs, h2_ioniz_rates, h2_ioniz_rates_tot, h2_ioniz_indexes, 
		h2_ioniz_min_nb, h2_nb, rate, enloss_deriv_ioniz);

	ydot_data[h2_nb] -= rate;
	ydot_data[h2p_nb] += rate;
	ydot_data[el_nb] += rate;

	rate /= y_data[h2_nb];
	for (i = 0; i < nb_lev_h2; i++) {
		ydot_data[h2eq_nb + i] -= y_data[h2eq_nb + i] * rate;
	}

	derivatives_ionization(y_data, ydot_data, he_ioniz_cs, he_ioniz_rates, he_ioniz_rates_tot, he_ioniz_indexes,
		he_ioniz_min_nb, he_nb, rate, enloss_deriv_ioniz);

	ydot_data[he_nb] -= rate;
	ydot_data[hep_nb] += rate;
	ydot_data[el_nb] += rate;

	rate /= y_data[he_nb];
	for (i = 0; i < nb_lev_hei; i++) {
		ydot_data[heieq_nb + i] -= y_data[heieq_nb + i] * rate;
	}

	derivatives_ionization(y_data, ydot_data, h_ioniz_cs, h_ioniz_rates, h_ioniz_rates_tot, h_ioniz_indexes,
		h_ioniz_min_nb, h_nb, rate, enloss_deriv_ioniz);

	ydot_data[h_nb] -= rate;
	ydot_data[hp_nb] += rate;
	ydot_data[el_nb] += rate;

	derivatives_ionization(y_data, ydot_data, h2p_ioniz_cs, h2p_ioniz_rates, h2p_ioniz_rates_tot, h2p_ioniz_indexes,
		h2p_ioniz_min_nb, h2p_nb, rate, enloss_deriv_ioniz);

	ydot_data[h2p_nb] -= rate;
	ydot_data[hp_nb] += 2.*rate;
	ydot_data[el_nb] += rate;

	derivatives_ionization(y_data, ydot_data, hep_ioniz_cs, hep_ioniz_rates, hep_ioniz_rates_tot, hep_ioniz_indexes,
		hep_ioniz_min_nb, hep_nb, rate, enloss_deriv_ioniz);

	ydot_data[hep_nb] -= rate;
	ydot_data[hepp_nb] += rate;
	ydot_data[el_nb] += rate;

	// H2 rotational excitation, lower level quantum nb j0 must be given,
	enloss_deriv_h2_rot = 0.;
	derivatives_h2_rotational_exc(y_data, ydot_data, h2_j02_rates, h2_j02_indexes, h2_j20_rates, h2_j20_indexes, 0, enloss_deriv_h2_rot);
	derivatives_h2_rotational_exc(y_data, ydot_data, h2_j13_rates, h2_j13_indexes, h2_j31_rates, h2_j31_indexes, 1, enloss_deriv_h2_rot);
	derivatives_h2_rotational_exc(y_data, ydot_data, h2_j24_rates, h2_j24_indexes, h2_j42_rates, h2_j42_indexes, 2, enloss_deriv_h2_rot);
	derivatives_h2_rotational_exc(y_data, ydot_data, h2_j35_rates, h2_j35_indexes, h2_j53_rates, h2_j53_indexes, 3, enloss_deriv_h2_rot);

	// H2 vibrational excitation,
	enloss_deriv_h2_vibr = 0.;
	derivatives_h2_vibrational_exc(y_data, ydot_data, h2_v01_rates, h2_v01_indexes, 1, enloss_deriv_h2_vibr);
	derivatives_h2_vibrational_exc(y_data, ydot_data, h2_v02_rates, h2_v02_indexes, 2, enloss_deriv_h2_vibr);

	// H2 electronic states excitation,
	// energy loss in [eV cm-3 s-1], rate in [cm-3 s-1]
	enloss_deriv_h2_electr = h2_solomon_diss_rate = h2_diss_exc_rate = thermal_energy = 0.;
	
	// S1u(X) -> S1u(B)
	derivatives_h2_electronic_exc(y_data, ydot_data, h2_di_b, h2_b_state_data, h2_v0_bstate_rates, h2_v0_bstate_indexes, 
		0, MAX_H2_VSTATES_B1SU, h2_b1su_min_nb, enloss_deriv_h2_electr, h2_solomon_diss_rate, thermal_energy);

	derivatives_h2_electronic_exc(y_data, ydot_data, h2_di_b, h2_b_state_data, h2_v1_bstate_rates, h2_v1_bstate_indexes, 
		1, MAX_H2_VSTATES_B1SU, h2_b1su_min_nb, enloss_deriv_h2_electr, h2_solomon_diss_rate, thermal_energy);

	// S1u(X) -> S1u(Bp)
	derivatives_h2_electronic_exc(y_data, ydot_data, h2_di_bp, h2_bp_state_data, h2_v0_bpstate_rates, h2_v0_bpstate_indexes,
		0, MAX_H2_VSTATES_BP1SU, h2_bp1su_min_nb, enloss_deriv_h2_electr, h2_solomon_diss_rate, thermal_energy);

	derivatives_h2_electronic_exc(y_data, ydot_data, h2_di_bp, h2_bp_state_data, h2_v1_bpstate_rates, h2_v1_bpstate_indexes,
		1, MAX_H2_VSTATES_BP1SU, h2_bp1su_min_nb, enloss_deriv_h2_electr, h2_solomon_diss_rate, thermal_energy);

	// S1u(X) -> P1u(C-/+)
	derivatives_h2_electronic_exc_2(y_data, ydot_data, h2_di_cplus, h2_di_cminus, h2_cplus_state_data, h2_cminus_state_data,
		h2_v0_cstate_rates, h2_v0_cstate_indexes, 0, MAX_H2_VSTATES_C1PU, h2_c1pu_min_nb, enloss_deriv_h2_electr, h2_solomon_diss_rate, thermal_energy);

	derivatives_h2_electronic_exc_2(y_data, ydot_data, h2_di_cplus, h2_di_cminus, h2_cplus_state_data, h2_cminus_state_data,
		h2_v1_cstate_rates, h2_v1_cstate_indexes, 1, MAX_H2_VSTATES_C1PU, h2_c1pu_min_nb, enloss_deriv_h2_electr, h2_solomon_diss_rate, thermal_energy);

	// S1u(X) -> P1u(D-/+)
	derivatives_h2_electronic_exc_2(y_data, ydot_data, h2_di_dplus, h2_di_dminus, h2_dplus_state_data, h2_dminus_state_data,
		h2_v0_dstate_rates, h2_v0_dstate_indexes, 0, MAX_H2_VSTATES_D1PU, h2_d1pu_min_nb, enloss_deriv_h2_electr, h2_solomon_diss_rate, thermal_energy);

	derivatives_h2_electronic_exc_2(y_data, ydot_data, h2_di_dplus, h2_di_dminus, h2_dplus_state_data, h2_dminus_state_data,
		h2_v1_dstate_rates, h2_v1_dstate_indexes, 1, MAX_H2_VSTATES_D1PU, h2_d1pu_min_nb, enloss_deriv_h2_electr, h2_solomon_diss_rate, thermal_energy);


	//  dissociative excitation through S1u(B) state
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v0_bstate_diss_rates, h2_v0_bstate_diss_indexes,
		0, h2_b1su_diss_min_nb, enloss_deriv_h2_electr, h2_diss_exc_rate);

	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v1_bstate_diss_rates, h2_v1_bstate_diss_indexes,
		1, h2_b1su_diss_min_nb, enloss_deriv_h2_electr, h2_diss_exc_rate);

	// through P1u(C) state
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v0_cstate_diss_rates, h2_v0_cstate_diss_indexes,
		0, h2_c1pu_diss_min_nb, enloss_deriv_h2_electr, h2_diss_exc_rate);

	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v1_cstate_diss_rates, h2_v1_cstate_diss_indexes,
		1, h2_c1pu_diss_min_nb, enloss_deriv_h2_electr, h2_diss_exc_rate);

	// through S1u(Bp) state
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v0_bpstate_diss_rates, h2_v0_bpstate_diss_indexes,
		0, h2_bp1su_diss_min_nb, enloss_deriv_h2_electr, h2_diss_exc_rate);

	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v1_bpstate_diss_rates, h2_v1_bpstate_diss_indexes,
		1, h2_bp1su_diss_min_nb, enloss_deriv_h2_electr, h2_diss_exc_rate);

	// through P1u(D) state
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v0_dstate_diss_rates, h2_v0_dstate_diss_indexes,
		0, h2_d1pu_diss_min_nb, enloss_deriv_h2_electr, h2_diss_exc_rate);

	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v1_dstate_diss_rates, h2_v1_dstate_diss_indexes,
		1, h2_d1pu_diss_min_nb, enloss_deriv_h2_electr, h2_diss_exc_rate);

	// through s3u(b) state (triplet)
	h2_diss_exc_triplet_rate = enloss_deriv_h2_electr_tr = 0.;
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v0_3bstate_diss_rates, h2_v0_3bstate_diss_indexes,
		0, h2_b3su_diss_min_nb, enloss_deriv_h2_electr_tr, h2_diss_exc_triplet_rate);

	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v1_3bstate_diss_rates, h2_v1_3bstate_diss_indexes,
		1, h2_b3su_diss_min_nb, enloss_deriv_h2_electr_tr, h2_diss_exc_triplet_rate);

	// the decrease of population densities of H2 is taken into account in the functions above,
	x = h2_solomon_diss_rate + h2_diss_exc_rate + h2_diss_exc_triplet_rate;
	ydot_data[h2_nb] -= x;
	ydot_data[h_nb] += 2. * x;

	energy_gain_n += thermal_energy;  // in erg,

	// Radiative transitions in H2 (only radiative decay)
	for (i = 1; i < nb_lev_h2; i++) {
		for (j = 0; j < i; j++) {
			if (h2_einst->arr[i][j] > 1.e-99)
			{
				x = h2_einst->arr[i][j] * y_data[h2eq_nb + i];
				ydot_data[h2eq_nb + j] += x;
				ydot_data[h2eq_nb + i] -= x;
			}
		}
	}

	// Excitation of HeI by electron impact
	enloss_deriv_hei = hei_exc_rate = 0.;
	derivatives_hei_exc(y_data, ydot_data, enloss_deriv_hei, hei_exc_rate);

	// Radiative transitions in HeI (only radiative decay),
	// the effect of He radiation is discussed in the book by B.T. Draine, Physics of the Interstellar and Intergalactic Medium (2011)
    // Transitions leading to 1s2s 1S will be followed by two-photon decay with A=51.0 s-1 (Drake 1986), 
	// and a total photon energy 20.62 eV; 56 % of these two photon decays produce a photon with hv > 13.6 eV.
	// optical depth to absorption, to H2/H ionization?
	for (i = 1; i < nb_lev_hei; i++) {
		for (j = 0; j < i; j++) {
			if (hei_einst->arr[i][j] > 1.e-99)
			{
				x = hei_einst->arr[i][j] * y_data[heieq_nb + i];
				ydot_data[heieq_nb + j] += x;
				ydot_data[heieq_nb + i] -= x;
			}
		}
	}
	
	// Momentum transfer between neutrals and ions,

	// Gas temperature:
	conc_n = y_data[h2_nb] + y_data[he_nb] + y_data[h_nb];
	conc_i = y_data[h2p_nb] + y_data[hep_nb] + y_data[hepp_nb] + y_data[hp_nb];
	
	nb_gain_n = ydot_data[h2_nb] + ydot_data[he_nb] + ydot_data[h_nb];
	nb_gain_i = ydot_data[h2p_nb] + ydot_data[hep_nb] + ydot_data[hepp_nb] + ydot_data[hp_nb];

	// energy gain must include gain of the particles,
	energy_gain_n += 1.5 * nb_gain_n * y_data[physeq_nb + 1] * BOLTZMANN_CONSTANT;  // [erg cm-3 s-1]
	energy_gain_i += 1.5 * nb_gain_i * y_data[physeq_nb] * BOLTZMANN_CONSTANT;

	ydot_data[physeq_nb] = (2. * energy_gain_n / (3. * BOLTZMANN_CONSTANT) - y_data[physeq_nb] * nb_gain_n) / conc_n;
	ydot_data[physeq_nb + 1] = (2. * energy_gain_i / (3. * BOLTZMANN_CONSTANT) - y_data[physeq_nb + 1] * nb_gain_i) / conc_i;
	
	// Dust grain charge
	// Electron escape length is 10 A (1 A = 1e-8 cm) for electron energies 10-40 eV (Weingartner & Draine, ApJSS 134, 263, 2001),
	if (dust_is_presented) {
		for (i = 0; i < nb_of_el_energies; i++) {
			x = y_data[physeq_nb + 2] / min_grain_charge;
			x = exp(-x * x);
			ydot_data[physeq_nb + 2] -= ELECTRON_DUST_SCATTERING_PROB *x* electron_velocities[i] * grain_cs * y_data[i];
			ydot_data[i] -= ELECTRON_DUST_SCATTERING_PROB *x* electron_velocities[i] * grain_cs * grain_nb_density * y_data[i];
		}
	}
	
	return 0;
}

void elspectra_evolution_data::derivatives_ionization(const realtype* y_data, realtype* ydot_data, 
	const electron_impact_ionization *cs, double ** rates, const double * rates_tot, spectra_data ** indexes, 
	int el_ion_min_nb, int target_nb, double& rate, double & enloss_deriv)
{
	int i, j, j0, l;
	double enl, x, y, r;
	
	enl = r = 0.;
#pragma omp parallel reduction(+: enl, r) private(i, j, j0, l, x, y)
	{
		double* arr_el = new double[nb_of_el_energies];
		memset(arr_el, 0, nb_of_el_energies * sizeof(double));

		// incident electron en_i -> scattered electron en_j
		// Note: more simple calculation of energy losses?
#pragma omp for schedule(dynamic, 1)
		for (i = el_ion_min_nb; i < nb_of_el_energies; i++) {
			x = 0.5 * (electron_energies[i] - cs->get_binding_energy());
			j0 = get_nb_electron_energy_array(x, 0);
			y = y_data[i] * y_data[target_nb];

			for (j = 0; j <= j0; j++) {
				if (rates[i][j] > 1.e-99) {
					const spectra_data& esd = indexes[i][j];

					x = rates[i][j] * y;  // [cm2 cm/s] *[cm-3 cm-3]

					arr_el[i] -= x;
					arr_el[j] += x;  // ejected electron

					enl += x * (electron_energies[j] - electron_energies[i]);

					// scattered electron:
					l = esd.i0;
					arr_el[l] += indexes[i][j].w1 * x;
					arr_el[l + 1] += indexes[i][j].w2 * x;

					enl += (esd.w1 * electron_energies[l] + esd.w2 * electron_energies[l + 1]) * x;  // [cm-3 s-1 eV],
				}
			}
			r += rates_tot[i] * y;
		}
#pragma omp critical
		{
			for (i = 0; i < nb_of_el_energies; i++) {
				ydot_data[i] += arr_el[i];
			}
		}
		delete[] arr_el;
	}
	enloss_deriv += enl;
	rate = r;
}

void elspectra_evolution_data::derivatives_h2_rotational_exc(const realtype* y_data, realtype* ydot_data, 
	const double* rates_up, const spectra_data* indexes_up, const double* rates_down, const spectra_data* indexes_down, int j0, double & enloss_deriv)
{
	int i;
	double x, t(0.);

	for (i = 0; i < nb_of_el_energies; i++) {
		// j0 -> j0 + 2
		const spectra_data& esd1 = indexes_up[i];
		 
		// [cm2 *cm/s] *[cm-3] *[cm-3]
		x = rates_up[i] * y_data[i] * y_data[h2eq_nb + j0];
		t += x;

		ydot_data[esd1.i0] += x * esd1.w1;
		ydot_data[esd1.i0 + 1] += x * esd1.w2;
		
		enloss_deriv += (esd1.w1 *electron_energies[esd1.i0] + esd1.w2 * electron_energies[esd1.i0 + 1] 
			- electron_energies[i]) *x;  // [cm-3 s-1 eV],
		
		ydot_data[i] -= x;

		// j0 + 2 -> j0
		const spectra_data& esd2 = indexes_down[i];

		// [cm2 *cm/s] *[cm-3] *[cm-3]
		x = rates_down[i] * y_data[i] * y_data[h2eq_nb + j0 + 2];
		t -= x;

		ydot_data[esd2.i0] += x * esd2.w1;
		ydot_data[esd2.i0 + 1] += x * esd2.w2;

		enloss_deriv += (esd2.w1 * electron_energies[esd2.i0] + esd2.w2 * electron_energies[esd2.i0 + 1]
			- electron_energies[i]) * x;  // [cm-3 s-1 eV],

		ydot_data[i] -= x;
	}
	ydot_data[h2eq_nb + j0] -= t;
	ydot_data[h2eq_nb + j0 + 2] += t;
}

void elspectra_evolution_data::derivatives_h2_vibrational_exc(const realtype* y_data, realtype* ydot_data, 
	const double* rates, spectra_data** indexes, int vibr_qnb_final, double & enloss_deriv)
{
	int i, j, dj, v, n, low, up, up1, up3;
	double enl, x, y, f;

	enl = 0.;
#pragma omp parallel reduction(+: enl) private(i, j, dj, v, n, low, up, up1, up3, x, y, f)
	{
		double* arr_h2 = new double[nb_lev_h2];
		memset(arr_h2, 0, nb_lev_h2 * sizeof(double));

		double* arr_el = new double[nb_of_el_energies];
		memset(arr_el, 0, nb_of_el_energies * sizeof(double));

#pragma omp for schedule(dynamic, 1)
		for (j = 0; j < MAX_J_H2_VIBR_INIT; j++) {
			low = h2_di->get_nb(v = 0, j);
			up1 = h2_di->get_nb(vibr_qnb_final, j - 2);
			up3 = h2_di->get_nb(vibr_qnb_final, j + 2);

			if (up1 != -1 && up3 != -1)
				x = h2_einst->arr[low][up3] / h2_einst->arr[low][up1];

			for (dj = -2; dj <= 2; dj += 2) {
				up = h2_di->get_nb(vibr_qnb_final, j + dj);

				if (up != -1) {
					n = 3 * j + (dj / 2 + 1);

					if (dj == 0) {
						if (up1 != -1 || up3 != -1)
							f = 0.5;
						else
							f = 1.;
					}
					else if (up1 == -1 || up3 == -1) {
						f = 0.5;
					}
					else {
						if (dj == -2)
							f = 0.5 / (x + 1.);
						else
							f = 0.5 * x / (x + 1.);
					}
					y = 0.;
					// the H2 excitation by collisions with electrons with the lowest energy are excluded, 
					for (i = 1; i < nb_of_el_energies; i++) {
						const spectra_data& esd = indexes[n][i];

						// [cm2 *cm/s] *[cm-3] *[cm-3]
						x = f * rates[i] * y_data[i] * y_data[h2eq_nb + low];
						
						arr_el[esd.i0] += x * esd.w1;
						arr_el[esd.i0 + 1] += x * esd.w2;

						enl += (esd.w1 * electron_energies[esd.i0] + esd.w2 * electron_energies[esd.i0 + 1]
							- electron_energies[i]) * x;  // [cm-3 s-1 eV],

						arr_el[i] -= x;
						y += x;
					}
					arr_h2[low] -= y;
					arr_h2[up] += y;
				}
			}
		}
#pragma omp critical
		{
			for (i = 0; i < nb_of_el_energies; i++) {
				ydot_data[i] += arr_el[i];
			}
			for (i = 0; i < nb_lev_h2; i++) {
				ydot_data[h2eq_nb + i] += arr_h2[i];
			}
		}
		delete[] arr_h2;
		delete[] arr_el;
	}
	enloss_deriv += enl;
}


// S1g(X) -> S1u(B,Bp), dJ = +/-1,  
// vi fixed, vf = [0,1,.., MAX_H2_VSTATES_B1SU-1]
void elspectra_evolution_data::derivatives_h2_electronic_exc(const realtype* y_data, realtype* ydot_data, 
	const energy_diagram* h2_di_exc, const vector<h2_energy_level_param> h2_exc_state_data, double** rates, spectra_data** indexes, 
	int vibr_qnb_init, int vibr_qnb_final_max, int el_min_nb, double& enloss_deriv, double& diss_rate, double& thermal_energy)
{
	int i, k, j, dj, vf, p, low, up;
	double dissr, enl, th_en, x, y, f;

	dissr = enl = th_en = 0.;
#pragma omp parallel reduction(+: dissr, enl, th_en) private(i, k, dj, j, vf, p, low, up, y, x, f)
	{
		double* arr_h2 = new double[nb_lev_h2];
		memset(arr_h2, 0, nb_lev_h2 * sizeof(double));

		double* arr_el = new double[nb_of_el_energies];
		memset(arr_el, 0, nb_of_el_energies * sizeof(double));

#pragma omp for schedule(dynamic, 1)
		for (vf = 0; vf < vibr_qnb_final_max; vf++) {
			for (j = 0; j < MAX_J_H2_VIBR_INIT; j++) 
			{
				low = h2_di->get_nb(vibr_qnb_init, j);
				for (dj = -1; dj <= 1; dj += 2) 
				{
					up = h2_di_exc->get_nb(vf, j + dj);
					if (up != -1) {
						p = 2 * (vf * MAX_J_H2_VIBR_INIT + j) + dj + 1;

						// the calculation of Honl-London factors
						// for S1g(X) -> S1u(B, Bp), j -> j + dj, upward
						f = hl_singlet_dl0(j, dj) * y_data[h2eq_nb + low];
						y = 0.;

						for (i = el_min_nb; i < nb_of_el_energies; i++) {
							const spectra_data& esd = indexes[p][i];

							// [cm2 *cm/s] *[cm-3] *[cm-3]
							x = f * rates[vf][i] * y_data[i];
							
							arr_el[esd.i0] += x * esd.w1;
							arr_el[esd.i0 + 1] += x * esd.w2;

							enl += (esd.w1 * electron_energies[esd.i0] + esd.w2 * electron_energies[esd.i0 + 1]
								- electron_energies[i]) * x;  // [cm-3 s-1 eV],

							arr_el[i] -= x;
							y += x;
						}
						
						// the upper level of the band transition in question (belongs to the excited electronic state):
						const h2_energy_level_param& lev_param_ref = h2_exc_state_data[up];

						arr_h2[low] -= y;
						for (k = 0; k < (int)(lev_param_ref.nb_of_decays); k++) {
							i = lev_param_ref.decay_level_nbs[k];
							arr_h2[i] += y * lev_param_ref.decay_probs[k];
						}
						// the downward transition (dissociative) to the vibrational continuum, [cm-3 s-1]
						dissr += y * lev_param_ref.diss_prob;
						th_en += y * lev_param_ref.diss_prob * lev_param_ref.kin_energy; // in erg
					}
				}
			}
		}
#pragma omp critical
		{
			for (i = 0; i < nb_of_el_energies; i++) {
				ydot_data[i] += arr_el[i];
			}
			for (i = 0; i < nb_lev_h2; i++) {
				ydot_data[h2eq_nb + i] += arr_h2[i];
			}
		}
		delete[] arr_h2;
		delete[] arr_el;
	}
	diss_rate += dissr;
	enloss_deriv += enl;
	thermal_energy += th_en;
}

// the collisions must not violate the ortho/para state!
// S1g(X) -> P1u(C,D), C- dJ = 0, C+ dJ = +/-1, Abgrall et al., A&A 253, 525 (1992)
void elspectra_evolution_data::derivatives_h2_electronic_exc_2(const realtype* y_data, realtype* ydot_data,
	const energy_diagram* h2_di_plus, const energy_diagram* h2_di_minus,
	const vector<h2_energy_level_param> h2_state_plus_data, const vector<h2_energy_level_param> h2_state_minus_data,
	double** rates, spectra_data** indexes, int vibr_qnb_init, int vibr_qnb_final_max, int el_min_nb, 
	double& enloss_deriv, double& diss_rate, double& thermal_energy)
{
	int i, k, j, dj, vf, p, low, up;
	double dissr, enl, th_en, x, y, f;

	dissr = enl = th_en = 0.;
#pragma omp parallel reduction(+: dissr, enl, th_en) private(i, k, j, dj, vf, p, low, up, x, y, f)
	{
		double* arr_h2 = new double[nb_lev_h2];
		memset(arr_h2, 0, nb_lev_h2 * sizeof(double));

		double* arr_el = new double[nb_of_el_energies];
		memset(arr_el, 0, nb_of_el_energies * sizeof(double));

#pragma omp for schedule(dynamic, 1)
		for (vf = 0; vf < vibr_qnb_final_max; vf++) {
			for (j = 0; j < MAX_J_H2_VIBR_INIT; j++) {
				low = h2_di->get_nb(vibr_qnb_init, j);

				for (dj = -1; dj <= 1; dj++) {
					if (dj == 0)
						up = h2_di_minus->get_nb(vf, j + dj);
					else
						up = h2_di_plus->get_nb(vf, j + dj);

					if (up != -1) {
						p = 3 * (vf * MAX_J_H2_VIBR_INIT + j) + dj + 1;

						// the calculation of Honl-London factors
						f = hl_singlet_dl1(j, dj)* y_data[h2eq_nb + low];
						y = 0.;

						for (i = el_min_nb; i < nb_of_el_energies; i++) {
							const spectra_data& esd = indexes[p][i];

							// [cm2 *cm/s] *[cm-3] *[cm-3]
							x = f * rates[vf][i] * y_data[i];
							
							arr_el[esd.i0] += x * esd.w1;
							arr_el[esd.i0 + 1] += x * esd.w2;

							enl += (esd.w1 * electron_energies[esd.i0] + esd.w2 * electron_energies[esd.i0 + 1]
								- electron_energies[i]) * x;  // [cm-3 s-1 eV],

							arr_el[i] -= x;
							y += x;
						}

						arr_h2[low] -= y;
						if (dj == 0) {
							// the upper level of the band transition in question (belongs to the excited electronic state):
							const h2_energy_level_param& lev_param_ref = h2_state_minus_data[up];

							for (k = 0; k < (int)(lev_param_ref.nb_of_decays); k++) {
								i = lev_param_ref.decay_level_nbs[k];
								arr_h2[i] += y * lev_param_ref.decay_probs[k];
							}
							// the downward transition (dissociative) to the vibrational continuum, [cm-3 s-1]
							dissr += y * lev_param_ref.diss_prob;
							th_en += y * lev_param_ref.diss_prob * lev_param_ref.kin_energy; // in erg
						}
						else {
							const h2_energy_level_param& lev_param_ref = h2_state_plus_data[up];

							for (k = 0; k < (int)(lev_param_ref.nb_of_decays); k++) {
								i = lev_param_ref.decay_level_nbs[k];
								arr_h2[i] += y * lev_param_ref.decay_probs[k];
							}
							// the downward transition (dissociative) to the vibrational continuum, [cm-3 s-1]
							dissr += y * lev_param_ref.diss_prob;
							th_en += y * lev_param_ref.diss_prob * lev_param_ref.kin_energy; // in erg
						}
					}
				}
			}
		}
#pragma omp critical
		{
			for (i = 0; i < nb_of_el_energies; i++) {
				ydot_data[i] += arr_el[i];
			}
			for (i = 0; i < nb_lev_h2; i++) {
				ydot_data[h2eq_nb + i] += arr_h2[i];
			}
		}
		delete[] arr_h2;
		delete[] arr_el;
	}
	diss_rate += dissr;
	enloss_deriv += enl;
	thermal_energy += th_en;
}

void elspectra_evolution_data::derivatives_h2_electronic_diss(const realtype* y_data, realtype* ydot_data, 
	const double* rates, spectra_data** indexes, int vibr_qnb_init, int el_min_nb, double& enloss_deriv, double& diss_rate)
{
	int i, j, k;
	double x, y, enl(0.), dissr(0.);

#pragma omp parallel reduction(+: dissr, enl) private(i, k, j, x, y)
	{
		double* arr_h2 = new double[nb_lev_h2];
		memset(arr_h2, 0, nb_lev_h2 * sizeof(double));

		double* arr_el = new double[nb_of_el_energies];
		memset(arr_el, 0, nb_of_el_energies * sizeof(double));

#pragma omp for schedule(dynamic, 1)
		for (j = 0; j < MAX_J_H2_VIBR_INIT; j++) {
			k = h2_di->get_nb(vibr_qnb_init, j);
			y = 0.;

			for (i = el_min_nb; i < nb_of_el_energies; i++) {
				const spectra_data& esd = indexes[j][i];

				// [cm2 *cm/s] *[cm-3] *[cm-3]
				x = rates[i] * y_data[i] * y_data[h2eq_nb + k];
				
				arr_el[esd.i0] += x * esd.w1;
				arr_el[esd.i0 + 1] += x * esd.w2;

				enl += (esd.w1 * electron_energies[esd.i0] + esd.w2 * electron_energies[esd.i0 + 1]
					- electron_energies[i]) * x;  // [cm-3 s-1 eV],

				arr_el[i] -= x;
				y += x;
			}
			arr_h2[k] -= y;
			dissr += y;
		}
#pragma omp critical
		{
			for (i = 0; i < nb_of_el_energies; i++) {
				ydot_data[i] += arr_el[i];
			}
			for (i = 0; i < nb_lev_h2; i++) {
				ydot_data[h2eq_nb + i] += arr_h2[i];
			}
		}
		delete[] arr_h2;
		delete[] arr_el;
	}
	enloss_deriv += enl;
	diss_rate += dissr;
}

void elspectra_evolution_data::derivatives_hei_exc(const realtype* y_data, realtype* ydot_data, double& enloss_deriv, double& exc_rate)
{
	int i, m, k;
	double x, y, enl(0.), excr(0.);

#pragma omp parallel reduction(+: excr, enl) private(i, k, m, x, y)
	{
		double* arr_hei = new double[nb_lev_hei];
		memset(arr_hei, 0, nb_lev_hei * sizeof(double));

		double* arr_el = new double[nb_of_el_energies];
		memset(arr_el, 0, nb_of_el_energies * sizeof(double));

#pragma omp for schedule(dynamic, 1)
		for (m = 0; m < nb_coll_trans_hei; m++) {
			k = hei_level_nbs[m];
			y = 0.;

			for (i = hei_min_nb; i < nb_of_el_energies; i++) {
				const spectra_data& esd = hei_indexes[m][i];

				// excitation from lower level of HeI is considered
				// [cm2 *cm/s] *[cm-3] *[cm-3]
				x = hei_rates[m][i] * y_data[i] * y_data[heieq_nb];
				
				arr_el[esd.i0] += x * esd.w1;
				arr_el[esd.i0 + 1] += x * esd.w2;

				enl += (esd.w1 * electron_energies[esd.i0] + esd.w2 * electron_energies[esd.i0 + 1]
					- electron_energies[i]) * x;  // [cm-3 s-1 eV],

				arr_el[i] -= x;
				y += x;
			}
			arr_hei[0] -= y;
			arr_hei[k] += y;
			excr += y;
		}
#pragma omp critical
		{
			for (i = 0; i < nb_of_el_energies; i++) {
				ydot_data[i] += arr_el[i];
			}
			for (i = 0; i < nb_lev_hei; i++) {
				ydot_data[heieq_nb + i] += arr_hei[i];
			}
		}
		delete[] arr_hei;
		delete[] arr_el;
	}
	enloss_deriv += enl;
	exc_rate += excr;
}


void elspectra_evolution_data::set_tolerances(N_Vector abs_tol)
{
	int i, i0;
	for (i = 0; i < nb_of_el_energies; i++) {
		NV_Ith_S(abs_tol, i) = ABS_ELSPECTRA_ERROR_SOLVER;
	}
	for (i = nb_of_el_energies; i < nb_of_el_energies + NB_OF_CHEM_SPECIES; i++) {
		NV_Ith_S(abs_tol, i) = ABS_CONCENTRATION_ERROR_SOLVER;
	}
	i0 = nb_of_el_energies + NB_OF_CHEM_SPECIES;
	for (i = i0; i < i0 + nb_lev_h2; i++) {
		NV_Ith_S(abs_tol, i) = ABS_POPULATION_H2_ERROR_SOLVER;
	}
	for (i = i0 + nb_lev_h2; i < nb_of_equat; i++) {
		NV_Ith_S(abs_tol, i) = ABS_PARAMETER_ERROR_SOLVER;
	}
}
