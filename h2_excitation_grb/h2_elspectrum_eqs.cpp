#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include<omp.h>
#include<cmath>
#include<limits>
#include<cstring>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<sstream>

#include "interpolation.h"
#include "h2_elspectrum_eqs.h"
#include "h2_microphysics.h"

#include "constants.h"
#include "spectroscopy.h"
#include "h2_parameters.h"

using namespace std;
#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "h2_elspectrum_eqs.cpp"


const double log_coulomb_constant = log(1.5 * pow(BOLTZMANN_CONSTANT, 1.5)
	/ (sqrt(M_PI) * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE)); // *T_e^1.5/n_e^0.5

const double el_ion_scatt_const = 8. * SQRT_PI * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE 
	/(ATOMIC_MASS_UNIT * EV_TO_ERGS);   // [eV cm^4 s^-2]


elspectra_evolution_data::elspectra_evolution_data(const string& data_path, double cht, const vector<double>& en_grid, int verb)
	: conc_h_tot(cht), dust_is_presented(false), grain_cs(0.), grain_nb_density(0.),
	enloss_rate_mt(0.), enloss_rate_h2_rot(0.), enloss_rate_h2_vibr(0.), enloss_rate_h2_singlet(0.), enloss_diss_h2_triplet(0.),
	enloss_rate_ioniz(0.), enloss_rate_coloumb_ions(0.), enloss_rate_coloumb_el(0.), enloss_rate_hei(0.), 
	conc_n(0.), conc_i(0.), energy_gain_n(0.), energy_gain_i(0.), nb_gain_n(0.), nb_gain_i(0.), 
	h2_solomon_diss_rate(0.), h2_diss_exc_singlet_rate(0.), h2_diss_exc_triplet_rate(0.), hei_exc_rate(0.)
{
	int i, j, dj, li, lf, vi, vf, isotope, electronic_state, verbosity = 1;
	double en, den, vel;
	
	string fname, name;
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
	
	// elel_scatt_rates = new double[nb_of_el_energies];
	// elion_scatt_rates = new double[nb_of_el_energies];

	el_enloss_h2_mt[0] = el_enloss_he_mt[0] = 0.; 
	// elel_scatt_rates[0] = elion_scatt_rates[0] = 0.;

	for (i = 1; i < nb_of_el_energies; i++) {
		// it is assumed that electron energy is the kinetic energy, in eV
		en = electron_energies[i];  // the centre of the interval,
		vel = electron_velocities[i];
		// 
		den = 0.5 * (electron_energy_bin_size[i] + electron_energy_bin_size[i - 1]);
		
		// energy losses are normalized on target (H2, H, He) concentration and energy bin size, [cm2 eV cm/s eV-1]
		// m_h2 = 2 a.m.u., m_he = 4 a.m.u.
		el_enloss_h2_mt[i] = (*elastic_h2_el_cs)(en) *en *vel * ELECTRON_MASS/(ATOMIC_MASS_UNIT * den); 
		el_enloss_he_mt[i] = (*elastic_he_el_cs)(en) *en * vel * 0.5* ELECTRON_MASS / (ATOMIC_MASS_UNIT *den); 

		// electron-electron collisions,
		// energy loss rate coefficient, [eV cm^3 s^-1]
		// elel_scatt_rates[i] = 4.*M_PI * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE 
		//	/ (ELECTRON_MASS * vel * EV_TO_ERGS);

		// electron-ion collisions,
		// rates are normalized on electron energy bin size, [eV cm^3 s^-1 eV-1]
		//elion_scatt_rates[i] = 4. * M_PI * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE
		//	/ (ATOMIC_MASS_UNIT * vel * EV_TO_ERGS * den);
	}

	// Electron impact ionization,
	h2_electron_ionization_data h2_ioniz_data;
	he_electron_ionization_data he_ioniz_data;
	h_electron_ionization_data  h_ioniz_data;
	h2p_electron_ionization_data h2p_ioniz_data;
	hep_electron_ionization_data hep_ioniz_data;

	h2_ioniz_cs
		= new electron_impact_ionization(h2_ioniz_data, verbosity);

	he_ioniz_cs
		= new electron_impact_ionization(he_ioniz_data, verbosity);

	h_ioniz_cs
		= new electron_impact_ionization(h_ioniz_data, verbosity);

	h2p_ioniz_cs
		= new electron_impact_ionization(h2p_ioniz_data, verbosity);

	hep_ioniz_cs
		= new electron_impact_ionization(hep_ioniz_data, verbosity);

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

	// H2 molecule data:
	nb_lev_h2 = 301;	
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

	// Pure rotational excitation of H2 levels,
	h2_rot_cs = new cross_section_table_mccc * [2 * NB_OF_H2_VSTATES_X1SU * MAX_J_H2];
	memset(h2_rot_cs, 0, 2 * NB_OF_H2_VSTATES_X1SU * MAX_J_H2 * sizeof(h2_rot_cs[0]));

	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		for (j = 0; j < MAX_J_H2; j++) {
			for (dj = -2; dj <= 2; dj += 4) 
			{
				li = h2_di->get_nb(vi, j);
				lf = h2_di->get_nb(vi, j + dj);

				if (li >= 0 && lf >= 0) {
					ss.clear();
					ss.str("");
					ss << "coll_h2/MCCC-el-H2-rot-X-X/vi=";
					ss << vi;
					ss << "/Ji=";
					ss << j;
					ss << "/MCCC-el-H2-X1Sg_vf=";
					ss << vi;
					ss << "_Jf=";
					ss << j + dj;
					ss << ".X1Sg_vi=";
					ss << vi;
					ss << "_Ji=";
					ss << j;
					ss << ".txt";
					
					name = data_path + ss.str();
					ifstream f(name.c_str());
					
					// check if file with data exists,  
					if (f.good()) {
						fname = ss.str();
						
						i = 2 * (vi * MAX_J_H2 + j) + (dj + 2) / 4;
						h2_rot_cs[i] = new cross_section_table_mccc(data_path, fname);
					}
					else if (verbosity) {
						cout << "There is no data for H2 transition (v,J)->(v,J): " << vi << " " << j << " " << vi << " " << j + dj << endl;
					}
				}
			}
		}
	}
	init_tables_h2_rotational_exc(h2_rot_cs, h2_rot_rates, verbosity);

	//
	// Ro-vibrational transitions of H2 within ground electronic state,
	// there are 3 transitions for each initial vibration nb, for each j in this state, to NB_OF_H2_VSTATES_X1SU-1 final vibration states,
	h2_rovibr_cs = new cross_section_table_mccc * [3 * NB_OF_H2_VSTATES_X1SU * (NB_OF_H2_VSTATES_X1SU-1) * MAX_J_H2];
	memset(h2_rovibr_cs, 0, 3 * NB_OF_H2_VSTATES_X1SU * (NB_OF_H2_VSTATES_X1SU - 1) * MAX_J_H2 * sizeof(h2_rovibr_cs[0]));
	
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		for (j = 0; j < MAX_J_H2; j++) {
			for (vf = 0; vf < NB_OF_H2_VSTATES_X1SU; vf++) {
				if (vi != vf) {
					for (dj = -2; dj <= 2; dj += 2)
					{
						li = h2_di->get_nb(vi, j);
						lf = h2_di->get_nb(vi, j + dj);

						if (li >= 0 && lf >= 0) {
							ss.clear();
							ss.str("");
							ss << "coll_h2/MCCC-el-H2-rot-X-X/vi=";
							ss << vi;
							ss << "/Ji=";
							ss << j;
							ss << "/MCCC-el-H2-X1Sg_vf=";
							ss << vf;
							ss << "_Jf=";
							ss << j + dj;
							ss << ".X1Sg_vi=";
							ss << vi;
							ss << "_Ji=";
							ss << j;
							ss << ".txt";

							name = data_path + ss.str();
							ifstream f(name.c_str());

							// check if file with data exists,  
							if (f.good()) {
								fname = ss.str();

								i = 3 * (vi * (NB_OF_H2_VSTATES_X1SU - 1) * MAX_J_H2 + vf * MAX_J_H2 + j) + (dj + 2) / 2;
								if (vf > vi)
									i -= 3 * MAX_J_H2;

								h2_rovibr_cs[i] = new cross_section_table_mccc(data_path, fname);
							}
							else if (verbosity) {
								cout << "There is no data for H2 transition (v,J)->(v,J): " << vi << " " << j << " " << vf << " " << j + dj << endl;
							}
						}
					}
				}
			}
		}
	}
	init_tables_h2_rovibr_exc(h2_rovibr_cs, h2_rovibr_rates, verbosity);

	//
	// Excitation of electronically excited states of H2,
	// S1g+(X) -> S1u+(B),  
	// Note, 0 =< vi < NB_OF_H2_VSTATES_X1SU,
	h2_bstate_cs = new cross_section_table_mccc *[NB_OF_H2_VSTATES_X1SU * MAX_H2_VSTATES_B1SU];
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		for (vf = 0; vf < MAX_H2_VSTATES_B1SU; vf++) {
			ss.clear();
			ss.str("");
			
			ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
			ss << vi;
			ss << "/MCCC-el-H2-B1Su_vf=";
			ss << vf;
			ss << ".X1Sg_vi=";
			ss << vi;
			ss << ".txt";
			
			// check for the existence of the file?..,
			fname = ss.str();
			i = vi * MAX_H2_VSTATES_B1SU + vf;
			h2_bstate_cs[i] = new cross_section_table_mccc(data_path, fname);
		}
	}
	init_tables_h2_electronic_exc(h2_bstate_cs, h2_bstate_rates, MAX_H2_VSTATES_B1SU, verbosity);
	h2_b1su_min_nb = calc_min_energy_of_transition(h2_di_b, MAX_H2_VSTATES_B1SU);

	// S1g+(X) -> P1u(C-/+)
	h2_cstate_cs = new cross_section_table_mccc * [NB_OF_H2_VSTATES_X1SU * MAX_H2_VSTATES_C1PU];
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		for (vf = 0; vf < MAX_H2_VSTATES_C1PU; vf++) {
			ss.clear();
			ss.str("");
			
			ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
			ss << vi;
			ss << "/MCCC-el-H2-C1Pu_vf=";
			ss << vf;
			ss << ".X1Sg_vi=";
			ss << vi;
			ss << ".txt";

			fname = ss.str();
			i = vi * MAX_H2_VSTATES_C1PU + vf;
			h2_cstate_cs[i] = new cross_section_table_mccc(data_path, fname);
		}
	}
	init_tables_h2_electronic_exc(h2_cstate_cs, h2_cstate_rates, MAX_H2_VSTATES_C1PU, verbosity);
	
	// the level energy difference of C- and C+ is neglected, C- state have higher number of levels,
	h2_c1pu_min_nb = calc_min_energy_of_transition(h2_di_cminus, MAX_H2_VSTATES_C1PU);

	// S1g+(X) -> S1u+(Bp)
	h2_bpstate_cs = new cross_section_table_mccc * [NB_OF_H2_VSTATES_X1SU * MAX_H2_VSTATES_BP1SU];
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		for (vf = 0; vf < MAX_H2_VSTATES_BP1SU; vf++) {
			ss.clear();
			ss.str("");
			
			ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
			ss << vi;
			ss << "/MCCC-el-H2-Bp1Su_vf=";
			ss << vf;
			ss << ".X1Sg_vi=";
			ss << vi;
			ss << ".txt";

			fname = ss.str();
			i = vi * MAX_H2_VSTATES_BP1SU + vf;
			h2_bpstate_cs[i] = new cross_section_table_mccc(data_path, fname);
		}
	}
	init_tables_h2_electronic_exc(h2_bpstate_cs, h2_bpstate_rates, MAX_H2_VSTATES_BP1SU, verbosity);
	h2_bp1su_min_nb = calc_min_energy_of_transition(h2_di_bp, MAX_H2_VSTATES_BP1SU);

	// S1g+(X) -> P1u(D-/+)
	h2_dstate_cs = new cross_section_table_mccc * [NB_OF_H2_VSTATES_X1SU * MAX_H2_VSTATES_D1PU];
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		for (vf = 0; vf < MAX_H2_VSTATES_D1PU; vf++) {
			ss.clear();
			ss.str("");
			ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
			ss << vi;
			ss << "/MCCC-el-H2-D1Pu_vf=";
			ss << vf;
			ss << ".X1Sg_vi=";
			ss << vi;
			ss << ".txt";

			fname = ss.str();
			i = vi * MAX_H2_VSTATES_D1PU + vf;
			h2_dstate_cs[i] = new cross_section_table_mccc(data_path, fname);
		}
	}
	init_tables_h2_electronic_exc(h2_dstate_cs, h2_dstate_rates, MAX_H2_VSTATES_D1PU, verbosity);
	
	// the level energy difference of D- and D+ is neglected, D- state have much higher number of levels,
	h2_d1pu_min_nb = calc_min_energy_of_transition(h2_di_dminus, MAX_H2_VSTATES_D1PU);

	// Electronic dissociative excitation of H2
	// S1g+(X) -> S1u+(B),
	// vi -> dissociation, 0 =< vi < NB_OF_H2_VSTATES_X1SU,
	h2_bstate_diss_cs = new cross_section_table_mccc * [NB_OF_H2_VSTATES_X1SU];
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		ss.clear();
		ss.str("");
		
		ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
		ss << vi;
		ss << "/MCCC-el-H2-B1Su_DE.X1Sg_vi=";
		ss << vi;
		ss << ".txt";

		fname = ss.str();
		h2_bstate_diss_cs[vi] = new cross_section_table_mccc(data_path, fname);
	}
	init_tables_h2_electronic_diss(h2_bstate_diss_cs, h2_bstate_diss_rates, h2_b1su_diss_min_nb, verbosity);

	// S1g+(X) -> P1u(C),
	h2_cstate_diss_cs = new cross_section_table_mccc * [NB_OF_H2_VSTATES_X1SU];
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		ss.clear();
		ss.str("");
		
		ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
		ss << vi;
		ss << "/MCCC-el-H2-C1Pu_DE.X1Sg_vi=";
		ss << vi;
		ss << ".txt";

		fname = ss.str();
		h2_cstate_diss_cs[vi] = new cross_section_table_mccc(data_path, fname);
	}
	init_tables_h2_electronic_diss(h2_cstate_diss_cs, h2_cstate_diss_rates, h2_c1pu_diss_min_nb, verbosity);

	// S1g+(X) -> S1u+(Bp),
	h2_bpstate_diss_cs = new cross_section_table_mccc * [NB_OF_H2_VSTATES_X1SU];
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		ss.clear();
		ss.str("");
		
		ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
		ss << vi;
		ss << "/MCCC-el-H2-Bp1Su_DE.X1Sg_vi=";
		ss << vi;
		ss << ".txt";

		fname = ss.str();
		h2_bpstate_diss_cs[vi] = new cross_section_table_mccc(data_path, fname);
	}
	init_tables_h2_electronic_diss(h2_bpstate_diss_cs, h2_bpstate_diss_rates, h2_bp1su_diss_min_nb, verbosity);

	// S1g+(X) -> P1u(D),
	h2_dstate_diss_cs = new cross_section_table_mccc * [NB_OF_H2_VSTATES_X1SU];
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		ss.clear();
		ss.str("");
		
		ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
		ss << vi;
		ss << "/MCCC-el-H2-D1Pu_DE.X1Sg_vi=";
		ss << vi;
		ss << ".txt";

		fname = ss.str();
		h2_dstate_diss_cs[vi] = new cross_section_table_mccc(data_path, fname);
	}
	init_tables_h2_electronic_diss(h2_dstate_diss_cs, h2_dstate_diss_rates, h2_d1pu_diss_min_nb, verbosity);

	// S1g(X) -> S3u+(b) (dissociative state),
	// An important impact of these cross sections in the 7-12 eV energy interval (Padovani et al., A&A 658, A 189, 2022, sections 2.1,2.2),
	h2_3bstate_diss_cs = new cross_section_table_mccc * [NB_OF_H2_VSTATES_X1SU];
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		ss.clear();
		ss.str("");
		
		ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
		ss << vi;
		ss << "/MCCC-el-H2-b3Su_DE.X1Sg_vi=";
		ss << vi;
		ss << ".txt";

		fname = ss.str();
		h2_3bstate_diss_cs[vi] = new cross_section_table_mccc(data_path, fname);
	}
	init_tables_h2_electronic_diss(h2_3bstate_diss_cs, h2_3bstate_diss_rates, h2_b3su_diss_min_nb, verbosity);

	// Dissociative excitation to other levels? e.g., a3Sg

	// Initialization of HeI excitation cross sections,
	init_tables_hei_electronic_exc(data_path);


	h2eq_nb = nb_of_el_energies + NB_OF_CHEM_SPECIES;
	heieq_nb = h2eq_nb + nb_lev_h2;
	physeq_nb = heieq_nb + nb_lev_hei;
	nb_of_equat = nb_of_el_energies + NB_OF_CHEM_SPECIES + nb_lev_h2 + nb_lev_hei + 3;

	el_nb = nb_of_el_energies + EL_NB;
	h_nb  = nb_of_el_energies + H_NB;
	hp_nb = nb_of_el_energies + H_P_NB;
	h2_nb = nb_of_el_energies + H2_NB;
	h2p_nb = nb_of_el_energies + H2_P_NB;
	he_nb = nb_of_el_energies + HE_NB;
	hep_nb = nb_of_el_energies + HE_P_NB;
	hepp_nb = nb_of_el_energies + HE_PP_NB;
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

	delete[] h2_bstate_cs;
	delete[] h2_bpstate_cs;
	delete[] h2_cstate_cs;
	delete[] h2_dstate_cs;
	
	delete[] h2_bstate_diss_cs;
	delete[] h2_cstate_diss_cs;
	delete[] h2_bpstate_diss_cs;
	delete[] h2_dstate_diss_cs;
	delete[] h2_3bstate_diss_cs;

	delete[] h2_rot_cs;
	delete[] h2_rovibr_cs;
	
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
	
	delete[] h2_ioniz_rates_tot;
	delete[] he_ioniz_rates_tot;
	delete[] h_ioniz_rates_tot;
	delete[] h2p_ioniz_rates_tot;
	delete[] hep_ioniz_rates_tot;

	free_2d_array(h2_ioniz_rates);
	free_2d_array(he_ioniz_rates);
	free_2d_array(h_ioniz_rates);
	free_2d_array(h2p_ioniz_rates);
	free_2d_array(hep_ioniz_rates);

	free_2d_array(h2_ioniz_indexes);
	free_2d_array(he_ioniz_indexes);
	free_2d_array(h_ioniz_indexes);
	free_2d_array(h2p_ioniz_indexes);
	free_2d_array(hep_ioniz_indexes);

	free_2d_array(h2_bstate_rates);
	free_2d_array(h2_cstate_rates);
	free_2d_array(h2_bpstate_rates);
	free_2d_array(h2_dstate_rates);

	free_2d_array(h2_bstate_diss_rates);
	free_2d_array(h2_cstate_diss_rates);
	free_2d_array(h2_bpstate_diss_rates);
	free_2d_array(h2_dstate_diss_rates);
	free_2d_array(h2_3bstate_diss_rates);

	free_2d_array(h2_rot_rates);
	free_2d_array(h2_rovibr_rates);
}

void elspectra_evolution_data::init_tables_ionization(electron_impact_ionization * cs,
	double *&rates_tot, double **&rates, spectra_data **& indexes, int &ioniz_min_nb, int verbosity)
{
	const double err = 1.e-2;
	int i, j, j_max, k;
	double en_max, y, fluct;

	if (verbosity) {
		cout << left << "Calculation of e-ionization rates, " << cs->get_name() << endl 
			<< "	(simple method, without averaging over energy intervals), " << endl;
	}

	// nb of the projectile energy interval which contains E_b - target electron binding energy,
	ioniz_min_nb = get_nb_electron_energy_array(cs->get_binding_energy());
	ioniz_min_nb++;  // the minimal nb of energy interval, that lies entirely above the ionization threshold, 

	rates_tot = new double[nb_of_el_energies];
	memset(rates_tot, 0, nb_of_el_energies * sizeof(double));

	rates = alloc_2d_array<double>(nb_of_el_energies, nb_of_el_energies);
	memset(*rates, 0, nb_of_el_energies * nb_of_el_energies * sizeof(double));

	indexes = alloc_2d_array<spectra_data>(nb_of_el_energies, nb_of_el_energies);

	// incident electron energy en_i, ejected electron energy en_j, 
	for (i = ioniz_min_nb; i < nb_of_el_energies; i++) {
		// upper bound of the energy of the ejected electron, (t - 1)/2 (limit in the integration)
		en_max = 0.5 * (electron_energies[i] - cs->get_binding_energy());
		// the number of the energy interval where the limiting energy is located,
		j_max = get_nb_electron_energy_array(en_max);

		for (j = 0; j <= j_max; j++) {
			y = electron_energies[i] - electron_energies[j] - cs->get_binding_energy();
			k = get_nb_electron_energy_array(y);

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
			
			if (j < j_max)
				rates[i][j] = cs->get_int_cs(electron_energies[i], electron_energies_grid[j], electron_energies_grid[j+1]);
			else
				rates[i][j] = cs->get_int_cs(electron_energies[i], electron_energies_grid[j], en_max);

			rates_tot[i] += rates[i][j];  // [cm2]	
		}

		// check the integrated cross section,
		fluct = fabs(cs->get_int_cs(electron_energies[i]) / (rates_tot[i] + 1.e-99));
		if (fabs(fluct - 1.) > err && verbosity) {
			cout << left << "    total cs (relative) difference for interval, accurate/in table " << i << ": " << fluct << endl;
		}
		
		for (j = 0; j <= j_max; j++) {
			rates[i][j] *= electron_velocities[i];  // [cm2 cm/s], 
		}
		rates_tot[i] *= electron_velocities[i];  // [cm2 cm/s],
	}
}


int elspectra_evolution_data::calc_min_energy_of_transition(const energy_diagram* h2_di_elexc, int vibr_qnb_final_max)
{
	int i, j, l, vf, dj, el_min_nb;
	double en, en_min;

	// eV, ionization energy of H2
	en_min = 15.43;

	for (i = 0; i < nb_lev_h2; i++) {
		for (vf = 0; vf < vibr_qnb_final_max; vf++) {
			for (dj = -1; dj <= 1; dj++)
			{
				j = rounding(h2_di->lev_array[i].j);
				l = h2_di_elexc->get_nb(vf, j + dj);

				if (l != -1) {
					en = (h2_di_elexc->lev_array[l].energy - h2_di->lev_array[i].energy) * CM_INVERSE_TO_EV;
					if (en < en_min)
						en_min = en;
				}
			}
		}
	}
	el_min_nb = get_nb_electron_energy_array(en_min);
	return el_min_nb;
}


void elspectra_evolution_data::init_tables_h2_electronic_exc(cross_section_table_mccc** cs, double**& rates, 
	int vibr_qnb_final_max, int verbosity)
{
	int i, l, vi, vf;
	
	// the cross sections are given for vi -> vf	
	rates = alloc_2d_array<double>(NB_OF_H2_VSTATES_X1SU * vibr_qnb_final_max, nb_of_el_energies);
	memset(*rates, 0, NB_OF_H2_VSTATES_X1SU * vibr_qnb_final_max * nb_of_el_energies * sizeof(rates[0][0]));

	// pumping by electron impact is taken into account for H2 vibrational quanta of the ground electronic state vi < NB_OF_H2_VSTATES_X1SU
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		for (vf = 0; vf < vibr_qnb_final_max; vf++) {
			for (i = 0; i < nb_of_el_energies; i++) 
			{
				l = vi * vibr_qnb_final_max + vf;
				// entire electron energy bin must be higher than the threshold energy,
				if ((*cs[vi])(electron_energies_grid[i]) > 0.)
					rates[l][i] = (*cs[l])(electron_energies[i]) * electron_velocities[i];  // [cm2 *cm/s]
			}
		}
	}
}

void elspectra_evolution_data::init_tables_h2_electronic_diss(cross_section_table_mccc** cs, double**& rates, int& el_min_nb, int verbosity)
{
	int i, vi;
	double en, en_min;

	// eV, ionization energy of H2
	en_min = 15.43;

	// the cross sections are given for vi -> dissociative continuum	
	rates = alloc_2d_array<double>(NB_OF_H2_VSTATES_X1SU, nb_of_el_energies);
	memset(*rates, 0, NB_OF_H2_VSTATES_X1SU * nb_of_el_energies * sizeof(rates[0][0]));

	// pumping by electron impact is taken into account for H2 vibrational quanta of the ground electronic state vi < NB_OF_H2_VSTATES_X1SU
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		for (i = 0; i < nb_of_el_energies; i++) 
		{
			// entire electron energy bin must be higher than the threshold energy,
			if ((*cs[vi])(electron_energies_grid[i]) > 0.)
				rates[vi][i] = (*cs[vi])(electron_energies[i]) * electron_velocities[i];  // [cm2 *cm/s]
		}
		
		// is subtracted energy of the highest level, v = 12, j = 10	
		en = cs[vi]->get_threshold_energy() - 36104.6 * CM_INVERSE_TO_EV;
		if (en < en_min)
			en_min = en;
	}
	el_min_nb = get_nb_electron_energy_array(en_min);
}

void elspectra_evolution_data::init_tables_h2_rotational_exc(cross_section_table_mccc** cs, double**& rates, int verbosity)
{
	int i, j, dj, n, vi, li, lf;

	// (vi, ji) -> (vi, jf), 
	// jf = ji-2, ji + 2
	rates = alloc_2d_array<double>(2 * NB_OF_H2_VSTATES_X1SU * MAX_J_H2, nb_of_el_energies);
	memset(*rates, 0, 2 * NB_OF_H2_VSTATES_X1SU * MAX_J_H2 * nb_of_el_energies * sizeof(rates[0][0]));

	// vi < NB_OF_H2_VSTATES_X1SU - nb of vibrational states taken into account in the simulations,
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		for (j = 0; j < MAX_J_H2; j++) {
			for (dj = -2; dj <= 2; dj += 4) 
			{
				n = 2 * (vi * MAX_J_H2 + j) + (dj + 2)/4;
				li = h2_di->get_nb(vi, j);
				lf = h2_di->get_nb(vi, j + dj);

				if (li >= 0 && lf >= 0 && cs[n] != 0) {
					for (i = 0; i < nb_of_el_energies; i++)
					{
						// entire electron energy bin must be higher than the threshold energy,
						if ((*cs[n])(electron_energies_grid[i]) > 0.)
							rates[n][i] = (*cs[n])(electron_energies[i]) * electron_velocities[i];  // [cm2 *cm/s]
					}
				}
			}
		}
	}
}

void elspectra_evolution_data::init_tables_h2_rovibr_exc(cross_section_table_mccc **cs, double**& rates, int verbosity)
{
	int i, j, dj, n, vi, vf, li, lf;

	// (vi, ji) -> (vf, jf), 
	// vi != vf, jf = ji-2, ji, ji + 2
	rates = alloc_2d_array<double>(3 * NB_OF_H2_VSTATES_X1SU * (NB_OF_H2_VSTATES_X1SU - 1) * MAX_J_H2, nb_of_el_energies);
	memset(*rates, 0, 3 * NB_OF_H2_VSTATES_X1SU * (NB_OF_H2_VSTATES_X1SU - 1) * MAX_J_H2 * nb_of_el_energies * sizeof(rates[0][0]));

	// vi < NB_OF_H2_VSTATES_X1SU - nb of vibrational states taken into account in the simulations,
	for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
		for (j = 0; j < MAX_J_H2; j++) {
			for (vf = 0; vf < NB_OF_H2_VSTATES_X1SU; vf++) {
				if (vi != vf) {
					for (dj = -2; dj <= 2; dj += 2)
					{
						n = 3 * (vi * (NB_OF_H2_VSTATES_X1SU - 1) * MAX_J_H2 + vf * MAX_J_H2 + j) + (dj + 2) / 2;
						if (vf > vi)
							n -= 3 * MAX_J_H2;
						
						li = h2_di->get_nb(vi, j);
						lf = h2_di->get_nb(vf, j + dj);

						if (li >= 0 && lf >= 0 && cs[n] != 0) {
							for (i = 0; i < nb_of_el_energies; i++)
							{
								// entire electron energy bin must be higher than the threshold energy,
								if ((*cs[n])(electron_energies_grid[i]) > 0.)
									rates[n][i] = (*cs[n])(electron_energies[i]) * electron_velocities[i];  // [cm2 *cm/s]
							}
						}
					}
				}
			}
		}
	}
}

// He, old
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

			k = get_nb_electron_energy_array(z);
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
	hei_min_nb = get_nb_electron_energy_array(en_min);
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

int elspectra_evolution_data::get_nb_electron_energy_array(double energy) const
{
	int i;
	locate_index(electron_energies_grid, nb_of_el_energies + 1, energy, i);
	
	if (i < 0)
		i = 0;
	else if (i > nb_of_el_energies - 1)
		i = nb_of_el_energies - 1;
	return i;
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

int elspectra_evolution_data::get_vibr_nb_h2(int nb) const
{
	if (nb < 0 || nb >= nb_lev_h2)
		return -1;
	return h2_di->lev_array[nb].v;
}

void elspectra_evolution_data::get_el_energy_losses(double& mt, double& h2_rot, double &h2_vibr, double &h2_electr_sin, double &h2_electr_tri,
	double& ion, double& col_ions, double& col_el, double& hei) const
{
	// [eV cm-3 s-1]
	mt = enloss_rate_mt;
	h2_rot = enloss_rate_h2_rot;
	h2_vibr = enloss_rate_h2_vibr;
	h2_electr_sin = enloss_rate_h2_singlet + enloss_diss_h2_singlet;  // via excitation of singlet H2 states (including dissociation)
	h2_electr_tri = enloss_diss_h2_triplet;							  // triplet H2 states (only dissociation yet)
	ion = enloss_rate_ioniz;
	col_ions = enloss_rate_coloumb_ions;
	col_el = enloss_rate_coloumb_el;
	hei = enloss_rate_hei;
}

void elspectra_evolution_data::get_h2_process_rates(double& sol_diss, double& diss_exc_sin, double& diss_exc_tri, double & hei_exc)
{
	// [cm-3 s-1]
	sol_diss = h2_solomon_diss_rate;  
	diss_exc_sin = h2_diss_exc_singlet_rate;
	diss_exc_tri = h2_diss_exc_triplet_rate;
	hei_exc = hei_exc_rate;
}

int elspectra_evolution_data::f(realtype t, N_Vector y, N_Vector ydot)
{
	int i, j;
	double x, rate, log_coulomb, th_energy_deriv;  // z, vel_i;

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
	enloss_rate_mt = 0.;
	for (i = 1; i < nb_of_el_energies; i++) {
		// [cm2 eV cm/s eV-1] * [cm-3] * [cm-3]
		x = (el_enloss_h2_mt[i]* y_data[h2_nb] + el_enloss_he_mt[i] * y_data[he_nb]) * y_data[i];

		ydot_data[i] -= x;
		ydot_data[i-1] += x;

		// [cm-3 s-1 eV], electrons lose energy, < 0
		enloss_rate_mt -= x * (electron_energies[i] - electron_energies[i - 1]);
	}
	energy_gain_n -= EV_TO_ERGS * enloss_rate_mt;  // [erg cm-3 s-1], neutrals gain energy in this process, must be > 0.

	// Coulomb losses
	// Electron scattering on free electrons,
	enloss_rate_coloumb_el = 0.;
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
			enloss_rate_coloumb_el -= x * (electron_energies[i] - electron_energies[i - 1]);
		}
		if (i < nb_of_el_energies - 1 && x > 0.) {
			x *= 2. * y_data[i] * log_coulomb / (electron_energy_bin_size[i] + electron_energy_bin_size[i + 1]);
			
			ydot_data[i] -= x;
			ydot_data[i + 1] += x;
			enloss_rate_coloumb_el += x * (electron_energies[i + 1] - electron_energies[i]);
		}
	}
#endif

	// Electron scattering on ions,
	/*z = log_coulomb * (0.5 * y_data[h2p_nb] + 0.25 * y_data[hep_nb] + y_data[hepp_nb] + y_data[hp_nb]);
	enloss_rate_coloumb_ions = 0.;
	for (i = 1; i < nb_of_el_energies; i++) {
		// rates are normalized on electron energy bin size,
		x = y_data[i] * elion_scatt_rates[i] * z;
		ydot_data[i] -= x;
		ydot_data[i - 1] += x;

		// [cm-3 s-1 eV], < 0.
		enloss_rate_coloumb_ions -= x * (electron_energies[i] - electron_energies[i - 1]);
	}*/

	/*
	enloss_rate_coloumb_ions = 0.;
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
			enloss_rate_coloumb_ions += z*(electron_energies[i] - electron_energies[i-1]);
		}
		if (i < nb_of_el_energies - 1 && z > 0.) {
			z /= 0.5 * (electron_energy_bin_size[i] + electron_energy_bin_size[i + 1]);
			ydot_data[i] -= z;
			ydot_data[i + 1] += z;

			// [cm-3 s-1 eV], > 0
			enloss_rate_coloumb_ions += z * (electron_energies[i + 1] - electron_energies[i]);
		}
	}
	energy_gain_i -= EV_TO_ERGS * enloss_rate_coloumb_ions;  // [erg cm-3 s-1], opposite sign (electrons loose energy, ions gain),
	*/
	//
	// Electron impact ionization of H2, He, H, H2+, He+
	enloss_rate_ioniz = rate = 0.;
	derivatives_ionization(y_data, ydot_data, h2_ioniz_cs, h2_ioniz_rates, h2_ioniz_rates_tot, h2_ioniz_indexes, 
		h2_ioniz_min_nb, h2_nb, rate, enloss_rate_ioniz);

	ydot_data[h2_nb] -= rate;
	ydot_data[h2p_nb] += rate;
	ydot_data[el_nb] += rate;

	rate /= y_data[h2_nb];
	for (i = 0; i < nb_lev_h2; i++) {
		ydot_data[h2eq_nb + i] -= y_data[h2eq_nb + i] * rate;
	}

	derivatives_ionization(y_data, ydot_data, he_ioniz_cs, he_ioniz_rates, he_ioniz_rates_tot, he_ioniz_indexes,
		he_ioniz_min_nb, he_nb, rate, enloss_rate_ioniz);

	ydot_data[he_nb] -= rate;
	ydot_data[hep_nb] += rate;
	ydot_data[el_nb] += rate;

	rate /= y_data[he_nb];
	for (i = 0; i < nb_lev_hei; i++) {
		ydot_data[heieq_nb + i] -= y_data[heieq_nb + i] * rate;
	}

	derivatives_ionization(y_data, ydot_data, h_ioniz_cs, h_ioniz_rates, h_ioniz_rates_tot, h_ioniz_indexes,
		h_ioniz_min_nb, h_nb, rate, enloss_rate_ioniz);

	ydot_data[h_nb] -= rate;
	ydot_data[hp_nb] += rate;
	ydot_data[el_nb] += rate;

	derivatives_ionization(y_data, ydot_data, h2p_ioniz_cs, h2p_ioniz_rates, h2p_ioniz_rates_tot, h2p_ioniz_indexes,
		h2p_ioniz_min_nb, h2p_nb, rate, enloss_rate_ioniz);

	ydot_data[h2p_nb] -= rate;
	ydot_data[hp_nb] += 2.*rate;
	ydot_data[el_nb] += rate;

	derivatives_ionization(y_data, ydot_data, hep_ioniz_cs, hep_ioniz_rates, hep_ioniz_rates_tot, hep_ioniz_indexes,
		hep_ioniz_min_nb, hep_nb, rate, enloss_rate_ioniz);

	ydot_data[hep_nb] -= rate;
	ydot_data[hepp_nb] += rate;
	ydot_data[el_nb] += rate;

	// H2 electronic states excitation,
	th_energy_deriv = 0.;  // [erg cm-3 s-1], thermal energy, gained by neutral gas per s,
	enloss_rate_h2_singlet = enloss_diss_h2_singlet = enloss_diss_h2_triplet = 0.;   // [eV cm-3 s-1], electron energy loss rate
	h2_solomon_diss_rate = h2_diss_exc_singlet_rate = h2_diss_exc_triplet_rate = 0.;  // [cm-3 s-1],   H2 dissociation rate

	// S1u+(X) -> S1u+(B)
	derivatives_h2_electronic_exc_dl0(y_data, ydot_data, h2_di_b, h2_b_state_data, h2_bstate_rates, MAX_H2_VSTATES_B1SU, 
		h2_b1su_min_nb, enloss_rate_h2_singlet, h2_solomon_diss_rate, th_energy_deriv);

	// S1u(X) -> S1u+(Bp)
	derivatives_h2_electronic_exc_dl0(y_data, ydot_data, h2_di_bp, h2_bp_state_data, h2_bpstate_rates, MAX_H2_VSTATES_BP1SU, 
		h2_bp1su_min_nb, enloss_rate_h2_singlet, h2_solomon_diss_rate, th_energy_deriv);

	// S1u(X) -> P1u(C-/+)
	derivatives_h2_electronic_exc_dl1(y_data, ydot_data, h2_di_cplus, h2_di_cminus, h2_cplus_state_data, h2_cminus_state_data,
		h2_cstate_rates, MAX_H2_VSTATES_C1PU, h2_c1pu_min_nb, enloss_rate_h2_singlet, h2_solomon_diss_rate, th_energy_deriv);

	// S1u(X) -> P1u(D-/+)
	derivatives_h2_electronic_exc_dl1(y_data, ydot_data, h2_di_dplus, h2_di_dminus, h2_dplus_state_data, h2_dminus_state_data,
		h2_dstate_rates, MAX_H2_VSTATES_D1PU, h2_d1pu_min_nb, enloss_rate_h2_singlet, h2_solomon_diss_rate, th_energy_deriv);

	// dissociative excitation,
	// S1u+(B)
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_bstate_diss_cs, h2_bstate_diss_rates, h2_b1su_diss_min_nb, 
		enloss_diss_h2_singlet, h2_diss_exc_singlet_rate);

	// S1u+(Bp)
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_bpstate_diss_cs, h2_bpstate_diss_rates, h2_bp1su_diss_min_nb,
		enloss_diss_h2_singlet, h2_diss_exc_singlet_rate);

	// P1u(C)
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_cstate_diss_cs, h2_cstate_diss_rates, h2_c1pu_diss_min_nb,
		enloss_diss_h2_singlet, h2_diss_exc_singlet_rate);

	// P1u(D)
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_dstate_diss_cs, h2_dstate_diss_rates, h2_d1pu_diss_min_nb,
		enloss_diss_h2_singlet, h2_diss_exc_singlet_rate);

	// s3u(b) (triplet)
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_3bstate_diss_cs, h2_3bstate_diss_rates, h2_b3su_diss_min_nb, 
		enloss_diss_h2_triplet, h2_diss_exc_triplet_rate);

	// the decrease of population densities of H2 is taken into account in the functions above,
	x = h2_solomon_diss_rate + h2_diss_exc_singlet_rate + h2_diss_exc_triplet_rate;
	ydot_data[h2_nb] -= x;
	ydot_data[h_nb] += 2. * x;

	energy_gain_n += th_energy_deriv;  // [erg cm-3 s-1]
	
	enloss_rate_h2_rot = enloss_rate_h2_vibr = 0.;  // [eV cm-3 s-1]
#if (ROVIBRATIONAL_EXC_LOSSES)
	// pure H2 rotational excitation,
	derivatives_h2_rotational_exc(y_data, ydot_data, h2_rot_rates, enloss_rate_h2_rot);

	// H2 ro-vibrational excitation, 
	derivatives_h2_rovibr_exc(y_data, ydot_data, h2_rovibr_rates, enloss_rate_h2_vibr);
	enloss_rate_h2_vibr += enloss_rate_h2_rot;  // including pure rotational transitions,
#endif

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

	/*
	// Excitation of HeI by electron impact
	enloss_rate_hei = hei_exc_rate = 0.;
	derivatives_hei_exc(y_data, ydot_data, enloss_rate_hei, hei_exc_rate);

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
	*/
	return 0;
}

void elspectra_evolution_data::derivatives_ionization(const realtype* y_data, realtype* ydot_data, const electron_impact_ionization *cs, 
	double ** rates, const double * rates_tot, spectra_data ** indexes, int el_ion_min_nb, int target_nb, double& rate, double & enloss_rate)
{
	int i, j, j_max, l;
	double enl, en_max, x, y, r;
	
	enl = r = 0.;
#pragma omp parallel reduction(+: enl, r) private(i, j, j_max, l, en_max, x, y)
	{
		double* arr_el = new double[nb_of_el_energies];
		memset(arr_el, 0, nb_of_el_energies * sizeof(double));

		// incident electron en_i -> scattered electron en_j
#pragma omp for schedule(dynamic, 1)
		for (i = el_ion_min_nb; i < nb_of_el_energies; i++) 
		{
			en_max = 0.5 * (electron_energies[i] - cs->get_binding_energy());
			j_max = get_nb_electron_energy_array(en_max);
			y = y_data[i] * y_data[target_nb];

			for (j = 0; j <= j_max; j++) {
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

					enl += (esd.w1 * electron_energies[l] + esd.w2 * electron_energies[l + 1]) * x;
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
	// check: energy loss = rate * binding_energy
	enloss_rate += enl;  // [cm-3 s-1 eV]
	rate = r;
}


void elspectra_evolution_data::derivatives_h2_electronic_exc_dl0(const realtype* y_data, realtype* ydot_data, const energy_diagram* h2_di_exc, 
	const vector<h2_energy_level_param> h2_exc_state_data, double** rates, 
	int vibr_qnb_final_max, int el_min_nb, double& enloss_rate, double& diss_rate, double& th_energy_deriv)
{
	int i, k, j, n, dj, vi, vf, low, up;
	double dissr, en, enl, th_en, x, y, f;

	dissr = enl = th_en = 0.;
#pragma omp parallel reduction(+: dissr, enl, th_en) private(en, i, k, n, dj, j, vi, vf, low, up, y, x, f)
	{
		int i0;
		double w1, w2;

		double* arr_h2 = new double[nb_lev_h2];
		memset(arr_h2, 0, nb_lev_h2 * sizeof(double));

		double* arr_el = new double[nb_of_el_energies];
		memset(arr_el, 0, nb_of_el_energies * sizeof(double));

#pragma omp for schedule(dynamic, 1)
		for (low = 0; low < nb_lev_h2; low++) 
		{
			const energy_level& low_lev_ref = h2_di->lev_array[low];
			vi = low_lev_ref.v;

			if (vi < NB_OF_H2_VSTATES_X1SU) {
				j = rounding(h2_di->lev_array[low].j);

				for (vf = 0; vf < vibr_qnb_final_max; vf++) {
					for (dj = -1; dj <= 1; dj += 2)
					{
						up = h2_di_exc->get_nb(vf, j + dj);
						if (up != -1) {
							const energy_level& up_lev_ref = h2_di_exc->lev_array[up];
							const h2_energy_level_param& lev_param_ref = h2_exc_state_data[up];
							
							// Honl-London factors, 
							// for S1g(X) -> S1u(B), S1u(Bp), j -> j + dj, 
							f = hl_singlet_dl0(j, dj) * y_data[h2eq_nb + low];
							y = 0.;
							n = vi * vibr_qnb_final_max + vf;

							for (i = el_min_nb; i < nb_of_el_energies; i++) 
							{
								en = electron_energies[i] - (up_lev_ref.energy - low_lev_ref.energy) * CM_INVERSE_TO_EV;
								if (en < 0.)
									en = 0.;

								k = get_nb_electron_energy_array(en);

								if (k == 0 && en < electron_energies[k]) {
									i0 = 0;
									w1 = 1.;
								}
								else {
									if (en < electron_energies[k])
										k--;

									i0 = k;
									w1 = (electron_energies[k + 1] - en) / (electron_energies[k + 1] - electron_energies[k]);
									if (w1 > 1.)
										w1 = 1.;
								}
								w2 = 1. - w1;

								// [cm-3] *[cm-3] *[cm2 *cm/s] = [cm-3 s-1]
								x = f * y_data[i] * rates[n][i];

								arr_el[i0] += x * w1;
								arr_el[i0 + 1] += x * w2;

								// must be equal to x * delta_E
								enl += (w1 * electron_energies[i0] + w2 * electron_energies[i0 + 1]
									- electron_energies[i]) * x;  // [cm-3 s-1 eV],

								arr_el[i] -= x;
								y += x;
							}

							arr_h2[low] -= y;
							for (k = 0; k < (int)(lev_param_ref.nb_of_decays); k++) 
							{
								i = lev_param_ref.decay_level_nbs[k];
								arr_h2[i] += y * lev_param_ref.decay_probs[k];
							}
							// the downward transition (dissociative) to the vibrational continuum, 
							dissr += y * lev_param_ref.diss_prob;  // [cm-3 s-1]

							// thermal energy gain by neutral gas per s,
							th_en += y * lev_param_ref.diss_prob * lev_param_ref.kin_energy;  // [erg cm-3 s-1]
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
	enloss_rate += enl;
	th_energy_deriv += th_en;
}


void elspectra_evolution_data::derivatives_h2_electronic_exc_dl1(const realtype* y_data, realtype* ydot_data,
	const energy_diagram* h2_di_plus, const energy_diagram* h2_di_minus,
	const vector<h2_energy_level_param> h2_state_plus_data, const vector<h2_energy_level_param> h2_state_minus_data,
	double** rates, int vibr_qnb_final_max, int el_min_nb, double& enloss_rate, double& diss_rate, double& th_energy_deriv)
{
	int i, k, j, n, dj, vi, vf, low, up;
	double dissr, en, enl, th_en, x, y, f;

	dissr = enl = th_en = 0.;
#pragma omp parallel reduction(+: dissr, enl, th_en) private(en, i, k, n, dj, j, vi, vf, low, up, y, x, f)
	{
		int i0;
		double w1, w2;

		double* arr_h2 = new double[nb_lev_h2];
		memset(arr_h2, 0, nb_lev_h2 * sizeof(double));

		double* arr_el = new double[nb_of_el_energies];
		memset(arr_el, 0, nb_of_el_energies * sizeof(double));

#pragma omp for schedule(dynamic, 1)
		for (low = 0; low < nb_lev_h2; low++)
		{
			const energy_level& low_lev_ref = h2_di->lev_array[low];
			vi = low_lev_ref.v;

			if (vi < NB_OF_H2_VSTATES_X1SU) {
				j = rounding(h2_di->lev_array[low].j);

				for (vf = 0; vf < vibr_qnb_final_max; vf++) {
					for (dj = -1; dj <= 1; dj += 2)
					{
						up = h2_di_plus->get_nb(vf, j + dj);
						if (up != -1) 
						{
							const energy_level& up_lev_ref = h2_di_plus->lev_array[up];
							const h2_energy_level_param& lev_param_ref = h2_state_plus_data[up];

							// Honl-London factors,
							// for S1g(X) -> P1u(C), P1u(D), j -> j + dj, 
							f = hl_singlet_dl1(j, dj) * y_data[h2eq_nb + low];
							y = 0.;
							n = vi * vibr_qnb_final_max + vf;

							for (i = el_min_nb; i < nb_of_el_energies; i++) 
							{
								en = electron_energies[i] - (up_lev_ref.energy - low_lev_ref.energy) * CM_INVERSE_TO_EV;
								if (en < 0.)
									en = 0.;

								k = get_nb_electron_energy_array(en);

								if (k == 0 && en < electron_energies[k]) {
									i0 = 0;
									w1 = 1.;
								}
								else {
									if (en < electron_energies[k])
										k--;

									i0 = k;
									w1 = (electron_energies[k + 1] - en) / (electron_energies[k + 1] - electron_energies[k]);
									if (w1 > 1.)
										w1 = 1.;
								}
								w2 = 1. - w1;

								// [cm-3] *[cm-3] *[cm2 *cm/s] = [cm-3 s-1]
								x = f * y_data[i] * rates[n][i];

								arr_el[i0] += x * w1;
								arr_el[i0 + 1] += x * w2;

								// must be equal to x * delta_E
								enl += (w1 * electron_energies[i0] + w2 * electron_energies[i0 + 1]
									- electron_energies[i]) * x;  // [cm-3 s-1 eV],

								arr_el[i] -= x;
								y += x;
							}

							arr_h2[low] -= y;
							for (k = 0; k < (int)(lev_param_ref.nb_of_decays); k++)
							{
								i = lev_param_ref.decay_level_nbs[k];
								arr_h2[i] += y * lev_param_ref.decay_probs[k];
							}
							// the downward transition (dissociative) to the vibrational continuum, 
							dissr += y * lev_param_ref.diss_prob;  // [cm-3 s-1]

							// thermal energy gain by neutral gas per s,
							th_en += y * lev_param_ref.diss_prob * lev_param_ref.kin_energy;  // [erg cm-3 s-1]
						}
					}

					dj = 0;
					up = h2_di_minus->get_nb(vf, j);

					if (up != -1)
					{
						const energy_level& up_lev_ref = h2_di_minus->lev_array[up];
						const h2_energy_level_param& lev_param_ref = h2_state_minus_data[up];	

						// Honl-London factors, 
						f = hl_singlet_dl1(j, 0) * y_data[h2eq_nb + low];
						y = 0.;
						n = vi * vibr_qnb_final_max + vf;

						for (i = el_min_nb; i < nb_of_el_energies; i++) 
						{
							en = electron_energies[i] - (up_lev_ref.energy - low_lev_ref.energy) * CM_INVERSE_TO_EV;
							if (en < 0.)
								en = 0.;

							k = get_nb_electron_energy_array(en);

							if (k == 0 && en < electron_energies[k]) {
								i0 = 0;
								w1 = 1.;
							}
							else {
								if (en < electron_energies[k])
									k--;

								i0 = k;
								w1 = (electron_energies[k + 1] - en) / (electron_energies[k + 1] - electron_energies[k]);
								if (w1 > 1.)
									w1 = 1.;
							}
							w2 = 1. - w1;

							// [cm-3] *[cm-3] *[cm2 *cm/s] = [cm-3 s-1]
							x = f * y_data[i] * rates[n][i];

							arr_el[i0] += x * w1;
							arr_el[i0 + 1] += x * w2;

							// must be equal to x * delta_E
							enl += (w1 * electron_energies[i0] + w2 * electron_energies[i0 + 1]
								- electron_energies[i]) * x;  // [cm-3 s-1 eV],

							arr_el[i] -= x;
							y += x;
						}

						arr_h2[low] -= y;
						for (k = 0; k < (int)(lev_param_ref.nb_of_decays); k++)
						{
							i = lev_param_ref.decay_level_nbs[k];
							arr_h2[i] += y * lev_param_ref.decay_probs[k];
						}

						dissr += y * lev_param_ref.diss_prob;
						th_en += y * lev_param_ref.diss_prob * lev_param_ref.kin_energy;  
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
	enloss_rate += enl;
	th_energy_deriv += th_en;
}


void elspectra_evolution_data::derivatives_h2_electronic_diss(const realtype* y_data, realtype* ydot_data, cross_section_table_mccc** cs, 
	double** rates, int el_min_nb, double& enloss_rate, double& diss_rate)
{
	int i, k, l, vi, low;
	double dissr, en, den, enl, th_en, x, y;

	dissr = enl = th_en = 0.;
#pragma omp parallel reduction(+: dissr, enl, th_en) private(en, den, i, k, l, vi, low, y, x)
	{
		int i0;
		double w1, w2;

		double* arr_h2 = new double[nb_lev_h2];
		memset(arr_h2, 0, nb_lev_h2 * sizeof(double));

		double* arr_el = new double[nb_of_el_energies];
		memset(arr_el, 0, nb_of_el_energies * sizeof(double));

#pragma omp for schedule(dynamic, 1)
		for (low = 0; low < nb_lev_h2; low++)
		{
			const energy_level& low_lev_ref = h2_di->lev_array[low];
			vi = low_lev_ref.v;

			if (vi < NB_OF_H2_VSTATES_X1SU) {
				l = h2_di->get_nb(vi, 0);
				
				// energy loss by the electron at the collision,
				den = cs[vi]->get_threshold_energy()
					- (low_lev_ref.energy - h2_di->lev_array[l].energy) * CM_INVERSE_TO_EV;

				y = 0.;
				for (i = el_min_nb; i < nb_of_el_energies; i++) 
				{	
					en = electron_energies[i] - den;
					if (en < 0.)
						en = 0.;

					k = get_nb_electron_energy_array(en);

					if (k == 0 && en < electron_energies[k]) {
						i0 = 0;
						w1 = 1.;
					}
					else {
						if (en < electron_energies[k])
							k--;

						i0 = k;
						w1 = (electron_energies[k + 1] - en) / (electron_energies[k + 1] - electron_energies[k]);
						if (w1 > 1.)
							w1 = 1.;
					}
					w2 = 1. - w1;

					// [cm-3] *[cm-3] *[cm2 *cm/s] = [cm-3 s-1]
					x = y_data[i] * y_data[h2eq_nb + low] * rates[vi][i];

					arr_el[i0] += x * w1;
					arr_el[i0 + 1] += x * w2;

					// must be equal to x * delta_E
					enl += (w1 * electron_energies[i0] + w2 * electron_energies[i0 + 1]
						- electron_energies[i]) * x;  // [cm-3 s-1 eV],

					arr_el[i] -= x;
					y += x;
				}
				arr_h2[low] -= y;
				dissr += y;
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
	enloss_rate += enl;
}


void elspectra_evolution_data::derivatives_h2_rotational_exc(const realtype* y_data, realtype* ydot_data, double**& rates, double& enloss_rate)
{
	int i, k, n, j, dj, vi, li, lf;
	double en, enl, x, y;

	enl = 0.;
#pragma omp parallel reduction(+: enl) private(en, i, k, n, dj, j, vi, li, lf, y, x)
	{
		int i0;
		double w1, w2;

		double* arr_h2 = new double[nb_lev_h2];
		memset(arr_h2, 0, nb_lev_h2 * sizeof(double));

		double* arr_el = new double[nb_of_el_energies];
		memset(arr_el, 0, nb_of_el_energies * sizeof(double));

#pragma omp for schedule(dynamic, 1)
		for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
			for (j = 0; j < MAX_J_H2; j++) {
				for (dj = -2; dj <= 2; dj += 4)
				{
					li = h2_di->get_nb(vi, j);
					lf = h2_di->get_nb(vi, j + dj);

					if (li >= 0 && lf >= 0) {
						y = 0.;
						n = 2 * (vi * MAX_J_H2 + j) + (dj + 2) / 4;
						
						const energy_level& li_ref = h2_di->lev_array[li];
						const energy_level& lf_ref = h2_di->lev_array[lf];
						
						for (i = 0; i < nb_of_el_energies; i++)
						{
							en = electron_energies[i] - (lf_ref.energy - li_ref.energy) * CM_INVERSE_TO_EV;
							if (en < 0.)
								en = 0.;

							k = get_nb_electron_energy_array(en);

							if (k == 0 && en < electron_energies[k]) {
								i0 = 0;
								w1 = 1.;
							}
							else {
								if (en < electron_energies[k])
									k--;

								i0 = k;
								w1 = (electron_energies[k + 1] - en) / (electron_energies[k + 1] - electron_energies[k]);
								if (w1 > 1.)
									w1 = 1.;
							}
							w2 = 1. - w1;

							// [cm-3] *[cm-3] *[cm2 *cm/s] = [cm-3 s-1]
							x = y_data[i] * y_data[h2eq_nb + li] * rates[n][i];

							arr_el[i0] += x * w1;
							arr_el[i0 + 1] += x * w2;

							// must be equal to x * delta_E
							enl += (w1 * electron_energies[i0] + w2 * electron_energies[i0 + 1]
								- electron_energies[i]) * x;  // [cm-3 s-1 eV],

							arr_el[i] -= x;
							y += x;
						}

						arr_h2[li] -= y;
						arr_h2[lf] += y;
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
	enloss_rate += enl;
}

void elspectra_evolution_data::derivatives_h2_rovibr_exc(const realtype* y_data, realtype* ydot_data, double**& rates, double& enloss_rate)
{
	int i, k, n, j, dj, vi, vf, li, lf;
	double en, enl, x, y;

	enl = 0.;
#pragma omp parallel reduction(+: enl) private(en, i, k, n, dj, j, vi, vf, li, lf, y, x)
	{
		int i0;
		double w1, w2;

		double* arr_h2 = new double[nb_lev_h2];
		memset(arr_h2, 0, nb_lev_h2 * sizeof(double));

		double* arr_el = new double[nb_of_el_energies];
		memset(arr_el, 0, nb_of_el_energies * sizeof(double));

#pragma omp for schedule(dynamic, 1)
		for (vi = 0; vi < NB_OF_H2_VSTATES_X1SU; vi++) {
			for (j = 0; j < MAX_J_H2; j++) {
				for (vf = 0; vf < NB_OF_H2_VSTATES_X1SU; vf++) {
					if (vi != vf) {
						for (dj = -2; dj <= 2; dj += 2)
						{
							li = h2_di->get_nb(vi, j);
							lf = h2_di->get_nb(vf, j + dj);

							if (li >= 0 && lf >= 0) {
								y = 0.;

								n = 3 * (vi * (NB_OF_H2_VSTATES_X1SU - 1) * MAX_J_H2 + vf * MAX_J_H2 + j) + (dj + 2) / 2;
								if (vf > vi)
									n -= 3 * MAX_J_H2;

								const energy_level& li_ref = h2_di->lev_array[li];
								const energy_level& lf_ref = h2_di->lev_array[lf];

								for (i = 0; i < nb_of_el_energies; i++)
								{
									en = electron_energies[i] - (lf_ref.energy - li_ref.energy) * CM_INVERSE_TO_EV;
									if (en < 0.)
										en = 0.;

									k = get_nb_electron_energy_array(en);

									if (k == 0 && en < electron_energies[k]) {
										i0 = 0;
										w1 = 1.;
									}
									else {
										if (en < electron_energies[k])
											k--;

										i0 = k;
										w1 = (electron_energies[k + 1] - en) / (electron_energies[k + 1] - electron_energies[k]);
										if (w1 > 1.)
											w1 = 1.;
									}
									w2 = 1. - w1;

									// [cm-3] *[cm-3] *[cm2 *cm/s] = [cm-3 s-1]
									x = y_data[i] * y_data[h2eq_nb + li] * rates[n][i];

									arr_el[i0] += x * w1;
									arr_el[i0 + 1] += x * w2;

									// must be equal to x * delta_E
									enl += (w1 * electron_energies[i0] + w2 * electron_energies[i0 + 1]
										- electron_energies[i]) * x;  // [cm-3 s-1 eV],

									arr_el[i] -= x;
									y += x;
								}

								arr_h2[li] -= y;
								arr_h2[lf] += y;
							}
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
	enloss_rate += enl;
}


void elspectra_evolution_data::derivatives_hei_exc(const realtype* y_data, realtype* ydot_data, double& enloss_rate, double& exc_rate)
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
	enloss_rate += enl;
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
		NV_Ith_S(abs_tol, i) = ABS_POPULATION_ERROR_SOLVER;
	}
	i0 += nb_lev_h2;
	for (i = i0; i < i0 + nb_lev_hei; i++) {
		NV_Ith_S(abs_tol, i) = ABS_POPULATION_ERROR_SOLVER;
	}
	i0 += nb_lev_hei;
	for (i = i0; i < nb_of_equat; i++) {
		NV_Ith_S(abs_tol, i) = ABS_PARAMETER_ERROR_SOLVER;
	}
}


/*
* 
	cross_section_integral_2d ion_cs_int_2d(cs, 1.e-6);
	cross_section_integral_weight ion_cs_int_weight(cs, 1.e-6);

void init_tables_coulomb(electron_impact_ionization*, double*&, double**&, spectra_data_adv**&, int verbosity = 1);
    
void elspectra_evolution_data::init_tables_coulomb(electron_impact_ionization* elel_cs,
	double*& el_ion_rates_tot, double**& el_ion_rates, spectra_data_adv**& el_ion_indexes, int verbosity)
{
}

	n = 0;
#ifdef _OPENMP
	n = omp_get_thread_num();
#endif

*/


//
// OLD
/*
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
	// the level energy difference of D- and D+ is neglected, D- state have much higher number of levels,
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
*/

/*
	//
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
		*/


		// OLD
		/*delete rot_h2_j02_cs;
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


		delete[] elel_scatt_rates;
		delete[] elion_scatt_rates;

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
		*/
/*
* 
	// H2 electronic states excitation,
	// energy loss in [eV cm-3 s-1], rate in [cm-3 s-1]
	enloss_rate_h2_singlet = h2_solomon_diss_rate = h2_diss_exc_singlet_rate = thermal_energy = 0.;
	
	// S1u(X) -> S1u(B)
	derivatives_h2_electronic_exc(y_data, ydot_data, h2_di_b, h2_b_state_data, h2_v0_bstate_rates, h2_v0_bstate_indexes, 
		0, MAX_H2_VSTATES_B1SU, h2_b1su_min_nb, enloss_rate_h2_singlet, h2_solomon_diss_rate, thermal_energy);

	derivatives_h2_electronic_exc(y_data, ydot_data, h2_di_b, h2_b_state_data, h2_v1_bstate_rates, h2_v1_bstate_indexes, 
		1, MAX_H2_VSTATES_B1SU, h2_b1su_min_nb, enloss_rate_h2_singlet, h2_solomon_diss_rate, thermal_energy);

	// S1u(X) -> S1u(Bp)
	derivatives_h2_electronic_exc(y_data, ydot_data, h2_di_bp, h2_bp_state_data, h2_v0_bpstate_rates, h2_v0_bpstate_indexes,
		0, MAX_H2_VSTATES_BP1SU, h2_bp1su_min_nb, enloss_rate_h2_singlet, h2_solomon_diss_rate, thermal_energy);

	derivatives_h2_electronic_exc(y_data, ydot_data, h2_di_bp, h2_bp_state_data, h2_v1_bpstate_rates, h2_v1_bpstate_indexes,
		1, MAX_H2_VSTATES_BP1SU, h2_bp1su_min_nb, enloss_rate_h2_singlet, h2_solomon_diss_rate, thermal_energy);

	// S1u(X) -> P1u(C-/+)
	derivatives_h2_electronic_exc_2(y_data, ydot_data, h2_di_cplus, h2_di_cminus, h2_cplus_state_data, h2_cminus_state_data,
		h2_v0_cstate_rates, h2_v0_cstate_indexes, 0, MAX_H2_VSTATES_C1PU, h2_c1pu_min_nb, enloss_rate_h2_singlet, h2_solomon_diss_rate, thermal_energy);

	derivatives_h2_electronic_exc_2(y_data, ydot_data, h2_di_cplus, h2_di_cminus, h2_cplus_state_data, h2_cminus_state_data,
		h2_v1_cstate_rates, h2_v1_cstate_indexes, 1, MAX_H2_VSTATES_C1PU, h2_c1pu_min_nb, enloss_rate_h2_singlet, h2_solomon_diss_rate, thermal_energy);

	// S1u(X) -> P1u(D-/+)
	derivatives_h2_electronic_exc_2(y_data, ydot_data, h2_di_dplus, h2_di_dminus, h2_dplus_state_data, h2_dminus_state_data,
		h2_v0_dstate_rates, h2_v0_dstate_indexes, 0, MAX_H2_VSTATES_D1PU, h2_d1pu_min_nb, enloss_rate_h2_singlet, h2_solomon_diss_rate, thermal_energy);

	derivatives_h2_electronic_exc_2(y_data, ydot_data, h2_di_dplus, h2_di_dminus, h2_dplus_state_data, h2_dminus_state_data,
		h2_v1_dstate_rates, h2_v1_dstate_indexes, 1, MAX_H2_VSTATES_D1PU, h2_d1pu_min_nb, enloss_rate_h2_singlet, h2_solomon_diss_rate, thermal_energy);


	//  dissociative excitation through S1u(B) state
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v0_bstate_diss_rates, h2_v0_bstate_diss_indexes,
		0, h2_b1su_diss_min_nb, enloss_rate_h2_singlet, h2_diss_exc_singlet_rate);

	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v1_bstate_diss_rates, h2_v1_bstate_diss_indexes,
		1, h2_b1su_diss_min_nb, enloss_rate_h2_singlet, h2_diss_exc_singlet_rate);

	// through P1u(C) state
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v0_cstate_diss_rates, h2_v0_cstate_diss_indexes,
		0, h2_c1pu_diss_min_nb, enloss_rate_h2_singlet, h2_diss_exc_singlet_rate);

	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v1_cstate_diss_rates, h2_v1_cstate_diss_indexes,
		1, h2_c1pu_diss_min_nb, enloss_rate_h2_singlet, h2_diss_exc_singlet_rate);

	// through S1u(Bp) state
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v0_bpstate_diss_rates, h2_v0_bpstate_diss_indexes,
		0, h2_bp1su_diss_min_nb, enloss_rate_h2_singlet, h2_diss_exc_singlet_rate);

	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v1_bpstate_diss_rates, h2_v1_bpstate_diss_indexes,
		1, h2_bp1su_diss_min_nb, enloss_rate_h2_singlet, h2_diss_exc_singlet_rate);

	// through P1u(D) state
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v0_dstate_diss_rates, h2_v0_dstate_diss_indexes,
		0, h2_d1pu_diss_min_nb, enloss_rate_h2_singlet, h2_diss_exc_singlet_rate);

	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v1_dstate_diss_rates, h2_v1_dstate_diss_indexes,
		1, h2_d1pu_diss_min_nb, enloss_rate_h2_singlet, h2_diss_exc_singlet_rate);

	// through s3u(b) state (triplet)
	h2_diss_exc_triplet_rate = enloss_rate_h2_triplet = 0.;
	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v0_3bstate_diss_rates, h2_v0_3bstate_diss_indexes,
		0, h2_b3su_diss_min_nb, enloss_rate_h2_triplet, h2_diss_exc_triplet_rate);

	derivatives_h2_electronic_diss(y_data, ydot_data, h2_v1_3bstate_diss_rates, h2_v1_3bstate_diss_indexes,
		1, h2_b3su_diss_min_nb, enloss_rate_h2_triplet, h2_diss_exc_triplet_rate);

*/

/*
// old
void elspectra_evolution_data::init_tables_h2_rotational_exc_old(const cross_section* cs, double *& rates_up, spectra_data*& indexes_up,
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
		k = get_nb_electron_energy_array(z);
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

// old
void elspectra_evolution_data::init_tables_h2_vibrational_exc_old(const cross_section* cs, double*& rates,
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
void elspectra_evolution_data::init_tables_h2_electronic_exc_old(const energy_diagram* h2_di_elexc, cross_section_table_mccc** cs,
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

// old
void elspectra_evolution_data::init_tables_h2_electronic_diss_old(cross_section_table_mccc* cs, double*& rates,
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

*/
/*
* 
// OLD
void elspectra_evolution_data::derivatives_h2_rotational_exc(const realtype* y_data, realtype* ydot_data, 
	const double* rates_up, const spectra_data* indexes_up, const double* rates_down, const spectra_data* indexes_down, int j0, double & enloss_rate)
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
		
		enloss_rate += (esd1.w1 *electron_energies[esd1.i0] + esd1.w2 * electron_energies[esd1.i0 + 1] 
			- electron_energies[i]) *x;  // [cm-3 s-1 eV],
		
		ydot_data[i] -= x;

		// j0 + 2 -> j0
		const spectra_data& esd2 = indexes_down[i];

		// [cm2 *cm/s] *[cm-3] *[cm-3]
		x = rates_down[i] * y_data[i] * y_data[h2eq_nb + j0 + 2];
		t -= x;

		ydot_data[esd2.i0] += x * esd2.w1;
		ydot_data[esd2.i0 + 1] += x * esd2.w2;

		enloss_rate += (esd2.w1 * electron_energies[esd2.i0] + esd2.w2 * electron_energies[esd2.i0 + 1]
			- electron_energies[i]) * x;  // [cm-3 s-1 eV],

		ydot_data[i] -= x;
	}
	ydot_data[h2eq_nb + j0] -= t;
	ydot_data[h2eq_nb + j0 + 2] += t;
}

void elspectra_evolution_data::derivatives_h2_vibrational_exc(const realtype* y_data, realtype* ydot_data, 
	const double* rates, spectra_data** indexes, int vibr_qnb_final, double & enloss_rate)
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
	enloss_rate += enl;
}

// S1g(X) -> S1u(B,Bp), dJ = +/-1,  
// vi fixed, vf = [0,1,.., MAX_H2_VSTATES_B1SU-1]
void elspectra_evolution_data::derivatives_h2_electronic_exc_old(const realtype* y_data, realtype* ydot_data, 
	const energy_diagram* h2_di_exc, const vector<h2_energy_level_param> h2_exc_state_data, double** rates, spectra_data** indexes, 
	int vibr_qnb_init, int vibr_qnb_final_max, int el_min_nb, double& enloss_rate, double& diss_rate, double& thermal_energy)
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
	enloss_rate += enl;
	thermal_energy += th_en;
}

// the collisions must not violate the ortho/para state!
// S1g(X) -> P1u(C,D), C- dJ = 0, C+ dJ = +/-1, Abgrall et al., A&A 253, 525 (1992)
void elspectra_evolution_data::derivatives_h2_electronic_exc_2(const realtype* y_data, realtype* ydot_data,
	const energy_diagram* h2_di_plus, const energy_diagram* h2_di_minus,
	const vector<h2_energy_level_param> h2_state_plus_data, const vector<h2_energy_level_param> h2_state_minus_data,
	double** rates, spectra_data** indexes, int vibr_qnb_init, int vibr_qnb_final_max, int el_min_nb, 
	double& enloss_rate, double& diss_rate, double& thermal_energy)
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
	enloss_rate += enl;
	thermal_energy += th_en;
}

void elspectra_evolution_data::derivatives_h2_electronic_diss_old(const realtype* y_data, realtype* ydot_data, 
	const double* rates, spectra_data** indexes, int vibr_qnb_init, int el_min_nb, double& enloss_rate, double& diss_rate)
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
	enloss_rate += enl;
	diss_rate += dissr;
}

*/

// Must be rewritten
// Excitation of H2 rotational levels,
/*
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
*/

/*
	enloss_rate_h2_rot = 0.;
	derivatives_h2_rotational_exc(y_data, ydot_data, h2_j02_rates, h2_j02_indexes, h2_j20_rates, h2_j20_indexes, 0, enloss_rate_h2_rot);
	derivatives_h2_rotational_exc(y_data, ydot_data, h2_j13_rates, h2_j13_indexes, h2_j31_rates, h2_j31_indexes, 1, enloss_rate_h2_rot);
	derivatives_h2_rotational_exc(y_data, ydot_data, h2_j24_rates, h2_j24_indexes, h2_j42_rates, h2_j42_indexes, 2, enloss_rate_h2_rot);
	derivatives_h2_rotational_exc(y_data, ydot_data, h2_j35_rates, h2_j35_indexes, h2_j53_rates, h2_j53_indexes, 3, enloss_rate_h2_rot);

	// H2 vibrational excitation,
	enloss_rate_h2_vibr = 0.;
	derivatives_h2_vibrational_exc(y_data, ydot_data, h2_v01_rates, h2_v01_indexes, 1, enloss_rate_h2_vibr);
	derivatives_h2_vibrational_exc(y_data, ydot_data, h2_v02_rates, h2_v02_indexes, 2, enloss_rate_h2_vibr);
	*/
