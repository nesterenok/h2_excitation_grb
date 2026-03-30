#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

#include "constants.h"
#include "spectroscopy.h"
#include "h2_microphysics.h"
#include "h2_uv_pump_data.h"
#include "h2_electron_cs_table.h"

using namespace std;

void save_cross_section_table(const string& output_path, const string& data_path, double ionization_fraction, double thermal_electron_temp, int ji,
	const string& str)
{
	const int vi = 0;
	const double energy_min = 0.01;  // eV
	const double energy_max = 1.e+7;
	const int nb_of_bins_per_order = 100;

	bool is_extrapolation_on;
	int i, k, nb_of_energies, lf, vf, dj, electronic_state, isotope, nb_lev_h2, nb_lev_h2_b, nb_lev_h2_bp, nb_lev_h2_ef,
		nb_lev_h2_cminus, nb_lev_h2_dminus, nb_lev_h2_a3, nb_lev_h2_c3, nb_lev_h2_d3, verbosity = 1;
	
	double t, enloss_mt_h2, enloss_ion_h2, enloss_ion_h2_rel, enloss_rot_h2, enloss_vibr_h2_v01, enloss_vibr_h2_v02,
		enloss_singlet_h2, enloss_singlet_h2_tot, enloss_diss_singlet_h2, enloss_diss_triplet_h2_b, enloss_triplet_h2_add, enloss_triplet_h2,
		enloss_mt_he, enloss_ion_he, enloss_ion_he_rel, enloss_coulomb;

	double cs_vibr_h2_v01, cs_vibr_h2_v02, cs_vibr_h2_v03, cs_singlet_h2, cs_singlet_h2_tot, cs_diss_singlet_h2, cs_diss_triplet_h2_b,
		cs_triplet_h2_add, cs_triplet_h2;
	double* electron_energies, * electron_velocities;

	string path, fname;
	ofstream output;
	stringstream ss;

	t = pow(10., 1. / nb_of_bins_per_order);

	nb_of_energies = (int)(nb_of_bins_per_order * log10(energy_max / energy_min)) + 1;
	electron_energies = new double[nb_of_energies];

	electron_energies[0] = energy_min;
	for (i = 1; i < nb_of_energies; i++) {
		electron_energies[i] = electron_energies[i - 1] * t;
	}

	electron_velocities = new double[nb_of_energies];
	for (i = 0; i < nb_of_energies; i++) {
		electron_velocities[i] = SPEED_OF_LIGHT * sqrt(2. * ELECTRON_MASS_EV * electron_energies[i] + electron_energies[i] * electron_energies[i])
			/ (ELECTRON_MASS_EV + electron_energies[i]);  // [cm/s]
	}
	// H2 energy levels
		// new H2 molecule data (Roueff et al. 2019), max number of levels 302,
	nb_lev_h2 = 302;
	molecule h2_mol("H2", isotope = 1, 2. * ATOMIC_MASS_UNIT);

	h2_diagram_roueff2019* h2_di
		= new h2_diagram_roueff2019(data_path, h2_mol, nb_lev_h2, verbosity);

	nb_lev_h2_b = 882;
	h2_diagram* h2_di_b
		= new h2_diagram(electronic_state = 1, data_path, h2_mol, nb_lev_h2_b, verbosity);

	nb_lev_h2_bp = 108;
	h2_diagram* h2_di_bp
		= new h2_diagram(electronic_state = 4, data_path, h2_mol, nb_lev_h2_bp, verbosity);

	// max angular momentum jmax = 26
	nb_lev_h2_ef = 793;
	h2_diagram* h2_di_ef
		= new h2_diagram(electronic_state = 10, data_path, h2_mol, nb_lev_h2_ef, verbosity);

	// we does not take into account the difference in energy levels of C- and C+ states
	nb_lev_h2_cminus = 251;
	h2_diagram* h2_di_cminus
		= new h2_diagram(electronic_state = 3, data_path, h2_mol, nb_lev_h2_cminus, verbosity);

	nb_lev_h2_dminus = 336;
	h2_diagram* h2_di_dminus
		= new h2_diagram(electronic_state = 6, data_path, h2_mol, nb_lev_h2_dminus, verbosity);

	nb_lev_h2_a3 = 384;
	h2_diagram* h2_di_a3
		= new h2_diagram(electronic_state = 101, data_path, h2_mol, nb_lev_h2_a3, verbosity);

	nb_lev_h2_c3 = 419;
	h2_diagram* h2_di_c3
		= new h2_diagram(electronic_state = 102, data_path, h2_mol, nb_lev_h2_c3, verbosity);

	nb_lev_h2_d3 = 412;
	h2_diagram* h2_di_d3
		= new h2_diagram(electronic_state = 103, data_path, h2_mol, nb_lev_h2_d3, verbosity);

	// H2 cross sections,
	// Initialization of cross section data for elastic scattering:
	cross_section_table_vs1* elastic_h2_el_cs
		= new cross_section_table_vs1(data_path, "elastic_scattering/e-h2_mt_cs.txt");

	// H2, electron impact ionization
	h2_electron_ionization_data h2_ioniz_data;

	electron_impact_ionization* h2_ioniz_cs
		= new electron_impact_ionization(h2_ioniz_data, verbosity);

	electron_impact_ionization_relativistic* h2_ioniz_cs_rel
		= new electron_impact_ionization_relativistic(h2_ioniz_data, verbosity);

	// H2 pure rotational excitation,
	// J = 0 -> 2, J = 1 -> 3
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-vi=0-rovibrational_excitation/Ni=";
	ss << ji;
	ss << "/MCCC-el-H2-X1Sg_vf=0_Nf=";
	ss << ji + 2;
	ss << ".X1Sg_vi=0_Ni=";
	ss << ji;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_rot_cs
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = false, verbosity);

	// H2 vibrational excitation, 
	// the sum of (v, J) = (0, ji) -> (1, ji) and (0, ji) -> (1, ji + 2)
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-rot-X-X/vi=0/Ji=";
	ss << ji;
	ss << "/MCCC-el-H2-X1Sg_vf=1_Jf=";
	ss << ji;
	ss << ".X1Sg_vi=0_Ji=";
	ss << ji;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_vibr_rot_cs_v01_dj0
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = false, verbosity);

	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-rot-X-X/vi=0/Ji=";
	ss << ji;
	ss << "/MCCC-el-H2-X1Sg_vf=1_Jf=";
	ss << (ji + 2);
	ss << ".X1Sg_vi=0_Ji=";
	ss << ji;
	ss << ".txt";
	fname = ss.str();

	cross_section_table_mccc* h2_vibr_rot_cs_v01_dj2
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = false, verbosity);

	// v = 0 -> 2,
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-rot-X-X/vi=0/Ji=";
	ss << ji;
	ss << "/MCCC-el-H2-X1Sg_vf=2_Jf=";
	ss << ji;
	ss << ".X1Sg_vi=0_Ji=";
	ss << ji;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_vibr_rot_cs_v02_dj0
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = false, verbosity);

	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-rot-X-X/vi=0/Ji=";
	ss << ji;
	ss << "/MCCC-el-H2-X1Sg_vf=2_Jf=";
	ss << (ji + 2);
	ss << ".X1Sg_vi=0_Ji=";
	ss << ji;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_vibr_rot_cs_v02_dj2
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = false, verbosity);

	// v = 0 -> 3,
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-rot-X-X/vi=0/Ji=";
	ss << ji;
	ss << "/MCCC-el-H2-X1Sg_vf=3_Jf=";
	ss << ji;
	ss << ".X1Sg_vi=0_Ji=";
	ss << ji;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_vibr_rot_cs_v03_dj0
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = false, verbosity);

	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-rot-X-X/vi=0/Ji=";
	ss << ji;
	ss << "/MCCC-el-H2-X1Sg_vf=3_Jf=";
	ss << (ji + 2);
	ss << ".X1Sg_vi=0_Ji=";
	ss << ji;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_vibr_rot_cs_v03_dj2
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = false, verbosity);

	// v = 0 -> 4,
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-rot-X-X/vi=0/Ji=";
	ss << ji;
	ss << "/MCCC-el-H2-X1Sg_vf=4_Jf=";
	ss << ji;
	ss << ".X1Sg_vi=0_Ji=";
	ss << ji;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_vibr_rot_cs_v04_dj0
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = false, verbosity);

	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-rot-X-X/vi=0/Ji=";
	ss << ji;
	ss << "/MCCC-el-H2-X1Sg_vf=4_Jf=";
	ss << (ji + 2);
	ss << ".X1Sg_vi=0_Ji=";
	ss << ji;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_vibr_rot_cs_v04_dj2
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = false, verbosity);

	// H2 electronic excitation (singlet)
	// 1Sg+(X) -> 1Su+(B), 
	cross_section_table_mccc** h2_bstate_cs
		= new cross_section_table_mccc * [6 * MAX_NB_H2_VSTATES_B1SU];

	for (vf = 0; vf < MAX_NB_H2_VSTATES_B1SU; vf++) {
		for (dj = -5; dj <= 5; dj += 2) {
			i = 6 * vf + (dj + 5) / 2;
			h2_bstate_cs[i] = 0;

			ss.clear();
			ss.str("");
			ss << "coll_h2/MCCC-el-H2-X1Sg-excitation-j-resolved/X1Sg-B1Su/vi=";
			ss << vi;
			ss << "/Ji=";
			ss << ji;
			ss << "/MCCC-el-H2-B1Su_vf=";
			ss << vf;
			ss << "_Jf=";
			ss << ji + dj;
			ss << ".X1Sg_vi=";
			ss << vi;
			ss << "_Ji=";
			ss << ji;
			ss << ".txt";

			fname = data_path + ss.str();
			ifstream f(fname.c_str());

			// check if file with data exists,  
			if (f.good()) {
				fname = ss.str();
				h2_bstate_cs[i] = new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, false);
			}
		}
	}

	// S1g+(X) -> S1u+(Bp)
	cross_section_table_mccc** h2_bpstate_cs
		= new cross_section_table_mccc * [6 * MAX_NB_H2_VSTATES_BP1SU];

	for (vf = 0; vf < MAX_NB_H2_VSTATES_BP1SU; vf++) {
		for (dj = -5; dj <= 5; dj += 2) {
			i = 6 * vf + (dj + 5) / 2;
			h2_bpstate_cs[i] = 0;

			ss.clear();
			ss.str("");
			ss << "coll_h2/MCCC-el-H2-X1Sg-excitation-j-resolved/X1Sg-Bp1Su/vi=";
			ss << vi;
			ss << "/Ji=";
			ss << ji;
			ss << "/MCCC-el-H2-Bp1Su_vf=";
			ss << vf;
			ss << "_Jf=";
			ss << ji + dj;
			ss << ".X1Sg_vi=";
			ss << vi;
			ss << "_Ji=";
			ss << ji;
			ss << ".txt";

			fname = data_path + ss.str();
			ifstream f(fname.c_str());

			// check if file with data exists,  
			if (f.good()) {
				fname = ss.str();
				h2_bpstate_cs[i] = new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, false);
			}
		}
	}

	// S1g+(X) -> P1u(C-/+)
	cross_section_table_mccc** h2_cstate_cs
		= new cross_section_table_mccc * [11 * MAX_NB_H2_VSTATES_C1PU];

	for (vf = 0; vf < MAX_NB_H2_VSTATES_C1PU; vf++) {
		for (dj = -5; dj <= 5; dj += 1) {
			i = 11 * vf + dj + 5;
			h2_cstate_cs[i] = 0;

			ss.clear();
			ss.str("");
			ss << "coll_h2/MCCC-el-H2-X1Sg-excitation-j-resolved/X1Sg-C1Pu/vi=";
			ss << vi;
			ss << "/Ji=";
			ss << ji;
			ss << "/MCCC-el-H2-C1Pu_vf=";
			ss << vf;
			ss << "_Jf=";
			ss << ji + dj;
			ss << ".X1Sg_vi=";
			ss << vi;
			ss << "_Ji=";
			ss << ji;
			ss << ".txt";

			fname = data_path + ss.str();
			ifstream f(fname.c_str());

			// check if file with data exists,  
			if (f.good()) {
				fname = ss.str();
				h2_cstate_cs[i] = new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, false);
			}
		}
	}

	// S1g+(X) -> P1u(D-/+)
	cross_section_table_mccc** h2_dstate_cs
		= new cross_section_table_mccc * [11 * MAX_NB_H2_VSTATES_D1PU];

	for (vf = 0; vf < MAX_NB_H2_VSTATES_D1PU; vf++) {
		for (dj = -5; dj <= 5; dj += 1) {
			i = 11 * vf + dj + 5;
			h2_dstate_cs[i] = 0;

			ss.clear();
			ss.str("");
			ss << "coll_h2/MCCC-el-H2-X1Sg-excitation-j-resolved/X1Sg-D1Pu/vi=";
			ss << vi;
			ss << "/Ji=";
			ss << ji;
			ss << "/MCCC-el-H2-D1Pu_vf=";
			ss << vf;
			ss << "_Jf=";
			ss << ji + dj;
			ss << ".X1Sg_vi=";
			ss << vi;
			ss << "_Ji=";
			ss << ji;
			ss << ".txt";

			fname = data_path + ss.str();
			ifstream f(fname.c_str());

			// check if file with data exists,  
			if (f.good()) {
				fname = ss.str();
				h2_dstate_cs[i] = new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, false);
			}
		}
	}

	// 1Sg+ (X) -> 1Sg (EF)
	cross_section_table_mccc** h2_efstate_cs
		= new cross_section_table_mccc * [5 * MAX_NB_H2_VSTATES_EF1SG];

	for (vf = 0; vf < MAX_NB_H2_VSTATES_EF1SG; vf++) {
		for (dj = -4; dj <= 4; dj += 2)
		{
			i = 5 * vf + (dj + 4) / 2;
			h2_efstate_cs[i] = 0;

			ss.clear();
			ss.str("");
			ss << "coll_h2/MCCC-el-H2-X1Sg-excitation-j-resolved/X1Sg-EF1Sg/vi=";
			ss << vi;
			ss << "/Ji=";
			ss << ji;
			ss << "/MCCC-el-H2-EF1Sg_vf=";
			ss << vf;
			ss << "_Jf=";
			ss << ji + dj;
			ss << ".X1Sg_vi=";
			ss << vi;
			ss << "_Ji=";
			ss << ji;
			ss << ".txt";

			fname = data_path + ss.str();
			ifstream f(fname.c_str());

			// check if file with data exists,  
			if (f.good()) {
				fname = ss.str();
				h2_efstate_cs[i] = new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, false);
			}
		}
	}

	//
	// Electronic dissociative excitation of H2
	// 1Sg+(X) -> 1Su+(B),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-B1Su_DE.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_bstate_diss_cs
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, verbosity);

	// 1Sg+(X) -> 1Su+(Bp),	
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-Bp1Su_DE.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_bpstate_diss_cs
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, verbosity);

	// 1Sg+(X) -> 1Sg+(EF), about 0.007	 of the total direct dissociation energy loss
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-EF1Sg_DE.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_efstate_diss_cs
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, verbosity);

	// 1Sg+(X) -> 1Pu(C),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-C1Pu_DE.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_cstate_diss_cs
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, verbosity);

	// 1Sg+(X) -> 1Pu(D),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-D1Pu_DE.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_dstate_diss_cs
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, verbosity);

	// dissociative excitation to triplet states,
	// S1g(X) -> S3u+(b) (dissociative state),
	ss.clear();
	ss.str("");
	
	if (vi == 0) {
		// several lines with data are commented in the file, in order to make energy threshold equal to 7.0 = 4.5 + 2.5 eV;
		// 2.5 is minimal kinetic energy of H atoms released for excitation from v = 0; Trevisan & Tennyson, Plasma Phys. Control. Fusion 44, 1263 (2002);
		ss << "coll_h2/MCCC-el-H2-b3Su_DE.X1Sg_vi=";
		ss << vi;
		ss << "_changed.txt";
	}
	else {
		ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
		ss << vi;
		ss << "/MCCC-el-H2-b3Su_DE.X1Sg_vi=";
		ss << vi;
		ss << ".txt";
	}

	fname = ss.str();
	cross_section_table_mccc* h2_3bstate_diss_cs
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, verbosity);

	//
	// Triplet states
	// 1Sg(X) -> 3Sg+(a)
	// Excitation to bound state, without excitation to dissociation continuum that is 0.1 per cent,
	cross_section_table_mccc** h2_3astate_cs
		= new cross_section_table_mccc * [5 * MAX_H2_VSTATES_a3Sg];

	for (vf = 0; vf < MAX_H2_VSTATES_a3Sg; vf++) {
		for (dj = -4; dj <= 4; dj += 2)
		{
			i = 5 * vf + (dj + 4) / 2;
			h2_3astate_cs[i] = 0;

			ss.clear();
			ss.str("");
			ss << "coll_h2/MCCC-el-H2-X1Sg-excitation-j-resolved/X1Sg-a3Sg/vi=";
			ss << vi;
			ss << "/Ji=";
			ss << ji;
			ss << "/MCCC-el-H2-a3Sg_vf=";
			ss << vf;
			ss << "_Jf=";
			ss << ji + dj;
			ss << ".X1Sg_vi=";
			ss << vi;
			ss << "_Ji=";
			ss << ji;
			ss << ".txt";

			fname = data_path + ss.str();
			ifstream f(fname.c_str());

			// check if file with data exists,  
			if (f.good()) {
				fname = ss.str();
				h2_3astate_cs[i] = new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, false);
			}
		}
	}

	// 1Sg(X) -> 3Pu (c) 
	// Excitation to bound state, without excitation to dissiociation continuum, that is less than 0.5 per cent,
	cross_section_table_mccc** h2_3cstate_cs
		= new cross_section_table_mccc * [11 * MAX_H2_VSTATES_c3Pu];

	for (vf = 0; vf < MAX_H2_VSTATES_c3Pu; vf++) {
		for (dj = -5; dj <= 5; dj++)
		{
			i = 11 * vf + dj + 5;
			h2_3cstate_cs[i] = 0;

			ss.clear();
			ss.str("");
			ss << "coll_h2/MCCC-el-H2-X1Sg-excitation-j-resolved/X1Sg-c3Pu/vi=";
			ss << vi;
			ss << "/Ji=";
			ss << ji;
			ss << "/MCCC-el-H2-c3Pu_vf=";
			ss << vf;
			ss << "_Jf=";
			ss << ji + dj;
			ss << ".X1Sg_vi=";
			ss << vi;
			ss << "_Ji=";
			ss << ji;
			ss << ".txt";

			fname = data_path + ss.str();
			ifstream f(fname.c_str());

			// check if file with data exists,  
			if (f.good()) {
				fname = ss.str();
				h2_3cstate_cs[i] = new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, false);
			}
		}
	}

	// 1Sg(X) -> 3Pu (d) 
	// Excitation to bound state, without excitation to dissociation continuum that is about 1 per cent,
	cross_section_table_mccc** h2_3dstate_cs
		= new cross_section_table_mccc * [11 * MAX_H2_VSTATES_d3Pu];

	for (vf = 0; vf < MAX_H2_VSTATES_d3Pu; vf++) {
		for (dj = -5; dj <= 5; dj++)
		{
			i = 11 * vf + dj + 5;
			h2_3dstate_cs[i] = 0;

			ss.clear();
			ss.str("");
			ss << "coll_h2/MCCC-el-H2-X1Sg-excitation-j-resolved/X1Sg-d3Pu/vi=";
			ss << vi;
			ss << "/Ji=";
			ss << ji;
			ss << "/MCCC-el-H2-d3Pu_vf=";
			ss << vf;
			ss << "_Jf=";
			ss << ji + dj;
			ss << ".X1Sg_vi=";
			ss << vi;
			ss << "_Ji=";
			ss << ji;
			ss << ".txt";

			fname = data_path + ss.str();
			ifstream f(fname.c_str());

			// check if file with data exists,  
			if (f.good()) {
				fname = ss.str();
				h2_3dstate_cs[i] = new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, false);
			}
		}
	}

	// There are not ro-vibrationally resolved data for states below,
	// States below are not dissociative, but are treated as dissociative in order to estimate their effect
	// 1Sg(X) -> 3Su+(e),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-e3Su_total.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_3estate_cs
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, verbosity);

	// 1Sg(X) -> 3Sg+(h),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-h3Sg_total.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_3hstate_cs
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, verbosity);

	// 1Sg(X) -> 3Sg+(g),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-g3Sg_total.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_3gstate_cs
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, verbosity);

	// 1Sg(X) -> 3Pg(i),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-i3Pg_total.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_3istate_cs
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, verbosity);

	// 1Sg(X) -> 3Dg(j),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-j3Dg_total.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_3jstate_cs
		= new cross_section_table_mccc(data_path, fname, is_extrapolation_on = true, verbosity);


	// Saving energy losses:
	fname = output_path + "energy_loss_electron_h2_ji";
	fname += to_string(ji);
	fname += "_x";
	fname += str;
	fname += ".txt";
	output.open(fname.c_str());

	output << left << "!Energy losses in [cm2 eV] of electron in pure H2 gas (normalized to one H2 molecule), initial level j = " << ji << endl
		<< "!for MCCC cross sections - at higher energies, the cross sections are extrapolated as a power law, cs = cs_0*(E/E_0)^(gamma);" << endl
		<< "!Coulomb losses - Swartz et al., J. of Geophys. Res. 76, p. 8425 (1971); ionization fraction: " << ionization_fraction << endl
		<< "!Vibrational excitation v = 0 -> 1: the sum of (v,j) = (0,j) -> (1,j) and (0,j) -> (1,j+2); the same for v = 0 -> 2; 0 -> 3 and etc.;" << endl
		<< "!tot(exc-singl) - total excitation (without direct dissociation) of singlet states, B, Bp, EF, C, D;" << endl;

	output << left << setw(14) << "!Energy(eV) "
		<< setw(12) << "m.t."
		<< setw(12) << "ioniz"
		<< setw(12) << "ioniz_rel"
		<< setw(15) << "ji-ji+2"
		<< setw(12) << "v=0-1"
		<< setw(15) << "v=0-2"
		<< setw(12) << "1Su+(B) "
		<< setw(12) << "1Su+(Bp) "
		<< setw(12) << "1Sg+(EF) "
		<< setw(12) << "1Pu(C) "
		<< setw(12) << "1Pu(D) "
		<< setw(15) << "tot(exc-singl)"
		<< setw(12) << "diss-singl "
		<< setw(12) << "3Su+(b) "
		<< setw(12) << "tripl(a,c,d) "
		<< setw(12) << "tripl(add) "
		<< setw(12) << "coulomb" << endl;

	for (i = 0; i < nb_of_energies; i++)
	{
		output.precision(6);
		output << left << setw(14) << electron_energies[i];
		output.precision(4);

		enloss_mt_h2 = (*elastic_h2_el_cs)(electron_energies[i]) * electron_energies[i] * ELECTRON_MASS / ATOMIC_MASS_UNIT;
		output << left << setw(12) << enloss_mt_h2;

		enloss_ion_h2 = h2_ioniz_cs->get_energy_loss(electron_energies[i]);
		enloss_ion_h2_rel = h2_ioniz_cs_rel->get_energy_loss(electron_energies[i]);

		output << left << setw(12) << enloss_ion_h2 << setw(15) << enloss_ion_h2_rel;

		lf = h2_di->get_nb(0, ji + 2);  // ji is small, 0 or 1
		enloss_rot_h2 = (*h2_rot_cs)(electron_energies[i]) * (h2_di->lev_array[lf].energy - h2_di->lev_array[ji].energy) * CM_INVERSE_TO_EV;
		output << left << setw(12) << enloss_rot_h2;

		// (v, J) = (0,j) -> (1,j) and (0,j) -> (1,j+2)
		lf = h2_di->get_nb(1, ji);
		enloss_vibr_h2_v01 = (*h2_vibr_rot_cs_v01_dj0)(electron_energies[i])
			* (h2_di->lev_array[lf].energy - h2_di->lev_array[ji].energy) * CM_INVERSE_TO_EV;

		lf = h2_di->get_nb(1, ji + 2);
		enloss_vibr_h2_v01 += (*h2_vibr_rot_cs_v01_dj2)(electron_energies[i])
			* (h2_di->lev_array[lf].energy - h2_di->lev_array[ji].energy) * CM_INVERSE_TO_EV;

		output << left << setw(12) << enloss_vibr_h2_v01;

		// v = 0 -> 2
		lf = h2_di->get_nb(2, ji);
		enloss_vibr_h2_v02 = (*h2_vibr_rot_cs_v02_dj0)(electron_energies[i])
			* (h2_di->lev_array[lf].energy - h2_di->lev_array[ji].energy) * CM_INVERSE_TO_EV;

		lf = h2_di->get_nb(2, ji + 2);
		enloss_vibr_h2_v02 += (*h2_vibr_rot_cs_v02_dj2)(electron_energies[i])
			* (h2_di->lev_array[lf].energy - h2_di->lev_array[ji].energy) * CM_INVERSE_TO_EV;

		output << left << setw(15) << enloss_vibr_h2_v02;

		// 1Su+(B)
		enloss_singlet_h2_tot = enloss_singlet_h2 = 0.;
		for (vf = 0; vf < MAX_NB_H2_VSTATES_B1SU; vf++) {
			for (dj = -5; dj <= 5; dj += 2)
			{
				k = 6 * vf + (dj + 5) / 2;
				lf = h2_di_b->get_nb(vf, ji + dj);

				if (h2_bstate_cs[k] != 0 && lf != -1) {
					enloss_singlet_h2 += (*h2_bstate_cs[k])(electron_energies[i])
						* (h2_di_b->lev_array[lf].energy - h2_di->lev_array[ji].energy);
				}
			}
		}
		enloss_singlet_h2 *= CM_INVERSE_TO_EV;  // to eV cm2

		enloss_singlet_h2_tot += enloss_singlet_h2;
		output << left << setw(12) << enloss_singlet_h2;

		// 1Su+(Bp)
		enloss_singlet_h2 = 0.;
		for (vf = 0; vf < MAX_NB_H2_VSTATES_BP1SU; vf++) {
			for (dj = -5; dj <= 5; dj += 2)
			{
				k = 6 * vf + (dj + 5) / 2;
				lf = h2_di_bp->get_nb(vf, ji + dj);

				if (h2_bpstate_cs[k] != 0 && lf != -1) {
					enloss_singlet_h2 += (*h2_bpstate_cs[k])(electron_energies[i])
						* (h2_di_bp->lev_array[lf].energy - h2_di->lev_array[ji].energy);
				}
			}
		}
		enloss_singlet_h2 *= CM_INVERSE_TO_EV;

		enloss_singlet_h2_tot += enloss_singlet_h2;
		output << left << setw(12) << enloss_singlet_h2;

		// 1Sg+(EF)
		enloss_singlet_h2 = 0.;
		for (vf = 0; vf < MAX_NB_H2_VSTATES_EF1SG; vf++) {
			for (dj = -4; dj <= 4; dj += 2)
			{
				k = 5 * vf + (dj + 4) / 2;
				lf = h2_di_ef->get_nb(vf, ji + dj);

				if (h2_efstate_cs[k] != 0 && lf != -1) {
					enloss_singlet_h2 += (*h2_efstate_cs[k])(electron_energies[i])
						* (h2_di_ef->lev_array[lf].energy - h2_di->lev_array[ji].energy);
				}
			}
		}
		enloss_singlet_h2 *= CM_INVERSE_TO_EV;

		enloss_singlet_h2_tot += enloss_singlet_h2;
		output << left << setw(12) << enloss_singlet_h2;

		// 1Pu(C)
		enloss_singlet_h2 = 0.;
		for (vf = 0; vf < MAX_NB_H2_VSTATES_C1PU; vf++) {
			for (dj = -5; dj <= 5; dj += 1)
			{
				k = 11 * vf + dj + 5;
				lf = h2_di_cminus->get_nb(vf, ji + dj);

				if (h2_cstate_cs[k] != 0 && lf != -1) {
					enloss_singlet_h2 += (*h2_cstate_cs[k])(electron_energies[i])
						* (h2_di_cminus->lev_array[lf].energy - h2_di->lev_array[ji].energy);
				}
			}
		}
		enloss_singlet_h2 *= CM_INVERSE_TO_EV;

		enloss_singlet_h2_tot += enloss_singlet_h2;
		output << left << setw(12) << enloss_singlet_h2;

		// 1Pu(D)
		enloss_singlet_h2 = 0.;
		for (vf = 0; vf < MAX_NB_H2_VSTATES_D1PU; vf++) {
			for (dj = -5; dj <= 5; dj += 1)
			{
				k = 11 * vf + dj + 5;
				lf = h2_di_dminus->get_nb(vf, ji + dj);

				if (h2_dstate_cs[k] != 0 && lf != -1) {
					enloss_singlet_h2 += (*h2_dstate_cs[k])(electron_energies[i])
						* (h2_di_dminus->lev_array[lf].energy - h2_di->lev_array[ji].energy);
				}
			}
		}
		enloss_singlet_h2 *= CM_INVERSE_TO_EV;

		enloss_singlet_h2_tot += enloss_singlet_h2;
		output << left << setw(12) << enloss_singlet_h2 << setw(15) << enloss_singlet_h2_tot;

		// direct dissociation through singlet states,
		enloss_diss_singlet_h2 = (*h2_bstate_diss_cs)(electron_energies[i]) * h2_bstate_diss_cs->get_threshold_energy();
		enloss_diss_singlet_h2 += (*h2_bpstate_diss_cs)(electron_energies[i]) * h2_bpstate_diss_cs->get_threshold_energy();
		enloss_diss_singlet_h2 += (*h2_efstate_diss_cs)(electron_energies[i]) * h2_efstate_diss_cs->get_threshold_energy();

		enloss_diss_singlet_h2 += (*h2_cstate_diss_cs)(electron_energies[i]) * h2_cstate_diss_cs->get_threshold_energy();
		enloss_diss_singlet_h2 += (*h2_dstate_diss_cs)(electron_energies[i]) * h2_dstate_diss_cs->get_threshold_energy();
		output << left << setw(12) << enloss_diss_singlet_h2;

		// dissociation through 3Su+(b)
		enloss_diss_triplet_h2_b = (*h2_3bstate_diss_cs)(electron_energies[i]) * h2_3bstate_diss_cs->get_threshold_energy();
		output << left << setw(12) << enloss_diss_triplet_h2_b;

		// excitation of triplet states (a, c, d)
		enloss_triplet_h2 = 0.;
		for (vf = 0; vf < MAX_H2_VSTATES_a3Sg; vf++) {
			for (dj = -4; dj <= 4; dj += 2)
			{
				k = 5 * vf + (dj + 4) / 2;
				lf = h2_di_a3->get_nb(vf, ji + dj);

				if (h2_3astate_cs[k] != 0 && lf != -1) {
					enloss_triplet_h2 += (*h2_3astate_cs[k])(electron_energies[i])
						* (h2_di_a3->lev_array[lf].energy - h2_di->lev_array[ji].energy);
				}
			}
		}
		for (vf = 0; vf < MAX_H2_VSTATES_c3Pu; vf++) {
			for (dj = -5; dj <= 5; dj++)
			{
				k = 11 * vf + dj + 5;
				lf = h2_di_c3->get_nb(vf, ji + dj);

				if (h2_3cstate_cs[k] != 0 && lf != -1) {
					enloss_triplet_h2 += (*h2_3cstate_cs[k])(electron_energies[i])
						* (h2_di_c3->lev_array[lf].energy - h2_di->lev_array[ji].energy);
				}
			}
		}
		for (vf = 0; vf < MAX_H2_VSTATES_d3Pu; vf++) {
			for (dj = -5; dj <= 5; dj++)
			{
				k = 11 * vf + dj + 5;
				lf = h2_di_d3->get_nb(vf, ji + dj);

				if (h2_3dstate_cs[k] != 0 && lf != -1) {
					enloss_triplet_h2 += (*h2_3dstate_cs[k])(electron_energies[i])
						* (h2_di_d3->lev_array[lf].energy - h2_di->lev_array[ji].energy);
				}
			}
		}
		enloss_triplet_h2 *= CM_INVERSE_TO_EV;
		output << left << setw(12) << enloss_triplet_h2;

		// dissociation through 3Su+(e), 3Sg+(h), 3Sg+(g), ...
		enloss_triplet_h2_add = (*h2_3estate_cs)(electron_energies[i]) * h2_3estate_cs->get_threshold_energy()
			+ (*h2_3hstate_cs)(electron_energies[i]) * h2_3hstate_cs->get_threshold_energy()
			+ (*h2_3gstate_cs)(electron_energies[i]) * h2_3gstate_cs->get_threshold_energy()
			+ (*h2_3istate_cs)(electron_energies[i]) * h2_3istate_cs->get_threshold_energy()
			+ (*h2_3jstate_cs)(electron_energies[i]) * h2_3jstate_cs->get_threshold_energy();

		output << left << setw(12) << enloss_triplet_h2_add;


		// energy loss in [eV cm2], per H2, number density of electrons ne = 2. * ionization_fraction * n_H2
		enloss_coulomb = 2. * calc_coulomb_losses_thermal_electrons(electron_energies[i], ionization_fraction, thermal_electron_temp)
			/ electron_velocities[i];

		output << left << setw(12) << enloss_coulomb;
		output << endl;
	}
	output.close();

	//
	// Saving cross sections,
	fname = output_path + "cs_electron_h2_ji";
	fname += to_string(ji);
	fname += ".txt";
	output.open(fname.c_str());

	output << left << "!Cross sections in [cm2] of electron interaction with the H2, initial level j = " << ji << endl
		<< "!for MCCC cross sections - at higher energies, the cross sections are extrapolated as a power law, cs = cs_0*(E/E_0)^(gamma)" << endl
		<< setw(14) << "!Energy(eV) "
		<< setw(12) << "m.t."
		<< setw(12) << "ion"
		<< setw(15) << "ion_rel"
		<< setw(12) << "J=0-2"
		<< setw(12) << "v=0-1"
		<< setw(12) << "v=0-2"
		<< setw(15) << "v=0-3"
		<< setw(12) << "1Su+(B) "
		<< setw(12) << "1Su+(Bp) "
		<< setw(12) << "1Sg+(EF) "
		<< setw(12) << "1Pu(C) "
		<< setw(12) << "1Pu(D) "
		<< setw(15) << "total "
		<< setw(12) << "diss-singl "
		<< setw(12) << "3Su+(b) "
		<< setw(12) << "tripl(a,c,d) "
		<< setw(12) << "tripl(add) " << endl;

	for (i = 0; i < nb_of_energies; i++)
	{
		output.precision(6);
		output << left << setw(14) << electron_energies[i];

		output.precision(4);
		output << left << setw(12) << (*elastic_h2_el_cs)(electron_energies[i])
			<< setw(12) << h2_ioniz_cs->get_int_cs(electron_energies[i])
			<< setw(15) << h2_ioniz_cs_rel->get_int_cs(electron_energies[i]);

		output << left << setw(12) << (*h2_rot_cs)(electron_energies[i]);

		// (v, J) = (0,j) -> (1,j) and (0,j) -> (1,j+2)
		cs_vibr_h2_v01 = (*h2_vibr_rot_cs_v01_dj0)(electron_energies[i]);
		cs_vibr_h2_v01 += (*h2_vibr_rot_cs_v01_dj2)(electron_energies[i]);
		output << left << setw(12) << cs_vibr_h2_v01;

		// v = 0 -> 2
		cs_vibr_h2_v02 = (*h2_vibr_rot_cs_v02_dj0)(electron_energies[i]);
		cs_vibr_h2_v02 += (*h2_vibr_rot_cs_v02_dj2)(electron_energies[i]);
		output << left << setw(12) << cs_vibr_h2_v02;

		// v = 0 -> 3
		cs_vibr_h2_v03 = (*h2_vibr_rot_cs_v03_dj0)(electron_energies[i]);
		cs_vibr_h2_v03 += (*h2_vibr_rot_cs_v03_dj2)(electron_energies[i]);
		output << left << setw(15) << cs_vibr_h2_v03;

		cs_singlet_h2_tot = cs_singlet_h2 = 0.;
		for (vf = 0; vf < MAX_NB_H2_VSTATES_B1SU; vf++) {
			for (dj = -5; dj <= 5; dj += 2)
			{
				k = 6 * vf + (dj + 5) / 2;
				if (h2_bstate_cs[k] != 0) {
					cs_singlet_h2 += (*h2_bstate_cs[k])(electron_energies[i]);
				}
			}
		}
		cs_singlet_h2_tot += cs_singlet_h2;
		output << left << setw(12) << cs_singlet_h2;

		cs_singlet_h2 = 0.;
		for (vf = 0; vf < MAX_NB_H2_VSTATES_BP1SU; vf++) {
			for (dj = -5; dj <= 5; dj += 2)
			{
				k = 6 * vf + (dj + 5) / 2;
				if (h2_bpstate_cs[k] != 0) {
					cs_singlet_h2 += (*h2_bpstate_cs[k])(electron_energies[i]);
				}
			}
		}
		cs_singlet_h2_tot += cs_singlet_h2;
		output << left << setw(12) << cs_singlet_h2;

		cs_singlet_h2 = 0.;
		for (vf = 0; vf < MAX_NB_H2_VSTATES_EF1SG; vf++) {
			for (dj = -4; dj <= 4; dj += 2)
			{
				k = 5 * vf + (dj + 4) / 2;
				if (h2_efstate_cs[k] != 0) {
					cs_singlet_h2 += (*h2_efstate_cs[k])(electron_energies[i]);
				}
			}
		}
		cs_singlet_h2_tot += cs_singlet_h2;
		output << left << setw(12) << cs_singlet_h2;

		cs_singlet_h2 = 0.;
		for (vf = 0; vf < MAX_NB_H2_VSTATES_C1PU; vf++) {
			for (dj = -5; dj <= 5; dj += 1)
			{
				k = 11 * vf + dj + 5;
				if (h2_cstate_cs[k] != 0) {
					cs_singlet_h2 += (*h2_cstate_cs[k])(electron_energies[i]);
				}
			}
		}
		cs_singlet_h2_tot += cs_singlet_h2;
		output << left << setw(12) << cs_singlet_h2;

		cs_singlet_h2 = 0.;
		for (vf = 0; vf < MAX_NB_H2_VSTATES_D1PU; vf++) {
			for (dj = -5; dj <= 5; dj += 1)
			{
				k = 11 * vf + dj + 5;
				if (h2_dstate_cs[k] != 0) {
					cs_singlet_h2 += (*h2_dstate_cs[k])(electron_energies[i]);
				}
			}
		}
		cs_singlet_h2_tot += cs_singlet_h2;
		output << left << setw(12) << cs_singlet_h2 << setw(15) << cs_singlet_h2_tot;

		//
		cs_diss_singlet_h2 = (*h2_bstate_diss_cs)(electron_energies[i]);
		cs_diss_singlet_h2 += (*h2_bpstate_diss_cs)(electron_energies[i]);
		cs_diss_singlet_h2 += (*h2_efstate_diss_cs)(electron_energies[i]);

		cs_diss_singlet_h2 += (*h2_cstate_diss_cs)(electron_energies[i]);
		cs_diss_singlet_h2 += (*h2_dstate_diss_cs)(electron_energies[i]);
		output << left << setw(12) << cs_diss_singlet_h2;

		// triplet b
		cs_diss_triplet_h2_b = (*h2_3bstate_diss_cs)(electron_energies[i]);
		output << left << setw(12) << cs_diss_triplet_h2_b;

		// a, c, d
		cs_triplet_h2 = 0.;
		for (vf = 0; vf < MAX_H2_VSTATES_a3Sg; vf++) {
			for (dj = -4; dj <= 4; dj += 2)
			{
				k = 5 * vf + (dj + 4) / 2;
				if (h2_3astate_cs[k] != 0) {
					cs_triplet_h2 += (*h2_3astate_cs[k])(electron_energies[i]);
				}
			}
		}
		for (vf = 0; vf < MAX_H2_VSTATES_c3Pu; vf++) {
			for (dj = -5; dj <= 5; dj++)
			{
				k = 11 * vf + dj + 5;
				if (h2_3cstate_cs[k] != 0) {
					cs_triplet_h2 += (*h2_3cstate_cs[k])(electron_energies[i]);
				}
			}
		}
		for (vf = 0; vf < MAX_H2_VSTATES_d3Pu; vf++) {
			for (dj = -5; dj <= 5; dj++)
			{
				k = 11 * vf + dj + 5;
				if (h2_3dstate_cs[k] != 0) {
					cs_triplet_h2 += (*h2_3dstate_cs[k])(electron_energies[i]);
				}
			}
		}
		output << left << setw(12) << cs_triplet_h2;

		// triplet e
		cs_triplet_h2_add = (*h2_3estate_cs)(electron_energies[i])
			+ (*h2_3hstate_cs)(electron_energies[i])
			+ (*h2_3gstate_cs)(electron_energies[i])
			+ (*h2_3istate_cs)(electron_energies[i])
			+ (*h2_3jstate_cs)(electron_energies[i]);

		output << left << setw(12) << cs_triplet_h2_add;

		output << endl;
	}
	output.close();

	// Saving cross sections,
	// only vibrational excitations,
	fname = output_path + "cs_electron_h2_rovibr_ji";
	fname += to_string(ji);
	fname += ".txt";
	output.open(fname.c_str());

	output << left << "!Cross sections in [cm2] of electron interaction with the H2, initial level j = " << ji << endl
		<< "!for MCCC cross sections - at higher energies, the cross sections are extrapolated as a power law, cs = cs_0*(E/E_0)^(gamma)" << endl
		<< "!(v_low - v_up) j_low - j_up" << endl
		<< setw(14) << "!Energy(eV) "
		<< setw(12) << "(0-1)j-j"
		<< setw(12) << "(0-1)j-j+2"
		<< setw(12) << "(0-2)j-j"
		<< setw(12) << "(0-2)j-j+2"
		<< setw(12) << "(0-3)j-j"
		<< setw(12) << "(0-3)j-j+2"
		<< setw(12) << "(0-4)j-j"
		<< setw(12) << "(0-4)j-j+2" << endl;

	for (i = 0; i < nb_of_energies / 2; i++)
	{
		output.precision(6);
		output << left << setw(14) << electron_energies[i];
		output.precision(3);

		// (v, J) = (0,j) -> (1,j) and (0,j) -> (1,j+2)
		output << left << setw(12) << (*h2_vibr_rot_cs_v01_dj0)(electron_energies[i])
			<< setw(12) << (*h2_vibr_rot_cs_v01_dj2)(electron_energies[i]);

		// v = 0 -> 2
		output << left << setw(12) << (*h2_vibr_rot_cs_v02_dj0)(electron_energies[i])
			<< setw(12) << (*h2_vibr_rot_cs_v02_dj2)(electron_energies[i]);

		// v = 0 -> 3
		output << left << setw(12) << (*h2_vibr_rot_cs_v03_dj0)(electron_energies[i])
			<< setw(12) << (*h2_vibr_rot_cs_v03_dj2)(electron_energies[i]);

		// v = 0 -> 4
		output << left << setw(12) << (*h2_vibr_rot_cs_v04_dj0)(electron_energies[i])
			<< setw(12) << (*h2_vibr_rot_cs_v04_dj2)(electron_energies[i]);

		output << endl;
	}
	output.close();

	//
	// Helium
	// He, momentum transfer cs
	cross_section_table_vs1* elastic_he_el_cs
		= new cross_section_table_vs1(data_path, "elastic_scattering/e-he_mt_cs.txt");

	// He, electron impact ionization
	he_electron_ionization_data he_ioniz_data;

	electron_impact_ionization* he_ioniz_cs
		= new electron_impact_ionization(he_ioniz_data, verbosity);

	electron_impact_ionization_relativistic* he_ioniz_cs_rel
		= new electron_impact_ionization_relativistic(he_ioniz_data, verbosity);

	fname = output_path + "energy_loss_electron_he.txt";
	output.open(fname.c_str());

	output << left << "!Energy losses in [cm2 eV] of electron in He," << endl
		<< setw(14) << "!Energy(eV) "
		<< setw(12) << "He(mt) "
		<< setw(12) << "He(ion) "
		<< setw(12) << "He(ion_rel) " << endl;

	for (i = 0; i < nb_of_energies; i++)
	{
		output.precision(6);
		output << left << setw(14) << electron_energies[i];
		output.precision(3);

		enloss_mt_he = (*elastic_he_el_cs)(electron_energies[i]) * electron_energies[i] * 0.5 * ELECTRON_MASS / ATOMIC_MASS_UNIT;
		output << left << setw(12) << enloss_mt_he;

		enloss_ion_he = he_ioniz_cs->get_energy_loss(electron_energies[i]);
		enloss_ion_he_rel = he_ioniz_cs_rel->get_energy_loss(electron_energies[i]);

		output << left << setw(12) << enloss_ion_he << setw(12) << enloss_ion_he_rel;
		output << endl;
	}
	output.close();

	fname = output_path + "cs_electron_he.txt";
	output.open(fname.c_str());

	output << left << "!Cross sections in [cm2] of electron interaction with the He," << endl
		<< setw(14) << "!Energy(eV) "
		<< setw(12) << "He(mt) "
		<< setw(12) << "He(ion) "
		<< setw(12) << "He(ion_rel) " << endl;

	for (i = 0; i < nb_of_energies; i++)
	{
		output.precision(6);
		output << left << setw(14) << electron_energies[i];

		output.precision(3);
		output << left << setw(12) << (*elastic_he_el_cs)(electron_energies[i])
			<< setw(12) << he_ioniz_cs->get_int_cs(electron_energies[i])
			<< setw(12) << he_ioniz_cs_rel->get_int_cs(electron_energies[i]) << endl;
	}
	output.close();
}
