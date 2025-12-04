#pragma once
#include <cstring>
#include <vector>
#include "spectroscopy.h"

// Number of excited electronic states taken into account,
// order - B, C, Bp, D, EF, Bpp, Dp, b, a, c, d
#define NB_EXC_ELECTRONIC_STATES 11

// Maximal number of rotational levels of the ground electronic state X1Sg (for vibrational state v = 0),
#define MAX_NB_J_H2 32  // 0,1,..,31

// Nb of vibrational states of electronic state taken into account,
// this nb is reconciled with CLOUDY data (reference in CLOUDY files, Abgrall et al., A&AS, 141, 297-300, 2000);
#define MAX_NB_H2_VSTATES_X1SG 15   // 0,1,..,14 (including) in CLOUDY, the same as in Dabrowski, Can. J. Phys. 62, p. 1639, 1984
#define MAX_NB_H2_VSTATES_B1SU 38   // 0,1,..,37 in CLOUDY
#define MAX_NB_H2_VSTATES_C1PU 14   // 0,1,..,13 in CLOUDY
#define MAX_NB_H2_VSTATES_BP1SU 10  // 0,1,..,9 in CLOUDY
#define MAX_NB_H2_VSTATES_D1PU 19   // 0,1,..,18 in CLOUDY - for D-; for D+ v=0,1,2; (for C-+, D-+ levels with j=0 were commented in original file);
#define MAX_NB_H2_VSTATES_EF1SG 33   // 0,1,..32 in data by Scarlett, Glass-Maujean et al. (1984)

// are not used
#define MAX_NB_H2_VSTATES_BPP1SU 69  // 0,1,..68  M. Glass-Maujean, private communication (2025),
#define MAX_NB_H2_VSTATES_DP1PU 18   // 0,1,..15 for D'-; 0,1,..17 for D'+ 

#define MAX_H2_VSTATES_a3Sg 21   // 0,1,...,20  Wolniewicz L., Molecular Physics 105, 1497 (2007); but 0,1,...,21 Scarlett et al. (2023), private comm.
#define MAX_H2_VSTATES_c3Pu 21   // 0,1,...,20  Scarlett et al. (2023), private comm.
#define MAX_H2_VSTATES_d3Pu 21   // 0,1,...,20  Scarlett et al. (2023), private comm.




// The data on Lyman and Werner band transitions is initialized by this function - vector<transition>,
// H2 energy levels of the excited electronic state and of the ground state must be given,
void h2_band_transitions(const std::string& path, std::vector<transition> & einstein_coeff_vector, 
	const energy_diagram* h2_excited_di, const energy_diagram* h2_di, int verbosity = 1);

// name of the H2 electronic band must be given (to include it in the file name),
// the distribution of transition energies
void print_statistics(std::string output_path, std::string name, const std::vector<transition> &einstein_coeff_vector);

// auxiliary class
class dissociation_data_unit {
public:
	int v, j;
	double prob, kin_energy;  // probability in the data file in s-1,
	
	dissociation_data_unit();
	dissociation_data_unit(const dissociation_data_unit &);
};


// auxiliary class, is used in the constructor of h2_excited_state_data class,
struct decay_channel {
	int nb;       // number of the level
	double rate;  // in s-1,
};

// the class contains data on dissociation probabilities of a given energy level of excited state,
struct h2_energy_level_param {
public:
	// nb of excited level in the set of energy levels,
	int nb, nb_of_decays;
	// dissociation probability is dimensionless, total decay rate in s-1, total decay includes dissociation, 
	// kinetic energy released in the dissociation, in erg,
	double diss_prob, tot_decay, kin_energy;
	
	// arrays containing decay channel information: the number of level, decay probability (dimensionless),
	int * decay_level_nbs;
	double * decay_probs;

	// the relation operators are needed to sort transition vector by energy:
	bool operator == (h2_energy_level_param obj) const { return (nb == obj.nb); }
	bool operator != (h2_energy_level_param obj) const { return (nb != obj.nb); }
	bool operator < (h2_energy_level_param obj) const { return (nb < obj.nb); }
	bool operator > (h2_energy_level_param obj) const { return (nb > obj.nb); }
	
	// the assignment operator,
	h2_energy_level_param& operator=(const h2_energy_level_param&);

	h2_energy_level_param(int level_nb, int nb_of_decays, double diss_prob, double tot_decay, double kin_energy);
	h2_energy_level_param(const h2_energy_level_param &);
	~h2_energy_level_param();
};

// A -> X
void h2_excited_state_data(const std::string& path, std::vector<h2_energy_level_param> & level_data, 
	const energy_diagram* h2_excited_di, const std::vector<transition> &h2_band, int diss_data_exists, int verbosity);

// A -> B -> X, the level data (calculated by previous routine) already exist for B,
// "up" is for A state, "low" is for B state,
void h2_excited_state_data(const std::string& path, std::vector<h2_energy_level_param>& level_data_up, std::vector<h2_energy_level_param> level_data_low,
	const energy_diagram* h2_excited_di_up, const energy_diagram* h2_excited_di_low, 
	const std::vector<transition>& h2_band_up, const std::vector<transition>& h2_band_low, int verbosity);


// name of the H2 electronic band must be given (to include it in the file name),
// prints the lifetimes of the levels,
void print_statistics(std::string output_path, std::string name, const energy_diagram* h2_excited_di, 
	const std::vector<h2_energy_level_param>& level_param_vector);


//
// Honl-London factors, normalized,
// Whiting & Nicholls, ApJ Suppl. Ser. 235, 27 (1974); 
// additional: Gay et al. ApJ 746, 78 (2012);
class honl_london {
public:
	// the rotational quantum nb j of the lower level,
	virtual double operator()(int j, int dj) { return 0.; }
};

// S1g(X) -> S1u(B), S1u(Bp), dL = 0, dJ = +/-1 (absorption)
class honl_london_singlet_dl0 : public honl_london {
public:
	double operator()(int j, int dj);
	honl_london_singlet_dl0() { ; }
};

// S1g(X) -> P1u(C), P1u(D), dL = 1, dJ = -1, 0, 1
class honl_london_singlet_dl1 : public honl_london {
public:
	double operator()(int j, int dj);
	honl_london_singlet_dl1() { ; }
};
