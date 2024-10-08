#pragma once
#include <cstring>
#include <vector>
#include "spectroscopy.h"

// The data on Lyman and Werner band transitions is initialized by this function - vector<transition>,
// H2 energy levels of the excited electronic state and of the ground state must be given,
void h2_band_transitions(const std::string& path, std::vector<transition> & einstein_coeff_vector, 
	const energy_diagram* h2_excited_di, const energy_diagram* h2_di, int verbosity = 1);

// name of the H2 electronic band must be given (to include it in the file name),
// the distribution of transition energies
void print_statistics(std::string output_path, std::string name, const std::vector<transition> &einstein_coeff_vector);


// auxiliary class, is used in the constructor of h2_excited_state_data class,
struct decay_channel {
	int nb;
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

// 
void h2_excited_state_data(const std::string& path, std::vector<h2_energy_level_param> & level_data, 
	const energy_diagram*, const std::vector<transition> &h2_band, int verbosity = 1);

// name of the H2 electronic band must be given (to include it in the file name),
// prints the lifetimes of the levels,
void print_statistics(std::string output_path, std::string name, const energy_diagram* h2_excited_di, 
	const std::vector<h2_energy_level_param>& level_param_vector);


class honl_london {
public:
	// the rotational quantum nb j of the lower level,
	virtual double operator()(int j, int dj) { return 0.; }
};

// for S1g(X) -> S1u(B), adsorption, j -> j + dj
class honl_london_singlet_dl0 : public honl_london {
public:
	double operator()(int j, int dj);
	honl_london_singlet_dl0() { ; }
};

class honl_london_singlet_dl1 : public honl_london {
public:
	double operator()(int j, int dj);
	honl_london_singlet_dl1() { ; }
};
