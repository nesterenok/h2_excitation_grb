#pragma once
#include <cstring>
#include <vector>
#include "spectroscopy.h"
#include "dynamic_array.h"


// pure virtual class:
class cross_section
{
public:
	double en_thr, en_max;  // energy threshold in eV, 
	// the function returns the cross section in cm2, energy in eV
	virtual double operator() (double energy) const { return 0.; }
	virtual double operator() (int lev_nb, double energy, int index) const { return 0.; }

	// returns threshold energy in eV,
	double get_threshold_energy() const { return en_thr; }

	cross_section();
	virtual ~cross_section() { ; }
};

// be careful with processes having energy threshold,
class cross_section_table : public cross_section
{
protected:
	int nb_cs;
	double gamma;
	double* en_arr, * cs_arr;  // energy in eV

public:
	// energy in eV, the cross section in cm2 is returned,
	// if the energy is < min energy, the cross section at minimal energy is returned,
	// if the energy is > max energy, the cross section approximation is used: sigma = sigma_0*(E/E_0)^(gamma)
	virtual double operator() (double energy) const;

	// parameter gamma must be initialized in the constructor:
	cross_section_table();
	virtual ~cross_section_table();
};


// Momentum transfer cross section,
// cross section values in the data file are in 1.e-16 cm2,
// H2 - e-: Yoon et al., J. of Phys Chem Data 37, p.913, 2008 (E <= 100 eV); Dalgarno et al. ApJ Suppl. Ser. 125, p.237, 1999 (E > 100 eV);
//          see also Pinto & Galli, A&A 484, 17 (2008); 
// He - e-: Crompton et al., Aust. J. Phys. 23, p.667, 1970 (E <= 6 eV); Milloy & Crompton, Physical Review A 15, p.1847, 1977 (E <= 12 eV);
//          Dalgarno et al. ApJ Suppl. Ser. 125, p.237, 1999 (E > 12 eV);
// H - e-:  is not implemented yet
class cross_section_table_vs1 : public cross_section_table
{
public:
	// path to the data folder, and file name with the path within the data folder must be given:
	cross_section_table_vs1(const std::string& data_path, const std::string& name);
	~cross_section_table_vs1() { ; }
};

// cross section values in the data file are in a0^2 (Bohr radius squared), a0 = 0.529e-8 cm
class cross_section_table_mccc : public cross_section_table
{
public:
	// path to the data folder, and file name with the path within the data folder must be given:
	cross_section_table_mccc(const std::string& data_path, const std::string& name);
	~cross_section_table_mccc() { ; }
};

//
// Excitation of H2 vibration states,
// Janev et al. "Collision Processes in Low-Temperature Hydrogen Plasmas", FZ-Juelich Report No. 4105 (2003)
class h2_excitation_vibr02_cs : public cross_section
{
protected:
	double alpha, gamma;

public:
	double operator() (double energy) const;
	h2_excitation_vibr02_cs();
};


// Electron impact ionization (only simple atoms and molecules: H, He, H2, He+, H2+,..)
// Kim & Rudd, Physical Review A 50, p.3954 (1994); 
// Kim et al. Physical Review A 62, 052710 (2000);
class oscillator_strength
{
protected: 
	// approximation parameters of oscillator strength, dimensionless
	double af, bf, cf, df, ef, ff; 

public:
	virtual double operator() (double w) const = 0;
	void assign_coeff(double a, double b, double c, double d, double e, double f);

	oscillator_strength();
	virtual ~oscillator_strength();
};

// df/dw
// for calculation of singly differential cross section,
class diff_oscillator_strength : public oscillator_strength
{
public:	
	double operator() (double w) const;
	diff_oscillator_strength();
};

// 1/(1+w) * df/dw
// for total cross section calculations (is used in integration),
class diff_oscillator_strength_aux : public oscillator_strength
{
public:
	double operator() (double w) const;
	diff_oscillator_strength_aux();
};

//
// Because the scattered and ejected electrons are indistinguishable, it is customary to call the faster one of the
// two (after a collision) the primary electron, and the slower one - the secondary electron,
// Non-relativistic cross sections,
class electron_impact_ionization
{
protected:
	int		ne;  // number of electrons,
	double  orbital_bind_en;  // orbital binding energy, in eV 
	double  orbital_kin_en;   // orbital kinetic energy of the target electron, normalized on binding energy, dimensionless
	double  s;   // normalization parameter, [cm2]
	double  ni;  // oscillator strength integral, dimensionless
	double  proj_en;  // projectile energy, in eV 
	double  err;      // relative error in the integration
	double  c1, c2;   // auxiliary parameters,
	std::string name;

	diff_oscillator_strength		diff_osc_str;
	diff_oscillator_strength_aux	diff_osc_str_aux;

public:
	// energy in eV, returns differential cross section, per unit of energy of ejected electron [cm2 eV-1]
	// ejected_energy - energy of the slowest electron after collision,
	virtual double operator() (double projectile_energy, double ejected_energy) const;
	virtual double operator() (double ejected_energy) const;  // projectile energy must be assigned before function call,

	// returns cross section (integrated on ejected electron energy), [cm2]
	virtual double get_int_cs(double projectile_energy);
	virtual double get_int_cs(double projectile_energy, double ej_energy_min, double ej_energy_max);

	// returns energy in eV,
	double get_binding_energy() const { return orbital_bind_en; }
	std::string get_name() const { return name; }

	// energy in eV;
	void set_projectile_energy(double en) { proj_en = en; }

	electron_impact_ionization();
	virtual ~electron_impact_ionization();
};

// Relativistic cross section for electron impact ionization,



// Binary-Encounter-Dipole (BED) model for all species
// dissociative recombination cross section is about 10 per cent - is not taken into account here (Yoon et al., J.Phys. Chem. Ref. Data 37, 913, 2008)
// H2 + e- -> H2+ + e- + e-
class h2_electron_impact_ionization : public electron_impact_ionization {
public:
	h2_electron_impact_ionization(int verbosity);
};

// only single ionization of helium is considered,
// He + e- -> He+ + e- + e-  (total cross section is about 3e-17 cm2 at the peak about 100 eV),
// double ionization of helium is about 1e-19 cm2 at the peak at energy 250 eV, Ralchenko et al., Atomic Data and Nuclear Data Tables 94, 603 (2008)
class he_electron_impact_ionization : public electron_impact_ionization {
public:
	he_electron_impact_ionization(int verbosity);
};

// H + e- -> H+ + e- + e-
class h_electron_impact_ionization : public electron_impact_ionization {
public:
	h_electron_impact_ionization(int verbosity);
};

// He+ + e- -> He++ + e- + e-
// Binary-Encounter-Dipole (BED) model, 
// for one-electron ions the differential oscillator strength as for H atom (Kim & Rudd, 1994),
class hep_electron_impact_ionization : public electron_impact_ionization {
public:
	hep_electron_impact_ionization(int verbosity);
};

// H2+ + e- -> H+ + H+ + e- + e-,
// is accompanied by ion destruction
class h2p_electron_impact_ionization : public electron_impact_ionization {
public:
	h2p_electron_impact_ionization(int verbosity);
};


class cross_section_integral_2d
{
protected:
	double err, x1, x2, z1, z2;
	electron_impact_ionization* el_ion;

public:
	double operator() (double x);
	// x - initial projectile energy, 
	// z - ejected electron energy (slowest electron after collision),
	// returns integral of cross section, [cm2 eV]
	double get(double xx1, double xx2, double zz1, double zz2);

	cross_section_integral_2d(electron_impact_ionization* ei, double err);
};

class cross_section_integral_weight
{
protected:
	double err, x1, x2, y1, y2, z1, z2;
	electron_impact_ionization* el_ion;

public:
	double operator() (double x);
	// x - initial projectile energy, 
	// y - final projectile energy (fastest electron after collision), 
	// z - energy of the ejected electron, z <= y,
	// returns integral of diff. cross section on [x1, x2] and [z1, z2] but the scattered energy is lying in the given interval [y1, y2], [cm2 eV],
	double get(double xx1, double xx2, double yy1, double yy2, double zz1, double zz2);

	cross_section_integral_weight(electron_impact_ionization* ei, double err);
};


// HeI excitation by electron impact,
// Ralchenko et al., Atomic Data and Nuclear Data Tables 94, 603 (2008)
// initial (lower) -> final (upper) 
class hei_electron_excitation : public cross_section
{
protected:
	double a1, a2, a3, a4, a5, a6, f;
public:
	virtual double operator() (double en) const { return 0.; }
	// threshold energy in eV, statistical weight of the initial (lower) state must be given,
	hei_electron_excitation(double aa1, double aa2, double aa3, double aa4, double aa5, double aa6, double en_thr, int stat_weight);
};

// dS = 0, dL = +-1
class hei_electron_excitation_dipole_allowed : public hei_electron_excitation
{
public:
	double operator() (double en) const;
	hei_electron_excitation_dipole_allowed(double aa1, double aa2, double aa3, double aa4, double aa5, double aa6, double en_thr, int stat_weight);
};

// dS = 0, dL <> +-1
class hei_electron_excitation_dipole_forbidden : public hei_electron_excitation
{
public:
	double operator() (double en) const;
	hei_electron_excitation_dipole_forbidden(double aa1, double aa2, double aa3, double aa4, double aa5, double aa6, double en_thr, int stat_weight);
};

// dS <> 0
class hei_electron_excitation_spin_forbidden : public hei_electron_excitation
{
public:
	double operator() (double en) const;
	hei_electron_excitation_spin_forbidden(double aa1, double aa2, double aa3, double aa4, double aa5, double aa6, double en_thr, int stat_weight);
};