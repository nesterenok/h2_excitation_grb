
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <limits>

#include "utils.h"
#include "constants.h"
#include "h2_microphysics.h"
#include "integration.h"

#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "h2_microphysics.cpp"
using namespace std;


cross_section::cross_section() : en_thr(0.), en_max(1.e+99)
{;}

cross_section_table::cross_section_table() : nb_cs(0), gamma(0.), en_arr(0), cs_arr(0)
{;}

cross_section_table::~cross_section_table() {
	delete[] en_arr;
	delete[] cs_arr;
}

// be careful with processes having energy threshold,
double cross_section_table::operator()(double energy) const
{
	int l;
	double answ;
	// if values in question are out of range of the array with a dimension dim, returned index < 0 or index >= dim-1
	locate_index(en_arr, nb_cs, energy, l);

	if (l < 0) {
		answ = cs_arr[0];  // no energy threshold
	}
	else if (l >= nb_cs - 1) {
		answ = cs_arr[nb_cs - 1] * pow(energy / en_arr[nb_cs - 1], gamma);
	}
	else answ = cs_arr[l] + (energy - en_arr[l]) * (cs_arr[l + 1] - cs_arr[l]) / (en_arr[l + 1] - en_arr[l]);

	return answ;
}


// the class for momentum transfer cross sections,
cross_section_table_vs1::cross_section_table_vs1(const std::string& data_path, const std::string& name)
	:cross_section_table()
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i;
	double en, cs, reduced_mass;

	string file_name;
	ifstream input;

	file_name = data_path + name;
	input.open(file_name.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with cross section data " << file_name << endl;
		exit(1);
	}
	// comment lines are read (line length must be < MAX_TEXT_LINE_WIDTH):
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	// reduced mass is necessary for conversion from energy units to velocity:
	input >> gamma >> reduced_mass >> nb_cs;

	en_arr = new double[nb_cs];
	cs_arr = new double[nb_cs];

	for (i = 0; i < nb_cs; i++) {
		input >> en >> cs;
		// the energy data in the file have dimension eV,
		en_arr[i] = en;
		// the cross section data in the file have dimension 1.e-16 cm2,
		cs_arr[i] = cs * 1.e-16;
	}
	input.close();
	en_thr = en_arr[0];
}

cross_section_table_mccc::cross_section_table_mccc(const std::string& data_path, const std::string& name)
	:cross_section_table()
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i;
	double en, cs;
	vector<double> env, csv;  // auxiliary vectors,

	string file_name;
	stringstream ss;
	ifstream input;

	file_name = data_path + name;
	input.open(file_name.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with cross section data " << file_name << endl;
		exit(1);
	}
	env.clear();
	csv.clear();

	while (!input.eof())
	{
		// comment lines are read:
		do {
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		} while (text_line[0] == '#');

		if (text_line[0] == '\0') // check for empty line at the file end;
			break;

		ss.clear();
		ss.str(text_line);

		ss >> en >> cs;  // energy in eV, cross section data in a0^2, 
		cs *= BOHR_RADIUS * BOHR_RADIUS;  // conversion to cm2

		env.push_back(en);
		csv.push_back(cs);
	}
	input.close();

	nb_cs = (int)env.size();
	en_arr = new double[nb_cs];
	cs_arr = new double[nb_cs];

	for (i = 0; i < nb_cs; i++) {
		en_arr[i] = env[i];
		cs_arr[i] = csv[i];
	}

	en_thr = en_arr[0];
	gamma = log(cs_arr[nb_cs - 1] / cs_arr[nb_cs - 2]) / log(en_arr[nb_cs - 1] / en_arr[nb_cs - 2]);
}


h2_excitation_vibr02_cs::h2_excitation_vibr02_cs()
{
	en_thr = 1.003;  // in eV,
	alpha = 1.e-16 * 5.78 * 0.628 / (en_thr * en_thr * en_thr * en_thr);  // in cm-2
	gamma = 6.11 / sqrt(en_thr);
}

double h2_excitation_vibr02_cs::operator()(double energy) const
{
	if (energy < en_thr)
		return 0.;

	double x = en_thr / energy;
	return alpha * x * x * pow(1. - x, gamma);
}


//
// Electron impact ionization classes
oscillator_strength::oscillator_strength() : af(0.), bf(0.), cf(0.), df(0.), ef(0.), ff(0.)
{;}

oscillator_strength::~oscillator_strength()
{;}

void oscillator_strength::assign_coeff(double a, double b, double c, double d, double e, double f)
{
	af = a;
	bf = b;
	cf = c;
	df = d;
	ef = e;
	ff = f;
}

diff_oscillator_strength::diff_oscillator_strength()
	: oscillator_strength()
{;}

double diff_oscillator_strength::operator()(double w) const {
	double x, x2;
	x = 1. / (1. + w);
	x2 = x * x;
	return (af + bf * x + cf * x2 + df * x2 * x + ef * x2 * x2 + ff * x2 * x2 * x) * x;
}

diff_oscillator_strength_aux::diff_oscillator_strength_aux()
	: oscillator_strength()
{;}

double diff_oscillator_strength_aux::operator()(double w) const {
	double x, x2;
	x = 1. / (1. + w);
	x2 = x * x;
	return (af + bf * x + cf * x2 + df * x2 * x + ef * x2 * x2 + ff * x2 * x2 * x) * x2;
}


electron_impact_ionization::electron_impact_ionization()
	: ne(0), s(0.), ni(0.), orbital_kin_en(0.), orbital_bind_en(0.), proj_en(0.), name(""), err(1.e-6), c1(0.), c2(0.)
{;}

electron_impact_ionization::~electron_impact_ionization()
{;}

// it is assumed in the formula, that ejected energy is the energy of the slowest electron,
double electron_impact_ionization::operator()(double projectile_energy, double ejected_energy) const
{
	if (ejected_energy < 0. || projectile_energy < 0.)
		return 0.;

	double cs, w, wp, t;

	t = projectile_energy / orbital_bind_en;
	if (t < 1.)
		return 0.;

	w = ejected_energy / orbital_bind_en;
	wp = t - w - 1.;
	if (wp < 0.)
		return 0.;

	// formula (44) in Kim & Rudd (1994),
	cs = c1 / (t + orbital_kin_en + 1.)
		* (1. / ((1. + w) * (1. + w)) - 1. / ((1. + w) * (1. + wp)) + 1. / ((1. + wp) * (1. + wp))
			+ log(t) * diff_osc_str(w) / (c2 * (1. + w)));
	return cs;
}

double electron_impact_ionization::operator()(double ejected_energy) const
{
	double cs;
	// saved projectile energy is in eV,
	cs = electron_impact_ionization::operator()(proj_en, ejected_energy);  
	return cs;
}

double electron_impact_ionization::get_int_cs(double projectile_energy)
{
	double cs, d, lnt, t;
	t = projectile_energy / orbital_bind_en;
	if (t < 1.)
		return 0.;

	lnt = log(t);
	d = qromb<diff_oscillator_strength_aux>(diff_osc_str_aux, 0., 0.5 * (t - 1.), err, false);
	
	// equations (55) and (56) in Kim & Rudd (1994),
	cs = s / (ne * (1. + t + orbital_kin_en)) * (d * lnt + c2 * (1. - 1. / t - lnt / (1. + t)));
	return cs;
}

double electron_impact_ionization::get_int_cs(double projectile_energy, double ej_energy_min, double ej_energy_max)
{
	double cs, d, lnt, z, w, t, e1, e2;
	t = projectile_energy / orbital_bind_en;
	if (t < 1.)
		return 0.;

	e1 = ej_energy_min / orbital_bind_en;
	e2 = ej_energy_max / orbital_bind_en;

	lnt = log(t);
	d = qromb<diff_oscillator_strength_aux>(diff_osc_str_aux, e1, e2, err, false);
	
	z = log((1. + e2) * (t - e1) / ((1. + e1) * (t - e2)));
	w = (e2 - e1) * (1. / ((1. + e2) * (1. + e1)) + 1. / ((t - e2) * (t - e1)));
	
	cs = s / (ne * (1. + t + orbital_kin_en)) * (d * lnt + c2 * (-z / (1. + t) + w));
	return cs;
}


h2_electron_impact_ionization::h2_electron_impact_ionization(int verbosity)
	: electron_impact_ionization()
{
	name = "H2 + e- -> H2+ + e- + e-";
	orbital_bind_en = 15.43;  // in eV,
	orbital_kin_en = 25.68 / orbital_bind_en;   // must be normalized,

	ne = 2;
	s = 4. * M_PI * BOHR_RADIUS * BOHR_RADIUS * ne
		* RYDBERG_ENERGY_EV * RYDBERG_ENERGY_EV / (orbital_bind_en * orbital_bind_en);

	// Kim & Rudd (1994)
	diff_osc_str.assign_coeff(0., 0., 1.1262, 6.3982, -7.8055, 2.144);
	diff_osc_str_aux.assign_coeff(0., 0., 1.1262, 6.3982, -7.8055, 2.144);

	// the values of N_i are provided by Kim & Rudd (1994) for simple atoms/molecules (H, He, H2, Ne),
	// this data may be used to check the integration:
	ni = qromb<diff_oscillator_strength>(diff_osc_str, 0., 100., err, false);  // upper integration limit is arbitrary

	if (verbosity) {
		cout << scientific;
		cout.precision(3);
		cout << left << "H2 electron impact ionization, " << endl
			<< "    parameter Ni from integration: " << ni << endl
			<< "    from table in Kim & Rudd (1994): 1.173" << endl;
	}
	c1 = s * (2. * ne - ni) / (orbital_bind_en * ne);
	c2 = 2. * ne - ni;
}

he_electron_impact_ionization::he_electron_impact_ionization(int verbosity)
	: electron_impact_ionization()
{
	name = "He + e- -> He+ + e- + e-";
	orbital_bind_en = 24.59;  // in eV,
	orbital_kin_en = 39.51 / orbital_bind_en;   // must be normalized,

	ne = 2;
	s = 4. * M_PI * BOHR_RADIUS * BOHR_RADIUS * ne
		* RYDBERG_ENERGY_EV * RYDBERG_ENERGY_EV / (orbital_bind_en * orbital_bind_en);

	// Kim & Rudd (1994)
	diff_osc_str.assign_coeff(0., 0., 12.178, -29.585, 31.251, -12.175);
	diff_osc_str_aux.assign_coeff(0., 0., 12.178, -29.585, 31.251, -12.175);

	ni = qromb<diff_oscillator_strength>(diff_osc_str, 0., 100., err, false);

	if (verbosity) {
		cout << scientific;
		cout.precision(3);
		cout << left << "He electron impact ionization, " << endl
			<< "    parameter Ni from integration: " << ni << endl
			<< "    from table in Kim & Rudd (1994): 1.605" << endl;
	}
	c1 = s * (2. * ne - ni) / (orbital_bind_en * ne);
	c2 = 2. * ne - ni;
}

h_electron_impact_ionization::h_electron_impact_ionization(int verbosity)
{
	name = "H + e- -> H+ + e- + e-";
	orbital_bind_en = 13.606;  // in eV,
	orbital_kin_en = 13.606;

	ne = 1;
	s = 4. * M_PI * BOHR_RADIUS * BOHR_RADIUS * ne
		* RYDBERG_ENERGY_EV * RYDBERG_ENERGY_EV / (orbital_bind_en * orbital_bind_en);

	// Kim & Rudd (1994),
	diff_osc_str.assign_coeff(0., -0.022473, 1.1775, -0.46264, 0.089064, 0.);
	diff_osc_str_aux.assign_coeff(0., -0.022473, 1.1775, -0.46264, 0.089064, 0.);

	ni = qromb<diff_oscillator_strength>(diff_osc_str, 0., 100., err, false);

	if (verbosity) {
		cout << scientific;
		cout.precision(3);
		cout << left << "H electron impact ionization, " << endl
			<< "    parameter Ni from integration: " << ni << endl
			<< "    from table in Kim & Rudd (1994): 0.4343" << endl;
	}
	c1 = s * (2. * ne - ni) / (orbital_bind_en * ne);
	c2 = 2. * ne - ni;
}

hep_electron_impact_ionization::hep_electron_impact_ionization(int verbosity)
{
	name = "He+ + e- -> He++ + e- + e-";
	orbital_bind_en = 54.416;  // in eV, Draine, "Physics of the interstellar and intergalactic medium" (2011)
	orbital_kin_en = 0.;     // for one-electron ions is set to zero,

	ne = 1;
	s = 4. * M_PI * BOHR_RADIUS * BOHR_RADIUS * ne
		* RYDBERG_ENERGY_EV * RYDBERG_ENERGY_EV / (orbital_bind_en * orbital_bind_en);

	// Kim & Rudd (1994), as for H atom
	diff_osc_str.assign_coeff(0., -0.022473, 1.1775, -0.46264, 0.089064, 0.);
	diff_osc_str_aux.assign_coeff(0., -0.022473, 1.1775, -0.46264, 0.089064, 0.);

	ni = qromb<diff_oscillator_strength>(diff_osc_str, 0., 100., err, false);

	if (verbosity) {
		cout << scientific;
		cout.precision(3);
		cout << left << "He+ electron impact ionization, " << endl
			<< "    parameter Ni from integration: " << ni << endl;
	}
	c1 = s * (2. * ne - ni) / (orbital_bind_en * ne);
	c2 = 2. * ne - ni;
}

h2p_electron_impact_ionization::h2p_electron_impact_ionization(int verbosity)
{
	name = "H2+ + e- -> H+ + H+ + e- + e-";
	orbital_bind_en = 16.25;  // in eV, Wikipedia
	orbital_kin_en = 0.;     // for one-electron ions is set to zero,

	ne = 1;
	s = 4. * M_PI * BOHR_RADIUS * BOHR_RADIUS * ne
		* RYDBERG_ENERGY_EV * RYDBERG_ENERGY_EV / (orbital_bind_en * orbital_bind_en);

	// Kim & Rudd (1994), as for H atom
	diff_osc_str.assign_coeff(0., -0.022473, 1.1775, -0.46264, 0.089064, 0.);
	diff_osc_str_aux.assign_coeff(0., -0.022473, 1.1775, -0.46264, 0.089064, 0.);

	ni = qromb<diff_oscillator_strength>(diff_osc_str, 0., 100., err, false);

	if (verbosity) {
		cout << scientific;
		cout.precision(3);
		cout << left << "H2+ electron impact ionization, " << endl
			<< "    parameter Ni from integration: " << ni << endl;
	}
	c1 = s * (2. * ne - ni) / (orbital_bind_en * ne);
	c2 = 2. * ne - ni;
}


cross_section_integral_2d::cross_section_integral_2d(electron_impact_ionization* ei, double e)
	:x1(0.), x2(0.), z1(0.), z2(0.), el_ion(ei), err(e)
{
	;
}

double cross_section_integral_2d::operator()(double x)
{
	if (x < el_ion->get_binding_energy())
		return 0.;

	double a;
	el_ion->set_projectile_energy(x);

	// ejected electron has the lowest energy,
	a = 0.5 * (x - el_ion->get_binding_energy());
	if (z2 < a)
		a = z2;

	if (a < z1)
		return 0.;

	a = qromb<electron_impact_ionization>(*el_ion, z1, a, 0.1 * err, false);
	return a;
}

double cross_section_integral_2d::get(double xx1, double xx2, double zz1, double zz2)
{
	x1 = xx1;
	x2 = xx2;
	z1 = zz1;
	z2 = zz2;

	return qromb<cross_section_integral_2d>(*this, x1, x2, err, false);
}

cross_section_integral_weight::cross_section_integral_weight(electron_impact_ionization* ei, double e)
	:x1(0.), x2(0.), y1(0.), y2(0.), z1(0.), z2(0.), el_ion(ei), err(e)
{
	;
}

double cross_section_integral_weight::operator()(double x)
{
	double a, b, c;
	el_ion->set_projectile_energy(x);

	a = x - y2 - el_ion->get_binding_energy();
	if (z1 > a)
		a = z1;

	b = x - y1 - el_ion->get_binding_energy();
	if (z2 < b)
		b = z2;

	// ejected electron has the lowest energy,
	c = 0.5 * (x - el_ion->get_binding_energy());
	if (c < b)
		b = c;

	if (a >= b)
		return 0.;

	a = qromb<electron_impact_ionization>(*el_ion, a, b, 0.1 * err, false);
	return a;
}

double cross_section_integral_weight::get(double xx1, double xx2, double yy1, double yy2, double zz1, double zz2)
{
	x1 = xx1;
	x2 = xx2;
	y1 = yy1;
	y2 = yy2;
	z1 = zz1;
	z2 = zz2;

	return qromb<cross_section_integral_weight>(*this, x1, x2, err, false);
}


// HeI excitation by electron impact,
hei_electron_excitation::hei_electron_excitation(double aa1, double aa2, double aa3, double aa4, double aa5, double aa6,
	double et, int stat_weight) : a1(aa1), a2(aa2), a3(aa3), a4(aa4), a5(aa5), a6(aa6) {
	en_thr = et;
	f = M_PI * BOHR_RADIUS * BOHR_RADIUS * RYDBERG_ENERGY_EV / stat_weight;
}

hei_electron_excitation_dipole_allowed::hei_electron_excitation_dipole_allowed(double aa1, double aa2, double aa3, double aa4, double aa5, double aa6,
	double en_thr, int stat_weight) : hei_electron_excitation(aa1, aa2, aa3, aa4, aa5, aa6, en_thr, stat_weight)
{
	;
}

double hei_electron_excitation_dipole_allowed::operator()(double en) const
{
	if (en < en_thr)
		return 0.;

	double answ, x;

	x = en / en_thr;
	answ = (a1 * log(x) + a2 + (a3 + a4 / x + a5 / (x * x)) / x) * (1. + x) / (a6 + x) * (f / en);

	if (answ < 0.)
		answ = 0.;

	return answ;
}

hei_electron_excitation_dipole_forbidden::hei_electron_excitation_dipole_forbidden(double aa1, double aa2, double aa3, double aa4, double aa5, double aa6,
	double en_thr, int stat_weight) : hei_electron_excitation(aa1, aa2, aa3, aa4, aa5, aa6, en_thr, stat_weight)
{
	;
}

double hei_electron_excitation_dipole_forbidden::operator()(double en) const
{
	if (en < en_thr)
		return 0.;

	double answ, x, x2;

	x = en / en_thr;
	x2 = x * x;
	answ = (a1 + (a2 + a3 / x + a4 / (x2)) / x) * x2 / (a5 + x2) * (f / en);

	if (answ < 0.)
		answ = 0.;

	return answ;
}

hei_electron_excitation_spin_forbidden::hei_electron_excitation_spin_forbidden(double aa1, double aa2, double aa3, double aa4, double aa5, double aa6,
	double en_thr, int stat_weight) : hei_electron_excitation(aa1, aa2, aa3, aa4, aa5, aa6, en_thr, stat_weight)
{
	;
}

double hei_electron_excitation_spin_forbidden::operator()(double en) const
{
	if (en < en_thr)
		return 0.;

	double answ, x;

	x = en / en_thr;
	answ = (a1 + (a2 + a3 / x + a4 / (x * x)) / x) / (a5 + x * x) * (f / en);

	if (answ < 0.)
		answ = 0.;

	return answ;
}
