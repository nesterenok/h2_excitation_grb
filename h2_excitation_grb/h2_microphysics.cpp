
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
#include "h2_uv_pump_data.h"
#include "integration.h"

#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "h2_microphysics.cpp"
using namespace std;


double calc_coulomb_losses_thermal_electrons(double energy, double ne, double te)
{
	if (energy < KELVIN_TO_EV * te)  // energy in eV, temperature in K
		return 0.;

	double losses = 2.e-4 * pow(ne, 0.97)
		* pow((energy - KELVIN_TO_EV * te) / (energy - 0.53 * KELVIN_TO_EV * te), 2.36) / pow(energy, 0.44);  // [eV s-1]
	
	return losses;
}

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


electron_ionization_data::electron_ionization_data()
	: ne(0), ni(0.), orbital_kin_en(0.), orbital_bind_en(0.), name("empty"), af(0.), bf(0.), cf(0.), df(0.), ef(0.), ff(0.)
{;}

h2_electron_ionization_data::h2_electron_ionization_data()
	: electron_ionization_data()
{
	name = "H2 + e- -> H2+ + e- + e-";
	orbital_bind_en = 15.43;  // in eV,
	orbital_kin_en = 25.68;   // in eV,

	// Kim & Rudd (1994)
	ne = 2;
	ni = 1.173;
	af = 0.;
	bf = 0.;
	cf = 1.1262;
	df = 6.3982;
	ef = -7.8055; 
	ff = 2.144;
}

he_electron_ionization_data::he_electron_ionization_data()
	: electron_ionization_data()
{
	name = "He + e- -> He+ + e- + e-";
	orbital_bind_en = 24.59;  // in eV,
	orbital_kin_en = 39.51;   // in eV,

	// Kim & Rudd (1994)
	ne = 2;
	ni = 1.605;
	af = 0.;
	bf = 0.;
	cf = 12.178;
	df = -29.585;
	ef = 31.251;
	ff = -12.175;
}

h_electron_ionization_data::h_electron_ionization_data()
	: electron_ionization_data()
{
	name = "H + e- -> H+ + e- + e-";
	orbital_bind_en = 13.606;  // in eV,
	orbital_kin_en = 13.606;

	// Kim & Rudd (1994),
	ne = 1;
	ni = 0.4343;
	af = 0.;
	bf = -0.022473;
	cf = 1.1775;
	df = -0.46264;
	ef = 0.089064; 
	ff = 0.;
}

hep_electron_ionization_data::hep_electron_ionization_data()
	: electron_ionization_data()

{
	name = "He+ + e- -> He++ + e- + e-";
	orbital_bind_en = 54.416;  // in eV, Draine, "Physics of the interstellar and intergalactic medium" (2011)
	orbital_kin_en = 0.;       // for one-electron ions is set to zero,

	// Kim & Rudd (1994), oscillation strength coefficients as for H atom (for one-electron ions)
	ne = 1;
	ni = -1.;
	af = 0.;
	bf = -0.022473;
	cf = 1.1775;
	df = -0.46264;
	ef = 0.089064; 
	ff = 0.;
}

// may be, BEB model must be used
h2p_electron_ionization_data::h2p_electron_ionization_data()
{
	name = "H2+ + e- -> H+ + H+ + e- + e-";
	orbital_bind_en = 16.25;  // in eV, Wikipedia
	orbital_kin_en = 0.;      // for one-electron ions is set to zero,

	// Kim & Rudd (1994), oscillation strength coefficients as for H atom (for one-electron ions)
	// BED model is used, possibly more simple BEQ or BEB model should be used
	ne = 1;
	ni = -1.;
	af = 0.;
	bf = -0.022473;
	cf = 1.1775;
	df = -0.46264;
	ef = 0.089064;
	ff = 0.;
}

electron_impact_ionization::electron_impact_ionization(const electron_ionization_data & data, int verbosity)
{
	double af, bf, cf, df, ef, ff;

	af = data.af;
	bf = data.bf;
	cf = data.cf;
	df = data.df;
	ef = data.ef;
	ff = data.ff;

	ne = data.ne;
	orbital_bind_en = data.orbital_bind_en;
	orbital_kin_en = data.orbital_kin_en;
	name = data.name;

	u = orbital_kin_en / orbital_bind_en;
	s = 4. * M_PI * BOHR_RADIUS * BOHR_RADIUS * ne
		* RYDBERG_ENERGY_EV * RYDBERG_ENERGY_EV / (orbital_bind_en * orbital_bind_en); 

	diff_osc_str.assign_coeff(af, bf, cf, df, ef, ff);
	diff_osc_str_aux.assign_coeff(af, bf, cf, df, ef, ff);

	// the values of N_i are provided by Kim & Rudd (1994) for simple atoms/molecules (H, He, H2, Ne),
	// this data may be used to check the integration:
	ni = qromb<diff_oscillator_strength>(diff_osc_str, 0., 100., err, false);  // upper integration limit is arbitrary

	if (verbosity) {
		cout << scientific;
		cout.precision(3);
		cout << left << "reaction: " << name << endl
			<< "    parameter Ni from integration: " << ni << endl
			<< "    from table in Kim & Rudd (1994): " << data.ni << endl;
	}

	c1 = 2. - ni / ne;	
	err = 1.e-6;
	proj_en = 0.;
}

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
	cs = s / (orbital_bind_en * (t + u + 1.))
		* ( c1 * (1. / ((1. + w) * (1. + w)) - 1. / ((1. + w) * (1. + wp)) + 1. / ((1. + wp) * (1. + wp)))
			+ log(t) * diff_osc_str(w) / (ne * (1. + w)) );
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
	d /= ne;

	// equations (55) and (56) in Kim & Rudd (1994),
	cs = s / (1. + t + u) * (d * lnt + c1 * (1. - 1. / t - lnt / (1. + t)));
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
	d /= ne;

	z = log((1. + e2) * (t - e1) / ((1. + e1) * (t - e2)));
	w = (e2 - e1) * (1. / ((1. + e2) * (1. + e1)) + 1. / ((t - e2) * (t - e1)));
	
	cs = s / (1. + t + u) * (d * lnt + c1 * (-z / (1. + t) + w));
	return cs;
}

double electron_impact_ionization::get_energy_loss(double projectile_energy)
{
	double loss, t, d1, d2, lnt, lnt2;

	t = projectile_energy / orbital_bind_en;
	if (t < 1.)
		return 0.;

	d1 = qromb<diff_oscillator_strength>(diff_osc_str, 0., 0.5 * (t - 1.), err, false);
	d2 = qromb<diff_oscillator_strength_aux>(diff_osc_str_aux, 0., 0.5 * (t - 1.), err, false);

	lnt = log(t);
	lnt2 = log(0.5 * (t + 1.));

	// two parts of energy losses - ejected electron energy and binding energy,
	loss = orbital_bind_en * s / (1. + t + u) * ( lnt * (d1 - d2) / ne + c1 * (3. * lnt2 - lnt * (2. - 1. / (1. + t))) );
	loss += orbital_bind_en * get_int_cs(projectile_energy); 
	
	return loss;
}


// Relativistic formulas for cross sections
electron_impact_ionization_relativistic::electron_impact_ionization_relativistic(const electron_ionization_data & data, int verbosity)
	: electron_impact_ionization(data, verbosity), srel(0.), bp(0.), up(0.), beta_b2(0.), beta_u2(0.)
{
	// constants that are used in relativistic formula of cross section,
	srel = 4.*M_PI * BOHR_RADIUS * BOHR_RADIUS 
		* FINE_STRUCTURE_CONSTANT * FINE_STRUCTURE_CONSTANT * FINE_STRUCTURE_CONSTANT * FINE_STRUCTURE_CONSTANT * 0.5 * ne;

	bp = orbital_bind_en / ELECTRON_MASS_EV;
	up = orbital_kin_en / ELECTRON_MASS_EV;

	beta_b2 = 1. - 1. / ((1. + bp) * (1. + bp));
	beta_u2 = 1. - 1. / ((1. + up) * (1. + up));
}

double electron_impact_ionization_relativistic::operator()(double projectile_energy, double ejected_energy) const
{
	if (ejected_energy < 0. || projectile_energy < 0.)
		return 0.;

	double t, tp, beta_t2;
	double cs, w, wp;

	t = projectile_energy / orbital_bind_en;
	if (t < 1.)
		return 0.;

	w = ejected_energy / orbital_bind_en;
	wp = t - w - 1.;
	if (wp < 0.)
		return 0.;

	tp = projectile_energy / ELECTRON_MASS_EV;
	beta_t2 = 1. - 1. / ((1. + tp) * (1. + tp));
	
	// formula (19) in Kim et al. (2000),
	// must be deleted by binding energy,
	cs = srel / ((beta_t2 + beta_u2 + beta_b2) * bp * orbital_bind_en) * (
		-c1 / ((1. + w) * (1. + wp)) * (1. + 2. * tp) / ((1. + 0.5 * tp) * (1. + 0.5 * tp))
		+ c1 * (1. / ((1. + w) * (1. + w)) + 1. / ((1. + wp) * (1. + wp)) + bp * bp / ((1. + 0.5 * tp) * (1. + 0.5 * tp)))
		+ diff_osc_str(w) / (ne * (1. + w)) * (log(beta_t2 / (1. - beta_t2)) - beta_t2 - log(2. * bp)));
	return cs;
}

double electron_impact_ionization_relativistic::operator()(double ejected_energy) const
{
	double cs;
	// the projectile energy is in eV, must be initialized before,
	cs = electron_impact_ionization_relativistic::operator()(proj_en, ejected_energy);
	return cs;
}

double electron_impact_ionization_relativistic::get_int_cs(double projectile_energy)
{
	double cs, x, d, t, tp, beta_t2;

	t = projectile_energy / orbital_bind_en;
	if (t < 1.)
		return 0.;

	tp = projectile_energy / ELECTRON_MASS_EV;
	beta_t2 = 1. - 1. / ((1. + tp) * (1. + tp));

	d = qromb<diff_oscillator_strength_aux>(diff_osc_str_aux, 0., 0.5 * (t - 1.), err, false);
	d /= ne;
	x = (1. + 0.5 * tp) * (1. + 0.5 * tp);

	cs = srel / ((beta_t2 + beta_u2 + beta_b2) * bp) * (
		d *(log(beta_t2 /(1. - beta_t2)) - beta_t2 - log(2.*bp))
		+ c1 *(-log(t) * (1. + 2. * tp) /((1. + t) * x) + 1. - 1./t  + 0.5 * bp *bp * (t - 1.) / x) );
	return cs;
}

double electron_impact_ionization_relativistic::get_int_cs(double projectile_energy, double ej_energy_min, double ej_energy_max)
{
	double cs, d, lnt, x, z, w, t, tp, beta_t2, e1, e2;
	t = projectile_energy / orbital_bind_en;
	if (t < 1.)
		return 0.;

	tp = projectile_energy / ELECTRON_MASS_EV;
	beta_t2 = 1. - 1. / ((1. + tp) * (1. + tp));

	e1 = ej_energy_min / orbital_bind_en;
	e2 = ej_energy_max / orbital_bind_en;

	lnt = log(t);
	d = qromb<diff_oscillator_strength_aux>(diff_osc_str_aux, e1, e2, err, false);
	d /= ne;

	x = (1. + 0.5 * tp) * (1. + 0.5 * tp);
	z = log((1. + e2) * (t - e1) / ((1. + e1) * (t - e2)));
	w = (e2 - e1) * ( 1. / ((1. + e2) * (1. + e1)) + 1. / ((t - e2) * (t - e1)) );

	cs = srel / ((beta_t2 + beta_u2 + beta_b2) * bp) * (
		d * (log(beta_t2 / (1. - beta_t2)) - beta_t2 - log(2. * bp)) 
		+ c1 * (-z * (1. + 2. * tp) / ((1. + t) * x)
		+ w + bp * bp *(e2 - e1) / x) );
	return cs;
}

double electron_impact_ionization_relativistic::get_energy_loss(double projectile_energy)
{
	double loss, x, d1, d2, t, tp, beta_t2, lnt, lnt2;

	t = projectile_energy / orbital_bind_en;
	if (t < 1.)
		return 0.;

	tp = projectile_energy / ELECTRON_MASS_EV;
	beta_t2 = 1. - 1. / ((1. + tp) * (1. + tp));

	d1 = qromb<diff_oscillator_strength>(diff_osc_str, 0., 0.5 * (t - 1.), err, false);
	d2 = qromb<diff_oscillator_strength_aux>(diff_osc_str_aux, 0., 0.5 * (t - 1.), err, false);
	
	lnt = log(t);
	lnt2 = log(0.5*(1. + t));
	x = (1. + 0.5 * tp) * (1. + 0.5 * tp);

	loss = srel * orbital_bind_en / ((beta_t2 + beta_u2 + beta_b2) * bp) * (
		(d1 - d2) / ne * (log(beta_t2 / (1. - beta_t2)) - beta_t2 - log(2. * bp))
		+ c1 * ((1. + t) * lnt2 - t * lnt) * (1. + 2. * tp) / ((1. + t) * x) 
		+ c1 * (2. * lnt2 - lnt)
		+ c1 * 0.5 * bp * bp * (t - 1.) / x );

	loss += orbital_bind_en * get_int_cs(projectile_energy);
	return loss;
}

electron_impact_ionization_relativistic::~electron_impact_ionization_relativistic()
{;}

//
// is not used
cross_section_integral_2d::cross_section_integral_2d(electron_impact_ionization* ei, double e)
	:x1(0.), x2(0.), z1(0.), z2(0.), ionization(ei), err(e)
{;}

double cross_section_integral_2d::operator()(double x)
{
	if (x < ionization->get_binding_energy())
		return 0.;

	double a;
	ionization->set_projectile_energy(x);

	// maximal energy of the ejected electron (ejected electron has the lowest energy among final electrons):
	a = 0.5 * (x - ionization->get_binding_energy());
	if (z2 < a)
		a = z2;

	if (a < z1)
		return 0.;

	a = qromb<electron_impact_ionization>(*ionization, z1, a, 0.1 * err, false);
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
	:x1(0.), x2(0.), y1(0.), y2(0.), z1(0.), z2(0.), ionization(ei), err(e)
{;}

double cross_section_integral_weight::operator()(double x)
{
	double a, b, max;
	ionization->set_projectile_energy(x);

	// lower limit of the integration:
	a = x - y2 - ionization->get_binding_energy();
	if (z1 > a)
		a = z1;

	// upper limit of the integration:
	b = x - y1 - ionization->get_binding_energy();
	if (z2 < b)
		b = z2;

	// ejected electron has the lowest energy,
	max = 0.5 * (x - ionization->get_binding_energy());
	if (max < b)
		b = max;

	if (a >= b)
		return 0.;

	// integral over the energy interval of ejected electron, 
	a = qromb<electron_impact_ionization>(*ionization, a, b, 0.1 * err, false);
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

	// integral over the energy interval of incident particle:
	return qromb<cross_section_integral_weight>(*this, x1, x2, err, false);
}


//
// Vibrationally-resolved electronic excitation of H2 ground state S1g+(X), 
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

	// By default, the abrupt cutoff of energy cross sections is used at high energies,
	// the cross section at this energy may not be zero - the reaction without threshold,
	en_thr = en_arr[0];
}

// be careful with processes having energy threshold,
double cross_section_table_mccc::operator()(double energy) const
{
	int l;
	double answ;
	// if values in question are out of range of the array with a dimension dim, returned index < 0 or index >= dim-1
	locate_index(en_arr, nb_cs, energy, l);

	if (l < 0) {
		answ = cs_arr[0];  // no energy threshold
	}
	else if (l >= nb_cs - 1) {
		answ = 0.;
	}
	else answ = cs_arr[l] + (energy - en_arr[l]) * (cs_arr[l + 1] - cs_arr[l]) / (en_arr[l + 1] - en_arr[l]);

	return answ;
}

//
// HeI excitation by electron impact,
//
hei_electron_excitation::hei_electron_excitation(double aa1, double aa2, double aa3, double aa4, double aa5, double aa6,
	double et, int stat_weight) : a1(aa1), a2(aa2), a3(aa3), a4(aa4), a5(aa5), a6(aa6) {
	en_thr = et;
	f = M_PI * BOHR_RADIUS * BOHR_RADIUS * RYDBERG_ENERGY_EV / stat_weight;
}

hei_electron_excitation_dipole_allowed::hei_electron_excitation_dipole_allowed(double aa1, double aa2, double aa3, double aa4, double aa5, double aa6,
	double en_thr, int stat_weight) : hei_electron_excitation(aa1, aa2, aa3, aa4, aa5, aa6, en_thr, stat_weight)
{;}

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
{;}

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
{;}

double hei_electron_excitation_spin_forbidden::operator()(double en) const
{
	if (en < en_thr)
		return 0.;

	double answ, x, x2;

	x = en / en_thr;
	x2 = x * x;
	answ = (a1 + (a2 + a3 / x + a4 / x2) / x) / (a5 + x2) * (f / en);

	if (answ < 0.)
		answ = 0.;

	return answ;
}

// 
void save_cross_section_table(const string& output_path, const string& data_path)
{
	const double ionization_fraction = 1.e-6;
	const double thermal_electron_temp = 1.e-3;  // eV

	const double energy_min = 0.01;  // eV
	const double energy_max = 1.e+7;
	const int nb_of_bins_per_order = 30;

	int i, nb_of_energies, vi, vf, ji, verbosity = 1;
	double t, enloss_mt_h2, enloss_ion_h2, enloss_ion_h2_rel, enloss_rot_h2, enloss_vibr_h2_01, enloss_vibr_h2_02,
		enloss_singlet_h2, enloss_diss_singlet_h2, enloss_diss_triplet_h2, enloss_mt_he, enloss_ion_he, enloss_ion_he_rel, enloss_coulomb,
		cs_vibr_h2_01, cs_vibr_h2_02, cs_vibr_h2_03, cs_singlet_h2, cs_diss_singlet_h2, cs_diss_triplet_h2;
	double *electron_energies, *electron_velocities;

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

	// H2 cross sections,
	// Initialization of cross section data for elastic scattering:
	cross_section_table_vs1 *elastic_h2_el_cs
		= new cross_section_table_vs1(data_path, "elastic_scattering/e-h2_mt_cs.txt");

	// H2, electron impact ionization
	h2_electron_ionization_data h2_ioniz_data;

	electron_impact_ionization* h2_ioniz_cs
		= new electron_impact_ionization(h2_ioniz_data, verbosity);

	electron_impact_ionization_relativistic* h2_ioniz_cs_rel
		= new electron_impact_ionization_relativistic(h2_ioniz_data, verbosity);

	// H2 pure rotational excitation,
	// J = 0 -> 2
	path = "coll_h2/MCCC-el-H2-rot-X-X/vi=0/Ji=0/";
	fname = path + "MCCC-el-H2-X1Sg_vf=0_Jf=2.X1Sg_vi=0_Ji=0.txt";
	
	cross_section_table_mccc * h2_rot_cs02 
		= new cross_section_table_mccc(data_path, fname);
	
	// H2 vibrational excitation, 
	// the sum of (v, J) = (0, ji) -> (1, ji) and (0, ji) -> (1, ji + 2)
	ji = 0;
	
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
	
	cross_section_table_mccc* h2_vibr_rot_cs1_0
		= new cross_section_table_mccc(data_path, fname);

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

	cross_section_table_mccc* h2_vibr_rot_cs1_2
		= new cross_section_table_mccc(data_path, fname);

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

	cross_section_table_mccc* h2_vibr_rot_cs2_0
		= new cross_section_table_mccc(data_path, fname);

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

	cross_section_table_mccc* h2_vibr_rot_cs2_2
		= new cross_section_table_mccc(data_path, fname);

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

	cross_section_table_mccc* h2_vibr_rot_cs3_0
		= new cross_section_table_mccc(data_path, fname);

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

	cross_section_table_mccc* h2_vibr_rot_cs3_2
		= new cross_section_table_mccc(data_path, fname);

	// H2 electronic excitation (singlet)
	vi = 0;

	// S1g+(X) -> S1u+(B), 
	cross_section_table_mccc **h2_bstate_cs 
		= new cross_section_table_mccc * [MAX_H2_VSTATES_B1SU];
	
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
		h2_bstate_cs[vf] = new cross_section_table_mccc(data_path, fname);
	}

	// S1g+(X) -> P1u(C-/+)
	cross_section_table_mccc** h2_cstate_cs 
		= new cross_section_table_mccc * [MAX_H2_VSTATES_C1PU];
	
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
		h2_cstate_cs[vf] = new cross_section_table_mccc(data_path, fname);
	}
	
	// S1g+(X) -> S1u+(Bp)
	cross_section_table_mccc ** h2_bpstate_cs 
		= new cross_section_table_mccc * [MAX_H2_VSTATES_BP1SU];
	
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
		h2_bpstate_cs[vf] = new cross_section_table_mccc(data_path, fname);
	}

	// S1g+(X) -> P1u(D-/+)
	cross_section_table_mccc ** h2_dstate_cs 
		= new cross_section_table_mccc * [MAX_H2_VSTATES_D1PU];
	
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
		h2_dstate_cs[vf] = new cross_section_table_mccc(data_path, fname);
	}
	
	// Electronic dissociative excitation of H2
	// S1g+(X) -> S1u+(B),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-B1Su_DE.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc *h2_bstate_diss_cs 
		= new cross_section_table_mccc(data_path, fname);
	
	// S1g+(X) -> P1u(C),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-C1Pu_DE.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc *h2_cstate_diss_cs 
		= new cross_section_table_mccc(data_path, fname);
	
	// S1g+(X) -> S1u+(Bp),	
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-Bp1Su_DE.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc *h2_bpstate_diss_cs 
		= new cross_section_table_mccc(data_path, fname);
	
	// S1g+(X) -> P1u(D),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-D1Pu_DE.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc * h2_dstate_diss_cs 
		= new cross_section_table_mccc(data_path, fname);
	
	// dissociative excitation to triplet states,
	// S1g(X) -> S3u+(b) (dissociative state),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-b3Su_DE.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc * h2_3bstate_diss_cs 
		= new cross_section_table_mccc(data_path, fname);
	
	// S1g(X) -> S3u+(h) (bound state),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-h3Sg_DE.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_3hstate_diss_cs
		= new cross_section_table_mccc(data_path, fname);

	// S1g(X) -> S3u+(e) (bound state),
	ss.clear();
	ss.str("");

	ss << "coll_h2/MCCC-el-H2-X1Sg-excitation/vi=";
	ss << vi;
	ss << "/MCCC-el-H2-e3Su_DE.X1Sg_vi=";
	ss << vi;
	ss << ".txt";

	fname = ss.str();
	cross_section_table_mccc* h2_3estate_diss_cs
		= new cross_section_table_mccc(data_path, fname);

	// Saving energy losses:
	fname = output_path + "energy_loss_electron_h2.txt";
	output.open(fname.c_str());

	output << left << "!Energy losses in [cm2 eV] of electron in the H2 (normalized to one H2 molecule)," << endl
		<< "!for MCCC cross sections - at higher energies, the cross sections are extrapolated as a power law, cs = cs_0*(E/E_0)^(gamma)" << endl
		<< "!Coulomb losses - Swartz et al., J. of Geophys. Res. 76, p. 8425 (1971);" << endl
		<< "!Vibrational excitation v = 0 -> 1: the sum of (v,j) = (0,0) -> (1,0) and (0,0) -> (1,2); the same for v = 0 -> 2;" << endl;

	output << left << setw(14) << "!Energy(eV) "
		<< setw(12) << "m.t."
		<< setw(12) << "ioniz"
		<< setw(12) << "ioniz_rel"
		<< setw(12) << "J=0-2" 
		<< setw(12) << "v=0-1"
		<< setw(12) << "v=0-2"
		<< setw(12) << "el-singlet " 
		<< setw(12) << "diss-singlet " 
		<< setw(12) << "diss-triplet " 
		<< setw(12) << "coulomb" << endl;

	for (i = 0; i < nb_of_energies; i++) 
	{
		output.precision(6);
		output << left << setw(14) << electron_energies[i];
		output.precision(3);
		
		enloss_mt_h2 = (*elastic_h2_el_cs)(electron_energies[i]) * electron_energies[i] * ELECTRON_MASS / ATOMIC_MASS_UNIT;
		output << left << setw(12) << enloss_mt_h2;

		enloss_ion_h2 = h2_ioniz_cs->get_energy_loss(electron_energies[i]);
		enloss_ion_h2_rel = h2_ioniz_cs_rel->get_energy_loss(electron_energies[i]);
		
		output << left << setw(12) << enloss_ion_h2 << setw(12) << enloss_ion_h2_rel;

		enloss_rot_h2 = (*h2_rot_cs02)(electron_energies[i]) * h2_rot_cs02->get_threshold_energy();
		output << left << setw(12) << enloss_rot_h2;

		// (v, J) = (0,0) -> (1,0) and (0,0) -> (1,2)
		enloss_vibr_h2_01 = (*h2_vibr_rot_cs1_0)(electron_energies[i]) * h2_vibr_rot_cs1_0->get_threshold_energy();
		enloss_vibr_h2_01 += (*h2_vibr_rot_cs1_2)(electron_energies[i]) * h2_vibr_rot_cs1_2->get_threshold_energy();
		output << left << setw(12) << enloss_vibr_h2_01;

		// v = 0 -> 2
		enloss_vibr_h2_02 = (*h2_vibr_rot_cs2_0)(electron_energies[i]) * h2_vibr_rot_cs2_0->get_threshold_energy();
		enloss_vibr_h2_02 += (*h2_vibr_rot_cs2_2)(electron_energies[i]) * h2_vibr_rot_cs2_2->get_threshold_energy();
		output << left << setw(12) << enloss_vibr_h2_02;

		enloss_singlet_h2 = 0.;
		for (vf = 0; vf < MAX_H2_VSTATES_B1SU; vf++) {
			enloss_singlet_h2 += (*h2_bstate_cs[vf])(electron_energies[i]) * h2_bstate_cs[vf]->get_threshold_energy();
		}

		for (vf = 0; vf < MAX_H2_VSTATES_C1PU; vf++) {
			enloss_singlet_h2 += (*h2_cstate_cs[vf])(electron_energies[i]) * h2_cstate_cs[vf]->get_threshold_energy();
		}

		for (vf = 0; vf < MAX_H2_VSTATES_BP1SU; vf++) {
			enloss_singlet_h2 += (*h2_bpstate_cs[vf])(electron_energies[i]) * h2_bpstate_cs[vf]->get_threshold_energy();
		}

		for (vf = 0; vf < MAX_H2_VSTATES_D1PU; vf++) {
			enloss_singlet_h2 += (*h2_dstate_cs[vf])(electron_energies[i]) * h2_dstate_cs[vf]->get_threshold_energy();
		}
		output << left << setw(12) << enloss_singlet_h2;
		
		enloss_diss_singlet_h2 = (*h2_bstate_diss_cs)(electron_energies[i]) * h2_bstate_diss_cs->get_threshold_energy();
		enloss_diss_singlet_h2 += (*h2_cstate_diss_cs)(electron_energies[i]) * h2_cstate_diss_cs->get_threshold_energy();
		enloss_diss_singlet_h2 += (*h2_bpstate_diss_cs)(electron_energies[i]) * h2_bpstate_diss_cs->get_threshold_energy();
		enloss_diss_singlet_h2 += (*h2_dstate_diss_cs)(electron_energies[i]) * h2_dstate_diss_cs->get_threshold_energy();

		output << left << setw(12) << enloss_diss_singlet_h2;

		enloss_diss_triplet_h2 = (*h2_3bstate_diss_cs)(electron_energies[i]) * h2_3bstate_diss_cs->get_threshold_energy();
		//enloss_diss_triplet_h2 += (*h2_3hstate_diss_cs)(electron_energies[i]) * h2_3hstate_diss_cs->get_threshold_energy();
		//enloss_diss_triplet_h2 += (*h2_3estate_diss_cs)(electron_energies[i]) * h2_3estate_diss_cs->get_threshold_energy();

		output << left << setw(12) << enloss_diss_triplet_h2;

		// energy loss in [eV cm2], per H2,
		enloss_coulomb = 2. * calc_coulomb_losses_thermal_electrons(electron_energies[i], ionization_fraction, thermal_electron_temp) 
			/ electron_velocities[i];

		output << left << setw(12) << enloss_coulomb;
		output << endl;
	}
	output.close();

	// Saving cross sections,
    fname = output_path + "cs_electron_h2.txt";
	output.open(fname.c_str());

	output << left << "!Cross sections in [cm2] of electron interaction with the H2," << endl
		<< "!for MCCC cross sections - at higher energies, the cross sections are extrapolated as a power law, cs = cs_0*(E/E_0)^(gamma)" << endl
		<< setw(14) << "!Energy(eV) "
		<< setw(12) << "m.t."
		<< setw(12) << "ion"
		<< setw(12) << "ion_rel"
		<< setw(12) << "J=0-2" 
		<< setw(12) << "v=0-1" 
		<< setw(12) << "v=0-2"
		<< setw(12) << "v=0-3"
		<< setw(12) << "el-singlet " 
		<< setw(12) << "diss-singlet "
		<< setw(12) << "diss-triplet " << endl;

	for (i = 0; i < nb_of_energies; i++)
	{
		output.precision(6);
		output << left << setw(14) << electron_energies[i];

		output.precision(3);
		output << left << setw(12) << (*elastic_h2_el_cs)(electron_energies[i])
			<< setw(12) << h2_ioniz_cs->get_int_cs(electron_energies[i])
			<< setw(12) << h2_ioniz_cs_rel->get_int_cs(electron_energies[i])
			<< setw(12) << (*h2_rot_cs02)(electron_energies[i]);

		// (v, J) = (0,0) -> (1,0) and (0,0) -> (1,2)
		cs_vibr_h2_01 = (*h2_vibr_rot_cs1_0)(electron_energies[i]);
		cs_vibr_h2_01 += (*h2_vibr_rot_cs1_2)(electron_energies[i]);
		output << left << setw(12) << cs_vibr_h2_01;

		// v = 0 -> 2
		cs_vibr_h2_02 = (*h2_vibr_rot_cs2_0)(electron_energies[i]);
		cs_vibr_h2_02 += (*h2_vibr_rot_cs2_2)(electron_energies[i]);
		output << left << setw(12) << cs_vibr_h2_02;

		// v = 0 -> 3
		cs_vibr_h2_03 = (*h2_vibr_rot_cs3_0)(electron_energies[i]);
		cs_vibr_h2_03 += (*h2_vibr_rot_cs3_2)(electron_energies[i]);
		output << left << setw(12) << cs_vibr_h2_03;

		cs_singlet_h2 = 0.;
		for (vf = 0; vf < MAX_H2_VSTATES_B1SU; vf++) {
			cs_singlet_h2 += (*h2_bstate_cs[vf])(electron_energies[i]);
		}

		for (vf = 0; vf < MAX_H2_VSTATES_C1PU; vf++) {
			cs_singlet_h2 += (*h2_cstate_cs[vf])(electron_energies[i]);
		}

		for (vf = 0; vf < MAX_H2_VSTATES_BP1SU; vf++) {
			cs_singlet_h2 += (*h2_bpstate_cs[vf])(electron_energies[i]);
		}

		for (vf = 0; vf < MAX_H2_VSTATES_D1PU; vf++) {
			cs_singlet_h2 += (*h2_dstate_cs[vf])(electron_energies[i]);
		}
		output << left << setw(12) << cs_singlet_h2;


		cs_diss_singlet_h2 = (*h2_bstate_diss_cs)(electron_energies[i]);
		cs_diss_singlet_h2 += (*h2_cstate_diss_cs)(electron_energies[i]);
		cs_diss_singlet_h2 += (*h2_bpstate_diss_cs)(electron_energies[i]);
		cs_diss_singlet_h2 += (*h2_dstate_diss_cs)(electron_energies[i]);

		output << left << setw(12) << cs_diss_singlet_h2;

		cs_diss_triplet_h2 = (*h2_3bstate_diss_cs)(electron_energies[i]);
		//cs_diss_triplet_h2 += (*h2_3hstate_diss_cs)(electron_energies[i]);
		//cs_diss_triplet_h2 += (*h2_3estate_diss_cs)(electron_energies[i]);

		output << left << setw(12) << cs_diss_triplet_h2;
		output << endl;
	}
	output.close();

	// Saving cross sections,
	// only vibrational excitations,
	fname = output_path + "cs_electron_h2_rovibr.txt";
	output.open(fname.c_str());

	output << left << "!Cross sections in [cm2] of electron interaction with the H2," << endl
		<< "!for MCCC cross sections - at higher energies, the cross sections are extrapolated as a power law, cs = cs_0*(E/E_0)^(gamma)" << endl
		<< "!(v_low - v_up) j_low - j_up" << endl
		<< setw(14) << "!Energy(eV) "
		<< setw(12) << "(0-1)0-0"
		<< setw(12) << "(0-1)0-2"
		<< setw(12) << "(0-2)0-0"
		<< setw(12) << "(0-2)0-2"
		<< setw(12) << "(0-3)0-0"
		<< setw(12) << "(0-3)0-2" << endl;

	for (i = 0; i < nb_of_energies / 2; i++)
	{
		output.precision(6);
		output << left << setw(14) << electron_energies[i];
		output.precision(3);

		// (v, J) = (0,0) -> (1,0) and (0,0) -> (1,2)
		output << left << setw(12) << (*h2_vibr_rot_cs1_0)(electron_energies[i]) 
			<< setw(12) << (*h2_vibr_rot_cs1_2)(electron_energies[i]);
		// v = 0 -> 2
		output << left << setw(12) << (*h2_vibr_rot_cs2_0)(electron_energies[i])
			<< setw(12) << (*h2_vibr_rot_cs2_2)(electron_energies[i]);
		
		// v = 0 -> 3
		output << left << setw(12) << (*h2_vibr_rot_cs3_0)(electron_energies[i])
			<< setw(12) << (*h2_vibr_rot_cs3_2)(electron_energies[i]);

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


void calc_helium_lifetimes(const string& output_path, const string& data_path)
{
	int verbosity = 1;
	int i, j, nb_lev_hei, isotope;
	double* life_times;

	string fname;
	ofstream output;

	nb_lev_hei = 31;
	molecule ion_hei("HeI", isotope = 1, 4. * ATOMIC_MASS_UNIT);

	ion_diagram* hei_di = new ion_diagram(data_path, ion_hei, nb_lev_hei, verbosity);
	ion_einstein_coeff* hei_einst = new ion_einstein_coeff(data_path, hei_di, verbosity);

	life_times = new double[nb_lev_hei];
	memset(life_times, 0, nb_lev_hei*sizeof(double));

	for (i = 1; i < nb_lev_hei; i++) {
		for (j = 0; j < i; j++) {
			life_times[i] += hei_einst->arr[i][j];
		}
	}

	fname = output_path + "helium_state_lifetimes.txt";
	output.open(fname.c_str());

	output << "!Lifetimes of Helium atom (HeI) states, in [s-1]" << endl;
	output.precision(3);

	for (i = 0; i < nb_lev_hei; i++) {
		output << left << setw(5) << i << setw(12) << hei_di->lev_array[i].name << setw(12) << life_times[i] << endl;
	}
	output.close();

	delete hei_di;
	delete hei_einst;
}

// Is not used,
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
