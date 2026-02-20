
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


h2_ioniz_cs_table_straub1996::h2_ioniz_cs_table_straub1996(const std::string& data_path, const std::string& name)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i;
	double en, cs1, cs2;

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
	
	input >> nb_cs;

	en_arr = new double[nb_cs];
	cs_arr = new double[nb_cs];

	for (i = 0; i < nb_cs; i++) {
		input >> en >> cs1 >> cs2;   // reading data on two channels,
		// the energy data in the file have dimension eV,
		en_arr[i] = en;
		// the cross section data have dimension 1.e-17 cm2 for H2+ channel and 1.e-18 cm2 for H+ channel
		cs_arr[i] = cs1 * 1.e-17 + cs2 * 1.e-18;
	}
	input.close();
	en_thr = en_arr[0];
	
	// used in extrapolation:
	// gamma = log(cs_arr[nb_cs - 1] / cs_arr[nb_cs - 3]) / log(en_arr[nb_cs - 1] / en_arr[nb_cs - 3]);
	gamma = -1.0808558;  // calculated explicitly based on three last points of the cross section data 
}


h2_dissioniz_cs_table_straub1996::h2_dissioniz_cs_table_straub1996(const std::string& data_path, const std::string& name)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i;
	double en, cs1, cs2;

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

	input >> nb_cs;

	en_arr = new double[nb_cs];
	cs_arr = new double[nb_cs];

	for (i = 0; i < nb_cs; i++) {
		input >> en >> cs1 >> cs2;   // reading data on two channels,
		// the energy data in the file have dimension eV,
		en_arr[i] = en;
		// the cross section data have dimension 1.e-17 cm2 for H2+ channel and 1.e-18 cm2 for H+ channel
		// only dissociative ionization is considered here,
		cs_arr[i] = cs2 * 1.e-18;
	}
	input.close();
	en_thr = en_arr[0];

	// used in extrapolation:
	// the same value as for ionization data, the ratio of H2+ and H+ channels will be the same at high energies
	gamma = -1.0808558;  
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

// Note, for accurate calculation, double ionization of He must be taken into account, He + e- -> He++ + e- + e- +e-
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
	: err(1.e-6)
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
// Electronic excitation of H2 ground state S1g+(X), 
cross_section_table_mccc::cross_section_table_mccc(const std::string& data_path, const std::string& name, bool is_extr_on, int verbosity)
	:cross_section_table(), is_extrapolation_on(is_extr_on)
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

	// the cross section at this energy may not be zero - the reaction without threshold,
	en_thr = en_arr[0];
	gamma = gamma2 = log(cs_arr[nb_cs - 1] / cs_arr[nb_cs - 2]) / log(en_arr[nb_cs - 1] / en_arr[nb_cs - 2]);
	
	// the value derived from ionization cross section data by Kim & Rudd, Physical Review A 50, p.3954 (1994); 
	// no relativistic effects, 
	// E = 500 - 1.e+3 eV  gamma = -0.7610825481
	// E = 1.e+3 - 1.e+4 eV  gamma2 = -0.8579127924
	if (gamma > -0.8579127924 && gamma < -0.7610825481) {
		gamma2 = -0.8579127924;
	}
	else if (gamma > -0.7610825481) {
		gamma = -0.7610825481;
		gamma2 = -0.8579127924;
	}

	en0 = 1.e+3;  // eV 
	if (en_arr[nb_cs - 1] < en0) {
		cs0 = cs_arr[nb_cs - 1] * pow(en0 / en_arr[nb_cs - 1], gamma);
	}
	else {
		en0 = en_arr[nb_cs - 1];
		cs0 = cs_arr[nb_cs - 1];
	}
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
		if (is_extrapolation_on) {
			if (energy < en0) {
				answ = cs_arr[nb_cs - 1] * pow(energy / en_arr[nb_cs - 1], gamma);
			}
			else {
				answ = cs0 * pow(energy / en0, gamma2);
			}
		}
		else {
			answ = 0.;
		}
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
//
// 
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


//
//
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

//
// is not used
cross_section_integral_2d::cross_section_integral_2d(electron_impact_ionization* ei, double e)
	:x1(0.), x2(0.), z1(0.), z2(0.), ionization(ei), err(e)
{
	;
}

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
{
	;
}

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
