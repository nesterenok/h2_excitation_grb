#pragma once

// SUNDIALS CVODE solver headers
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */

#include "special_functions.h"
#include "h2_microphysics.h"
#include "spectroscopy.h"
#include "coll_rates_h2.h"
#include "h2_uv_pump_data.h"
#include <cstring>

//
// Electronic H2 states, taken into account: 
//      ground state   - S1g+(X), 
//      singlet states - S1u+(B,Bp), P1u(C+-,D+-)
//      triplet states - 
//      triplet states (dissociative) - S3u+(b)
// Transitions:
//      S1g+(X) -> S1u(B,Bp),    singlet transition, dLambda = 0, Lambda = S = 0 
//      S1g+(X) -> P1u(C+-,D+-), singlet transition, dLambda = 1, Lambda = P = 1 
//

// Contribution of triplet states to the dissociation (cross section at maximum):
// S3u+(b) - 1.64*a0^2, S3u+(e) - 0.024*a0^2, S3u+(h) - 0.045*a0^2,

//
// Structure of the data array:
// electron_spectrum [cm-3] (nb of electrons per cm-3 in energy interval), 
// concentration [cm-3] e-, H, H+, H2, H2+, He, He+, He++,
// H2 level population densities, [cm-3] 
// HeI level populations, [cm-3]
// gas temperature (of neutrals and ions) [K], dust grain charge


// auxiliary class to contain data on electron spectrum evolution
struct spectra_data {
    int i0;
    double w1, w2;
    spectra_data() : i0(0), w1(0), w2(0) { ; }
};

class elspectra_evolution_data
{
protected:
    bool dust_is_presented;
    int el_nb, h_nb, hp_nb, h2_nb, h2p_nb, he_nb, hep_nb, hepp_nb, 
        nb_lev_h2, nb_lev_h2_b, nb_lev_h2_cplus, nb_lev_h2_cminus, nb_lev_h2_bp, nb_lev_h2_dplus, nb_lev_h2_dminus, nb_lev_hei,
        h2_ioniz_min_nb, he_ioniz_min_nb, h_ioniz_min_nb, h2p_ioniz_min_nb, hep_ioniz_min_nb, 
        h2_b1su_min_nb, h2_c1pu_min_nb, h2_bp1su_min_nb, h2_d1pu_min_nb, 
        h2_b1su_diss_min_nb, h2_c1pu_diss_min_nb, h2_bp1su_diss_min_nb, h2_d1pu_diss_min_nb, h2_b3su_diss_min_nb,
        hei_min_nb, nb_coll_trans_hei, nb_of_el_energies, nb_of_equat, h2eq_nb, heieq_nb, physeq_nb, min_grain_charge;
    int* indices;

    double conc_h_tot, ioniz_fract, grain_cs, grain_nb_density,
        enloss_rate_mt, enloss_rate_h2_rot, enloss_rate_h2_rot_pos, enloss_rate_h2_vibr, enloss_rate_h2_vibr_pos,
        enloss_rate_h2_singlet, enloss_diss_h2_singlet, enloss_diss_h2_triplet, enloss_rate_ioniz,
        enloss_rate_hei, enloss_rate_coulomb_el, conc_n, conc_i, energy_gain_n, energy_gain_e, nb_gain_n, nb_gain_i, 
        diss_decay_heating_rate, neutral_coll_heating_rate;

    double h2_solomon_diss_rate, h2_diss_exc_singlet_rate, h2_diss_exc_triplet_rate, hei_exc_rate, h2_excit_electr_rate, 
        h2_excit_electr_bs_rate, h2_excit_electr_cp_rate, h2_excit_vibr_rate, h2_excit_vibr_v1_rate, h2_excit_vibr_v2_rate, h2_excit_vibr_v3_rate, 
        h2_excit_rot_rate, h2_vstate_el_exc_rate;

    // in eV, in electron_energies - the centre of energy interval is provided,
    double* electron_energies_grid, * electron_energy_bin_size, *electron_energies, *electron_velocities;      
    
    double* el_enloss_h2_mt, * el_enloss_he_mt, * el_enloss_elthermal;
    double* h2_ioniz_rates_tot, * he_ioniz_rates_tot, * h_ioniz_rates_tot, * h2p_ioniz_rates_tot, * hep_ioniz_rates_tot;
    double** h2_ioniz_rates, ** he_ioniz_rates, ** h_ioniz_rates, ** h2p_ioniz_rates, ** hep_ioniz_rates;

    double** h2_bstate_rates, ** h2_cstate_rates, ** h2_bpstate_rates, ** h2_dstate_rates;
    double** h2_bstate_diss_rates, ** h2_cstate_diss_rates, ** h2_bpstate_diss_rates, ** h2_dstate_diss_rates, ** h2_3bstate_diss_rates;
    double** h2_rot_rates, ** h2_rovibr_rates;
    double** hei_rates;
    double*  coll_partn_conc;

    const energy_diagram* hei_di, *h2_di, *h2_di_b, * h2_di_cminus, * h2_di_cplus, *h2_di_bp, * h2_di_dminus, * h2_di_dplus;
    const einstein_coeff* h2_einst, *hei_einst;
    
    const h2_collisions*  h2_coll;

    std::vector<transition> lyman_band_h2, werner_plus_band_h2, werner_minus_band_h2, bp_band_h2, dplus_band_h2, dminus_band_h2;
    std::vector<h2_energy_level_param> h2_b_state_data, h2_cplus_state_data, h2_cminus_state_data, 
        h2_bp_state_data, h2_dplus_state_data, h2_dminus_state_data;

    const cross_section_table_vs1 *elastic_h2_el_cs, * elastic_he_el_cs;
   
    cross_section_table_mccc** h2_bstate_cs, ** h2_cstate_cs, ** h2_bpstate_cs, ** h2_dstate_cs;
    cross_section_table_mccc** h2_bstate_diss_cs, ** h2_cstate_diss_cs, ** h2_bpstate_diss_cs, ** h2_dstate_diss_cs, ** h2_3bstate_diss_cs;
    cross_section_table_mccc** h2_rot_cs, **h2_rovibr_cs;

    h2_ioniz_cs_table_straub1996 * h2_ioniz_cs_straub;
    h2_dissioniz_cs_table_straub1996 * h2_dissioniz_cs_straub;

    electron_impact_ionization* h2_ioniz_cs, * he_ioniz_cs, * h_ioniz_cs, * h2p_ioniz_cs, * hep_ioniz_cs;
    
    spectra_data** h2_ioniz_indexes, ** he_ioniz_indexes, ** h_ioniz_indexes, ** h2p_ioniz_indexes, ** hep_ioniz_indexes;

    honl_london_singlet_dl0 hl_singlet_dl0;
    honl_london_singlet_dl1 hl_singlet_dl1;
    err_func error_function;

    // Helium
    // the arrays of the initial and final level nb, which are included in the array of cross sections,
    std::vector<int> hei_init_level_nbs, hei_fin_level_nbs;
    std::vector<const hei_electron_excitation*> hei_cs;

    // Ionization,
    // version without averaging over energy intervals,
    void init_tables_ionization(electron_impact_ionization*, double *& rates_tot, double**& rates, spectra_data**& indexes, 
        int &el_ion_min_nb, int verbosity = 1);
    
    // Note: parameter 'rate' stores the rate of the current reaction, while energy losses are summed
    void derivatives_ionization(const realtype* y_data, realtype* ydot_data, const electron_impact_ionization*, 
        double** rates, const double* rates_tot, spectra_data** indexes, int el_ion_min_nb, int target_nb, double& rate, double& enloss_rate);

    // H2 electronic excitation,
    void init_tables_h2_electronic_exc(cross_section_table_mccc**, double**& rates, int vibr_qnb_final_max, int verbosity = 1);

    // calculation the nb of electron energy interval that corresponds to the collisional transition with the minimal energy, 
    int calc_min_energy_of_transition(const energy_diagram* h2_di_elexc, int vibr_qnb_final_max);
 
    // all parameters with & are summed
    // thermal_energy - kinetic energy released in the dissociation, in erg cm-3 s-1,
    // S1g(X) -> S1u(B), S1u(Bp), dL = 0, dJ = +/-1,
    void derivatives_h2_electronic_exc_dl0(const realtype* y_data, realtype* ydot_data, const energy_diagram* h2_di_exc, 
        const std::vector<h2_energy_level_param> h2_exc_state_data, double** rates, int vibr_qnb_final_max, int el_min_nb, 
        double& enloss_rate, double& diss_rate, double & excit_rate, double& th_energy_deriv, double& vstate_exc_rate);

    // S1g(X) -> P1u(C), P1u(D), dJ = -1, 0, 1
    void derivatives_h2_electronic_exc_dl1(const realtype* y_data, realtype* ydot_data, 
        const energy_diagram* h2_di_plus, const energy_diagram* h2_di_minus,
        const std::vector<h2_energy_level_param> h2_state_plus_data, const std::vector<h2_energy_level_param> h2_state_minus_data, double** rates, 
        int vibr_qnb_final_max, int el_min_nb, double& enloss_rate, double& diss_rate, double& excit_rate, double& th_energy_deriv, double& vstate_exc_rate);

    // H2 dissociation 
    // through electronic excitation (direct),
    // the energy lost by the electron is equal to the threshold value given in the cross section data file 
    // (minus the energy of the initial state of H2 with given J),
    void init_tables_h2_electronic_diss(cross_section_table_mccc**, double**& rates, int& el_min_nb, int verbosity = 1);

    // all parameters with & (energy loss and dissociation rate) are summed,
    void derivatives_h2_electronic_diss(const realtype* y_data, realtype* ydot_data, cross_section_table_mccc** cs, double** rates, int el_min_nb,
        double& enloss_rate, double& diss_rate);

    // pure H2 rotational excitation,
    void init_tables_h2_rotational_exc(cross_section_table_mccc**, double**& rates, int verbosity = 1);

    void derivatives_h2_rotational_exc(const realtype* y_data, realtype* ydot_data, double**& rates, double& enloss_rate, double& enloss_rate_pos, 
        double & excit_rate);

    // H2 ro-vibrational excitation, vibration nb initial != final,
    void init_tables_h2_rovibr_exc(cross_section_table_mccc**, double**& rates, int verbosity = 1);

    void derivatives_h2_rovibr_exc(const realtype* y_data, realtype* ydot_data, double**& rates, double& enloss_rate, double& enloss_rate_pos, 
        double& excit_rate, double & excit_v1_rate, double & excit_v2_rate, double & excit_v3_rate);

    // All data relevant to HeI excitation are initialized in this subroutine, no input parameters, class members are used,
    // data by Ralchenko et al. Atomic Data and Nuclear Data Tables 94, 603 (2008),
    void init_hei_cross_sections(const std::string& data_path);
    void init_tables_hei_electronic_exc(const std::string& data_path);

    // all parameters with & (energy loss and excitation rate) are summed,
    void derivatives_hei_exc(const realtype* y_data, realtype* ydot_data, double& enloss_rate, double& exc_rate);

public:
    // function defining the ODE system for the evolution of the parameters
    int f(realtype t, N_Vector y, N_Vector ydot);
    
    void set_dust_parameters(bool dust_is_presented, double radius, double grain_nb_density);

    int get_nb_of_equat() const { return nb_of_equat; }
    int get_nb_of_el_en() const { return nb_of_el_energies; }
    int get_nb_of_h2_lev() const { return nb_lev_h2; }
    int get_nb_of_hei_lev() const { return nb_lev_hei; }
    int get_h2eq_nb() const { return h2eq_nb; }
    int get_heieq_nb() const { return heieq_nb; }
    int get_physeq_nb() const { return physeq_nb; }
    
    int get_vibr_nb_h2(int level_nb) const;   // returns the vibration quantum nb of the ground electronic state of H2,
    int get_level_nb_h2(int v, int j) const;  // returns the nb of the H2 level belonging to the ground electronic state, 

    // energy losses are < 0 if electrons loses energy, [eV cm-3 s-1]
    // energy loss via excitation of singlet (h2_electr_sin) and triplet states (h2_electr_tri), including excitation and dissociative excitation for singlet states,
    // neutral heating rate due to collisions H2-H2, H2-He, [eV cm-3 s-1]
    void get_el_energy_losses(double & mt, double & h2_rot, double& h2_rot_p, double &h2_vibr, double& h2_vibr_p, 
        double & h2_electr_sin, double & h2_electr_tri, double & ion, double & hei, double& el_coul, 
        double & diss_decay_heat, double & neut_heat_coll) const;
    
    // dissociation rates [cm-3 s-1],
    // through excitation of singlet states of H2 (diss_exc) and triplet (diss_exc_tr),
    void get_h2_diss_rates(double& sol_diss, double& diss_excit_sin, double& diss_excit_tri);

    // excitation rates [cm-3 s-1],
    void get_h2_process_rates(double & excit_electr_rate, double & excit_vibr_rate, double & excit_rot_rate, 
        double & excit_electr_bs_rate, double & excit_electr_cp_rate, double & excit_vibr_v1_rate, double & excit_vibr_v2_rate, 
        double& excit_vibr_v3_rate, double & el_excit_vstate_rate, double & hei_excit);

    double get_electron_energy(int i) const;      // returns the centre of the interval [i, i+1], energy in eV
    double get_electron_energy_bin(int i) const;  // returns the energy bin size, in eV,

    int get_nb_electron_energy_array(double energy) const;
    int get_nb_electron_energy_array(double energy, int i) const;  // given energy in eV, seed of nb value have to be given

    const energy_diagram* get_pointer_h2_levels() const { return h2_di; }
    
    // saving h2 level quantum numbers and energies of the ground electronic state, for which data are saved,
    void save_h2_levels_used(const std::string& output_path);

    void set_tolerances(N_Vector abs_tol);

    // concentration in cm-3, electron energy in eV,
    elspectra_evolution_data(const std::string& data_path, const std::string& output_path, double conc_h_tot, double ioniz_fract,
        const std::vector<double> &electron_energies_grid, int verb);

    ~elspectra_evolution_data();
};
