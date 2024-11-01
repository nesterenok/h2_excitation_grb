#pragma once

// SUNDIALS CVODE solver headers
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */

#include "special_functions.h"
#include "h2_microphysics.h"
#include "spectroscopy.h"
#include "h2_uv_pump_data.h"
#include <cstring>

//
// Electronic H2 states, taken into account: 
//      ground state - S1g+(X), 
//      singlet states - S1u+(B,Bp), P1u(C+-,D+-)
//      triplet states - 
//      triplet states (dissociative) - S3u+(b)
// Transitions:
//      S1g+(X) -> S1u(B,Bp),    singlet transition, dLambda = 0, Lambda = S = 0 
//      S1g+(X) -> P1u(C+-,D+-), singlet transition, dLambda = 1, Lambda = P = 1 
//

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

struct spectra_data_adv {
    int i0, nb;
    double* weights;

    spectra_data_adv() : i0(0), nb(0), weights(0) { ; }
    spectra_data_adv(const spectra_data_adv& obj);
    spectra_data_adv& operator = (const spectra_data_adv& obj);
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
    
    double conc_h_tot, grain_cs, grain_nb_density,
        enloss_rate_mt, enloss_rate_h2_rot, enloss_rate_h2_vibr, enloss_rate_h2_singlet, enloss_rate_h2_triplet, enloss_rate_ioniz,
        enloss_rate_coloumb_ions, enloss_rate_coloumb_el, enloss_rate_hei, conc_n, conc_i, energy_gain_n, energy_gain_i, nb_gain_n, nb_gain_i, 
        h2_solomon_diss_rate, h2_diss_exc_singlet_rate, h2_diss_exc_triplet_rate, hei_exc_rate;

    // in eV, in electron_energies - the centre of energy interval is provided,
    double* electron_energies_grid, * electron_energy_bin_size, *electron_energies, *electron_velocities;      
    
    // energy loss values are positive in the array,
    double* el_enloss_h2_mt, * el_enloss_he_mt; 
    double* h2_ioniz_rates_tot, * he_ioniz_rates_tot, * h_ioniz_rates_tot, * h2p_ioniz_rates_tot, * hep_ioniz_rates_tot;
    double** h2_ioniz_rates, ** he_ioniz_rates, ** h_ioniz_rates, ** h2p_ioniz_rates, ** hep_ioniz_rates;

    const energy_diagram* hei_di, *h2_di, *h2_di_b, * h2_di_cminus, * h2_di_cplus, *h2_di_bp, * h2_di_dminus, * h2_di_dplus;
    const einstein_coeff* h2_einst, *hei_einst;

    std::vector<transition> lyman_band_h2, werner_plus_band_h2, werner_minus_band_h2, bp_band_h2, dplus_band_h2, dminus_band_h2;
    std::vector<h2_energy_level_param> h2_b_state_data, h2_cplus_state_data, h2_cminus_state_data, 
        h2_bp_state_data, h2_dplus_state_data, h2_dminus_state_data;

    const cross_section_table_vs1 *elastic_h2_el_cs, * elastic_he_el_cs;
   
    cross_section_table_mccc** h2_bstate_cs, ** h2_cstate_cs, ** h2_bpstate_cs, ** h2_dstate_cs;
    cross_section_table_mccc** h2_bstate_diss_cs, ** h2_cstate_diss_cs, ** h2_bpstate_diss_cs, ** h2_dstate_diss_cs, ** h2_3bstate_diss_cs;

    electron_impact_ionization* h2_ioniz_cs, * he_ioniz_cs, * h_ioniz_cs, * h2p_ioniz_cs, * hep_ioniz_cs;

    double** h2_bstate_rates;

    honl_london_singlet_dl0 hl_singlet_dl0;
    honl_london_singlet_dl1 hl_singlet_dl1;

    err_func error_function;

    // OLD
    const cross_section_table_vs1 *rot_h2_j02_cs, * rot_h2_j13_cs, * rot_h2_j24_cs, * rot_h2_j35_cs, *vibr_h2_v01_cs;
    const h2_excitation_vibr02_cs* vibr_h2_v02_cs;
    
    cross_section_table_mccc* h2_v0_bstate_diss_cs, *h2_v1_bstate_diss_cs, *h2_v0_cstate_diss_cs, *h2_v1_cstate_diss_cs, 
        * h2_v0_bpstate_diss_cs, * h2_v1_bpstate_diss_cs, * h2_v0_dstate_diss_cs, * h2_v1_dstate_diss_cs, 
        * h2_v0_3bstate_diss_cs, * h2_v1_3bstate_diss_cs;
    
    cross_section_table_mccc** h2_v0_bstate_cs, ** h2_v1_bstate_cs, ** h2_v0_cstate_cs, ** h2_v1_cstate_cs, 
        ** h2_v0_bpstate_cs, ** h2_v1_bpstate_cs, ** h2_v0_dstate_cs, ** h2_v1_dstate_cs;
    
    std::vector<int> hei_level_nbs;
    std::vector<const hei_electron_excitation*> hei_cs;


    double* h2_j02_rates, * h2_j24_rates, * h2_j13_rates, * h2_j35_rates, * h2_j20_rates, * h2_j42_rates, * h2_j31_rates, * h2_j53_rates,
        * h2_v01_rates, * h2_v02_rates, * elel_scatt_rates, * elion_scatt_rates,
        * h2_v0_bstate_diss_rates, * h2_v1_bstate_diss_rates, * h2_v0_cstate_diss_rates, * h2_v1_cstate_diss_rates,
        * h2_v0_bpstate_diss_rates, * h2_v1_bpstate_diss_rates, * h2_v0_dstate_diss_rates, * h2_v1_dstate_diss_rates,
        * h2_v0_3bstate_diss_rates, * h2_v1_3bstate_diss_rates;


    double** h2_v0_bstate_rates, ** h2_v1_bstate_rates, ** h2_v0_cstate_rates, ** h2_v1_cstate_rates, ** h2_v0_bpstate_rates, ** h2_v1_bpstate_rates,
        ** h2_v0_dstate_rates, ** h2_v1_dstate_rates, ** coloumb_rates, ** hei_rates;

    // weights are dimensionless, the sum of weights is equal to 1:
    spectra_data* h2_j02_indexes, * h2_j13_indexes, * h2_j24_indexes, * h2_j35_indexes,
        * h2_j20_indexes, * h2_j31_indexes, * h2_j42_indexes, * h2_j53_indexes;
    spectra_data** h2_v01_indexes, ** h2_v02_indexes;
    spectra_data** h2_ioniz_indexes, ** he_ioniz_indexes, ** h_ioniz_indexes, ** h2p_ioniz_indexes, ** hep_ioniz_indexes;
    spectra_data** h2_v0_bstate_indexes, ** h2_v1_bstate_indexes, ** h2_v0_cstate_indexes, ** h2_v1_cstate_indexes,
        ** h2_v0_bpstate_indexes, ** h2_v1_bpstate_indexes, ** h2_v0_dstate_indexes, ** h2_v1_dstate_indexes,
        ** h2_v0_bstate_diss_indexes, ** h2_v1_bstate_diss_indexes, ** h2_v0_cstate_diss_indexes, ** h2_v1_cstate_diss_indexes,
        ** h2_v0_bpstate_diss_indexes, ** h2_v1_bpstate_diss_indexes, ** h2_v0_dstate_diss_indexes, ** h2_v1_dstate_diss_indexes,
        ** h2_v0_3bstate_diss_indexes, ** h2_v1_3bstate_diss_indexes,
        ** coloumb_indexes, ** hei_indexes;
    //

    // version without averaging over energy intervals,
    void init_tables_ionization(electron_impact_ionization*, double *& rates_tot, double**& rates, spectra_data**& indexes, 
        int &el_ion_min_nb, int verbosity = 1);
    
    // parameter 'rate' stores the rate of the current reaction, while energy losses are summed
    void derivatives_ionization(const realtype* y_data, realtype* ydot_data, const electron_impact_ionization*, 
        double** rates, const double* rates_tot, spectra_data** indexes, int el_ion_min_nb, int target_nb, double& rate, double& enloss_rate);

    // 
    void init_tables_h2_electronic_exc(const energy_diagram* h2_di_elexc, cross_section_table_mccc**, double**& rates,
        int vibr_qnb_final_max, int& el_min_nb, int verbosity = 1);

   // S1g(X) -> S1u(B), S1u(Bp) dJ = +/-1, electronic states
   // S1g(X) -> P1u(C), P1u(D), dJ = -1, 0, 1
   // all parameters with & are summed
   // thermal_energy - kinetic energy released in the dissociation, in erg cm-3 s-1,
    void derivatives_h2_electronic_exc(const realtype* y_data, realtype* ydot_data, const energy_diagram* h2_di_exc, 
        const std::vector<h2_energy_level_param> h2_exc_state_data, cross_section_table_mccc** , double** rates, 
        int vibr_qnb_final_max, int el_min_nb, double& enloss_rate, double& diss_rate, double& thermal_energy);


    // OLD
    // rotational excitation and de-excitation are taken into account,
    void init_tables_h2_rotational_exc(const cross_section*, double*& rates_up, spectra_data*& indexes_up, 
        double*& rates_down, spectra_data*& indexes_down, int j0, int verbosity = 1);
    void init_tables_h2_vibrational_exc(const cross_section*, double*& rates, spectra_data**& indexes, int vibr_qnb_final, int verbosity = 1);

    // the nb of final electronic sate of H2 is used in the routine (determines the nb of branches),
    void init_tables_h2_electronic_exc_old(const energy_diagram * h2_di_elexc, cross_section_table_mccc**, double**& rates, spectra_data**& indexes,
        int vibr_qnb_init, int max_h2_vstates, int & el_min_nb, int verbosity = 1);

    // the energy lost by the electron is equal to the threshold value given in the cross section data file 
    // (minus the energy of the initial state of H2 with given J),
    void init_tables_h2_electronic_diss(cross_section_table_mccc*, double*& rates, spectra_data**& indexes,
        int vibr_qnb_init, int& el_min_nb, int verbosity = 1);

    // All data relevant to HeI excitation are initialized in this subroutine, no input parameters, class members are used,
    void init_tables_hei_electronic_exc(const std::string& data_path);


    // j0 - the angular momentum of the lower H2 level, energy losses are summed
    void derivatives_h2_rotational_exc(const realtype* y_data, realtype* ydot_data, const double* rates_up, const spectra_data* indexes_up, 
        const double* rates_down, const spectra_data* indexes_down, int j0, double & enloss_rate);

    void derivatives_h2_vibrational_exc(const realtype* y_data, realtype* ydot_data,
        const double* rates, spectra_data** indexes, int vibr_qnb_final, double & enloss_rate);

    // S1g(X) -> S1u(B, Bp), dJ = +/-1, 
    // thermal_energy - kinetic energy released in the dissociation, in erg cm-3 s-1,
    // final vibration quantum nb is 0,1,..,vibr_qnb_final_max-1
    // all parameters with & are summed
    void derivatives_h2_electronic_exc_old(const realtype* y_data, realtype* ydot_data, const energy_diagram* h2_di_exc, 
        const std::vector<h2_energy_level_param> h2_exc_state_data, double** rates, spectra_data** indexes, 
        int vibr_qnb_init, int vibr_qnb_final_max, int el_min_nb, double& enloss_rate, double &diss_rate, double &thermal_energy);

    // S1g(X) -> P1u(C), dJ = -1, 0, 1,
    void derivatives_h2_electronic_exc_2(const realtype* y_data, realtype* ydot_data, 
        const energy_diagram* h2_di_plus, const energy_diagram* h2_di_minus, 
        const std::vector<h2_energy_level_param> h2_state_plus_data, const std::vector<h2_energy_level_param> h2_state_minus_data,
        double** rates, spectra_data** indexes, int vibr_qnb_init, int vibr_qnb_final_max, int el_min_nb, 
        double& enloss_rate, double& diss_rate, double& thermal_energy);

    // all parameters with & (energy loss and dissociation rate) are summed,
    void derivatives_h2_electronic_diss(const realtype* y_data, realtype* ydot_data,
        const double* rates, spectra_data** indexes, int vibr_qnb_init, int el_min_nb, double& enloss_rate, double& diss_rate);

    // all parameters with & (energy loss and dissociation rate) are summed,
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
    int get_vibr_nb_h2(int level_nb) const;  // returns the vibration quantum nb of the ground electronic state of H2,

    // energy losses are < 0 if electrons loses energy, [eV cm-3 s-1]
    // energy loss via excitation of singlet (h2_electr_sin) and triplet states (h2_electr_tri), including excitation and dissociative excitation for singlet states,
    void get_el_energy_losses(double & mt, double & h2_rot, double &h2_vibr, double & h2_electr_sin, double & h2_electr_tri, double & ion, 
        double & col_ions, double & col_el, double & hei) const;
    
    // rate of dissociation [cm-3 s-1],
    // through excitation of singlet states of H2 (diss_exc) and triplet (diss_exc_tr),
    void get_h2_process_rates(double & sol_diss, double & diss_exc, double & diss_exc_tr, double & hei_exc);

    const energy_diagram* get_pointer_h2_levels() const { return h2_di; }
    const energy_diagram* get_pointer_hei_levels() const { return hei_di; }

    double get_electron_energy(int i) const;      // returns the centre of the interval [i, i+1], energy in eV
    double get_electron_energy_bin(int i) const;  // returns the energy bin size, in eV,
    int get_nb_electron_energy_array(double energy, int i) const;  // given energy in eV, seed of nb value have to be given

    void set_tolerances(N_Vector abs_tol);

    // concentration in cm-3, electron energy in eV,
    // parameters of electron sepctra are given,
    elspectra_evolution_data(const std::string& data_path, double conc_h_tot, const std::vector<double> &electron_energies_grid, int verb);
    ~elspectra_evolution_data();
};
