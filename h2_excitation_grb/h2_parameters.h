#pragma once

const std::string chemical_species[] = { "e-", "H", "H+", "H2", "H2+", "He", "He+", "He++",
    "CI",  "CII",  "CIII",  "CIV",  "CV",  "CVI",  "CVII",
    "NI",  "NII",  "NIII",  "NIV",  "NV",  "NVI",  "NVII",  "NVIII",
    "OI",  "OII",  "OIII",  "OIV",  "OV",  "OVI",  "OVII",  "OVIII",  "OIX",
    "NeI", "NeII", "NeIII", "NeIV", "NeV", "NeVI", "NeVII", "NeVIII", "NeIX", "NeX", "NeXI",
    "MgI", "MgII", "MgIII", "MgIV", "MgV", "MgVI", "MgVII", "MgVIII", "MgIX", "MgX", "MgXI", "MgXII", "MgXIII",
    "SiI", "SiII", "SiIII", "SiIV", "SiV", "SiVI", "SiVII", "SiVIII", "SiIX", "SiX", "SiXI", "SiXII", "SiXIII", "SiXIV", "SiXV",
    "SI",  "SII",  "SIII",  "SIV",  "SV",  "SVI",  "SVII",  "SVIII",  "SIX",  "SX",  "SXI",  "SXII",  "SXIII",  "SXIV",  "SXV",  "SXVI",  "SXVII",
    "FeI", "FeII", "FeIII", "FeIV", "FeV", "FeVI", "FeVII", "FeVIII", "FeIX", "FeX", "FeXI", "FeXII", "FeXIII", "FeXIV", "FeXV", "FeXVI", "FeXVII",
    "FeXVIII", "FeXIX", "FeXX", "FeXXI", "FeXXII", "FeXXIII", "FeXXIV", "FeXXV", "FeXXVI", "FeXXVII"
};

#define NB_OF_CHEM_SPECIES 8

// nb of specimen in the array of variables,
#define EL_NB 0
#define H_NB 1
#define H_P_NB 2
#define H2_NB 3
#define H2_P_NB 4
#define HE_NB 5
#define HE_P_NB 6
#define HE_PP_NB 7


//-------------------------------------------------
// Numerical parameters
//-------------------------------------------------

#define REL_ERROR_SOLVER 1.e-6
#define ABS_CONCENTRATION_ERROR_SOLVER 1.e-14  // specimen concentration (cm-3)
#define ABS_POPULATION_ERROR_SOLVER 1.e-14     // specimen level population density (cm-3)
#define ABS_ELSPECTRA_ERROR_SOLVER 1.e-14
#define ABS_PARAMETER_ERROR_SOLVER 1.e-11  // e.g., grain radius (cm), temperature (K)
#define MAX_CONV_FAILS_SOLVER 100  // default value is 10;
#define MAX_ERR_TEST_FAILS_SOLVER 14  // default value is 7;


#define MAX_NB_STEPS 15
#define NB_OF_BINS_PER_ORDER_EL 100
#define MINIMAL_ABUNDANCE 1.e-99  // for saving in file

// time nb per order is used to initiate the time grid, 
// this number can not be very large (otherwise the data arrays with electron spectra will be large);
#define NB_OF_BINS_PER_ORDER_TIME 30
#define NB_OF_TIME_STEPS 1  // for saving electron spectra, the spectra is saved with this step in time point numbers, 


#define ELECTRON_ENERGY_FIXED 1.   // eV, the energy grid intervals are equal at E < E_fixed,
#define MAX_ELECTRON_ENERGY 1.e+6  // eV, may be lower than in the input file with initial spectrum,
#define MIN_MODEL_TIME 1.e+2   // in s, minimal model time at which data are saved,
#define MAX_MODEL_TIME 1.e+8   // maximal model time in s (the simulations are done up to this time),

// the number of vibrational states of the ground electronic state of H2 taken into account in electron-impact excitation,
// must be less or equal than maximal, MAX_H2_VSTATES_X1SU,
// v qnb = 0, 1,.., NB_OF_H2_VSTATES_X1SU-1
// very small difference between nb = 5 and 15;
#define NB_OF_H2_VSTATES_X1SU 10

// ----------------------------------------------------
// Electron spectra evolution
//-----------------------------------------------------

// this constant is used in the initialization of the electron spectrum data,
#define ELECTRON_SPECTRUM_TEST 0
#define ELECTRON_CONCENTRATION_TEST  0.1  // cm-3

// 1 - initialization of H2 population densities from file,
// 0 - H2 molecules are in lowest energy levels, J = 0 and 1,
// Note, there is no check yet about the consistency of energy levels used in this code and those used in the simulations before, 
#define H2_POP_DENS_INIT 0
// this parameter is used if H2 population densities are not read from file,  
#define ORTHO_TO_PARA_RATIO 0.1

// Parameters of the electron energy losses
// 1 - switch on, 0 - switch off,      
#define ROVIBRATIONAL_EXC_LOSSES 1
#define HELIUM_EXC_LOSSES 1

// if switch off, ionization of species is not taken into account, the number of electrons must be conserved:
#define IONIZATION_LOSSES 1

//
// The energy loss of fast electron on thermal electrons, approximate formula,
// Swartz et al. "Analytic expression for the energy-transfer from photoelectrons to thermal electrons", J. of Geophys. Res. 76, p. 8425, 1971;
#define CALC_EL_LOSSES_THERMAL_EL 1

// in K, T[eV] = 8.618e-5 * T[K]
#define THERMAL_EL_TEMPERATURE 100.  // in K 

#define IONIZATION_FRACTION_THERMAL_EL 1.e-6

// collisions of H2 with H2, He
#define H2_COLLISIONS_WITH_H2_HE 1
#define THERMAL_NEUTRAL_TEMPERATURE 15.  // in K