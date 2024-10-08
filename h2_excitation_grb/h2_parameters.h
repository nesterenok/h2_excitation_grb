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
#define ABS_POPULATION_H2_ERROR_SOLVER 1.e-14  // specimen level population density (cm-3)
#define ABS_ELSPECTRA_ERROR_SOLVER 1.e-14
#define ABS_PARAMETER_ERROR_SOLVER 1.e-11  // e.g., grain radius (cm), temperature (K)
#define MAX_CONV_FAILS_SOLVER 100  // default value is 10;
#define MAX_ERR_TEST_FAILS_SOLVER 14  // default value is 7;

#define NB_OF_BINS_PER_ORDER_EL 100


// ----------------------------------------------------
// Electron spectra evolution
//-----------------------------------------------------
// Weingartner & Draine, ApJSS 134, 263 (2001)
#define ELECTRON_DUST_SCATTERING_PROB 0.5
#define CALC_COLOUMB_EL_LOSSES 0
