
cmake_minimum_required(VERSION 3.1)

# Specify project name
PROJECT(gamma_ray_burst CXX)

# add_compile_options(-std=c11)
# Optimization
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Set names of source files
SET(source_files grb_propagation_eqs.cpp spectroscopy.cpp integration.cpp interpolation.cpp grb_microphysics.cpp grb_saving_data.cpp grb_spectra.cpp h2_uv_pump_data.cpp spectral_line_profile.cpp grb_abundances.cpp dynamic_array.cpp special_functions.cpp)

SET(SUNDIALS_INC_DIR /disk4/nester/cvode-5.7.0-install/include)

# Add path to SUNDIALS header files
INCLUDE_DIRECTORIES(${SUNDIALS_INC_DIR})

# Set search path for SUNDIALS libraries 
SET(SUNDIALS_LIB_DIR /disk4/nester/cvode-5.7.0-install/lib)

# Find the libraries
FIND_LIBRARY(SUNDIALS_SOLVER_LIB sundials_cvode ${SUNDIALS_LIB_DIR})
FIND_LIBRARY(SUNDIALS_NVEC_LIB sundials_nvecserial ${SUNDIALS_LIB_DIR})
FIND_LIBRARY(SUNDIALS_SUNMATRIXDENSE_LIB sundials_sunmatrixdense ${SUNDIALS_LIB_DIR})

#message(${SUNDIALS_SOLVER_LIB})

SET(SUNDIALS_EXTRA_LIB  -lm /usr/lib/x86_64-linux-gnu/librt.so)
SET(SUNDIALS_LIBS ${SUNDIALS_SOLVER_LIB} ${SUNDIALS_NVEC_LIB} ${SUNDIALS_SUNMATRIXDENSE_LIB} ${SUNDIALS_EXTRA_LIB})

find_package(OpenMP)
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

ADD_EXECUTABLE(gamma_ray_burst gamma_ray_burst.cpp ${source_files})
TARGET_LINK_LIBRARIES(gamma_ray_burst ${SUNDIALS_LIBS})

