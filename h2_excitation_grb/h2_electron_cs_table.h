#pragma once
#include <string>

// Saving cross section table
// ji = 0 or 1; initial H2 energy level,
void save_cross_section_table(const std::string& output_path, const std::string& data_path,
	double ionization_fraction, double thermal_electron_temp, int ji, const std::string & str);
