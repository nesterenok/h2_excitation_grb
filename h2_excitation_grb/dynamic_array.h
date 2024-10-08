#pragma once
// the dynamic array for floating point values;
struct dynamic_array
{
	int dim;
	double* arr;

	dynamic_array& operator = (const dynamic_array& obj);
	dynamic_array(int);
	dynamic_array(int, double*);
	~dynamic_array();

	dynamic_array(const dynamic_array& d_a);
};
