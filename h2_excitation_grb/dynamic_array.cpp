
#include "dynamic_array.h"
#include <cstring>
using namespace std;

dynamic_array::dynamic_array(int d)
{
	dim = d;
	arr = new double[dim];
	memset(arr, 0, dim * sizeof(double));
}

dynamic_array::dynamic_array(int d, double* arr2)
{
	dim = d;
	arr = new double[dim];
	memcpy(arr, arr2, dim * sizeof(double));
}

dynamic_array::~dynamic_array() {
	delete[] arr;
}

dynamic_array::dynamic_array(const dynamic_array& obj)
{
	dim = obj.dim;
	arr = new double[dim];
	memcpy(arr, obj.arr, dim * sizeof(double));
}

dynamic_array& dynamic_array::operator = (const dynamic_array& obj)
{
	dim = obj.dim;
	delete[] arr;
	arr = new double[dim];

	memcpy(arr, obj.arr, dim * sizeof(double));
	return *this;
}
