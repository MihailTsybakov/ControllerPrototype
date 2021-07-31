#ifndef MATRIX_SOLVERS_TESTS
#define MATRIX_SOLVERS_TESTS

#include <mkl.h>
#include <fstream>
#include <vector>

using namespace std;

namespace TESTIO {
template<class T>
void read_vector(ifstream& file, MKL_INT size, vector<T>& result)
{
  result.resize(size);
  for (MKL_INT i = 0; i < size; ++i)
    file >> result[i];
}

template<class T>
void read_csr_matrix_file(ifstream& file, MKL_INT& matrix_size, MKL_INT& non_zero,
  vector<MKL_INT>& ia, vector<MKL_INT>& ja, vector<T>& A)
{
  file >> matrix_size >> non_zero;
  read_vector(file, matrix_size + 1, ia);
  read_vector(file, non_zero, ja);
  read_vector(file, non_zero, A);
}
}

#endif
