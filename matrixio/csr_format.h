#ifndef CSR_FORMAT_H
#define CSR_FORMAT_H
#include "vec_format.h"

namespace MATRIXIO
{
template<class T>
void read_matrix_file_csr(std::string& file_name,
  std::vector<MKL_INT>& ia, std::vector<MKL_INT>& ja, std::vector<T>& A,
  MKL_INT& matrix_size, MKL_INT& non_zero,
  bool binary, bool matInt32 = false)
{
  std::ifstream file(file_name.c_str(), std::ifstream::in | std::ifstream::binary);
  if (!file.is_open())
  {
    ia.resize(1);
    ia[0] = 0;
    ja.resize(0);
    A.resize(0);
    std::cerr << "   ERROR: Cannot open file " << file_name << ".\n";
    return;
  }

  //if (matInt32)				
  int temp = static_cast<int>(matrix_size);
  if (binary)
  {
    if (matInt32) {
      file.read(reinterpret_cast<char*>(&temp), sizeof(temp));
      matrix_size = temp;
      file.read(reinterpret_cast<char*>(&temp), sizeof(temp));
      non_zero = temp;
    }
    else {
      file.read(reinterpret_cast<char*>(&matrix_size), sizeof(matrix_size));
      file.read(reinterpret_cast<char*>(&non_zero), sizeof(non_zero));
    }
  }
  else
  {
    file >> matrix_size >> non_zero;
  }
  if (matInt32) {
    std::vector<int> ia_temp;
    std::vector<int> ja_temp;
    read_vector(file, binary, matrix_size + 1, ia_temp, matInt32);
    read_vector(file, binary, non_zero, ja_temp, matInt32);
    read_vector(file, binary, non_zero, A, matInt32);
    ia.resize(matrix_size + 1);
    ja.resize(non_zero);
    for (MKL_INT i = 0; i < matrix_size + 1; ++i)
      ia[i] = ia_temp[i];
    for (MKL_INT i = 0; i < non_zero; ++i)
      ja[i] = ja_temp[i];
  }
  else {
    read_vector(file, binary, matrix_size + 1, ia, matInt32);
    read_vector(file, binary, non_zero, ja, matInt32);
    read_vector(file, binary, non_zero, A, matInt32);
  }

}

template<class T>
void save_crs_matrix_file(std::ofstream& file, bool binary,
  const std::vector<MKL_INT>& ia,
  const std::vector<MKL_INT>& ja,
  const std::vector<T>& A)
{
  const MKL_INT matrix_size = (MKL_INT)ia.size() - 1;
  const MKL_INT non_zero = (MKL_INT)ja.size();
  if (binary)
  {
    file.write(reinterpret_cast<const char*>(&matrix_size), sizeof(matrix_size));
    file.write(reinterpret_cast<const char*>(&non_zero), sizeof(non_zero));
  }
  else
  {
    file << matrix_size << " " << non_zero << "\n";
  }
  save_vector(file, binary, ia);
  save_vector(file, binary, ja);
  save_vector(file, binary, A);
  file.flush();
}

}

#endif //CSR_FORMAT_H
