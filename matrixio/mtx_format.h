#ifndef MM_FORMAT_H
#define MM_FORMAT_H

#include <stdio.h>
#include <string>
#include <vector>
#include <mkl.h>

namespace MATRIXIO
{

#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64

typedef char MM_typecode[4];

/********************* MM_typecode query functions ***************************/

#define mm_is_matrix(typecode)       ((typecode)[0]=='M')

#define mm_is_sparse(typecode)       ((typecode)[1]=='C')
#define mm_is_coordinate(typecode)   ((typecode)[1]=='C')
#define mm_is_dense(typecode)        ((typecode)[1]=='A')
#define mm_is_array(typecode)        ((typecode)[1]=='A')

#define mm_is_complex(typecode)      ((typecode)[2]=='C')
#define mm_is_real(typecode)         ((typecode)[2]=='R')
#define mm_is_pattern(typecode)	     ((typecode)[2]=='P')
#define mm_is_integer(typecode)      ((typecode)[2]=='I')

#define mm_is_symmetric(typecode)    ((typecode)[3]=='S')
#define mm_is_general(typecode)      ((typecode)[3]=='G')
#define mm_is_skew(typecode)         ((typecode)[3]=='K')
#define mm_is_hermitian(typecode)    ((typecode)[3]=='H')


/********************* MM_typecode modify functions ***************************/

#define mm_set_matrix(typecode)      ((*typecode)[0]='M')
#define mm_set_coordinate(typecode)  ((*typecode)[1]='C')
#define mm_set_array(typecode)       ((*typecode)[1]='A')
#define mm_set_dense(typecode)	     mm_set_array(typecode)
#define mm_set_sparse(typecode)	     mm_set_coordinate(typecode)

#define mm_set_complex(typecode)     ((*typecode)[2]='C')
#define mm_set_real(typecode)        ((*typecode)[2]='R')
#define mm_set_pattern(typecode)     ((*typecode)[2]='P')
#define mm_set_integer(typecode)     ((*typecode)[2]='I')


#define mm_set_symmetric(typecode)   ((*typecode)[3]='S')
#define mm_set_general(typecode)     ((*typecode)[3]='G')
#define mm_set_skew(typecode)        ((*typecode)[3]='K')
#define mm_set_hermitian(typecode)   ((*typecode)[3]='H')

#define mm_clear_typecode(typecode)  ((*typecode)[0]=(*typecode)[1]= (*typecode)[2]=' ',(*typecode)[3]='G')

#define mm_initialize_typecode(typecode) mm_clear_typecode(typecode)


/********************* Matrix Market error codes ***************************/


#define MM_COULD_NOT_READ_FILE       11
#define MM_PREMATURE_EOF             12
#define MM_NOT_MTX                   13
#define MM_NO_HEADER                 14
#define MM_UNSUPPORTED_TYPE          15
#define MM_LINE_TOO_LONG             16
#define MM_COULD_NOT_WRITE_FILE      17


/******************** Matrix Market internal definitions ********************

   MM_matrix_typecode: 4-character sequence

            ojbect 		sparse/   	data        storage
                  dense     	type        scheme

   string position:	 [0]        [1]			[2]         [3]

   Matrix typecode:  M(atrix)  C(oord)		R(eal)   	G(eneral)
                    A(array)	C(omplex)   H(ermitian)
                      P(attern)   S(ymmetric)
                        I(nteger)	K(kew)

 ***********************************************************************/

#define MM_MTX_STR        "matrix"
#define MM_ARRAY_STR      "array"
#define MM_DENSE_STR      "array"
#define MM_COORDINATE_STR "coordinate" 
#define MM_SPARSE_STR     "coordinate"
#define MM_COMPLEX_STR    "complex"
#define MM_REAL_STR       "real"
#define MM_INT_STR        "integer"
#define MM_GENERAL_STR    "general"
#define MM_SYMM_STR       "symmetric"
#define MM_HERM_STR       "hermitian"
#define MM_SKEW_STR       "skew-symmetric"
#define MM_PATTERN_STR    "pattern"

MKL_INT mm_read_banner(FILE* f, MM_typecode* matcode);

MKL_INT mm_read_mtx_crd_size(FILE* f, MKL_INT& M, MKL_INT& nz);

void convert_mtx_to_csr(const MKL_INT& matrix_size, std::vector<MKL_INT>& ia);

template<class T>
void read_matrix_file_mtx(const std::string& file_name,
  std::vector<MKL_INT>& ia, std::vector<MKL_INT>& ja, std::vector<T>& A,
  MKL_INT& matrix_size, MKL_INT& non_zero,
  bool binary, bool matInt32 = false)
{
  if (binary) {
    std::cerr << "   ERROR: MTX format couldn't be binary\n";
    return;
  }

  FILE* f = fopen(file_name.c_str(), "r");
  MM_typecode matcode;

  if (f == NULL)
  {
    std::cerr << "   ERROR: Cannot open or process matrix market file " << file_name << "\n";
    return;
  }

  mm_read_banner(f, &matcode);
  if (!mm_is_matrix(matcode) || !(mm_is_real(matcode) || mm_is_integer(matcode) || mm_is_pattern(matcode)))
  {
    std::cerr << "   ERROR: Matrix market file has bad typecode " << matcode << "\n";
    return;
  }

  mm_read_mtx_crd_size(f, matrix_size, non_zero);

  ia.resize(non_zero);
  ja.resize(non_zero);
  A.resize(non_zero);

  MKL_INT i, j;
  double val;

  for (MKL_INT el = 0; el < non_zero; el++) {
    if (mm_is_pattern(matcode)) {
      fscanf(f, "%d %d\n", &i, &j);
      val = 1.;
    }
    else
      fscanf(f, "%d %d %lg\n", &i, &j, &val);

    ia[el] = i - 1;
    ja[el] = j - 1;
    A[el] = val;

  }
}
}

#endif
