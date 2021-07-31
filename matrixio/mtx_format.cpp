#include "mtx_format.h"

namespace MATRIXIO
{

MKL_INT mm_is_valid(MM_typecode matcode)
{
  if (!mm_is_matrix(matcode)) return 0;
  if (mm_is_dense(matcode) && mm_is_pattern(matcode)) return 0;
  if (mm_is_real(matcode) && mm_is_hermitian(matcode)) return 0;
  if (mm_is_pattern(matcode) && (mm_is_hermitian(matcode) ||
    mm_is_skew(matcode))) return 0;
  return 1;
}

MKL_INT mm_read_banner(FILE* f, MM_typecode* matcode)
{
  char line[MM_MAX_LINE_LENGTH];
  char banner[MM_MAX_TOKEN_LENGTH];
  char mtx[MM_MAX_TOKEN_LENGTH];
  char crd[MM_MAX_TOKEN_LENGTH];
  char data_type[MM_MAX_TOKEN_LENGTH];
  char storage_scheme[MM_MAX_TOKEN_LENGTH];
  char* p;


  mm_clear_typecode(matcode);

  if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
    return MM_PREMATURE_EOF;

  if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type,
    storage_scheme) != 5)
    return MM_PREMATURE_EOF;

  for (p = mtx; *p != '\0'; *p = tolower(*p), p++);  /* convert to lower case */
  for (p = crd; *p != '\0'; *p = tolower(*p), p++);
  for (p = data_type; *p != '\0'; *p = tolower(*p), p++);
  for (p = storage_scheme; *p != '\0'; *p = tolower(*p), p++);

  /* check for banner */
  if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
    return MM_NO_HEADER;

  /* first field should be "mtx" */
  if (strcmp(mtx, MM_MTX_STR) != 0)
    return  MM_UNSUPPORTED_TYPE;
  mm_set_matrix(matcode);


  /* second field describes whether this is a sparse matrix (in coordinate
          storgae) or a dense array */

  if (strcmp(crd, MM_SPARSE_STR) == 0)
    mm_set_sparse(matcode);
  else
    if (strcmp(crd, MM_DENSE_STR) == 0)
      mm_set_dense(matcode);
    else
      return MM_UNSUPPORTED_TYPE;


  /* third field */

  if (strcmp(data_type, MM_REAL_STR) == 0)
    mm_set_real(matcode);
  else
    if (strcmp(data_type, MM_COMPLEX_STR) == 0)
      mm_set_complex(matcode);
    else
      if (strcmp(data_type, MM_PATTERN_STR) == 0)
        mm_set_pattern(matcode);
      else
        if (strcmp(data_type, MM_INT_STR) == 0)
          mm_set_integer(matcode);
        else
          return MM_UNSUPPORTED_TYPE;


  /* fourth field */

  if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
    mm_set_general(matcode);
  else
    if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
      mm_set_symmetric(matcode);
    else
      if (strcmp(storage_scheme, MM_HERM_STR) == 0)
        mm_set_hermitian(matcode);
      else
        if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
          mm_set_skew(matcode);
        else
          return MM_UNSUPPORTED_TYPE;


  return 0;
}

MKL_INT mm_read_mtx_crd_size(FILE* f, MKL_INT& M, MKL_INT& nz)
{
  char line[MM_MAX_LINE_LENGTH];
  MKL_INT num_items_read;

  /* set return null parameter values, in case we exit with errors */
  M = nz = 0;

  /* now continue scanning until you reach the end-of-comments */
  do
  {
    if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
      return MM_PREMATURE_EOF;
  } while (line[0] == '%');

  MKL_INT tmp;
  /* line[] is either blank or has M,N, nz */
  if (sscanf(line, "%d %d %d", &tmp, &M, &nz) == 3)
    return 0;

  else
    do
    {
      num_items_read = fscanf(f, "%d %d %d", tmp, &M, &nz);
      if (num_items_read == EOF) return MM_PREMATURE_EOF;
    } while (num_items_read != 3);

    return 0;
}

// https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatCreateMPIAIJWithArrays.html
void convert_mtx_to_csr(const MKL_INT& matrix_size, std::vector<MKL_INT>& ia)
{
  std::vector<MKL_INT> mtxia(ia);
  ia.resize(static_cast<size_t>(matrix_size) + 1);
  ia[0] = 0;

  size_t cur_row = 0;
  size_t num = 0;
  // new last element matrix_size==max_row+1 hence it'll not be added to
  // the new ia just will be used as terminator to run last if (row > cur_row)
  // and fulfill last rows with final num
  mtxia.push_back(matrix_size);
  for (const auto& row : mtxia) {
    if (row > cur_row)
    {
      ia[++cur_row] = ia[cur_row - 1] + num;
      // empty rows case
      for (auto i = cur_row; i < row; ++i) ia[i + 1] = ia[i];
      // row == cur_row => row isn't empty
      cur_row = row;
      num = 1;
    }
    else
    {
      ++num;
    }
  }
}

}
