#include "solvertest.h"

#include "mtx_format.h"
#include "csr_simple_format.h"

using namespace std;

void SolverTest::read_matrix_file(const std::string& file_name,
	const MatrixFileFormat file_format,
	std::vector<MKL_INT>& ia, std::vector<MKL_INT>& ja,
	std::vector<double>& A, MKL_INT& matrixSize, MKL_INT& nnz)
{
	// already loaded
	if (matrixSize > 0) return;

	switch (file_format) {
	case MatrixFileFormat::CSR: {
		ifstream A_file(file_name, ifstream::in | ifstream::binary);
		assert(("Can't open matrix file.", A_file.is_open()));
		TESTIO::read_csr_matrix_file(A_file, matrixSize, nnz, ia, ja, A);
		break;
	}
	case MatrixFileFormat::MTX: {
		MATRIXIO::read_matrix_file_mtx(file_name, ia, ja, A, matrixSize, nnz, false);
		MATRIXIO::convert_mtx_to_csr(matrixSize, ia);
		break;
	}
	default: {
		throw std::runtime_error("Error: Unsupported MatrixFileFormat in "
			"SolverTest::read_matrix_file");
	}
	}
	cout << "Matrix size = " << matrixSize << "\nNon_zero = " << nnz << "\n";
}

int LAESTest::test()
{
	read_matrix_file(A_name, file_format, ia, ja, a, sizea, nnza);

	int err = testSpecific();

	const int time = cr::duration_cast<cr::milliseconds>(end - start).count();
	cout << "Time == " << time << " ms" << endl;
	cout << "Error code == " << err << std::endl;

	return err;
}