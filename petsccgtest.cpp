#include "solvertest.h"

#include <process.h>
#include <petscksp.h>

using namespace std;

PetscErrorCode MyKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void* dummy)
{
	if (!(n % 50)) PetscPrintf(PETSC_COMM_WORLD,
		"Iteration %D KSP Residual norm %14.12e \n", n, rnorm);
	return 0;
}


int PETScCGTest::test()
{
	read_matrix_file(string(A_name), file_format, ia, ja, a, sizea, nnza);
	if (sizea == 0) return -1;

	int err = testSpecific();

	const int scatter_time = cr::duration_cast<cr::milliseconds>(scatter_end - scatter_start).count();
	const int time = cr::duration_cast<cr::milliseconds>(end - start).count();
	//cout << "Scatter time == " << scatter_time << " ms" << endl; // ~~~
	//cout << "Solve time   == " << time << " ms" << endl;
	cout << "Error code   == " << err << std::endl;
}

int localSize(int global_size, int mpi_rank, int mpi_size)
{
	return global_size / mpi_size + ((mpi_rank < (global_size % mpi_size)) ? 1 : 0);
}

void KSPSolveTask::syncMSize()
{
	MPI_Bcast(&matrix_size, 1, MPI_INT, 0, communicator);
}

void KSPSolveTask::prepareMarkup()
{
	loc_ia_sizes.resize(size);
	loc_ia_starts.resize(size);
	loc_ja_sizes.resize(size);
	loc_ja_starts.resize(size);
}

void KSPSolveTask::syncMarkup()
{
	MPI_Bcast(&(loc_ia_sizes.at(0)), size, MPI_INT, 0, communicator);
	MPI_Bcast(&(loc_ia_starts.at(0)), size, MPI_INT, 0, communicator);
	MPI_Bcast(&(loc_ja_sizes.at(0)), size, MPI_INT, 0, communicator);
	MPI_Bcast(&(loc_ja_starts.at(0)), size, MPI_INT, 0, communicator);
}

void KSPSolveTask::printContext() const
{
	std::cout << " Rank [" << rank << "]: Matrix size: " << matrix_size << std::endl;
	std::cout << " Rank [" << rank << "]: loc_num: " << loc_num << std::endl;
	std::cout << " Rank [" << rank << "]: loc_rows: " << loc_rows << std::endl;

	std::cout << " Rank [" << rank << "]: loc_ia_sizes: ";
	for (int i = 0; i < size; ++i) std::cout << loc_ia_sizes[i] << " ";
	std::cout << std::endl;

	std::cout << " Rank [" << rank << "]: loc_ja_sizes: ";
	for (int i = 0; i < size; ++i) std::cout << loc_ja_sizes[i] << " ";
	std::cout << std::endl;

	std::cout << " Rank [" << rank << "]: loc_ja: ";
	for (int i = 0; i < loc_num; ++i) std::cout << loc_ja[i] << " ";
	std::cout << std::endl;

	std::cout << " Rank [" << rank << "]: loc_ia: ";
	for (int i = 0; i < loc_rows; ++i) std::cout << loc_ia[i] << " ";
	std::cout << std::endl;
}

void KSPSolveTask::task()
{
	MPI_Barrier(communicator);
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &size);

	// Synchronizing matrix size
	syncMSize();

	// Resizing markup vectors
	if (rank) prepareMarkup();

	// Synchronizing markup
	syncMarkup();
	
	MPI_Barrier(communicator);

	// Scattering Ia:
	loc_rows = localSize(matrix_size, rank, size);
	loc_ia.resize(loc_rows + 1);
	loc_ia[0] = 0;
	
	MPI_Barrier(communicator);
	auto scatter_start = cr::system_clock::now();

	MPI_Scatterv(ia.data(), loc_ia_sizes.data(), loc_ia_starts.data(), MPI_INT,
		loc_ia.data() + 1, loc_rows, MPI_INT, 0, communicator);

	// Scattering Ja:
	loc_num = loc_ia[loc_ia.size() - 1];
	loc_ja.resize(loc_num);
	MPI_Scatterv(ja.data(), loc_ja_sizes.data(), loc_ja_starts.data(), MPI_INT,
		loc_ja.data(), loc_num, MPI_INT, 0, communicator);
	
	// Scattering A:
	loc_a.resize(loc_num);
	MPI_Scatterv(a.data(), loc_ja_sizes.data(), loc_ja_starts.data(), MPI_DOUBLE,
		loc_a.data(), loc_num, MPI_DOUBLE, 0, communicator);

	MPI_Barrier(communicator);
	auto scatter_end = cr::system_clock::now();

	int time = cr::duration_cast<cr::milliseconds>(scatter_end - scatter_start).count();
	if (!rank) std::cout << " Scatter elapsed: ~" << time << " ms" << std::endl;

	// Converting data
	PetscMalloc(loc_ia.size() * sizeof(PetscInt), &petsc_loc_ia);
	PetscMalloc(loc_ja.size() * sizeof(PetscInt), &petsc_loc_ja);
	PetscMalloc(loc_a.size() * sizeof(PetscScalar), &petsc_loc_a);
	for (size_t i = 0; i < loc_ia.size(); ++i) petsc_loc_ia[i] = loc_ia[i];
	for (size_t i = 0; i < loc_ja.size(); ++i) petsc_loc_ja[i] = loc_ja[i];
	for (size_t i = 0; i < loc_a.size(); ++i) petsc_loc_a[i] = loc_a[i];

	MatCreateMPIAIJWithArrays(communicator, loc_rows, PETSC_DECIDE,
		PETSC_DETERMINE, matrix_size, petsc_loc_ia, petsc_loc_ja,
		petsc_loc_a, &A);

	PetscFree(petsc_loc_a);
	PetscFree(petsc_loc_ia);
	PetscFree(petsc_loc_ja);

	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	MatCreateVecs(A, &b, NULL);
	VecDuplicate(b, &result);
	VecDuplicate(b, &ref_result);
	VecSet(result, 0.0);
	VecSet(ref_result, 1.0);
	MatMult(A, ref_result, b);

	KSPCreate(communicator, &ksp);
	KSPGetPC(ksp, &pc);
	PCSetType(pc, PCJACOBI);
	KSPMonitorSet(ksp, MyKSPMonitor, NULL, 0);
	KSPSetOperators(ksp, A, A);
	//KSPSetType(ksp, KSPCG);
	KSPSetType(ksp, KSPGMRES);
	KSPSetTolerances(ksp, 1e-8, PETSC_DEFAULT,
		PETSC_DEFAULT, 1'000);
	KSPSetFromOptions(ksp);

	// These calls are optional enable more precise profiling, 
	// since both will be called within KSPSolve() if they 
	// haven't been called already.
	KSPSetUp(ksp);

	auto solve_start = cr::system_clock::now();
	KSPSolve(ksp, b, result);
	auto solve_end = cr::system_clock::now();
	time = cr::duration_cast<cr::milliseconds>(solve_end - solve_start).count();

	if (!rank) std::cout << " Solve elapsed: ~" << time << " ms" << std::endl;

	// Convergence check
	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp, &reason);
	if (reason < 0) {
		std::cout << "Divergence detected: " << reason << std::endl;
	}
	else
	{
		PetscInt its;
		KSPGetIterationNumber(ksp, &its);
		PetscScalar res_norm_;
		KSPGetResidualNorm(ksp, &res_norm_);
		iterations = its;
		res_norm = res_norm_;
	}

	KSPDestroy(&ksp);
	MatDestroy(&A);
	VecDestroy(&ref_result);
	VecDestroy(&result);
	VecDestroy(&b);
}

int PETScCGTest::testSpecific()
{
	ProcessController* pc = ProcessController::getInstance();

	KSPSolveTask* solveTask = dynamic_cast<KSPSolveTask*>(pc->getTask(Task::KSPSolve));
	
	// Configuring settings
	solveTask->a = a;
	solveTask->ia = ia;
	solveTask->ja = ja;
	solveTask->matrix_size = sizea;

	solveTask->loc_ia_sizes.resize(pc->size());
	solveTask->loc_ja_sizes.resize(pc->size());
	solveTask->loc_ia_starts.resize(pc->size());
	solveTask->loc_ja_starts.resize(pc->size());

	// Calculating IA markup
	solveTask->loc_ia_sizes[0] = localSize(sizea, 0, pc->size());
	solveTask->loc_ia_starts[0] = 1;

	for (int rank = 1; rank < pc->size(); ++rank)
	{
		solveTask->loc_ia_sizes[rank] = localSize(sizea, rank, pc->size());
		solveTask->loc_ia_starts[rank] = solveTask->loc_ia_starts[rank - 1] +solveTask->loc_ia_sizes[rank - 1];
	}
	for (int rank = pc->size() - 1; rank >= 0; rank--)
	{
		for (int i = 0; i < solveTask->loc_ia_sizes[rank]; ++i)
		{
			solveTask->ia[solveTask->loc_ia_starts[rank] + i] -= solveTask->ia[solveTask->loc_ia_starts[rank] - 1];
		}
	}

	// Calculating JA markup
	int sum = 0;
	for (int rank = 0; rank < pc->size(); ++rank)
	{
		// Number of elements is in the last local ia
		int last_loc_ia_idx = solveTask->loc_ia_starts[rank] + solveTask->loc_ia_sizes[rank] - 1;
		solveTask->loc_ja_sizes[rank] = solveTask->ia[last_loc_ia_idx];
		solveTask->loc_ja_starts[rank] = sum;
		sum +=solveTask->loc_ja_sizes[rank];
	}

	// Setting and Invoking solver:
	start = cr::system_clock::now();
	pc->evaluateTask(Task::KSPSolve);
	end = cr::system_clock::now();

	std::cout << " Iterations taken: " << solveTask->iterations << std::endl;
	std::cout << " Residual norm: " << solveTask->res_norm << std::endl;

	pc->evaluateTask(Task::Shutdown);

	return 0;
}