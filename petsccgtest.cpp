#include "solvertest.h"
#include "petsccgtest.h"

#include <petscksp.h>
#include <csignal>

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
	//    Now time is measured in task...
	//cout << "Scatter time == " << scatter_time << " ms" << endl; // ~~~
	//cout << "Solve time   == " << time << " ms" << endl;
	cout << "Error code   == " << err << std::endl;
}

int localSize(int global_size, int mpi_rank, int mpi_size)
{
	return global_size / mpi_size + ((mpi_rank < (global_size % mpi_size)) ? 1 : 0);
}

void KSPSolveTask::syncMatrixSize()
{
	MPI_Bcast(&matrixSize, 1, MPI_INT, 0, communicator);
}

void KSPSolveTask::printContext() const
{
	std::cout << " Rank [" << MPI_rank << "]: Matrix size: " << matrixSize << std::endl;
	std::cout << " Rank [" << MPI_rank << "]: locNum: " << loc_num << std::endl;
	std::cout << " Rank [" << MPI_rank << "]: locRows: " << loc_rows << std::endl;

	std::cout << " Rank [" << MPI_rank << "]: locIASizes: ";
	for (int i = 0; i < MPI_size; ++i) std::cout << loc_ia_sizes[i] << " ";
	std::cout << std::endl;

	std::cout << " Rank [" << MPI_rank << "]: locJASizes: ";
	for (int i = 0; i < MPI_size; ++i) std::cout << loc_ja_sizes[i] << " ";
	std::cout << std::endl;

	std::cout << " Rank [" << MPI_rank << "]: locJA: ";
	for (int i = 0; i < loc_num; ++i) std::cout << loc_ja[i] << " ";
	std::cout << std::endl;

	std::cout << " Rank [" << MPI_rank << "]: locIA: ";
	for (int i = 0; i < loc_rows; ++i) std::cout << loc_ia[i] << " ";
	std::cout << std::endl;
}

void KSPSolveTask::task()
{
	PetscErrorCode ierr = _task();
	if (ierr != 0) throw ProcessTaskException("Error occured in _task().", ierr);
}

int KSPSolveTask::_task()
{

	auto mpiController = ProcessController::getInstance();

	PetscErrorCode ierr;
	
	ierr = PetscInitialize(0, nullptr, (char*)0, (char*)0); CHKERRQ(ierr);
	ierr = PetscPopSignalHandler();
	// Known Issue: PetscPopSignalHandler drops signal handlers to SIG_DEF
  // setting up signal handlers again is needed
	ProcessControllerSignalHandler::pushSignalHandler(false);

	MPI_rank = mpiController->MPIRank();
	MPI_size = mpiController->MPISize();		     
	communicator = mpiController->getCommunicator();

	syncMatrixSize();
	
	ierr = MPI_Barrier(communicator);

	// Scattering Ia:
	loc_rows = localSize(matrixSize, MPI_rank, MPI_size);
	loc_ia.resize(loc_rows + 1);
	loc_ia[0] = 0;
	
	ierr = MPI_Barrier(communicator);
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

	ierr = MPI_Barrier(communicator);
	auto scatter_end = cr::system_clock::now();

	int time = cr::duration_cast<cr::milliseconds>(scatter_end - scatter_start).count();
	if (!MPI_rank) std::cout << " Scatter elapsed: ~" << time << " ms" << std::endl;

	// Converting data
	ierr = PetscMalloc(loc_ia.size() * sizeof(PetscInt), &PETSc_loc_ia);	CHKERRQ(ierr);
	ierr = PetscMalloc(loc_ja.size() * sizeof(PetscInt), &PETSc_loc_ja);	CHKERRQ(ierr);
	ierr = PetscMalloc(loc_a.size() * sizeof(PetscScalar), &PETSc_loc_a);	CHKERRQ(ierr);

	for (size_t i = 0; i < loc_ia.size(); ++i) PETSc_loc_ia[i] = loc_ia[i];
	for (size_t i = 0; i < loc_ja.size(); ++i) PETSc_loc_ja[i] = loc_ja[i];
	for (size_t i = 0; i < loc_a.size(); ++i) PETSc_loc_a[i] = loc_a[i];

	ierr = MatCreateMPIAIJWithArrays(communicator, loc_rows, PETSC_DECIDE,
		PETSC_DETERMINE, matrixSize, PETSc_loc_ia, PETSc_loc_ja,
		PETSc_loc_a, &A);			CHKERRQ(ierr);

	ierr = PetscFree(PETSc_loc_a);	CHKERRQ(ierr);
	ierr = PetscFree(PETSc_loc_ia); CHKERRQ(ierr);
	ierr = PetscFree(PETSc_loc_ja); CHKERRQ(ierr);

	ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);		CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);		CHKERRQ(ierr);

	ierr = MatCreateVecs(A, &b, NULL);					CHKERRQ(ierr);
	ierr = VecDuplicate(b, &result);					CHKERRQ(ierr);
	ierr = VecDuplicate(b, &ref_result);				CHKERRQ(ierr);
	ierr = VecSet(result, 0.0);							CHKERRQ(ierr);
	ierr = VecSet(ref_result, 1.0);						CHKERRQ(ierr);
	ierr = MatMult(A, ref_result, b);					CHKERRQ(ierr);

	ierr = KSPCreate(communicator, &ksp);				CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &pc);							CHKERRQ(ierr);
	ierr = PCSetType(pc, PCJACOBI);						CHKERRQ(ierr);
	ierr = KSPMonitorSet(ksp, MyKSPMonitor, NULL, 0);   CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp, A, A);					CHKERRQ(ierr);
	ierr = KSPSetType(ksp, KSPCG);						CHKERRQ(ierr);
//	KSPSetType(ksp, KSPGMRES);
	ierr = KSPSetTolerances(ksp, 1e-8, PETSC_DEFAULT,
		PETSC_DEFAULT, 1'000);							CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);						CHKERRQ(ierr);

	// These calls are optional enable more precise profiling, 
	// since both will be called within KSPSolve() if they 
	// haven't been called already.
	ierr = KSPSetUp(ksp); CHKERRQ(ierr);

	auto solve_start = cr::system_clock::now();
	ierr = KSPSolve(ksp, b, result); CHKERRQ(ierr);
	auto solve_end = cr::system_clock::now();
	time = cr::duration_cast<cr::milliseconds>(solve_end - solve_start).count();

	if (!MPI_rank) std::cout << " Solve elapsed: ~" << time << " ms" << std::endl;

	// Convergence check
	KSPConvergedReason reason;
	ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);
	if (reason < 0) {
		std::cout << "Divergence detected: " << reason << std::endl;
	}
	else
	{
		PetscInt its;
		ierr = KSPGetIterationNumber(ksp, &its);	CHKERRQ(ierr);
		PetscScalar res_norm_;
		ierr = KSPGetResidualNorm(ksp, &res_norm_); CHKERRQ(ierr);
		iterations = its;
		res_norm = res_norm_;
	}

	ierr = KSPDestroy(&ksp);		CHKERRQ(ierr);
	ierr = MatDestroy(&A);			CHKERRQ(ierr);
	ierr = VecDestroy(&ref_result); CHKERRQ(ierr);
	ierr = VecDestroy(&result);		CHKERRQ(ierr);
	ierr = VecDestroy(&b);			CHKERRQ(ierr);
	
	ierr = PetscFinalize();			CHKERRQ(ierr);
	return ierr;
}

int PETScCGTest::testSpecific()
{
	auto pc = ProcessController::getInstance();

	auto solveTask = dynamic_pointer_cast<KSPSolveTask>(pc->getTask(Task::KSPSolve));
	
	// Configuring settings
	solveTask->a = a;
	solveTask->ia = ia;
	solveTask->ja = ja;
	solveTask->matrixSize = sizea;

	solveTask->loc_ia_sizes.resize(pc->MPISize());
	solveTask->loc_ja_sizes.resize(pc->MPISize());
	solveTask->loc_ia_starts.resize(pc->MPISize());
	solveTask->loc_ja_starts.resize(pc->MPISize());

	// Calculating IA markup
	solveTask->loc_ia_sizes[0] = localSize(sizea, 0, pc->MPISize());
	solveTask->loc_ia_starts[0] = 1;

	for (int MPIRank = 1; MPIRank < pc->MPISize(); ++MPIRank)
	{
		solveTask->loc_ia_sizes[MPIRank] = localSize(sizea, MPIRank, pc->MPISize());
		solveTask->loc_ia_starts[MPIRank] = solveTask->loc_ia_starts[MPIRank - 1] +solveTask->loc_ia_sizes[MPIRank - 1];
	}
	for (int MPIRank = pc->MPISize() - 1; MPIRank >= 0; MPIRank--)
	{
		for (int i = 0; i < solveTask->loc_ia_sizes[MPIRank]; ++i)
		{
			solveTask->ia[solveTask->loc_ia_starts[MPIRank] + i] -= solveTask->ia[solveTask->loc_ia_starts[MPIRank] - 1];
		}
	}

	// Calculating JA markup
	int sum = 0;
	for (int MPIRank = 0; MPIRank < pc->MPISize(); ++MPIRank)
	{
		// Number of elements is in the last local ia
		int last_loc_ia_idx = solveTask->loc_ia_starts[MPIRank] + solveTask->loc_ia_sizes[MPIRank] - 1;
		solveTask->loc_ja_sizes[MPIRank] = solveTask->ia[last_loc_ia_idx];
		solveTask->loc_ja_starts[MPIRank] = sum;
		sum +=solveTask->loc_ja_sizes[MPIRank];
	}

	// Setting and Invoking solver:

	pc->evaluateTask(Task::KSPSolve);
	

	std::cout << " Iterations taken: " << solveTask->iterations << std::endl;
	std::cout << " Residual norm: " << solveTask->res_norm << std::endl;

	return 0;
}
