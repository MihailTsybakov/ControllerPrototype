#include "solvertest.h"
#include "petsccgtest.h"

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
	std::cout << " Rank [" << MPIRank << "]: Matrix size: " << matrixSize << std::endl;
	std::cout << " Rank [" << MPIRank << "]: locNum: " << locNum << std::endl;
	std::cout << " Rank [" << MPIRank << "]: locRows: " << locRows << std::endl;

	std::cout << " Rank [" << MPIRank << "]: locIASizes: ";
	for (int i = 0; i < MPISize; ++i) std::cout << locIASizes[i] << " ";
	std::cout << std::endl;

	std::cout << " Rank [" << MPIRank << "]: locJASizes: ";
	for (int i = 0; i < MPISize; ++i) std::cout << locJASizes[i] << " ";
	std::cout << std::endl;

	std::cout << " Rank [" << MPIRank << "]: locJA: ";
	for (int i = 0; i < locNum; ++i) std::cout << locJA[i] << " ";
	std::cout << std::endl;

	std::cout << " Rank [" << MPIRank << "]: locIA: ";
	for (int i = 0; i < locRows; ++i) std::cout << locIA[i] << " ";
	std::cout << std::endl;
}

void KSPSolveTask::task()
{
	PetscErrorCode ierr;
	ierr = MPI_Barrier(communicator);
  
	ProcessController* mpiController = ProcessController::getInstance();
	MPIRank = mpiController->MPIRank();
	MPISize = mpiController->MPISize();
	communicator = mpiController->getCommunicator();

	syncMatrixSize();
	
	ierr = MPI_Barrier(communicator);

	// Scattering Ia:
	locRows = localSize(matrixSize, MPIRank, MPISize);
	locIA.resize(locRows + 1);
	locIA[0] = 0;
	
	ierr = MPI_Barrier(communicator);
	auto scatter_start = cr::system_clock::now();

	MPI_Scatterv(ia.data(), locIASizes.data(), locIAStarts.data(), MPI_INT,
		locIA.data() + 1, locRows, MPI_INT, 0, communicator);

	// Scattering Ja:
	locNum = locIA[locIA.size() - 1];
	locJA.resize(locNum);
	MPI_Scatterv(ja.data(), locJASizes.data(), locJAStarts.data(), MPI_INT,
		locJA.data(), locNum, MPI_INT, 0, communicator);
	
	// Scattering A:
	locA.resize(locNum);
	MPI_Scatterv(a.data(), locJASizes.data(), locJAStarts.data(), MPI_DOUBLE,
		locA.data(), locNum, MPI_DOUBLE, 0, communicator);

	ierr = MPI_Barrier(communicator);
	auto scatter_end = cr::system_clock::now();

	int time = cr::duration_cast<cr::milliseconds>(scatter_end - scatter_start).count();
	if (!MPIRank) std::cout << " Scatter elapsed: ~" << time << " ms" << std::endl;

	// Converting data
	ierr = PetscMalloc(locIA.size() * sizeof(PetscInt), &petscLocIA);
	ierr = PetscMalloc(locJA.size() * sizeof(PetscInt), &petscLocJA);
	ierr = PetscMalloc(locA.size() * sizeof(PetscScalar), &petscLocA);

	if (ierr != 0) throw ProcessTaskException("Error occured in PetscMalloc.", ierr);

	for (size_t i = 0; i < locIA.size(); ++i) petscLocIA[i] = locIA[i];
	for (size_t i = 0; i < locJA.size(); ++i) petscLocJA[i] = locJA[i];
	for (size_t i = 0; i < locA.size(); ++i) petscLocA[i] = locA[i];

	ierr = MatCreateMPIAIJWithArrays(communicator, locRows, PETSC_DECIDE,
		PETSC_DETERMINE, matrixSize, petscLocIA, petscLocJA,
		petscLocA, &A);

	if (ierr != 0) throw ProcessTaskException("Error occured in MatCreateMPIAIWithArrays.", ierr);

	ierr = PetscFree(petscLocA);
	ierr = PetscFree(petscLocIA);
	ierr = PetscFree(petscLocJA);

	if (ierr != 0) throw ProcessTaskException("Error occured while tried to release memory.", ierr);

	ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	if (ierr != 0) throw ProcessTaskException("Error: MatAssembly failed.", ierr);

	ierr = MatCreateVecs(A, &b, NULL);
	ierr = VecDuplicate(b, &result);
	ierr = VecDuplicate(b, &refResult);
	ierr = VecSet(result, 0.0);
	ierr = VecSet(refResult, 1.0);
	ierr = MatMult(A, refResult, b);

	if (ierr != 0) throw ProcessTaskException("Error occured while tried to form vectors.", ierr);

	ierr = KSPCreate(communicator, &ksp);
	ierr = KSPGetPC(ksp, &pc);
	ierr = PCSetType(pc, PCJACOBI);
	ierr = KSPMonitorSet(ksp, MyKSPMonitor, NULL, 0);
	ierr = KSPSetOperators(ksp, A, A);
	ierr = KSPSetType(ksp, KSPCG);
//	KSPSetType(ksp, KSPGMRES);
	ierr = KSPSetTolerances(ksp, 1e-8, PETSC_DEFAULT,
		PETSC_DEFAULT, 1'000);
	ierr = KSPSetFromOptions(ksp);

	// These calls are optional enable more precise profiling, 
	// since both will be called within KSPSolve() if they 
	// haven't been called already.
	ierr = KSPSetUp(ksp);

	if (ierr != 0) throw ProcessTaskException("Error: Solver setup failed.", ierr);

	auto solve_start = cr::system_clock::now();
	ierr = KSPSolve(ksp, b, result);
	auto solve_end = cr::system_clock::now();
	time = cr::duration_cast<cr::milliseconds>(solve_end - solve_start).count();

	if (ierr != 0) throw ProcessTaskException("Error: solving failed.", ierr);
	if (!MPIRank) std::cout << " Solve elapsed: ~" << time << " ms" << std::endl;

	// Convergence check
	KSPConvergedReason reason;
	ierr = KSPGetConvergedReason(ksp, &reason);
	if (reason < 0) {
		std::cout << "Divergence detected: " << reason << std::endl;
	}
	else
	{
		PetscInt its;
		ierr = KSPGetIterationNumber(ksp, &its);
		PetscScalar res_norm_;
		ierr = KSPGetResidualNorm(ksp, &res_norm_);
		iterations = its;
		resNorm = res_norm_;
	}

	if (ierr != 0) throw ProcessTaskException("Error: convergence check failed.", ierr);

	ierr = KSPDestroy(&ksp);
	ierr = MatDestroy(&A);
	ierr = VecDestroy(&refResult);
	ierr = VecDestroy(&result);
	ierr = VecDestroy(&b);

	if (ierr != 0) throw ProcessTaskException("Error: Releasing petsc memory failed.", ierr);
}

int PETScCGTest::testSpecific()
{
	ProcessController* pc = ProcessController::getInstance();

	KSPSolveTask* solveTask = dynamic_cast<KSPSolveTask*>(pc->getTask(Task::KSPSolve));
	
	// Configuring settings
	solveTask->a = a;
	solveTask->ia = ia;
	solveTask->ja = ja;
	solveTask->matrixSize = sizea;

	solveTask->locIASizes.resize(pc->MPISize());
	solveTask->locJASizes.resize(pc->MPISize());
	solveTask->locIAStarts.resize(pc->MPISize());
	solveTask->locJAStarts.resize(pc->MPISize());

	// Calculating IA markup
	solveTask->locIASizes[0] = localSize(sizea, 0, pc->MPISize());
	solveTask->locIAStarts[0] = 1;

	for (int MPIRank = 1; MPIRank < pc->MPISize(); ++MPIRank)
	{
		solveTask->locIASizes[MPIRank] = localSize(sizea, MPIRank, pc->MPISize());
		solveTask->locIAStarts[MPIRank] = solveTask->locIAStarts[MPIRank - 1] +solveTask->locIASizes[MPIRank - 1];
	}
	for (int MPIRank = pc->MPISize() - 1; MPIRank >= 0; MPIRank--)
	{
		for (int i = 0; i < solveTask->locIASizes[MPIRank]; ++i)
		{
			solveTask->ia[solveTask->locIAStarts[MPIRank] + i] -= solveTask->ia[solveTask->locIAStarts[MPIRank] - 1];
		}
	}

	// Calculating JA markup
	int sum = 0;
	for (int MPIRank = 0; MPIRank < pc->MPISize(); ++MPIRank)
	{
		// Number of elements is in the last local ia
		int last_loc_ia_idx = solveTask->locIAStarts[MPIRank] + solveTask->locIASizes[MPIRank] - 1;
		solveTask->locJASizes[MPIRank] = solveTask->ia[last_loc_ia_idx];
		solveTask->locJAStarts[MPIRank] = sum;
		sum +=solveTask->locJASizes[MPIRank];
	}

	// Setting and Invoking solver:
	try
	{
		pc->evaluateTask(Task::KSPSolve);
	}
	catch (const std::exception& exc)
	{
		std::cerr << " While tried to solve system, exception occured: \"" << exc.what() << "\"." << std::endl;
	}

	std::cout << " Iterations taken: " << solveTask->iterations << std::endl;
	std::cout << " Residual norm: " << solveTask->resNorm << std::endl;

	return 0;
}
