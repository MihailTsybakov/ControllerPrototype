#ifndef PROCESS_CONTROLLER
#define PROCESS_CONTROLLER

#include "mpi.h"
#include "ProcessTask.h"

#include <petscksp.h>
#include <map>
#include <iostream>
#include <vector>

// Custom class tasks
// ??
// ==================

class ProcessController
{
private:
	std::map<Task, ProcessTask*> task_map;
	int mpiRank, mpiSize;
	MPI_Comm communicator;

	ProcessController();
	static ProcessController* _instance;

public:

	ProcessController(ProcessController& to_copy) = delete;
	void operator=(const ProcessController& equals_to) = delete;

	int rank() const;
	int size() const;
	MPI_Comm getComm() const;

	void waitForTask();
	void evaluateTask(Task taskID);
	ProcessTask* getTask(Task taskID);

	static ProcessController* getInstance();
	~ProcessController();
};


/// ==================================================================

class KSPSolveTask : public ProcessTask
{
public:
	MPI_Comm communicator = PETSC_COMM_WORLD;
	int rank, size;

	void task() override;
	/// Synchronizes Matrix sizes between processes
	void syncMSize();
	/// Allocates memory for ia/ja/a markup before synchronizing
	void prepareMarkup();
	/// Synchronizes markup
	void syncMarkup();
	/// Check print
	void printContext() const;

	// === Task context ===
	
	//  -- Synchronized context --
	int matrix_size;

	std::vector<int> loc_ia_sizes;
	std::vector<int> loc_ja_sizes;
	std::vector<int> loc_ia_starts;
	std::vector<int> loc_ja_starts;

	//  -- Personal context --
	std::vector<int> ia;   // Root only
	std::vector<int> ja;   // Root only
	std::vector<double> a; // Root only

	Mat A;
	Vec result, ref_result, b;
	KSP ksp; // Krylov solver
	PC pc;   // Preconditioner
	PetscInt iterations;
	PetscScalar res_norm;

	std::vector<int> loc_ia;
	std::vector<int> loc_ja;
	std::vector<double> loc_a;

	PetscInt* petsc_loc_ia;
	PetscInt* petsc_loc_ja;
	PetscScalar* petsc_loc_a;

	int loc_num;
	int loc_rows;
};

class ShutdownTask : public ProcessTask
{
public:
	void task() override
	{
		PetscFinalize();
	}
};

#endif//PROCESS_CONTROLLER
