#ifndef PROCESS_CONTROLLER
#define PROCESS_CONTROLLER

#include "mpi.h"
#include "ProcessTask.h"

#include <petscksp.h>
#include <map>
#include <iostream>
#include <vector>

class ProcessController
{
private:
	std::map<Task, ProcessTask*> tasks_map;
	int MPI_rank, MPI_Size;
	MPI_Comm communicator;

	ProcessController();
	static ProcessController* _instance;

public:

	ProcessController(ProcessController&) = delete;
	void operator=(const ProcessController&) = delete;

	int MPIRank() const;
	int MPISize() const;
	MPI_Comm getCommunicator() const;

	void waitForTask();
	void evaluateTask(Task taskID);
	ProcessTask* getTask(Task taskID);

	static ProcessController* getInstance();
	~ProcessController();
};

#endif//PROCESS_CONTROLLER
