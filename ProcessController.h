#ifndef PROCESS_CONTROLLER
#define PROCESS_CONTROLLER

#include "mpi.h"
#include "ProcessTask.h"

#include <petscksp.h>
#include <map>
#include <iostream>
#include <vector>
#include <memory>

class ProcessController
{
private:
	std::map<Task, std::shared_ptr<ProcessTask>> tasks_map;
	int MPI_rank, MPI_size;
	MPI_Comm communicator;

	ProcessController();
	static std::shared_ptr<ProcessController> _instance;
public:
	void initialize();
	void finalize();

	ProcessController(ProcessController&) = delete;
	void operator=(const ProcessController&) = delete;

	int MPIRank() const;
	int MPISize() const;
	MPI_Comm getCommunicator() const;

	void waitForTask();
	void evaluateTask(Task taskID);
	std::shared_ptr<ProcessTask> getTask(Task taskID);

	static std::shared_ptr<ProcessController> getInstance();
	~ProcessController();
};

#endif//PROCESS_CONTROLLER
