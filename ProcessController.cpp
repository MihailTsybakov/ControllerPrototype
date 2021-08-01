#include "processcontroller.h"
#include "petsccgtest.h"

int ProcessController::MPIRank() const
{
	return mpiRank;
}

int ProcessController::MPISize() const
{
	return mpiSize;
}

MPI_Comm ProcessController::getCommunicator() const
{
	return communicator;
}

ProcessController::ProcessController()
{
	communicator = MPI_COMM_WORLD;
	PetscInitialize(0, nullptr, (char*)0, (char*)0);
	MPI_Comm_rank(communicator, &mpiRank);
	MPI_Comm_size(communicator, &mpiSize);

	tasksMap[Task::KSPSolve] = new KSPSolveTask;
	tasksMap[Task::Shutdown] = new ShutdownTask;
}

ProcessController::~ProcessController()
{
	if (_instance != nullptr) delete _instance;
	for (auto pair : tasksMap) delete pair.second;
}

void ProcessController::evaluateTask(Task taskID)
{
	// Sending command to helper processes
	if (!mpiRank) MPI_Bcast(&taskID, 1, MPI_INT, 0, communicator);

	// Checking taskID
	if (tasksMap.find(taskID) == tasksMap.end())
	{
		std::cout << " Rank [" << mpiRank << "]: Get unknown task (Integer Id = " 
              << static_cast<int>(taskID) << "). Shutting down..." << std::endl;
		tasksMap[Task::Shutdown]->task();
	}
	//std::cout << " Rank [" << mpiRank << "]: evaluating task " << static_cast<int>(taskID) << std::endl;
	tasksMap[taskID]->task();
}

void ProcessController::waitForTask()
{
	bool finishFlag = false;
	while (!finishFlag)
	{
		int taskID;
		MPI_Bcast(&taskID, 1, MPI_INT, 0, communicator);

		if (taskID == -1) finishFlag = true;
		evaluateTask(static_cast<Task>(taskID));
	}
}

ProcessTask* ProcessController::getTask(Task taskID)
{
	return tasksMap[taskID];
}

ProcessController* ProcessController::_instance = nullptr;

ProcessController* ProcessController::getInstance()
{
	if (_instance == nullptr)
	{
		_instance = new ProcessController;
	}
	return _instance;
}