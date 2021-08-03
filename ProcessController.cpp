#include "processcontroller.h"
#include "petsccgtest.h"
#include "process.h"

int ProcessController::MPIRank() const
{
	return MPI_rank;
}

int ProcessController::MPISize() const
{
	return MPI_Size;
}

MPI_Comm ProcessController::getCommunicator() const
{
	return communicator;
}

ProcessController::ProcessController()
{
	communicator = MPI_COMM_WORLD;
	//PetscInitialize(0, nullptr, (char*)0, (char*)0);
	MPI_Init(0, nullptr);
	MPI_Comm_rank(communicator, &MPI_rank);
	MPI_Comm_size(communicator, &MPI_Size);

#ifdef _DEBUG
	int my_pid = _getpid();
	std::cout << " Rank [" << MPI_rank << "]: My pid = " << my_pid << std::endl;
	if (!MPI_rank)
	{
		std::cout << "Press any key to continue..." << std::endl;
		getchar();
	}
	MPI_Barrier(communicator);
#endif

	tasks_map[Task::KSPSolve] = new KSPSolveTask;
	tasks_map[Task::Shutdown] = new ShutdownTask;
}

ProcessController::~ProcessController()
{
	if (_instance != nullptr) delete _instance;
	for (auto pair : tasks_map) delete pair.second;
}

void ProcessController::evaluateTask(Task taskID)
{
	// Sending command to helper processes
	if (!MPI_rank) MPI_Bcast(&taskID, 1, MPI_INT, 0, communicator);

	// Checking taskID
	if (tasks_map.find(taskID) == tasks_map.end())
	{
		std::cout << " Rank [" << MPI_rank << "]: Got unknown task (Integer Id = " 
              << static_cast<int>(taskID) << "). Shutting down..." << std::endl;
		tasks_map[Task::Shutdown]->task();
	}
	try
	{
		tasks_map[taskID]->task();
	}
	catch (const std::exception& exc)
	{
		std::cout << " Rank [" << MPI_rank << "]: caught an exception: \"" << exc.what() << "\"." << std::endl;
		std::cout << " Rank [" << MPI_rank << "]: Calling MPI_Abort()." << std::endl;
		MPI_Abort(communicator, MPI_ERR_OTHER);
	}
}

void ProcessController::waitForTask()
{
	bool finish_flag = false;
	while (!finish_flag)
	{
		int taskID;
		MPI_Bcast(&taskID, 1, MPI_INT, 0, communicator);

		if (taskID == -1) finish_flag = true;
		evaluateTask(static_cast<Task>(taskID));
	}
}

ProcessTask* ProcessController::getTask(Task taskID)
{
	return tasks_map[taskID];
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