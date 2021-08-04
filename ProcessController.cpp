#include "processcontroller.h"
#include "petsccgtest.h"
#include "process.h"

int ProcessController::MPIRank() const
{
	return MPI_rank;
}

int ProcessController::MPISize() const
{
	return MPI_size;
}

MPI_Comm ProcessController::getCommunicator() const
{
	return communicator;
}

ProcessController::ProcessController()
{
	communicator = MPI_COMM_WORLD;
	MPI_Init(0, nullptr);                   
	MPI_Comm_rank(communicator, &MPI_rank); 
	MPI_Comm_size(communicator, &MPI_size); 

#ifdef _DEBUG
	int my_pid = _getpid();
	std::cout << " Rank [" << MPI_rank << "]: My pid = " << my_pid << std::endl;
	if (!MPI_rank)
	{
		std::cout << " Press any key to continue..." << std::endl;
		getchar();
	}
	MPI_Barrier(communicator);
#endif

	tasks_map[Task::KSPSolve] = std::make_shared<KSPSolveTask>();
	tasks_map[Task::Shutdown] = std::make_shared<ShutdownTask>();
}

ProcessController::~ProcessController(){}

void ProcessController::evaluateTask(Task taskID)
{
	// Sending command to helper processes
	if (!MPI_rank) MPI_Bcast(&taskID, 1, MPI_INT, 0, communicator);

	// Checking taskID
	try
	{
		if (tasks_map.find(taskID) == tasks_map.end())
		{
			std::cout << " Rank [" << MPI_rank << "]: Got unknown task (Integer Id = "
				<< static_cast<int>(taskID) << ")." << std::endl;
			throw ProcessTaskException("Error: unknown taskID.");
		}

    namespace pcsh = ProcessControllerSignalHandler;
		pcsh::error_handler_MPI_rank = MPI_rank;
		pcsh::pushSignalHandler();

		tasks_map[taskID]->task();

    pcsh::popSignalHandler();
	}
	catch (const std::exception& exc)
	{
		std::cout << " Rank [" << MPI_rank << "]: caught an exception: \"" << exc.what() << "\"." << std::endl;
		std::cout << " Rank [" << MPI_rank << "]: Calling MPI_Abort()." << std::endl;
		MPI_Abort(communicator, MPI_ERR_OTHER);
	}
	catch (...)
	{
		std::cout << " Rank [" << MPI_rank << "]: caught an exception." << std::endl;
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
	exit(0);
}

std::shared_ptr<ProcessTask> ProcessController::getTask(Task taskID)
{
	return tasks_map[taskID];
}

std::shared_ptr<ProcessController> ProcessController::_instance = std::shared_ptr<ProcessController>(nullptr);

std::shared_ptr<ProcessController> ProcessController::getInstance()
{
	if (_instance.get() == nullptr)
	{
		_instance = std::shared_ptr<ProcessController>(new ProcessController);
	}
	return _instance;
}


void ProcessController::initialize()
{
	if (MPI_rank != 0) waitForTask();
}

void ProcessController::finalize()
{
	evaluateTask(Task::Shutdown);
}
