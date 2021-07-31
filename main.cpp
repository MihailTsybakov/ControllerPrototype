#include <iostream>
#include "mpi.h"

#include <memory>
#include <map>
#include <thread>
#include <chrono>

class ProcessTask;
class CalcPITask;
class ShutdownTask;

enum class Task : int
{
	CalcPI = 1,
	Shutdown = -1
};

class ProcessController
{
private:

	ProcessController();
	static ProcessController* _instance;
	int mpiRank, mpiSize;
public:
	std::map<Task, ProcessTask*> tasks_map;
	void waitForTask();
	ProcessController(ProcessController& to_copy) = delete;
	void operator=(const ProcessController& equals_to) = delete;

	void evaluateTask(Task taskID);
	int getRank() { return mpiRank; }
	int getSize() { return mpiSize; }

	static ProcessController* getInstance();
	~ProcessController()
	{
		if (_instance != nullptr) delete _instance;
		for (auto pair : tasks_map)
		{
			delete pair.second;
		}
	}
};

ProcessController* ProcessController::_instance = nullptr;
ProcessController* ProcessController::getInstance()
{
	if (_instance == nullptr)
	{
		_instance = new ProcessController();
	}
	return _instance;
}
/// =================================================================

class ProcessTask
{
private:

public:
	ProcessTask(){}
	virtual void task() = 0;
	virtual ~ProcessTask(){}
};


class CalcPITask : public ProcessTask
{
private:
	int rank;
	int size; 
public:
	int N;
	double d;
	double PI;

	void sync()
	{
		MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&d, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	void task() override
	{
		sync();
		ProcessController* pc_task = ProcessController::getInstance();
		rank = pc_task->getRank();
		size = pc_task->getSize();
		std::cout << " > rank [" << rank << "] working on PI" << std::endl;
		std::cout << " > rank [" << rank << "] N = " << N << std::endl;
		std::cout << " > rank [" << rank << "] d = " << d << std::endl;
		double I = 0.0, res = 0.0;
		for (int i = rank; i < N; i += size)
		{
			I += d / (1 + d * d * i * i);
		}
		MPI_Reduce(&I, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		PI = 4 * res;
	}
};

class ShutdownTask : public ProcessTask
{
public:
	void task() override
	{
		MPI_Finalize();
	}
};

ProcessController::ProcessController()
{
	MPI_Init(0, nullptr);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

	tasks_map[Task::CalcPI] = new CalcPITask;
	tasks_map[Task::Shutdown] = new ShutdownTask;
}

void ProcessController::waitForTask()
{
	bool finish_flag = false;
	while (!finish_flag)
	{
		int taskID;
		//syncContext();
		MPI_Bcast(&taskID, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (taskID == -1)
		{
			finish_flag = true;
		}
		evaluateTask(static_cast<Task>(taskID));
	}
}

void ProcessController::evaluateTask(Task taskID)
{
	double service_arg;
	if (!mpiRank)
	{
		//syncContext();
		MPI_Bcast(&taskID, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	tasks_map[taskID]->task();
	std::cout << " <logs> Rank [" << mpiRank << "]: Task " << static_cast<int>(taskID) << " done." << std::endl;
}

int main(int argc, char* argv[])
{
	ProcessController* pc = ProcessController::getInstance();

	if (pc->getRank())
	{
		pc->waitForTask();
		//MPI_Finalize();
		return 0;
	}
	int N = 1e7;
	double d = 1e-7;

	CalcPITask* task = dynamic_cast<CalcPITask*>(pc->tasks_map[Task::CalcPI]);
	task->N = N;
	task->d = d;

	pc->evaluateTask(Task::CalcPI);

	std::cout << " Pi = " << task->PI << std::endl;

	std::this_thread::sleep_for(std::chrono::milliseconds(1'000));
	pc->evaluateTask(Task::Shutdown);
	//MPI_Finalize();

	return 0;
}
