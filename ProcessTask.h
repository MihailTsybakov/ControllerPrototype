#ifndef PROCESS_TASK
#define PROCESS_TASK

enum class Task : int
{
	KSPSolve = 1,
	Shutdown = -1
};

class ProcessTask
{
public:
	ProcessTask() {};
	virtual void task() = 0;
	virtual ~ProcessTask(){};
};

#endif//PROCESS_TASK
