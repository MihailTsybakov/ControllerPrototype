#ifndef PROCESS_TASK
#define PROCESS_TASK

#include <petscksp.h>
#include "exceptions.h"

enum class Task : int
{
	KSPSolve = 1,
	Shutdown = -1
};

class ProcessTask
{
public:
	ProcessTask() {};
	virtual int task() = 0;
	virtual ~ProcessTask(){};
};

class ShutdownTask : public ProcessTask
{
public:
  int task() override
  {
	  PetscErrorCode ierr = PetscFinalize(); 
	  if (ierr != 0) throw ShutdownException("Shutdown Error", ierr);
  }
};


#endif//PROCESS_TASK
