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
	virtual void task() = 0;
	virtual ~ProcessTask(){};
};

class ShutdownTask : public ProcessTask
{
public:
  void task() override
  {
	  //PetscErrorCode ierr = PetscFinalize(); //CHKERRQ(ierr);
	  PetscErrorCode ierr = -1;
	  if (ierr != 0) throw ShutdownException("Shutdown Error", ierr);
  }
};


#endif//PROCESS_TASK
