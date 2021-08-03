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
	virtual void task() noexcept(false) = 0;
	virtual ~ProcessTask(){};
};

class ShutdownTask : public ProcessTask
{
public:
  void task() override
  {
	  PetscErrorCode ierr = MPI_Finalize(); 
	  if (ierr != 0) throw ProcessTaskException("PetscFinalize() error", ierr);
  }
};


#endif//PROCESS_TASK
