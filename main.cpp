#include <iostream>
#include <memory>
#include <map>
#include <thread>
#include <chrono>

#include "solvertest.h"
#include "petsccgtest.h"

int main(int argc, char* argv[])
{
	ProcessController* pc = ProcessController::getInstance();

	if (pc->MPIRank())
	{
		pc->waitForTask();
		return 0;
	}

  std::unique_ptr<SolverTest> test4k1 = std::make_unique<PETScCGTest>();
  test4k1->A_name = "TestHomoStress4k.mtx";
  test4k1->file_format = MatrixFileFormat::MTX;
  test4k1->test();

  std::unique_ptr<SolverTest> test4k2 = std::make_unique<PETScCGTest>();
  test4k2->A_name = "TestHomoStress4k.mtx";
  test4k2->file_format = MatrixFileFormat::MTX;
  test4k2->test();

/*
  std::unique_ptr<SolverTest> test400k = std::make_unique<PETScCGTest>();
  test400k->A_name = "TestHomoStress400k.mtx";
  test400k->file_format = MatrixFileFormat::MTX;
  test400k->test();
*/

  pc->evaluateTask(Task::Shutdown);

	return 0;
}
