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


	pc->evaluateTask(Task::Shutdown);

	return 0;
}
