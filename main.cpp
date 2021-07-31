#include <iostream>
#include <memory>
#include <map>
#include <thread>
#include <chrono>

#include "ProcessController.h"
#include "solvertest.h"

ProcessController* ProcessController::_instance = nullptr;

int main(int argc, char* argv[])
{
	ProcessController* pc = ProcessController::getInstance();

	if (pc->rank())
	{
		pc->waitForTask();
		return 0;
	}

	std::unique_ptr<SolverTest> test = std::make_unique<PETScCGTest>();
	test->argc = argc;
	test->argv = argv;
	test->A_name = "C:\\tcybakov\\MPI_Ksp_Solver\\repository\\Testing\\TestHomoStress4k.mtx";
	test->file_format = MatrixFileFormat::MTX;
	test->test();

	return 0;
}
