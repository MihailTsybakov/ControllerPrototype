#pragma once
#ifndef PETSC_SOLVERS_TESTS
#define PETSC_SOLVERS_TESTS

#include "solvertest.h"

class PETScCGTest : public LAESTest {
protected:
  string petscfilename;
  cr::time_point<cr::system_clock> scatter_start, scatter_end;

  int test();
  int testSpecific();
};

class KSPSolveTask : public ProcessTask
{
public:
  MPI_Comm communicator;
  int MPIRank, MPISize;

  void task() override;
  /// Synchronizes Matrix sizes between processes
  void syncMatrixSize();
  /// Allocates memory for ia/ja/a markup before synchronizing
  void prepareMarkup();
  /// Synchronizes markup
  void syncMarkup();
  /// Check print
  void printContext() const;

  // === Task context ===

  //  -- Synchronized context --
  int matrixSize;

  std::vector<int> locIASizes;
  std::vector<int> locJASizes;
  std::vector<int> locIAStarts;
  std::vector<int> locJAStarts;

  //  -- Personal context --
  std::vector<int> ia;   // Root only
  std::vector<int> ja;   // Root only
  std::vector<double> a; // Root only

  Mat A;
  Vec result, refResult, b;
  KSP ksp; // Krylov solver
  PC pc;   // Preconditioner
  PetscInt iterations;
  PetscScalar resNorm;

  std::vector<int> locIA;
  std::vector<int> locJA;
  std::vector<double> locA;

  PetscInt* petscLocIA;
  PetscInt* petscLocJA;
  PetscScalar* petscLocA;

  int locNum;
  int locRows;
};

#endif // PETSC_SOLVERS_TESTS