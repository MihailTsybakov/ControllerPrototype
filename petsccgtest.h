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
  int MPI_rank, MPI_size;

  int _task();
  void task() override;
  /// Synchronizes Matrix sizes between processes
  void syncMatrixSize();
  /// Check print
  void printContext() const;

  // === Task context ===

  //  -- Synchronized context --
  int matrixSize;

  std::vector<int> loc_ia_sizes;
  std::vector<int> loc_ja_sizes;
  std::vector<int> loc_ia_starts;
  std::vector<int> loc_ja_starts;

  //  -- Personal context --
  std::vector<int> ia;   // Root only
  std::vector<int> ja;   // Root only
  std::vector<double> a; // Root only

  Mat A;
  Vec result, ref_result, b;
  KSP ksp; // Krylov solver
  PC pc;   // Preconditioner
  PetscInt iterations;
  PetscScalar res_norm;

  std::vector<int> loc_ia;
  std::vector<int> loc_ja;
  std::vector<double> loc_a;

  PetscInt* PETSc_loc_ia;
  PetscInt* PETSc_loc_ja;
  PetscScalar* PETSc_loc_a;

  int loc_num;
  int loc_rows;
};

#endif // PETSC_SOLVERS_TESTS