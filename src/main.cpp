#include <iostream>

#include <vector>

#include "mpi.h"
#include "input.h"
#include "PFDD.h"

using namespace PFDD_NS;

int main(int argc, char **argv)

{
  MPI_Init(&argc, &argv);

  PFDD_C *pfdd_p = new PFDD_C(argc, argv, MPI_COMM_WORLD);
  pfdd_p->input->file();
  delete pfdd_p;
  //VPFFT::UnitTester::BasicTestMain();

  MPI_Finalize();
  return 0;
}
