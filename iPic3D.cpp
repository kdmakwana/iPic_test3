
#include <iomanip>
#include "iPic3D.h"
#include "AMRVAC_coupling.h"
#include "mpi.h"

using namespace iPic3D;

int main(int argc, char **argv) {

  iPic3D::c_Solver KCode;
  bool b_err = false;

  //#########################the hello world part (later AMRVAC)#################
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  printf("Hello world from processor %s, rank %d"
         " out of %d processors\n",
         processor_name, world_rank, world_size);

  //kdm addition 28-01-2016
  //cout << "global_int = " << global_int << endl;
  //global_int = 489; 
  

  //#################start of the iPic3D part#####################################
  //Get the master group of all processes available
  MPI_Group MPI_master_group;
  MPI_Comm_group(MPI_COMM_WORLD, &MPI_master_group);

  //Create a new sub-group of first 4 processors with ranks 0,1,2,3
  int num_sub_procs = 4;
  const int sub_ranks [4] = {0, 1, 2, 3};
  MPI_Group MPI_ipic_group; //the sub-group for iPic3D processes
  MPI_Group_incl(MPI_master_group, num_sub_procs, sub_ranks, &MPI_ipic_group);

  //Create a new communicator for the iPic sub-group, sort of like 
  //MPI_COMM_WORLD but only for the ipic processes
  //MPI_Comm MPI_IPIC_COMM;
  MPI_Comm_create_group(MPI_COMM_WORLD, MPI_ipic_group, 0, &MPI_COMM_IPIC3D);

  //check if I am on the Ipic communicator before starting ipic
  if (MPI_COMM_NULL != MPI_COMM_IPIC3D) {

    int group_size;
    MPI_Comm_size(MPI_COMM_IPIC3D, &group_size);
    printf("group size = %d \n", group_size);

    KCode.Init(argc, argv);
    KCode.InjectBoundaryParticles();
    KCode.GatherMoments();

    /////* ------------ */
    /////* 1- Main loop */
    /////* ------------ */

    for (int i = KCode.FirstCycle(); i <= KCode.LastCycle(); i++) {

      if (KCode.get_myrank() == 0) cout << " ======= Cycle " << i << " ======= " << endl;

      /* ----------------------------------------------------- */
      /* 2- Calculate fields and move particles                */
      /*    Exit if there is a memory issue with the particles */
      /* ----------------------------------------------------- */

      KCode.UpdateCycleInfo(i);
      KCode.CalculateField();

      b_err = KCode.ParticlesMover();

      if (!b_err) KCode.CalculateBField();
      if (!b_err) KCode.GatherMoments();
      if ( b_err) i = KCode.LastCycle() + 1;

      /* --------------- */
      /* 3- Output files */
      /* --------------- */

      KCode.WriteOutput(i);
      KCode.WriteConserved(i);
      KCode.WriteRestart(i);

    }

    //KCode.Finalize();
  
  }

  //MPI_Group_free(&MPI_ipic_group);
  //MPI_Comm_free(&MPI_COMM_IPIC3D);

  // Finalize the MPI environment.
  MPI_Finalize();

  return 0;
}
