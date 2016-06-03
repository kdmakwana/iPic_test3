#include <iomanip>
#include "iPic3D.h"
#include "runpic.h"
#include "AMRVAC_coupling.h"
using namespace iPic3D;

void myfunc(double nrg, MPI_Fint f_handle) {
  
  MPI_COMM_IPIC3D = MPI_Comm_f2c(f_handle);

  iPic3D::c_Solver KCode;
  bool b_err = false; 

  char  arg0[] = "./hello";
  char  arg1[] = "./GEM.inp";
  char* kargv[] = { &arg0[0], &arg1[0], NULL };
  int kargc = (int)(sizeof(kargv)/sizeof(kargv[0]))-1;
  cout << " kargv testing " << kargv[0] << kargv[1] << endl;

  //MPI_Comm IPIC_COMM; 
  //IPIC_COMM = MPI_Comm_f2c(f_handle);
  //MPI_COMM_IPIC3D = IPIC_COMM;

  KCode.Init(kargc, kargv);
  KCode.InjectBoundaryParticles();
  KCode.GatherMoments();
 
  //the main time advance loop
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

}



