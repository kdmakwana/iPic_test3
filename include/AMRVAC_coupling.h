/***************************************************************************
  AMRVAC_coupling.h
  -------------------
 ***************************************************************************/

#ifndef AMRVAC_COUPLING_H
#define AMRVAC_COUPLING_H

#include "mpi.h"
#include <math.h>
#include "Alloc.h"
#include <iostream>
#include <fstream>
using namespace std;

extern MPI_Comm MPI_COMM_IPIC3D; 

class amrvac_interf
{
  public :

    int global_int, nn_dim, ncomponents, nx_mhd, ny_mhd, nz_mhd ;
    double dx_mhd, dy_mhd, dz_mhd, xl_mhd, xh_mhd, yl_mhd, yh_mhd, zl_mhd, zh_mhd; 
    double ***rho_mhd;
    double ***prs_mhd;
    double ***Bx_mhd; //allocatable arrays for the fields
    double ***By_mhd;
    double ***Bz_mhd;
    double ***Ex_mhd; 
    double ***Ey_mhd;
    double ***Ez_mhd;
    double ***ux_mhd;
    double ***uy_mhd;
    double ***uz_mhd;
    double ***jx_mhd;
    double ***jy_mhd;
    double ***jz_mhd;
    double ***xval_mhd;
    double ***yval_mhd; 

    //####################### initializing sim box ###################################
    int init_sim_box()
    {

      ifstream myfile;
      myfile.open("pic.property", ios::in|ios::binary);
      if (myfile.is_open())
      {
        myfile.read( (char*) &nn_dim, sizeof(int));
        myfile.read( (char*) &ncomponents, sizeof(int));
        myfile.read( (char*) &nx_mhd, sizeof(int));
        myfile.read( (char*) &ny_mhd, sizeof(int));
        myfile.read( (char*) &dx_mhd, sizeof(double));
        myfile.read( (char*) &dy_mhd, sizeof(double));
        myfile.read( (char*) &xl_mhd, sizeof(double));
        myfile.read( (char*) &xh_mhd, sizeof(double));
        myfile.read( (char*) &yl_mhd, sizeof(double));
        myfile.read( (char*) &yh_mhd, sizeof(double));
        myfile.close();
        //adjusting for the guard cells and node values
        nx_mhd = nx_mhd+3;
        ny_mhd = ny_mhd+3;
        nz_mhd = nz_mhd+3;
        if (nn_dim == 2)
        {
          nz_mhd = 1;
          zl_mhd = 0.0;
          zh_mhd = 1.0;
          dz_mhd = 1.0;
        }
      }
      myfile.close();
      return 0;
    };
    //#####################end of initializing simulation box#######################


    //###########################initializing the field arays#######################
    int fields_from_amrvac()
    {
      double read_buf[16];
      Bx_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      By_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      Bz_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      Ex_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      Ey_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      Ez_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      xval_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      yval_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      ifstream myfile;
      myfile.open("pic.dat", ios::in|ios::binary);
      if (myfile.is_open())
      {
        for (int k = 0; k < nz_mhd; k++) {
          for (int j = 0; j < ny_mhd; j++) {
            for (int i = 0; i < nx_mhd; i++) {
              myfile.read( (char*) &read_buf, sizeof(read_buf));
              Bx_mhd[i][j][k] = read_buf[7];
              By_mhd[i][j][k] = read_buf[8];
              Bz_mhd[i][j][k] = read_buf[9];
              Ex_mhd[i][j][k] = read_buf[10];
              Ey_mhd[i][j][k] = read_buf[11];
              Ez_mhd[i][j][k] = read_buf[12];
              xval_mhd[i][j][k] = read_buf[0];
              yval_mhd[i][j][k] = read_buf[1];
            }
          }
        }         
      }
      myfile.close();

      //double value = linear_interpolation_1(xval_mhd, 113.0, 161.73, -4.3, dx_mhd, dy_mhd, dz_mhd,
      //          nx_mhd, ny_mhd, nz_mhd, true);
      //cout << "inside fields_from_amrvac xval (113.0,161.73,-4.3)  " << value << endl ; 
      //value = linear_interpolation_1(yval_mhd, 113.0, 161.73, -4.3, dx_mhd, dy_mhd, dz_mhd,
      //          nx_mhd, ny_mhd, nz_mhd, true);
      //cout << "inside fields_from_amrvac yval (113.0,161.73,-4.3)  " << value << endl ; 
      //value = linear_interpolation_1(xval_mhd, -113.0, -161.73, -4.3, dx_mhd, dy_mhd, dz_mhd,
      //          nx_mhd, ny_mhd, nz_mhd, true);
      //cout << "inside fields_from_amrvac xval (-113.0,-161.73,-4.3)  " << value << endl ; 
      //value = linear_interpolation_1(yval_mhd, -113.0, -161.73, -4.3, dx_mhd, dy_mhd, dz_mhd,
      //          nx_mhd, ny_mhd, nz_mhd, true);
      //cout << "inside fields_from_amrvac yval (-113.0,-161.73,-4.3)  " << value << endl ; 


      return 0;
    }
    //############################end of initializing fields########################
    
    
    //###########################initializing the velocity arays#######################
    int vel_from_amrvac()
    {
      double read_buf[16];
      rho_mhd= newArr3(double,nx_mhd,ny_mhd,nz_mhd);
      prs_mhd= newArr3(double,nx_mhd,ny_mhd,nz_mhd);
      ux_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      uy_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      uz_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      jx_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      jy_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      jz_mhd = newArr3(double,nx_mhd,ny_mhd,nz_mhd);    
      ifstream myfile;
      myfile.open("pic.dat", ios::in|ios::binary);
      if (myfile.is_open())
      {
        for (int k = 0; k < nz_mhd; k++) {
          for (int j = 0; j < ny_mhd; j++) {
            for (int i = 0; i < nx_mhd; i++) {
              myfile.read( (char*) &read_buf, sizeof(read_buf));
              rho_mhd[i][j][k]= read_buf[2];
              prs_mhd[i][j][k]= read_buf[3];
              ux_mhd[i][j][k] = read_buf[4];
              uy_mhd[i][j][k] = read_buf[5];
              uz_mhd[i][j][k] = read_buf[6];
              jx_mhd[i][j][k] = read_buf[13];
              jy_mhd[i][j][k] = read_buf[14];
              jz_mhd[i][j][k] = read_buf[15];
            } 
          } 
        }
      }
      myfile.close();
      return 0;
    }
    //############################end of initializing fields########################

    

    double getVix_mhd(double x, double y, double z){
      // NOTE: the axes in the input file are assumed to be identical to iPIC3D
      double value = linear_interpolation(ux_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      // The velocity in Km/s
      //value = value *1e3 / code_V;
      return(value);
    } 
    
    double getViy_mhd(double x, double y, double z){
      // NOTE: the axes in the input file are assumed to be identical to iPIC3D
      double value = linear_interpolation(uy_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      // The velocity in Km/s
      //value = value *1e3 / code_V;
      return(value);
    } 
    
    double getViz_mhd(double x, double y, double z){
      // NOTE: the axes in the input file are assumed to be identical to iPIC3D
      double value = linear_interpolation(uz_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      // The velocity in Km/s
      //value = value *1e3 / code_V;
      return(value);
    } 

    double getrho_mhd(double x, double y, double z){
      // NOTE: the axes in the input file are assumed to be identical to iPIC3D
      double value = linear_interpolation_1(rho_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      // The velocity in Km/s
      //value = value *1e3 / code_V;
      return(value);
    } 
    
    double getprs_mhd(double x, double y, double z){
      // NOTE: the axes in the input file are assumed to be identical to iPIC3D
      double value = linear_interpolation_1(prs_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      // The velocity in Km/s
      //value = value *1e3 / code_V;
      return(value);
    } 
    
    double getufx_mhd(double x, double y, double z){
      double value = linear_interpolation_1(ux_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      return(value);
    } 
    
    double getufy_mhd(double x, double y, double z){
      double value = linear_interpolation_1(uy_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      return(value);
    } 

    double getufz_mhd(double x, double y, double z){
      double value = linear_interpolation_1(uz_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      return(value);
    } 
    
    double getjfx_mhd(double x, double y, double z){
      double value = linear_interpolation_1(jx_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      return(value);
    } 
    
    double getjfy_mhd(double x, double y, double z){
      double value = linear_interpolation_1(jy_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      return(value);
    } 

    double getjfz_mhd(double x, double y, double z){
      double value = linear_interpolation_1(jz_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      return(value);
    } 

    double getbx_mhd(double x, double y, double z){
      double value = linear_interpolation_1(Bx_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      return(value);
    } 
    
    double getby_mhd(double x, double y, double z){
      double value = linear_interpolation_1(By_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      return(value);
    } 

    double getbz_mhd(double x, double y, double z){
      double value = linear_interpolation_1(Bz_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      return(value);
    } 

    double getex_mhd(double x, double y, double z){
      double value = linear_interpolation_1(Ex_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      return(value);
    } 
    
    double getey_mhd(double x, double y, double z){
      double value = linear_interpolation_1(Ey_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      return(value);
    } 
    
    double getez_mhd(double x, double y, double z){
      double value = linear_interpolation_1(Ez_mhd, x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true);
      return(value);
    }

    double getex_direct(double x, double y, double z){
      double value = getval_direct_file(x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true, 10);
      return(value);
    }
 
    double getey_direct(double x, double y, double z){
      double value = getval_direct_file(x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true, 11);
      return(value);
    }

    double getez_direct(double x, double y, double z){
      double value = getval_direct_file(x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true, 12);
      return(value);
    }
    
    double getbx_direct(double x, double y, double z){
      double value = getval_direct_file(x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true, 7);
      return(value);
    }
    
    double getby_direct(double x, double y, double z){
      double value = getval_direct_file(x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true, 8);
      return(value);
    }

    double getbz_direct(double x, double y, double z){
      double value = getval_direct_file(x, y, z, dx_mhd, dy_mhd, dz_mhd,
                nx_mhd, ny_mhd, nz_mhd, true, 9);
      return(value);
    }


    //###################linear interpolation from Gianni's iPic-FromFluid############################
    inline double linear_interpolation(double*** A, double xp, double yp, double zp,
                double ext_dx, double ext_dy, double ext_dz,
                int ext_nx, int ext_ny, int ext_nz, bool flag2d)
    {
      // NOTE: teh axis are switched between MHD and iPic3D!!! y and z need to be invereted
      // NGP approximation
      double xi[2]; double eta[2]; double zeta[2];
      double weight[2][2][2];

      double xstart=0;
      double ystart=0;
      double zstart=0;

      int ix = floor((xp - xstart) / ext_dx);
      int iy = floor((yp - ystart) / ext_dy);
      int iz = floor((zp - zstart) / ext_dz);
      if (ix < 0) ix=0;
      if (ix > ext_nx-2) ix = ext_nx-2;
      if (iy < 0) iy=0;
      if (iy > ext_ny-2) iy = ext_ny-2;
      if (iz < 0) iz=0;
      if (iz > ext_nz-2) iz = ext_nz-2;
      //iy = floor(ext_ny/2);

      xi[1]   = xp - ix * ext_dx;
      eta[1]  = yp - iy * ext_dy;
      zeta[1] = zp - iz * ext_dz;
      xi[0]   = -xp + (ix+1) * ext_dx;
      eta[0]  = -yp + (iy+1) * ext_dy;
      zeta[0] = -zp + (iz+1) * ext_dz;
      if(flag2d)
      {
        iz = floor(ext_nz/2);
        zeta[1]  = 0.0;
        zeta[0]  = ext_dz;
      }
      for (int ii=0; ii < 2; ii++)
        for (int jj=0; jj < 2; jj++)
          for(int kk=0; kk < 2; kk++){
              weight[ii][jj][kk] = xi[ii]*eta[jj]*zeta[kk]/ext_dx/ext_dy/ext_dz;
      }
      double value = 0.0;

      for (int ii=0; ii < 2; ii++)
        for (int jj=0; jj < 2; jj++)
          for(int kk=0; kk < 2; kk++){
            value += weight[ii][jj][kk] * A[ix+ii][iy+jj][iz+kk];
      }
      return(value);
    }
    //#################end of linear interpolation##########################################


    //###################kdms modified linear interpolation routine#########################
    inline double linear_interpolation_1(double*** A, double xp, double yp, double zp,
                double ext_dx, double ext_dy, double ext_dz,
                int ext_nx, int ext_ny, int ext_nz, bool flag2d)
    {
      double xi[2]; double eta[2]; double zeta[2];
      double weight[2][2][2];

      //starting from the guard cell corner
      double xstart=0.0-ext_dx;
      double ystart=0.0-ext_dy;
      double zstart=0.0-ext_dy;

      int ix = floor((xp - xstart) / ext_dx);
      int iy = floor((yp - ystart) / ext_dy);
      int iz = floor((zp - zstart) / ext_dz);
      if (ix < 0) ix=0;
      if (ix > ext_nx-2) ix = ext_nx-2;
      if (iy < 0) iy=0;
      if (iy > ext_ny-2) iy = ext_ny-2;
      if (iz < 0) iz=0;
      if (iz > ext_nz-2) iz = ext_nz-2;
      //iy = floor(ext_ny/2);

      xi[1]   = xp - xstart - ix * ext_dx;
      eta[1]  = yp - ystart - iy * ext_dy;
      zeta[1] = zp - zstart - iz * ext_dz;
      xi[0]   = (ix+1) * ext_dx - (xp-xstart);
      eta[0]  = (iy+1) * ext_dy - (yp-ystart);
      zeta[0] = (iz+1) * ext_dz - (zp-zstart);
      if(flag2d)
      {
        iz = floor(ext_nz/2);
        zeta[1]  = 0.0;
        zeta[0]  = ext_dz;
      }
      for (int ii=0; ii < 2; ii++)
        for (int jj=0; jj < 2; jj++)
          for(int kk=0; kk < 2; kk++){
              weight[ii][jj][kk] = xi[ii]*eta[jj]*zeta[kk]/ext_dx/ext_dy/ext_dz;
      }
      double value = 0.0;

      for (int ii=0; ii < 2; ii++)
        for (int jj=0; jj < 2; jj++)
          for(int kk=0; kk < 2; kk++){
            value += weight[ii][jj][kk] * A[ix+ii][iy+jj][iz+kk];
      }
      return(value);
    }
    //#################end of kdm's linear interpolation######################################

    
    //###################kdms modified linear interpolation routine#########################
    inline double getval_direct_file(double xp, double yp, double zp,
                double ext_dx, double ext_dy, double ext_dz,
                int ext_nx, int ext_ny, int ext_nz, bool flag2d, int which_field)
    {
      double xi[2]; double eta[2]; double zeta[2];
      double weight[2][2][2];
      double read_val;
      int offset;

      //starting from the guard cell corner
      double xstart=0.0-ext_dx;
      double ystart=0.0-ext_dy;
      double zstart=0.0-ext_dy;

      int ix = floor((xp - xstart) / ext_dx);
      int iy = floor((yp - ystart) / ext_dy);
      int iz = floor((zp - zstart) / ext_dz);
      if (ix < 0) ix=0;
      if (ix > ext_nx-2) ix = ext_nx-2;
      if (iy < 0) iy=0;
      if (iy > ext_ny-2) iy = ext_ny-2;
      if (iz < 0) iz=0;
      if (iz > ext_nz-2) iz = ext_nz-2;
      //iy = floor(ext_ny/2);

      xi[1]   = xp - xstart - ix * ext_dx;
      eta[1]  = yp - ystart - iy * ext_dy;
      zeta[1] = zp - zstart - iz * ext_dz;
      xi[0]   = (ix+1) * ext_dx - (xp-xstart);
      eta[0]  = (iy+1) * ext_dy - (yp-ystart);
      zeta[0] = (iz+1) * ext_dz - (zp-zstart);
      if(flag2d)
      {
        iz = floor(ext_nz/2);
        zeta[1]  = 0.0;
        zeta[0]  = ext_dz;
      }
      for (int ii=0; ii < 2; ii++)
        for (int jj=0; jj < 2; jj++)
          for(int kk=0; kk < 2; kk++){
              weight[ii][jj][kk] = xi[ii]*eta[jj]*zeta[kk]/ext_dx/ext_dy/ext_dz;
      }
      double value = 0.0;
      
      ifstream myfile;
      myfile.open("pic.dat", ios::in|ios::binary);

      for (int ii=0; ii < 2; ii++)
        for (int jj=0; jj < 2; jj++)
          for (int kk=0; kk < 2; kk++){
            offset = (iy+jj)*ext_nx+(ix+ii);
            myfile.seekg((offset*16+which_field)*8, ios::beg); 
            myfile.read( (char*) &read_val, sizeof(read_val));
            //cout << " ix, iy, readval = " << ix+ii << " " << iy+jj << " " << read_val << endl ;
            value += weight[ii][jj][kk] * read_val ;
      }

      myfile.close();

      return(value);
    }

};


#endif



