/* .-----------------------------------------------------------------------.
 * | SFcollapse1D                                                          |
 * | Gravitational collapse of scalar fields in spherical symmetry         |
 * |                                                                       |
 * | Copyright (c) 2020, Leonardo Werneck                                  |
 * |                                                                       |
 * | This program is free software: you can redistribute it and/or modify  |
 * | it under the terms of the GNU General Public License as published by  |
 * | the Free Software Foundation, either version 3 of the License, or     |
 * | (at your option) any later version.                                   |
 * |                                                                       |
 * | This program is distributed in the hope that it will be useful,       |
 * | but WITHOUT ANY WARRANTY; without even the implied warranty of        |
 * | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
 * | GNU General Public License for more details.                          |
 * |                                                                       |
 * | You should have received a copy of the GNU General Public License     |
 * | along with this program.  If not, see <https://www.gnu.org/licenses/>.|
 * .-----------------------------------------------------------------------.
 */

/* Basic includes */
#include <fstream>
#include <cmath>
#include <vector>
#include "macros.hpp"
#include "grid.hpp"
#include "gridfunction.hpp"
#include "evolution.hpp"

using namespace std;

/* Function prototypes needed by this file */
REAL pointwise_Newton_method( const int j, grid::parameters grid, const std::vector<REAL> Phi, const std::vector<REAL> Pi, const std::vector<REAL> a );
void rescaling_of_the_lapse( grid::parameters grid, const std::vector<REAL> a, std::vector<REAL> &alpha );

/* First time step driver function */
void evolution::first_time_step( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha ) {

#if( COORD_SYSTEM == SPHERICAL )
  evolution::first_time_step_Spherical( grid, phi, Phi, Pi, a, alpha );
#elif( COORD_SYSTEM == SINH_SPHERICAL )
  evolution::first_time_step_SinhSpherical( grid, phi, Phi, Pi, a, alpha );
#endif

}

/* Generical time step driver function */
void evolution::time_step( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha ) {

#if( COORD_SYSTEM == SPHERICAL )
  evolution::time_step_Spherical( grid, phi, Phi, Pi, a, alpha );
#elif( COORD_SYSTEM == SINH_SPHERICAL )
  evolution::time_step_SinhSpherical( grid, phi, Phi, Pi, a, alpha );
#endif

}

/* First time step in Spherical coordinates */
void evolution::first_time_step_Spherical( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha ) {

  DECLARE_GRID_PARAMETERS;

  /* .---------------------------------------.
   * | Step 1: Integrating from n=0 to n=1/2 |
   * .---------------------------------------.
   */
  const REAL half_dt      = 0.5 * dt;
  const REAL half_inv_dx0 = 0.5 * inv_dx0; // This is 1.0/(2.0*dx0)
  const int J             = Nx0Total-1;
  const int Jm1           = J-1;
  const int Jm2           = J-2;

  /* .----------------------------------------------------------------.
   * | Step 1.a: Apply inner boundary conditions to Phi, a, and alpha |
   * .----------------------------------------------------------------.
   */
  Phi.level_n[0]   = 0.0;
  a.level_n[0]     = 1.0;
  alpha.level_n[0] = 1.0;
  
  /* .--------------------------------------.
   * | Step 1.b: Integrate phi, Phi, and Pi |
   * .--------------------------------------.
   */
  LOOP(1,Nx0Total-1) {
    /* Step 1.b.i: Compute useful auxiliary variables */
    const REAL r_sqr_jm1               = SQR(x[0][j-1]);
    const REAL r_sqr_jp1               = SQR(x[0][j+1]);
    const REAL r_cbd_jm1               = r_sqr_jm1 * x[0][j-1];
    const REAL r_cbd_jp1               = r_sqr_jp1 * x[0][j+1];
    const REAL alpha_over_a_jm1        = alpha.level_nm1[j-1]/a.level_nm1[j-1];
    const REAL alpha_over_a_j          = alpha.level_nm1[j]  /a.level_nm1[j];
    const REAL alpha_over_a_jp1        = alpha.level_nm1[j+1]/a.level_nm1[j+1];
    const REAL alpha_Pi_over_a_jm1     = alpha_over_a_jm1 * Pi.level_nm1[j-1];
    const REAL alpha_Pi_over_a_j       = alpha_over_a_j   * Pi.level_nm1[j];
    const REAL alpha_Pi_over_a_jp1     = alpha_over_a_jp1 * Pi.level_nm1[j+1];
    const REAL alpha_Phi_r2_over_a_jm1 = r_sqr_jm1 * alpha_over_a_jm1 * Phi.level_nm1[j-1];
    const REAL alpha_Phi_r2_over_a_jp1 = r_sqr_jp1 * alpha_over_a_jp1 * Phi.level_nm1[j+1];

    /* Step 1.b.ii: Compute the RHSs of Phi */
    const REAL rhs_Phi = half_inv_dx0 * ( alpha_Pi_over_a_jp1 - alpha_Pi_over_a_jm1 );
    const REAL rhs_Pi  = 3.0 * ( alpha_Phi_r2_over_a_jp1 - alpha_Phi_r2_over_a_jm1 )/( r_cbd_jp1 - r_cbd_jm1 );

    /* Step 1.b.iii: Evolve Phi, Pi, and phi */
    if( j==1 ) phi.level_np1[0] = phi.level_nm1[0] + half_dt * alpha_Pi_over_a_jm1;
    phi.level_n[j] = phi.level_nm1[j] + half_dt * alpha_Pi_over_a_j;
    Phi.level_n[j] = Phi.level_nm1[j] + half_dt * rhs_Phi;
    Pi.level_n[j]  = Pi.level_nm1[j]  + half_dt * rhs_Pi;
  }
  /* Step 1.b.iv: Apply boundary conditions */
  {
    /* phi: Apply no-incoming radiation boundary condition */
    const REAL coeff1 = 1.0 - half_dt / x[0][J];
    const REAL coeff2 = half_dt * half_inv_dx0;
    phi.level_n[J]    = coeff1 * phi.level_nm1[J] - coeff2 * ( 3.0 * phi.level_nm1[J] - 4.0 * phi.level_nm1[Jm1] + phi.level_nm1[Jm2] );
    /* Phi: Compute at outer boundary from phi */
    Phi.level_n[J] = half_inv_dx0 * ( 3.0 * phi.level_n[J] - 4.0 * phi.level_n[Jm1] + phi.level_n[Jm2] );
    /* Pi: Apply Neumann boundary conditions at the inner boundary */
    Pi.level_n[0] = Pi.level_n[1];
    /* Pi: Compute from the evolution equation at the outer boundary */
    const REAL r_sqr_Jm2               = SQR(x[0][Jm2]);
    const REAL r_cbd_Jm2               = r_sqr_Jm2 * x[0][Jm2];
    const REAL r_sqr_Jm1               = SQR(x[0][Jm1]);
    const REAL r_sqr_J                 = SQR(x[0][J]);
    const REAL r_cbd_J                 = r_sqr_J   * x[0][J];
    const REAL alpha_over_a_Jm2        = alpha.level_nm1[Jm2]/a.level_nm1[Jm2];
    const REAL alpha_over_a_Jm1        = alpha.level_nm1[Jm1]/a.level_nm1[Jm1];
    const REAL alpha_over_a_J          = alpha.level_nm1[J]  /a.level_nm1[J];
    const REAL alpha_Phi_r2_over_a_Jm2 = r_sqr_Jm2 * alpha_over_a_Jm2 * Phi.level_nm1[Jm2];
    const REAL alpha_Phi_r2_over_a_Jm1 = r_sqr_Jm1 * alpha_over_a_Jm1 * Phi.level_nm1[Jm1];
    const REAL alpha_Phi_r2_over_a_J   = r_sqr_J   * alpha_over_a_J   * Phi.level_nm1[J];
    const REAL rhs_Pi                  = 3.0 * ( 3.0 * alpha_Phi_r2_over_a_J - 4.0 * alpha_Phi_r2_over_a_Jm1 + alpha_Phi_r2_over_a_Jm2 )/( r_cbd_J - r_cbd_Jm2 );
    Pi.level_n[J]                      = Pi.level_nm1[J]  + half_dt * rhs_Pi;
  }

  /* .---------------------------------.
   * | Step 1.c: Integrate a and alpha |
   * .---------------------------------.
   */
  LOOP(1,Nx0Total) {
    /* Step 1.c.i: Compute a */
    a.level_n[j] = pointwise_Newton_method(j,grid,Phi.level_n,Pi.level_n,a.level_n);

    /* Step 1.c.ii: Compute auxiliary quantities */
    const REAL b = a.level_n[j] + a.level_n[j-1];
    const REAL c = a.level_n[j] - a.level_n[j-1];
    const REAL midway_r = x[0][j] + x[0][j-1];
    const REAL d = ( 1.0 - 0.25 * SQR(b) )/midway_r - inv_dx0 * c / b;

    /* Step 1.c.iii: Compute alpha */
    alpha.level_n[j] = alpha.level_n[j-1]*( 1.0 - d*dx0 )/( 1.0 + d*dx0 );
  }
  /* Step 1.c.iv: Now rescale alpha */
  rescaling_of_the_lapse(grid,a.level_n,alpha.level_n);

  /* .---------------------------.
   * | Step 1.d: Update the time |
   * .---------------------------.
   */
  grid.t = half_dt;
  
  /* .----------------------------.
   * | Step 2: Integrating to n=1 |
   * .----------------------------.
   *
   * .----------------------------------------------------------------.
   * | Step 2.a: Apply inner boundary conditions to Phi, a, and alpha |
   * .----------------------------------------------------------------.
   */
  Phi.level_np1[0]   = 0.0;
  a.level_np1[0]     = 1.0;
  alpha.level_np1[0] = 1.0;
  
  /* .--------------------------------------.
   * | Step 2.b: Integrate phi, Phi, and Pi |
   * .--------------------------------------.
   */
  LOOP(1,Nx0Total-1) {
    /* Step 2.b.i: Compute useful auxiliary variables */
    const REAL r_sqr_jm1               = SQR(x[0][j-1]);
    const REAL r_cbd_jm1               = r_sqr_jm1 * x[0][j-1];
    const REAL r_sqr_jp1               = SQR(x[0][j+1]);
    const REAL r_cbd_jp1               = r_sqr_jp1 * x[0][j+1];
    const REAL alpha_over_a_jm1        = alpha.level_n[j-1]/a.level_n[j-1];
    const REAL alpha_over_a_j          = alpha.level_n[j]  /a.level_n[j];
    const REAL alpha_over_a_jp1        = alpha.level_n[j+1]/a.level_n[j+1];
    const REAL alpha_Pi_over_a_jm1     = alpha_over_a_jm1 * Pi.level_n[j-1];
    const REAL alpha_Pi_over_a_j       = alpha_over_a_j   * Pi.level_n[j];
    const REAL alpha_Pi_over_a_jp1     = alpha_over_a_jp1 * Pi.level_n[j+1];
    const REAL alpha_Phi_r2_over_a_jm1 = r_sqr_jm1 * alpha_over_a_jm1 * Phi.level_n[j-1];
    const REAL alpha_Phi_r2_over_a_jp1 = r_sqr_jp1 * alpha_over_a_jp1 * Phi.level_n[j+1];

    /* Step 2.b.ii: Compute the RHSs of Phi */
    const REAL rhs_Phi = half_inv_dx0 * ( alpha_Pi_over_a_jp1 - alpha_Pi_over_a_jm1 );
    const REAL rhs_Pi  = 3.0 * ( alpha_Phi_r2_over_a_jp1 - alpha_Phi_r2_over_a_jm1 )/( r_cbd_jp1 - r_cbd_jm1 );

    /* Step 2.b.iii: Evolve Phi and Pi */
    if( j==1 ) phi.level_np1[0] = phi.level_nm1[0] + dt * alpha_Pi_over_a_jm1;
    phi.level_np1[j] = phi.level_n[j]   + dt * alpha_Pi_over_a_j;
    Phi.level_np1[j] = Phi.level_nm1[j] + dt * rhs_Phi;
    Pi.level_np1[j]  = Pi.level_nm1[j]  + dt * rhs_Pi;
  }
  /* Step 2.b.iv: Apply boundary conditions */
  {
    /* phi: Apply no-incoming radiation boundary condition */
    const REAL coeff1 = dt / x[0][J];
    const REAL coeff2 = dt * half_inv_dx0;
    phi.level_np1[J]  = phi.level_nm1[J] - coeff1 * phi.level_n[J] - coeff2 * ( 3.0 * phi.level_n[J] - 4.0 * phi.level_n[Jm1] + phi.level_n[Jm2] );
    /* Phi: Compute at outer boundary from phi */
    Phi.level_np1[J] = half_inv_dx0 * ( 3.0 * phi.level_np1[J] - 4.0 * phi.level_np1[Jm1] + phi.level_np1[Jm2] );
    /* Pi: Apply Neumann boundary conditions at the inner boundary */
    Pi.level_np1[0] = -Pi.level_nm1[0] + Pi.level_nm1[1] + Pi.level_np1[1];
    /* Pi: Compute from the evolution equation at the outer boundary */
    const REAL r_sqr_Jm2               = SQR(x[0][Jm2]);
    const REAL r_cbd_Jm2               = r_sqr_Jm2 * x[0][Jm2];
    const REAL r_sqr_Jm1               = SQR(x[0][Jm1]);
    const REAL r_sqr_J                 = SQR(x[0][J]);
    const REAL r_cbd_J                 = r_sqr_J   * x[0][J];
    const REAL alpha_over_a_Jm2        = alpha.level_n[Jm2]/a.level_n[Jm2];
    const REAL alpha_over_a_Jm1        = alpha.level_n[Jm1]/a.level_n[Jm1];
    const REAL alpha_over_a_J          = alpha.level_n[J]  /a.level_n[J];
    const REAL alpha_Phi_r2_over_a_Jm2 = r_sqr_Jm2 * alpha_over_a_Jm2 * Phi.level_n[Jm2];
    const REAL alpha_Phi_r2_over_a_Jm1 = r_sqr_Jm1 * alpha_over_a_Jm1 * Phi.level_n[Jm1];
    const REAL alpha_Phi_r2_over_a_J   = r_sqr_J   * alpha_over_a_J   * Phi.level_n[J];
    const REAL rhs_Pi                  = 3.0 * ( 3.0 * alpha_Phi_r2_over_a_J - 4.0 * alpha_Phi_r2_over_a_Jm1 + alpha_Phi_r2_over_a_Jm2 )/( r_cbd_J - r_cbd_Jm2 );
    Pi.level_np1[J]                    = Pi.level_nm1[J] + dt * rhs_Pi;
  }

  /* .---------------------------------.
   * | Step 2.c: Integrate a and alpha |
   * .---------------------------------.
   */
  LOOP(1,Nx0Total) {
    /* Step 2.c.i: Compute a */
    a.level_np1[j] = pointwise_Newton_method(j,grid,Phi.level_np1,Pi.level_np1,a.level_np1);

    /* Step 2.c.ii: Compute auxiliary quantities */
    const REAL b = a.level_np1[j] + a.level_np1[j-1];
    const REAL c = a.level_np1[j] - a.level_np1[j-1];
    const REAL midway_r = x[0][j] + x[0][j-1];
    const REAL d = ( 1.0 - 0.25 * SQR(b) )/midway_r - inv_dx0 * c / b;

    /* Step 2.c.iii: Compute alpha */
    alpha.level_np1[j] = alpha.level_np1[j-1]*( 1.0 - d*dx0 )/( 1.0 + d*dx0 );
  }
  /* Step 2.c.iv: Now rescale alpha */
  rescaling_of_the_lapse(grid,a.level_np1,alpha.level_np1);

  /* .---------------------------.
   * | Step 2.d: Update the time |
   * .---------------------------.
   */
  grid.t = dt;
  
}

/* General time step in Spherical coordinates */
void evolution::time_step_Spherical( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha ) {

  DECLARE_GRID_PARAMETERS;

  /* .---------------------------------------------------------.
   * | This function performs an integration to the time level |
   * | n+1 given the gridfunctions at the timelevels n and n-1 |
   * .---------------------------------------------------------.
   */
  const int J             = Nx0Total-1;
  const int Jm1           = J-1;
  const int Jm2           = J-2;
  const REAL dt_times_2   = 2.0 * dt;
  const REAL half_inv_dx0 = 0.5 * inv_dx0; // This is 1.0/(2.0*dx0)

  /* .--------------------------------------------------------------.
   * | Step 1: Apply inner boundary conditions to Phi, a, and alpha |
   * .--------------------------------------------------------------.
   */
  Phi.level_np1[0]   = 0.0;
  a.level_np1[0]     = 1.0;
  alpha.level_np1[0] = 1.0;

  /* .------------------------------.
   * | Step 2: Integrate Phi and Pi |
   * .------------------------------.
   */
  LOOP(1,Nx0Total-1) {
    /* Step 2.a: Compute useful auxiliary variables */
    const REAL r_sqr_jm1               = SQR(x[0][j-1]);
    const REAL r_cbd_jm1               = r_sqr_jm1 * x[0][j-1];
    const REAL r_sqr_jp1               = SQR(x[0][j+1]);
    const REAL r_cbd_jp1               = r_sqr_jp1 * x[0][j+1];
    const REAL alpha_over_a_jm1        = alpha.level_n[j-1]/a.level_n[j-1];
    const REAL alpha_over_a_j_nm1      = alpha.level_nm1[j]/a.level_nm1[j];
    const REAL alpha_over_a_j_n        = alpha.level_n[j]  /a.level_n[j];
    const REAL alpha_over_a_jp1        = alpha.level_n[j+1]/a.level_n[j+1];
    const REAL alpha_Pi_over_a_jm1     = alpha_over_a_jm1 * Pi.level_n[j-1];
    const REAL alpha_Pi_over_a_j_nm1   = alpha_over_a_j_nm1 * Pi.level_nm1[j];
    const REAL alpha_Pi_over_a_j_n     = alpha_over_a_j_n * Pi.level_n[j];
    const REAL alpha_Pi_over_a_jp1     = alpha_over_a_jp1 * Pi.level_n[j+1];
    const REAL alpha_Phi_r2_over_a_jm1 = r_sqr_jm1 * alpha_over_a_jm1 * Phi.level_n[j-1];
    const REAL alpha_Phi_r2_over_a_jp1 = r_sqr_jp1 * alpha_over_a_jp1 * Phi.level_n[j+1];

    /* Step 2.b: Compute the RHSs of Phi */
    const REAL rhs_phi = 1.5 * alpha_Pi_over_a_j_n - 0.5 * alpha_Pi_over_a_j_nm1;
    const REAL rhs_Phi = half_inv_dx0 * ( alpha_Pi_over_a_jp1 - alpha_Pi_over_a_jm1 );
    const REAL rhs_Pi  = 3.0 * ( alpha_Phi_r2_over_a_jp1 - alpha_Phi_r2_over_a_jm1 )/( r_cbd_jp1 - r_cbd_jm1 );

    /* Step 2.c: Evolve Phi and Pi */
    if( j==1 ) {
      const REAL tmp0 = 3.0 * alpha.level_n[0] * Pi.level_n[0] / a.level_n[0];
      const REAL tmp1 = alpha.level_nm1[0] * Pi.level_nm1[0] / a.level_nm1[0];
      phi.level_np1[0] = phi.level_n[0] + 0.5 * dt * ( tmp0 - tmp1 );
    }
    phi.level_np1[j] = phi.level_n[j]   + dt * rhs_phi;
    Phi.level_np1[j] = Phi.level_nm1[j] + dt_times_2 * rhs_Phi;
    Pi.level_np1[j]  = Pi.level_nm1[j]  + dt_times_2 * rhs_Pi;
  }
  /* Step 2.d: Apply boundary conditions */
  {
    /* phi: Apply no-incoming radiation boundary condition */
    const REAL inv_dt = 1.0/dt;
    const REAL coeff  = 1.0 / ( 3.0*inv_dt + 2.0/x[0][J] + 3.0*inv_dx0 );
    phi.level_np1[J]  = coeff * ( inv_dt * ( 4.0*phi.level_n[J] - phi.level_nm1[J] ) + inv_dx0 * ( 4.0*phi.level_np1[J-1] - phi.level_np1[J-2] ) );
    /* Phi: Compute at outer boundary from phi */
    Phi.level_np1[J] = half_inv_dx0 * ( 3.0 * phi.level_np1[J] - 4.0 * phi.level_np1[Jm1] + phi.level_np1[Jm2] );
    /* Pi: Apply Neumann boundary conditions at the inner boundary */
    Pi.level_np1[0] = -Pi.level_nm1[0] + Pi.level_nm1[1] + Pi.level_np1[1];
    /* Pi: Compute from the evolution equation at the outer boundary */
    const REAL r_sqr_Jm2               = SQR(x[0][Jm2]);
    const REAL r_cbd_Jm2               = r_sqr_Jm2 * x[0][Jm2];
    const REAL r_sqr_Jm1               = SQR(x[0][Jm1]);
    const REAL r_sqr_J                 = SQR(x[0][J]);
    const REAL r_cbd_J                 = r_sqr_J   * x[0][J];
    const REAL alpha_over_a_Jm2        = alpha.level_n[Jm2]/a.level_n[Jm2];
    const REAL alpha_over_a_Jm1        = alpha.level_n[Jm1]/a.level_n[Jm1];
    const REAL alpha_over_a_J          = alpha.level_n[J]  /a.level_n[J];
    const REAL alpha_Phi_r2_over_a_Jm2 = r_sqr_Jm2 * alpha_over_a_Jm2 * Phi.level_n[Jm2];
    const REAL alpha_Phi_r2_over_a_Jm1 = r_sqr_Jm1 * alpha_over_a_Jm1 * Phi.level_n[Jm1];
    const REAL alpha_Phi_r2_over_a_J   = r_sqr_J   * alpha_over_a_J   * Phi.level_n[J];
    const REAL rhs_Pi                  = 3.0 * ( 3.0 * alpha_Phi_r2_over_a_J - 4.0 * alpha_Phi_r2_over_a_Jm1 + alpha_Phi_r2_over_a_Jm2 )/( r_cbd_J - r_cbd_Jm2 );
    Pi.level_np1[J]                    = Pi.level_nm1[J] + dt_times_2 * rhs_Pi;
  }

  /* .-------------------------------.
   * | Step 3: Integrate a and alpha |
   * .-------------------------------.
   */
  LOOP(1,Nx0Total) {
    /* Step 3.a: Compute a */
    a.level_np1[j] = pointwise_Newton_method(j,grid,Phi.level_np1,Pi.level_np1,a.level_np1);

    /* Step 3.b: Compute auxiliary quantities */
    const REAL b = a.level_np1[j] + a.level_np1[j-1];
    const REAL c = a.level_np1[j] - a.level_np1[j-1];
    const REAL midway_r = x[0][j] + x[0][j-1];
    const REAL d = ( 1.0 - 0.25 * SQR(b) )/midway_r - inv_dx0 * c / b;

    /* Step 3.c: Compute alpha */
    alpha.level_np1[j] = alpha.level_np1[j-1]*( 1.0 - d*dx0 )/( 1.0 + d*dx0 );
  }
  /* Step 3.d: Now rescale alpha */
  rescaling_of_the_lapse(grid,a.level_np1,alpha.level_np1);

  /* .-------------------------.
   * | Step 4: Update the time |
   * .-------------------------.
   */
  grid.t += dt;
  
}

/* First time step in SinhSpherical coordinates */
void evolution::first_time_step_SinhSpherical( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha ) {

  DECLARE_GRID_PARAMETERS;

  /* .----------------------------------------------------.
   * | This function performs an integration to the time  |
   * | level 1 given the gridfunctions at the timelevel 0 |
   * .----------------------------------------------------.
   */
  const int J             = Nx0Total-1;
  const int Jm1           = J-1;
  const int Jm2           = J-2;
  const REAL half_dt      = 0.5 * dt;
  const REAL half_inv_dx0 = 0.5 * inv_dx0; // This is 1.0/(2.0*dx0)
  const REAL W_times_sinh_invW_over_A = sinhW/A_over_sinh_inv_W;

  /* .-------------------------------------.
   * | Step 1: Integrate from n=0 to n=1/2 |
   * .-------------------------------------.
   *
   * .----------------------------------------------------------------.
   * | Step 1.a: Apply inner boundary conditions to Phi, a, and alpha |
   * .----------------------------------------------------------------.
   */
  Phi.level_n[0]   = 0.0;
  a.level_n[0]     = 1.0;
  alpha.level_n[0] = 1.0;

  /* .--------------------------------.
   * | Step 1.b: Integrate Phi and Pi |
   * .--------------------------------.
   */
  LOOP(1,Nx0Total-1) {
    /* Step 1.b.i: Compute useful auxiliary variables */
    const REAL sinh_x0_inv_W_jm1          = sinh( x[0][j-1] * inv_sinhW );
    const REAL sinh_x0_inv_W_j            = sinh( x[0][j]   * inv_sinhW );
    const REAL sinh_x0_inv_W_jp1          = sinh( x[0][j+1] * inv_sinhW );
    const REAL cosh_x0_inv_W_j            = cosh( x[0][j]   * inv_sinhW );
    const REAL alpha_over_a_jm1           = alpha.level_nm1[j-1]/a.level_nm1[j-1];
    const REAL alpha_over_a_j             = alpha.level_nm1[j]  /a.level_nm1[j];
    const REAL alpha_over_a_jp1           = alpha.level_nm1[j+1]/a.level_nm1[j+1];
    const REAL alpha_Pi_over_a_jm1        = alpha_over_a_jm1 * Pi.level_nm1[j-1];
    const REAL alpha_Pi_over_a_j          = alpha_over_a_j   * Pi.level_nm1[j];
    const REAL alpha_Pi_over_a_jp1        = alpha_over_a_jp1 * Pi.level_nm1[j+1];
    const REAL alpha_Phi_shx02_over_a_jm1 = SQR(sinh_x0_inv_W_jm1) * alpha_over_a_jm1 * Phi.level_nm1[j-1];
    const REAL alpha_Phi_shx02_over_a_jp1 = SQR(sinh_x0_inv_W_jp1) * alpha_over_a_jp1 * Phi.level_nm1[j+1];

    /* Step 1.b.ii: Compute the RHSs of Phi */
    const REAL rhs_phi = alpha_Pi_over_a_j;
    const REAL tmp0    = half_inv_dx0 * ( W_times_sinh_invW_over_A / cosh_x0_inv_W_j );
    const REAL rhs_Phi = tmp0 * ( alpha_Pi_over_a_jp1 - alpha_Pi_over_a_jm1 );
    const REAL rhs_Pi  = tmp0 / SQR(sinh_x0_inv_W_j) * ( alpha_Phi_shx02_over_a_jp1 - alpha_Phi_shx02_over_a_jm1 );

    /* Step 1.b.iii: Evolve Phi and Pi */
    if( j==1 ) {
      const REAL tmp0 = 3.0 * alpha.level_nm1[0] * Pi.level_nm1[0] / a.level_nm1[0];
      const REAL tmp1 = alpha.level_nm1[0] * Pi.level_nm1[0] / a.level_nm1[0];
      phi.level_n[0]  = phi.level_nm1[0] + 0.5 * dt * ( tmp0 - tmp1 );
    }
    phi.level_n[j] = phi.level_nm1[j] + half_dt * rhs_phi;
    Phi.level_n[j] = Phi.level_nm1[j] + half_dt * rhs_Phi;
    Pi.level_n[j]  = Pi.level_nm1[j]  + half_dt * rhs_Pi;
  }
  /* Step 1.b.iv: Apply boundary conditions */
  {
    /* phi: Apply no-incoming radiation boundary condition */
    const REAL tmp0  = 1.0 - dt / ( A_over_sinh_inv_W * sinh( x[0][J] * inv_sinhW ) );
    const REAL tmp1  = sinhW / ( A_over_sinh_inv_W * cosh( x[0][J] * inv_sinhW ) );
    const REAL tmp2  = half_inv_dx0 * ( 3.0 * phi.level_nm1[J] - 4.0 * phi.level_nm1[Jm1] + phi.level_nm1[Jm2] );
    phi.level_n[J] = tmp0 * phi.level_nm1[J] - tmp1 * tmp2;
    /* Phi: Compute at outer boundary from phi */
    Phi.level_n[J] = tmp1 / SQR(sinh(x[0][J]*inv_sinhW)) * half_inv_dx0 * ( 3.0 * phi.level_n[J] - 4.0 * phi.level_n[Jm1] + phi.level_n[Jm2] );
    /* Pi: Apply Neumann boundary conditions at the inner boundary */
    Pi.level_n[0] = Pi.level_n[1];
    /* Pi: Compute from the evolution equation at the outer boundary */
    const REAL sinh_x0_inv_W_Jm2          = sinh( x[0][J-2] * inv_sinhW );
    const REAL sinh_x0_inv_W_Jm1          = sinh( x[0][J-1] * inv_sinhW );
    const REAL sinh_x0_inv_W_J            = sinh( x[0][J]   * inv_sinhW );
    const REAL cosh_x0_inv_W_J            = cosh( x[0][J]   * inv_sinhW );
    const REAL alpha_over_a_Jm2           = alpha.level_nm1[Jm2]/a.level_nm1[Jm2];
    const REAL alpha_over_a_Jm1           = alpha.level_nm1[Jm1]/a.level_nm1[Jm1];
    const REAL alpha_over_a_J             = alpha.level_nm1[J]  /a.level_nm1[J];
    const REAL alpha_Phi_shx02_over_a_Jm2 = SQR(sinh_x0_inv_W_Jm2) * alpha_over_a_Jm2 * Phi.level_nm1[Jm2];
    const REAL alpha_Phi_shx02_over_a_Jm1 = SQR(sinh_x0_inv_W_Jm1) * alpha_over_a_Jm1 * Phi.level_nm1[Jm1];
    const REAL alpha_Phi_shx02_over_a_J   = SQR(sinh_x0_inv_W_J  ) * alpha_over_a_J   * Phi.level_nm1[J];
    const REAL coeff_Pi                   = half_inv_dx0 * ( W_times_sinh_invW_over_A / ( SQR(sinh_x0_inv_W_J) * cosh_x0_inv_W_J ) ) ;
    const REAL rhs_Pi                     = coeff_Pi * ( 3.0 * alpha_Phi_shx02_over_a_J - 4.0 * alpha_Phi_shx02_over_a_Jm1 + alpha_Phi_shx02_over_a_Jm2 );
    Pi.level_n[J]                         = Pi.level_nm1[J] + half_dt * rhs_Pi;
  }

  /* .---------------------------------.
   * | Step 1.c: Integrate a and alpha |
   * .---------------------------------.
   */
  LOOP(1,Nx0Total) {
    /* Step 1.c.i: Compute a */
    a.level_n[j] = pointwise_Newton_method(j,grid,Phi.level_n,Pi.level_n,a.level_n);

    /* Step 1.c.ii: Compute auxiliary quantities */
    const REAL b = a.level_n[j] + a.level_n[j-1];
    const REAL c = a.level_n[j] - a.level_n[j-1];
    const REAL midway_r = sinhW * tanh( inv_sinhW * 0.5 * (x[0][j] + x[0][j-1]) );
    const REAL d = ( 1.0 - 0.25 * SQR(b) )/( 2.0 * midway_r ) - inv_dx0 * c / b;

    /* Step 1.c.iii: Compute alpha */
    alpha.level_n[j] = alpha.level_n[j-1]*( 1.0 - d*dx0 )/( 1.0 + d*dx0 );
  }
  /* Step 1.c.iv: Now rescale alpha */
  rescaling_of_the_lapse(grid,a.level_n,alpha.level_n);

  /* .---------------------------.
   * | Step 1.d: Update the time |
   * .---------------------------.
   */
  grid.t += half_dt;

  /* .---------------------------------------------.
   * | Step 2: Integrate from n=0 and n=1/2 to n=1 |
   * .---------------------------------------------.
   *
   * .----------------------------------------------------------------.
   * | Step 2.a: Apply inner boundary conditions to Phi, a, and alpha |
   * .----------------------------------------------------------------.
   */
  Phi.level_np1[0]   = 0.0;
  a.level_np1[0]     = 1.0;
  alpha.level_np1[0] = 1.0;

  /* .--------------------------------.
   * | Step 2.b: Integrate Phi and Pi |
   * .--------------------------------.
   */
  LOOP(1,Nx0Total-1) {
    /* Step 2.b.i: Compute useful auxiliary variables */
    const REAL sinh_x0_inv_W_jm1          = sinh( x[0][j-1] * inv_sinhW );
    const REAL sinh_x0_inv_W_j            = sinh( x[0][j]   * inv_sinhW );
    const REAL sinh_x0_inv_W_jp1          = sinh( x[0][j+1] * inv_sinhW );
    const REAL cosh_x0_inv_W_j            = cosh( x[0][j]   * inv_sinhW );
    const REAL alpha_over_a_jm1           = alpha.level_n[j-1]/a.level_n[j-1];
    const REAL alpha_over_a_j_nm1         = alpha.level_nm1[j]/a.level_nm1[j];
    const REAL alpha_over_a_j_n           = alpha.level_n[j]  /a.level_n[j];
    const REAL alpha_over_a_jp1           = alpha.level_n[j+1]/a.level_n[j+1];
    const REAL alpha_Pi_over_a_jm1        = alpha_over_a_jm1 * Pi.level_n[j-1];
    const REAL alpha_Pi_over_a_j_nm1      = alpha_over_a_j_nm1 * Pi.level_nm1[j];
    const REAL alpha_Pi_over_a_j_n        = alpha_over_a_j_n * Pi.level_n[j];
    const REAL alpha_Pi_over_a_jp1        = alpha_over_a_jp1 * Pi.level_n[j+1];
    const REAL alpha_Phi_shx02_over_a_jm1 = SQR(sinh_x0_inv_W_jm1) * alpha_over_a_jm1 * Phi.level_n[j-1];
    const REAL alpha_Phi_shx02_over_a_jp1 = SQR(sinh_x0_inv_W_jp1) * alpha_over_a_jp1 * Phi.level_n[j+1];

    /* Step 2.b.ii: Compute the RHSs of Phi */
    const REAL rhs_phi = 1.5 * alpha_Pi_over_a_j_n - 0.5 * alpha_Pi_over_a_j_nm1;
    const REAL tmp0    = half_inv_dx0 * ( W_times_sinh_invW_over_A / cosh_x0_inv_W_j );
    const REAL rhs_Phi = tmp0 * ( alpha_Pi_over_a_jp1 - alpha_Pi_over_a_jm1 );
    const REAL rhs_Pi  = tmp0 / SQR(sinh_x0_inv_W_j) * ( alpha_Phi_shx02_over_a_jp1 - alpha_Phi_shx02_over_a_jm1 );

    /* Step 2.b.iii: Evolve Phi and Pi */
    if( j==1 ) {
      const REAL tmp0 = 3.0 * alpha.level_n[0] * Pi.level_n[0] / a.level_n[0];
      const REAL tmp1 = alpha.level_nm1[0] * Pi.level_nm1[0] / a.level_nm1[0];
      phi.level_np1[0] = phi.level_n[0] + 0.5 * dt * ( tmp0 - tmp1 );
    }
    phi.level_np1[j] = phi.level_n[j]   + dt * rhs_phi;
    Phi.level_np1[j] = Phi.level_nm1[j] + dt * rhs_Phi;
    Pi.level_np1[j]  = Pi.level_nm1[j]  + dt * rhs_Pi;
  }
  /* Step 2.b.iv: Apply boundary conditions */
  {
    /* phi: Apply no-incoming radiation boundary condition */
    const REAL inv_dt  = 1.0/dt;
    const REAL tmp0    = sinh( x[0][J] * inv_sinhW );
    const REAL tmp1    = cosh( x[0][J] * inv_sinhW );
    const REAL tmp2    = 3.0 * inv_dt;
    const REAL tmp3    = 2.0 / ( A_over_sinh_inv_W * tmp0 );
    const REAL tmp4    = 3.0 * sinhW * inv_dx0 / ( A_over_sinh_inv_W * tmp1 );
    const REAL tmp5    = 1.0 / ( tmp2 + tmp3 + tmp4 );
    const REAL tmp6    = inv_dt * ( 4.0 * phi.level_n[J] - phi.level_nm1[J] );
    const REAL tmp7    = sinhW / ( A_over_sinh_inv_W * tmp1 );
    const REAL tmp8    = inv_dx0 * ( 4.0 * phi.level_np1[Jm1] - phi.level_np1[Jm2] );
    const REAL tmp9    = tmp7 * tmp8;
    const REAL rhs_phi = tmp5 * ( tmp6 + tmp9 );
    phi.level_np1[J]   = rhs_phi;
    /* Phi: Compute at outer boundary from phi */
    Phi.level_np1[J] = half_inv_dx0 * tmp7 * ( 3.0 * phi.level_np1[J] - 4.0 * phi.level_np1[Jm1] + phi.level_np1[Jm2] );
    /* Pi: Apply Neumann boundary conditions at the inner boundary */
    Pi.level_np1[0] = -Pi.level_nm1[0] + Pi.level_nm1[1] + Pi.level_np1[1];
    /* Pi: Compute from the evolution equation at the outer boundary */
    const REAL sinh_x0_inv_W_Jm2          = sinh( x[0][J-2] * inv_sinhW );
    const REAL sinh_x0_inv_W_Jm1          = sinh( x[0][J-1] * inv_sinhW );
    const REAL sinh_x0_inv_W_J            = sinh( x[0][J]   * inv_sinhW );
    const REAL cosh_x0_inv_W_J            = cosh( x[0][J]   * inv_sinhW );
    const REAL alpha_over_a_Jm2           = alpha.level_n[Jm2]/a.level_n[Jm2];
    const REAL alpha_over_a_Jm1           = alpha.level_n[Jm1]/a.level_n[Jm1];
    const REAL alpha_over_a_J             = alpha.level_n[J]  /a.level_n[J];
    const REAL alpha_Phi_shx02_over_a_Jm2 = SQR(sinh_x0_inv_W_Jm2) * alpha_over_a_Jm2 * Phi.level_n[Jm2];
    const REAL alpha_Phi_shx02_over_a_Jm1 = SQR(sinh_x0_inv_W_Jm1) * alpha_over_a_Jm1 * Phi.level_n[Jm1];
    const REAL alpha_Phi_shx02_over_a_J   = SQR(sinh_x0_inv_W_J  ) * alpha_over_a_J   * Phi.level_n[J];
    const REAL coeff_Pi                   = half_inv_dx0 * ( W_times_sinh_invW_over_A / ( SQR(sinh_x0_inv_W_J) * cosh_x0_inv_W_J ) ) ;
    const REAL rhs_Pi                     = coeff_Pi * ( 3.0 * alpha_Phi_shx02_over_a_J - 4.0 * alpha_Phi_shx02_over_a_Jm1 + alpha_Phi_shx02_over_a_Jm2 );
    Pi.level_np1[J]                       = Pi.level_nm1[J] + dt * rhs_Pi;
  }

  /* .---------------------------------.
   * | Step 2.c: Integrate a and alpha |
   * .---------------------------------.
   */
  LOOP(1,Nx0Total) {
    /* Step 2.c.i: Compute a */
    a.level_np1[j] = pointwise_Newton_method(j,grid,Phi.level_np1,Pi.level_np1,a.level_np1);

    /* Step 2.c.ii: Compute auxiliary quantities */
    const REAL b = a.level_np1[j] + a.level_np1[j-1];
    const REAL c = a.level_np1[j] - a.level_np1[j-1];
    const REAL midway_r = sinhW * tanh( inv_sinhW * 0.5 * (x[0][j] + x[0][j-1]) );
    const REAL d = ( 1.0 - 0.25 * SQR(b) )/( 2.0 * midway_r ) - inv_dx0 * c / b;

    /* Step 2.c.iii: Compute alpha */
    alpha.level_np1[j] = alpha.level_np1[j-1]*( 1.0 - d*dx0 )/( 1.0 + d*dx0 );
  }
  /* Step 2.c.iv: Now rescale alpha */
  rescaling_of_the_lapse(grid,a.level_np1,alpha.level_np1);

  /* .---------------------------.
   * | Step 2.d: Update the time |
   * .---------------------------.
   */
  grid.t += half_dt;

}

/* Generic time step in SinhSpherical coordinates */
void evolution::time_step_SinhSpherical( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha ) {

  DECLARE_GRID_PARAMETERS;

  /* .---------------------------------------------------------.
   * | This function performs an integration to the time level |
   * | n+1 given the gridfunctions at the timelevels n and n-1 |
   * .---------------------------------------------------------.
   */
  const int J             = Nx0Total-1;
  const int Jm1           = J-1;
  const int Jm2           = J-2;
  const REAL dt_times_2   = 2.0 * dt;
  const REAL half_inv_dx0 = 0.5 * inv_dx0; // This is 1.0/(2.0*dx0)
  const REAL W_times_sinh_invW_over_A = sinhW/A_over_sinh_inv_W;

  /* .--------------------------------------------------------------.
   * | Step 1: Apply inner boundary conditions to Phi, a, and alpha |
   * .--------------------------------------------------------------.
   */
  Phi.level_np1[0]   = 0.0;
  a.level_np1[0]     = 1.0;
  alpha.level_np1[0] = 1.0;

  /* .------------------------------.
   * | Step 2: Integrate Phi and Pi |
   * .------------------------------.
   */
  LOOP(1,Nx0Total-1) {
    /* Step 2.a: Compute useful auxiliary variables */
    const REAL sinh_x0_inv_W_jm1          = sinh( x[0][j-1] * inv_sinhW );
    const REAL sinh_x0_inv_W_j            = sinh( x[0][j]   * inv_sinhW );
    const REAL sinh_x0_inv_W_jp1          = sinh( x[0][j+1] * inv_sinhW );
    const REAL cosh_x0_inv_W_j            = cosh( x[0][j]   * inv_sinhW );
    const REAL alpha_over_a_jm1           = alpha.level_n[j-1]/a.level_n[j-1];
    const REAL alpha_over_a_j_nm1         = alpha.level_nm1[j]/a.level_nm1[j];
    const REAL alpha_over_a_j_n           = alpha.level_n[j]  /a.level_n[j];
    const REAL alpha_over_a_jp1           = alpha.level_n[j+1]/a.level_n[j+1];
    const REAL alpha_Pi_over_a_jm1        = alpha_over_a_jm1 * Pi.level_n[j-1];
    const REAL alpha_Pi_over_a_j_nm1      = alpha_over_a_j_nm1 * Pi.level_nm1[j];
    const REAL alpha_Pi_over_a_j_n        = alpha_over_a_j_n * Pi.level_n[j];
    const REAL alpha_Pi_over_a_jp1        = alpha_over_a_jp1 * Pi.level_n[j+1];
    const REAL alpha_Phi_shx02_over_a_jm1 = SQR(sinh_x0_inv_W_jm1) * alpha_over_a_jm1 * Phi.level_n[j-1];
    const REAL alpha_Phi_shx02_over_a_jp1 = SQR(sinh_x0_inv_W_jp1) * alpha_over_a_jp1 * Phi.level_n[j+1];

    /* Step 2.b: Compute the RHSs of Phi */
    const REAL rhs_phi = 1.5 * alpha_Pi_over_a_j_n - 0.5 * alpha_Pi_over_a_j_nm1;
    const REAL tmp0    = half_inv_dx0 * ( W_times_sinh_invW_over_A / cosh_x0_inv_W_j );
    const REAL rhs_Phi = tmp0 * ( alpha_Pi_over_a_jp1 - alpha_Pi_over_a_jm1 );
    const REAL rhs_Pi  = tmp0 / SQR(sinh_x0_inv_W_j) * ( alpha_Phi_shx02_over_a_jp1 - alpha_Phi_shx02_over_a_jm1 );

    /* Step 2.c: Evolve Phi and Pi */
    if( j==1 ) {
      const REAL tmp0 = 3.0 * alpha.level_n[0] * Pi.level_n[0] / a.level_n[0];
      const REAL tmp1 = alpha.level_nm1[0] * Pi.level_nm1[0] / a.level_nm1[0];
      phi.level_np1[0] = phi.level_n[0] + 0.5 * dt * ( tmp0 - tmp1 );
    }
    phi.level_np1[j] = phi.level_n[j]   + dt * rhs_phi;
    Phi.level_np1[j] = Phi.level_nm1[j] + dt_times_2 * rhs_Phi;
    Pi.level_np1[j]  = Pi.level_nm1[j]  + dt_times_2 * rhs_Pi;
  }
  /* Step 2.d: Apply boundary conditions */
  {
    /* phi: Apply no-incoming radiation boundary condition */
    const REAL inv_dt  = 1.0/dt;
    const REAL tmp0    = sinh( x[0][J] * inv_sinhW );
    const REAL tmp1    = cosh( x[0][J] * inv_sinhW );
    const REAL tmp2    = 3.0 * inv_dt;
    const REAL tmp3    = 2.0 / ( A_over_sinh_inv_W * tmp0 );
    const REAL tmp4    = 3.0 * sinhW * inv_dx0 / ( A_over_sinh_inv_W * tmp1 );
    const REAL tmp5    = 1.0 / ( tmp2 + tmp3 + tmp4 );
    const REAL tmp6    = inv_dt * ( 4.0 * phi.level_n[J] - phi.level_nm1[J] );
    const REAL tmp7    = sinhW / ( A_over_sinh_inv_W * tmp1 );
    const REAL tmp8    = inv_dx0 * ( 4.0 * phi.level_np1[Jm1] - phi.level_np1[Jm2] );
    const REAL tmp9    = tmp7 * tmp8;
    const REAL rhs_phi = tmp5 * ( tmp6 + tmp9 );
    phi.level_np1[J]   = rhs_phi;
    /* Phi: Compute at outer boundary from phi */
    Phi.level_np1[J] = half_inv_dx0 * tmp7 * ( 3.0 * phi.level_np1[J] - 4.0 * phi.level_np1[Jm1] + phi.level_np1[Jm2] );
    /* Pi: Apply Neumann boundary conditions at the inner boundary */
    Pi.level_np1[0] = -Pi.level_nm1[0] + Pi.level_nm1[1] + Pi.level_np1[1];
    /* Pi: Compute from the evolution equation at the outer boundary */
    const REAL sinh_x0_inv_W_Jm2          = sinh( x[0][J-2] * inv_sinhW );
    const REAL sinh_x0_inv_W_Jm1          = sinh( x[0][J-1] * inv_sinhW );
    const REAL sinh_x0_inv_W_J            = sinh( x[0][J]   * inv_sinhW );
    const REAL cosh_x0_inv_W_J            = cosh( x[0][J]   * inv_sinhW );
    const REAL alpha_over_a_Jm2           = alpha.level_n[Jm2]/a.level_n[Jm2];
    const REAL alpha_over_a_Jm1           = alpha.level_n[Jm1]/a.level_n[Jm1];
    const REAL alpha_over_a_J             = alpha.level_n[J]  /a.level_n[J];
    const REAL alpha_Phi_shx02_over_a_Jm2 = SQR(sinh_x0_inv_W_Jm2) * alpha_over_a_Jm2 * Phi.level_n[Jm2];
    const REAL alpha_Phi_shx02_over_a_Jm1 = SQR(sinh_x0_inv_W_Jm1) * alpha_over_a_Jm1 * Phi.level_n[Jm1];
    const REAL alpha_Phi_shx02_over_a_J   = SQR(sinh_x0_inv_W_J  ) * alpha_over_a_J   * Phi.level_n[J];
    const REAL coeff_Pi                   = half_inv_dx0 * ( W_times_sinh_invW_over_A / ( SQR(sinh_x0_inv_W_J) * cosh_x0_inv_W_J ) ) ;
    const REAL rhs_Pi                     = coeff_Pi * ( 3.0 * alpha_Phi_shx02_over_a_J - 4.0 * alpha_Phi_shx02_over_a_Jm1 + alpha_Phi_shx02_over_a_Jm2 );
    Pi.level_np1[J]                       = Pi.level_nm1[J] + dt_times_2 * rhs_Pi;
  }

  /* .-------------------------------.
   * | Step 3: Integrate a and alpha |
   * .-------------------------------.
   */
  LOOP(1,Nx0Total) {
    /* Step 3.a: Compute a */
    a.level_np1[j] = pointwise_Newton_method(j,grid,Phi.level_np1,Pi.level_np1,a.level_np1);

    /* Step 3.b: Compute auxiliary quantities */
    const REAL b = a.level_np1[j] + a.level_np1[j-1];
    const REAL c = a.level_np1[j] - a.level_np1[j-1];
    const REAL midway_r = sinhW * tanh( inv_sinhW * 0.5 * (x[0][j] + x[0][j-1]) );
    const REAL d = ( 1.0 - 0.25 * SQR(b) )/( 2.0 * midway_r ) - inv_dx0 * c / b;

    /* Step 3.c: Compute alpha */
    alpha.level_np1[j] = alpha.level_np1[j-1]*( 1.0 - d*dx0 )/( 1.0 + d*dx0 );
  }
  /* Step 3.d: Now rescale alpha */
  rescaling_of_the_lapse(grid,a.level_np1,alpha.level_np1);

  /* .-------------------------.
   * | Step 4: Update the time |
   * .-------------------------.
   */
  grid.t += dt;
  
}

void evolution::output_mass_aspect_ratio( const int which_level, const int n, const grid::parameters grid, const gridfunction a ) {

  DECLARE_GRID_PARAMETERS;

  ofstream outfile;
  outfile.open("out/mass_"+string(to_string(Nt).length() - to_string(n).length(),'0')+to_string(n)+".dat");
  outfile.precision(15);
  LOOP(0,Nx0Total) {
    const REAL r_local = r_ito_x0[j];
    const REAL a_local = (which_level == -1) * a.level_nm1[j] + (which_level == 0) * a.level_n[j] + (which_level == 1) * a.level_np1[j];
    outfile << scientific << x[0][j] << " " << r_ito_x0[j]  << " " << 0.5 * r_local * ( 1.0 - 1.0/SQR(a_local) ) << endl;
    
  }
  outfile.close();
}
