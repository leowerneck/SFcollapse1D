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
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "macros.hpp"
#include "grid.hpp"
#include "gridfunction.hpp"
#include "evolution.hpp"
#include "utilities.hpp"

using namespace std;

/* Function to set the initial condition for all gridfunctions: phi, Phi, Pi, a, and alpha */
void evolution::initial_condition( grid::parameters grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha ) {

  DECLARE_GRID_PARAMETERS;

  LOOP(0,Nx0Total) {
    /* Set some useful auxiliary variables */
    const REAL r = r_ito_x0[j];
    const REAL factor = (r-R0)/SQR(DELTA);
    const REAL expfactor = (r-R0)*factor;
    const REAL exp_rmr0_over_deltasqrd = exp(-expfactor);

    /* Set the initial condition for phi and Pi */
    phi.level_nm1[j] = PHI0*exp_rmr0_over_deltasqrd;
    Pi.level_nm1[j]  = 0.0;

    if( j>0 ) {
      /* Set the initial condition for Phi */
      Phi.level_nm1[j] = -2.0*factor*PHI0*exp_rmr0_over_deltasqrd;
      
      /* Compute a */
      a.level_nm1[j] = evolution::pointwise_solution_of_the_Hamiltonian_constraint( j, grid, Phi.level_nm1, Pi.level_nm1, a.level_nm1 );

      /* Compute alpha */
      alpha.level_nm1[j] = evolution::pointwise_solution_of_the_polar_slicing_condition( j, grid, a.level_nm1, alpha.level_nm1 );

    }
    else {
      /* If at inner boundary, impose Phi = 0 */
      Phi.level_nm1[j] = 0.0;
      /* If at inner boundary, impose a = 1 */
      a.level_nm1[j] = 1.0;
      /* If at inner boundary, impose alpha = 1 */
      alpha.level_nm1[j] = 1.0;
    }
  }

  /* Now rescale alpha */
  evolution::rescaling_of_the_lapse(grid,a.level_nm1,alpha.level_nm1);

}

/* Function to step phi, Phi, and Pi forward in time */
void evolution::time_step_scalarfield_gridfunctions( const int n, const grid::parameters grid,
						     const vector<REAL> phi_n, const vector<REAL> Phi_n   , const vector<REAL> Pi_n   , const vector<REAL> a_n    , const vector<REAL> alpha_n  ,
						     const vector<REAL> Phi_nm1 , const vector<REAL> Pi_nm1 , const vector<REAL> a_nm1  , const vector<REAL> alpha_nm1, 
						           vector<REAL> &Phi_np1,       vector<REAL> &Pi_np1,       vector<REAL> &phi_np1 ) {

  DECLARE_GRID_PARAMETERS;

  const REAL Phi_Pi_dt_coeff = ( n==0 ? 0.5 * dt : 0.0 ) + ( n==1 ? dt : 0.0 ) + ( n>1 ? 2.0 * dt : 0.0 );
  const REAL phi_dt_coeff    = ( n==0 ? 0.5 * dt : 0.0 ) + ( n>=1 ? dt : 0.0 );

  //#pragma omp parallel for
  LOOP(1,Nx0Total-1) {
    /* Compute useful auxiliary variables to all RHSs */
    const REAL alpha_over_a_jm1      = alpha_n[j-1]/a_n[j-1];
    const REAL alpha_over_a_j_nm1    = alpha_nm1[j]/a_nm1[j];
    const REAL alpha_over_a_j_n      = alpha_n[j]  /a_n[j];
    const REAL alpha_over_a_jp1      = alpha_n[j+1]/a_n[j+1];
    const REAL alpha_Pi_over_a_jm1   = alpha_over_a_jm1 * Pi_n[j-1];
    const REAL alpha_Pi_over_a_j_nm1 = alpha_over_a_j_nm1 * Pi_nm1[j];
    const REAL alpha_Pi_over_a_j_n   = alpha_over_a_j_n * Pi_n[j];
    const REAL alpha_Pi_over_a_jp1   = alpha_over_a_jp1 * Pi_n[j+1];

    /* Compute the RHS of phi */
    if( j==1 ) {
      const REAL tmp0    = 1.5 * alpha_n[0]   * Pi_n[0]   / a_n[0];
      const REAL tmp1    = 0.5 * alpha_nm1[0] * Pi_nm1[0] / a_nm1[0];
      const REAL rhs_phi = tmp0 - tmp1;
      phi_np1[0]         = phi_n[0] + phi_dt_coeff * rhs_phi;
    }
    const REAL rhs_phi = 1.5 * alpha_Pi_over_a_j_n - 0.5 * alpha_Pi_over_a_j_nm1;
    
    phi_np1[j] = phi_n[j] + phi_dt_coeff * rhs_phi;

#if( COORD_SYSTEM == SPHERICAL )

    /* Compute useful auxiliary variables to the RHSs of Phi and Pi */
    const REAL r_sqr_jm1               = SQR(x[0][j-1]);
    const REAL r_cbd_jm1               = r_sqr_jm1 * x[0][j-1];
    const REAL r_sqr_jp1               = SQR(x[0][j+1]);
    const REAL r_cbd_jp1               = r_sqr_jp1 * x[0][j+1];
    const REAL alpha_Phi_r2_over_a_jm1 = r_sqr_jm1 * alpha_over_a_jm1 * Phi_n[j-1];
    const REAL alpha_Phi_r2_over_a_jp1 = r_sqr_jp1 * alpha_over_a_jp1 * Phi_n[j+1];

    /* Compute the RHSs of Phi and Pi */
    const REAL rhs_Phi = 0.5*inv_dx0 * ( alpha_Pi_over_a_jp1 - alpha_Pi_over_a_jm1 );
    const REAL rhs_Pi  = 3.0 * ( alpha_Phi_r2_over_a_jp1 - alpha_Phi_r2_over_a_jm1 )/( r_cbd_jp1 - r_cbd_jm1 );

#elif( COORD_SYSTEM == SINH_SPHERICAL )

    /* Compute useful auxiliary variables to the RHSs of Phi and Pi */
    const REAL sinh_x0_inv_W_jm1          = sinh( x[0][j-1] * inv_sinhW );
    const REAL sinh_x0_inv_W_j            = sinh( x[0][j]   * inv_sinhW );
    const REAL sinh_x0_inv_W_jp1          = sinh( x[0][j+1] * inv_sinhW );
    const REAL cosh_x0_inv_W_j            = cosh( x[0][j]   * inv_sinhW );
    const REAL alpha_Phi_shx02_over_a_jm1 = SQR(sinh_x0_inv_W_jm1) * alpha_over_a_jm1 * Phi_n[j-1];
    const REAL alpha_Phi_shx02_over_a_jp1 = SQR(sinh_x0_inv_W_jp1) * alpha_over_a_jp1 * Phi_n[j+1];
    const REAL coefficient = 0.5 * inv_dx0 * ( sinhW * sinh_inv_W / sinhA ) / cosh_x0_inv_W_j;

    /* Compute the RHSs of Phi and Pi */
    const REAL rhs_Phi = coefficient * ( alpha_Pi_over_a_jp1 - alpha_Pi_over_a_jm1 );
    const REAL rhs_Pi  = coefficient / SQR(sinh_x0_inv_W_j) * ( alpha_Phi_shx02_over_a_jp1 - alpha_Phi_shx02_over_a_jm1 );

#else
    utilities::SFcollapse1D_error(COORD_SYSTEM_ERROR);
#endif

    Phi_np1[j] = Phi_nm1[j] + Phi_Pi_dt_coeff * rhs_Phi;
    Pi_np1[j]  = Pi_nm1[j]  + Phi_Pi_dt_coeff * rhs_Pi;

  }

}

/* Function to apply outgoing radiation boundary conditions to phi, Phi, and Pi */
void evolution::apply_outgoing_radiation_bdry_cond( const int n, grid::parameters grid,
						    const vector<REAL> phi_nm1, const vector<REAL> Pi_nm1,
						    const vector<REAL> phi_n  , const vector<REAL> Phi_n, const vector<REAL> a_n, const vector<REAL> alpha_n,
						    vector<REAL> &phi_np1, vector<REAL> &Phi_np1, vector<REAL> &Pi_np1 ) {

  DECLARE_GRID_PARAMETERS;

  /* Declare useful auxiliary variables */
  const int J     = Nx0Total;
  const REAL tmp0 = - phi_n[J] / r_ito_x0[J];
  const REAL tmp1 = 0.5 * inv_dx0;
  const REAL tmp2 = - tmp1 * ( 3.0*phi_n[J] - 4.0*phi_n[J-1] + phi_n[J-2] );
  
#if( COORD_SYSTEM == SPHERICAL )

  /* RHS of phi in Spherical coordinates */
  const REAL rhs_phi = tmp0 + tmp2;
  
#elif( COORD_SYSTEM == SINH_SPHERICAL )

  /* RHS of phi in SinhSpherical coordinates */
  const REAL tmp3 = ( sinhW / A_over_sinh_inv_W ) / cosh( x[0][J] * inv_sinhW );
  const REAL rhs_phi = tmp0 + tmp3 * tmp2;
  
#else
  utilities::SFcollapse1D_error(COORD_SYSTEM_ERROR);
#endif

  /* Compute phi at the outer boundary */
  const REAL phi_Pi_coeff = ( n==0 ? 0.5 * dt : 0.0 ) + ( n==1 ? dt : 0.0 ) + ( n>1 ? 2.0 * dt : 0.0 );
  phi_np1[J] = phi_nm1[J] + phi_Pi_coeff * rhs_phi;

#if( COORD_SYSTEM == SPHERICAL )

  /* RHS of Phi in Spherical coordinates */
  const REAL rhs_Phi = tmp1 * ( 3.0*phi_np1[J] - 4.0*phi_np1[J-1] + phi_np1[J-2] );
  
#elif( COORD_SYSTEM == SINH_SPHERICAL )

  /* RHS of Phi in SinhSpherical coordinates */
  const REAL rhs_Phi = tmp1 * tmp3 * ( 3.0*phi_np1[J] - 4.0*phi_np1[J-1] + phi_np1[J-2] );
  
#else
  utilities::SFcollapse1D_error(COORD_SYSTEM_ERROR);
#endif

  /* Compute Phi at the outer boundary */
  Phi_np1[J] = rhs_Phi;

#if( COORD_SYSTEM == SPHERICAL )

  /* RHS of Pi in Spherical coordinates */
  const REAL r_sqd_Jm2 = SQR(x[0][J-2]);
  const REAL r_sqd_Jm1 = SQR(x[0][J-1]);
  const REAL r_sqd_J   = SQR(x[0][J]);
  const REAL r_cbd_Jm2 = r_sqd_Jm2 * x[0][J-2];
  const REAL r_cbd_J   = r_sqd_J   * x[0][J];
  const REAL coeff     = 3.0 / ( r_cbd_J - r_cbd_Jm2 );
  const REAL term1     = 3.0 * r_sqd_J   * ( alpha_n[J]   / a_n[J]   ) * Phi_n[J];
  const REAL term2     = 4.0 * r_sqd_Jm1 * ( alpha_n[J-1] / a_n[J-1] ) * Phi_n[J-1];
  const REAL term3     =       r_sqd_Jm2 * ( alpha_n[J-2] / a_n[J-2] ) * Phi_n[J-2];
  const REAL rhs_Pi    = coeff * ( term1 - term2 + term3 );
  
#elif( COORD_SYSTEM == SINH_SPHERICAL )
  
  const REAL sinh_x0_inv_W_Jm2          = sinh( x[0][J-2] * inv_sinhW );
  const REAL sinh_x0_inv_W_Jm1          = sinh( x[0][J-1] * inv_sinhW );
  const REAL sinh_x0_inv_W_J            = sinh( x[0][J]   * inv_sinhW );
  const REAL cosh_x0_inv_W_J            = cosh( x[0][J]   * inv_sinhW );
  const REAL alpha_over_a_Jm2           = alpha_n[J-2]/a_n[J-2];
  const REAL alpha_over_a_Jm1           = alpha_n[J-1]/a_n[J-1];
  const REAL alpha_over_a_J             = alpha_n[J]  /a_n[J];
  const REAL alpha_Phi_shx02_over_a_Jm2 = SQR(sinh_x0_inv_W_Jm2) * alpha_over_a_Jm2 * Phi_n[J-2];
  const REAL alpha_Phi_shx02_over_a_Jm1 = SQR(sinh_x0_inv_W_Jm1) * alpha_over_a_Jm1 * Phi_n[J-1];
  const REAL alpha_Phi_shx02_over_a_J   = SQR(sinh_x0_inv_W_J  ) * alpha_over_a_J   * Phi_n[J];
  const REAL coeff_Pi                   = tmp1 * ( sinhW * sinh_inv_W / sinhA ) / ( SQR(sinh_x0_inv_W_J) * cosh_x0_inv_W_J ) ;
  const REAL rhs_Pi                     = coeff_Pi * ( 3.0 * alpha_Phi_shx02_over_a_J - 4.0 * alpha_Phi_shx02_over_a_Jm1 + alpha_Phi_shx02_over_a_Jm2 );
  
#else
  utilities::SFcollapse1D_error(COORD_SYSTEM_ERROR);
#endif

  /* Compute Pi at the outer boundary */
  Pi_np1[J] = Pi_nm1[J] + phi_Pi_coeff * rhs_Pi;
    
}

/* Function to solve the Hamiltonian constraint */
REAL evolution::pointwise_solution_of_the_Hamiltonian_constraint( const int j, grid::parameters grid, const vector<REAL> Phi, const vector<REAL> Pi, const vector<REAL> a ) {

  DECLARE_GRID_PARAMETERS;

  /* Set auxiliary variables */
  const REAL A         = log(a[j-1]);               // A^{n+1}_{j}
  const REAL avgPhi    = 0.5*( Phi[j] + Phi[j-1] ); // 0.5*( Phi^{n+1}_{j+1} + Phi^{n+1}_{j} )
  const REAL avgPi     = 0.5*( Pi[j]  + Pi[j-1]  ); // 0.5*(  Pi^{n+1}_{j+1} +  Pi^{n+1}_{j} )
  const REAL PhiSqr    = SQR(avgPhi);
  const REAL PiSqr     = SQR(avgPi);
  const REAL midx0     = 0.5 * ( x[0][j] + x[0][j-1] );
#if( COORD_SYSTEM == SPHERICAL )
  const REAL PhiPiTerm = 2.0 * M_PI * midx0 * ( PhiSqr + PiSqr );
  const REAL half_invr = 0.5 / midx0;
#elif( COORD_SYSTEM == SINH_SPHERICAL )
  const REAL PhiPiTerm = 2.0 * M_PI * (SQR(sinhA) * inv_sinhW) * sinh(midx0*inv_sinhW) * cosh(midx0*inv_sinhW) / SQR(sinh_inv_W) * ( PhiSqr + PiSqr );
  const REAL half_invr = 0.5/( sinhW * tanh(midx0*inv_sinhW) );
#else
  utilities::SFcollapse1D_error(COORD_SYSTEM_ERROR);
#endif

  /* Set Newton's guess */
  REAL A_old = log(a[j-1]);

  /* Set variable for output */
  REAL A_new = A_old;

  /* Set a counter */
  REAL iter = 0;
  
  /* Perform Newton's method */
  do{

    /* Update A_old */
    A_old = A_new;
    
    /* Compute f and df */
    const REAL tmp0 = half_invr * exp(A_old+A);
    const REAL f  = inv_dx0 * (A_old - A) + tmp0 - half_invr - PhiPiTerm;
    const REAL df = inv_dx0 + tmp0;

    /* Update A_new */
    A_new = A_old - f/df;

    /* Increment the iterator */
    iter++;

  } while( ( fabs(A_new - A_old) > NEWTON_TOL ) && ( iter <= NEWTON_MAX_ITER ) );

  /* Check for convergence */
  if( iter > NEWTON_MAX_ITER ) cerr << "\n(Newton's method WARNING) Newton's method did not converge to a root! j = " << j << " | iter = " << iter << endl;

  /* Return the value of a */
  return( exp(A_new) );
  
}

/* Function to solve the polar slicing condition */
REAL evolution::pointwise_solution_of_the_polar_slicing_condition( const int j, grid::parameters grid, const vector<REAL> a, const vector<REAL> alpha ) {

  DECLARE_GRID_PARAMETERS;

  /* Step 3.b: Compute auxiliary quantities */
  const REAL b = a[j] + a[j-1];
  const REAL c = a[j] - a[j-1];
#if( COORD_SYSTEM == SPHERICAL )
  const REAL midway_r = 0.5 * ( x[0][j] + x[0][j-1] );
#elif( COORD_SYSTEM == SINH_SPHERICAL )
  const REAL midway_r = sinhW * tanh( inv_sinhW * 0.5 * (x[0][j] + x[0][j-1]) );
#else
  utilities::SFcollapse1D_error(COORD_SYSTEM_ERROR);
#endif
  const REAL d = ( 1.0 - 0.25 * SQR(b) )/( 2.0 * midway_r ) - inv_dx0 * c / b;

  /* Step 3.c: Compute alpha */
  return( alpha[j-1]*( 1.0 - d*dx0 )/( 1.0 + d*dx0 ) );
  
}

/* Function to perform the rescaling of the lapse function */
void evolution::rescaling_of_the_lapse( grid::parameters grid, const std::vector<REAL> a, std::vector<REAL> &alpha ) {

  DECLARE_GRID_PARAMETERS;

  /* Set the initial value of kappa */
  REAL kappa = a[0]/alpha[0];

  /* Loop over the grid, updating kappa if needed */
  LOOP(1,Nx0Total) {
    REAL kappa_new = a[j]/alpha[j];
    if( kappa_new < kappa )
      kappa = kappa_new;
  }

  /* Rescale the lapse */
  LOOP(0,Nx0Total) alpha[j] *= kappa;

}

/* Function to compute and output the mass aspect ratio */
void evolution::compute_and_output_mass_aspect_ratio( const int which_level, const int n, const grid::parameters grid, const gridfunction a ) {

  DECLARE_GRID_PARAMETERS;

  ofstream outfile;
  const int number_of_digits = 8;
  outfile.open("out/mass_"+string(number_of_digits - to_string(n).length(),'0')+to_string(n)+".dat");
  outfile.precision(15);
  LOOP(0,Nx0Total) {
    const REAL r_local = r_ito_x0[j];
    const REAL a_local = (which_level == -1) * a.level_nm1[j] + (which_level == 0) * a.level_n[j] + (which_level == 1) * a.level_np1[j];
    const REAL m_local = 0.5 * r_local * ( 1.0 - 1.0/(SQR(a_local)) );
    outfile << scientific << x[0][j] << " " << r_ito_x0[j]  << " " << m_local << endl;
    
  }
  outfile.close();

}
