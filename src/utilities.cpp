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
#include <iomanip>
#include <fstream>
#include <cmath>
#include "macros.hpp"
#include "grid.hpp"
#include "utilities.hpp"
#include "evolution.hpp"
#include "gridfunction.hpp"

using namespace std;

/* Check whether or not a regrid is necessary */
bool utilities::check_regrid_criterion( const grid::parameters grid, const realvec phi ) {

  DECLARE_GRID_PARAMETERS;

  /* We will adopt the same regrid criterion from Baumgarte, 2018 [1], which
   * is
   *
   * l := 1.0 / sqrt( phi_{rr} ) < 25 dr_{min} ,
   *
   * where l is a length scale we have introduced and dr_{min} is the minimum
   * radial step size on our current numerical grid (in our code this is given
   * by grid.ds_min). Due to the square root, we will use | phi_{rr} | to avoid
   * running into trouble.
   *
   * [1] T.W. Baumgarte, Aspherical deformations of the Choptuik spacetime,
   *     Phys. Rev. D 98 084012, 2018. ArXiV: https://arxiv.org/abs/1807.10342
   */
  LOOP(1,Nx0Total-1) {
    /* First, compute phi_{x} and phi_{xx} */
    const real phi_x  = inv_dx0*0.5  * ( phi[j+1] - phi[j-1] );
    const real phi_xx = inv_dx0_sqrd * ( phi[j+1] - 2.0*phi[j] + phi[j-1] );
    /* Then, compute phi_{rr}. We know that
     *
     * r = A sinh( x/w ) / sinh( 1/w )
     *
     * which implies
     *
     * => dr = [ (A/w) cosh( x/w ) / sinh( 1/w ) ] dx
     *
     * => pd_{r} = [ (w/A) sinh( 1/w ) / cosh( x/w ) ] pd_{x}
     *
     * => pd_{rr} = pd_{r} ( pd_{r} )
     *
     * => pd_{rr} = ( (w/A) sinh(1/w) )^{2} ( 1/cosh(x/w) ) pd_{x} [ ( 1/cosh(x/w) ) pd_{x} )
     *
     * => pd_{rr} = ( (w/A) sinh(1/w) )^{2} ( 1/cosh(x/w) )^{2} [ pd_{xx} - tanh(x/w)/w pd_{x} )
     */
    const real tmp0   = cosh( x[0][j] * inv_sinhW );
    const real tmp1   = tanh( x[0][j] * inv_sinhW ) * inv_sinhW;
    const real tmp2   = 1.0 / SQR(tmp0);
    const real tmp3   = SQR(sinhW) / SQR(A_over_sinh_inv_W);
    const real phi_rr = tmp3 * tmp2 * ( phi_xx - tmp1 * phi_x );
    const real l      = 1.0 / sqrt( abs(phi_rr) );

    /* If the criterion is met, return true */
    if( l < 25*ds_min ) return true;
  }

  /* If we reach this point, then the criterion has not been met */
  return false;
  
}

/* Regridding function */
void utilities::regrid( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha ) {

  DECLARE_GRID_PARAMETERS;

  /* This function is the main driver when performing a regrid.
   * Regrids can be performed when one the following is desired:
   *
   * .-----------------------------------------------------------.
   * | Spherical or SinhSpherical coordinates:                   |
   * .-----------------------------------------------------------.
   * |   1) Increase/decrease the number of radial points        |
   * |   2) Increase/decrease the position of the outer boundary |
   * .-----------------------------------------------------------.
   * | SinhSpherical coordinates only:                           |
   * .-----------------------------------------------------------.
   * |   3) Increase/decrease the point density near the origin  |
   * .-----------------------------------------------------------.
   *
   * Although we have listed all possible scenarios above, this initial
   * release of the code won't make it possible to regrid in order to
   * *decrease* the resolution of the grid, only to *increase* it.
   * Thus, in reality, the only available options are:
   *
   * 1) Increase the number of radial points
   * 2) Decrease the position of the outer boundary
   * 3) Increase the point density near the origin (SinhSpherical only)
   *
   * If we keep the number of radial points constant, the effect of regrid
   * option 2) above in Spherical coordinates is virtually the same as
   * option 3) in SinhSpherical coordinates: we increase the resolution
   * on a smaller portion of the computational domain. However, option 3)
   * has the advantage that the outer boundary remains unchanged, so we
   * don't need to worry about it getting into causal contact with the
   * center of the simulation after a regrid is performed (assuming it has
   * been placed out of causal contact in the first place).
   *
   * To choose which option above this algorithm should use, please set
   * the REGRID_OPTION macro appropriately in the macros.hpp file.
   */

#if( REGRID_OPTION == REGRID_RADIAL_POINTS )

  /* .------------------------------------------------.
   * | Option 1: Increase the number of radial points |
   * .------------------------------------------------.
   *
   * This option sets Nx0 -> Nx0_new, with Nx0_new > Nx0.
   */
  const int  Nx0_new    = (Nx0Total-1) * REGRID_FACTOR + 1;
  const real x0_max_new = x0_max;
#if( COORD_SYSTEM == SINH_SPHERICAL )
  const real sinhA_new  = sinhA;
  const real sinhW_new  = sinhW;
#endif

#elif( REGRID_OPTION == REGRID_OUTER_BOUNDARY )

  /* .-------------------------------------------------------.
   * | Option 2: Decrease the position of the outer boundary |
   * .-------------------------------------------------------.
   *
   * This option sets x0_max -> x0_max_new, with x0_max_new < x0_max.
   */
  const int  Nx0_new    = Nx0Total;
  const real x0_max_new = x0_max * REGRID_FACTOR;
#if( COORD_SYSTEM == SINH_SPHERICAL )
  const real sinhA_new  = sinhA * REGRID_FACTOR;
  const real sinhW_new  = sinhW;
#endif

#elif( (REGRID_OPTION == REGRID_POINT_DENSITY) && (COORD_SYSTEM == SINH_SPHERICAL) )

  /* .------------------------------------------------------.
   * | Option 3: Increase the point density near the origin |
   * .------------------------------------------------------.
   *
   * This option sets sinhW -> sinhW_new, with sinhW_new < sinhW.
   */
  const int  Nx0_new    = Nx0Total;
  const real x0_max_new = x0_max;
  const real sinhA_new  = sinhA;
  const real sinhW_new  = sinhW * REGRID_FACTOR;

#else
  utilities::SFcollapse1D_error( REGRIDDING_OPTION_ERROR );
#endif

  /* Set arrays for the values of r_star and x_star */
  realvec x_new(Nx0_new);
  realvec r_star(Nx0_new);
  realvec x_star(Nx0_new);

  /* Compute the values of r_star from the new parameters */
  const real dx0_new = (x0_max_new - x0_min)/((real)Nx0_new-1.0);
  LOOP(0,Nx0_new) {
    x_new[j]  = (j-Ngz0) * dx0_new;
#if( COORD_SYSTEM == SPHERICAL )
    r_star[j] = grid::compute_r_of_x( x_new[j] );
    x_star[j] = grid::compute_x_of_r( r_star[j] );
#elif( COORD_SYSTEM == SINH_SPHERICAL )
    r_star[j] = grid::compute_r_of_x( x_new[j] , sinhA_new, sinhW_new );
    x_star[j] = grid::compute_x_of_r( r_star[j], sinhA    , sinhW     );
#endif
  }

  cout << scientific << setprecision(15) << "Before: " << alpha.level_np1[0] << endl;

  // Set the gridfunctions at the origin and outer boundary
  phi.level_nm1[0]   = phi.level_np1[0];
  Phi.level_nm1[0]   = Phi.level_np1[0];
  Pi.level_nm1[0]    = Pi.level_np1[0];
  a.level_nm1[0]     = a.level_np1[0];
  alpha.level_nm1[0] = alpha.level_np1[0];

  // phi.level_nm1[max_idx] = phi.level_np1[max_idx];
  // Phi.level_nm1[max_idx] = Phi.level_np1[max_idx];
  // Pi.level_nm1[max_idx]  = Pi.level_np1[max_idx];

  /* Now perform the interpolations, one point at a time */
  LOOP(1,Nx0Total) {
    utilities::Lagrange_interpolator( j, REGRID_INTERP_STENCIL_SIZE, grid.x[0],
				      phi.level_np1, Phi.level_np1, Pi.level_np1, a.level_np1, alpha.level_np1,
				      phi.level_nm1, Phi.level_nm1, Pi.level_nm1, a.level_nm1, alpha.level_nm1,
				      x_star[j] );
  }

  phi.level_nm1[0] = - 3.0 * phi.level_nm1[2] + 4.0 * phi.level_nm1[1];

  /* Then perform a time step, similar to the initial time step */
  /* .--------------------------------------------------------------.
   * | Step 1: Apply inner boundary conditions to Phi, a, and alpha |
   * .--------------------------------------------------------------.
   */
  Phi.level_n[0]   = 0.0;
  a.level_n[0]     = 1.0;
  alpha.level_n[0] = 1.0;

  /* .-----------------------------------.
   * | Step 2: Integrate phi, Phi and Pi |
   * .-----------------------------------.
   */
  evolution::time_step_scalarfield_gridfunctions( 0, grid, 
						  phi.level_nm1, Phi.level_nm1, Pi.level_nm1, a.level_nm1, alpha.level_nm1, 
						  Phi.level_nm1, Pi.level_nm1, a.level_nm1, alpha.level_nm1,
						  Phi.level_n  , Pi.level_n  , phi.level_n );
  /* Step 2.d: Apply boundary conditions */
  evolution::apply_outgoing_radiation_bdry_cond( 0, grid,
						 phi.level_nm1, Pi.level_nm1,
						 phi.level_nm1, Phi.level_nm1, a.level_nm1, alpha.level_nm1,
						 phi.level_n  , Phi.level_n  , Pi.level_n );

  Pi.level_n[0] = - Pi.level_nm1[0] + Pi.level_n[1] + Pi.level_nm1[1];

  /* .-------------------------------.
   * | Step 3: Integrate a and alpha |
   * .-------------------------------.
   */
  // TODO: can this be parallelized in some way? Maybe SIMD?
  LOOP(1,grid.Nx0Total) {
    /* Step 3.a: Compute a */
    a.level_n[j] = evolution::pointwise_solution_of_the_Hamiltonian_constraint(j,grid,Phi.level_n,Pi.level_n,a.level_n);
    
    /* Step 3.b: Compute alpha */
    alpha.level_n[j] = evolution::pointwise_solution_of_the_polar_slicing_condition( j, grid, a.level_n, alpha.level_n );
  }
  /* Step 3.d: Now rescale alpha */
  evolution::rescaling_of_the_lapse(grid,a.level_n,alpha.level_n);
  
  /* .-------------------------.
   * | Step 4: Update the time |
   * .-------------------------.
   */
  grid.t += 0.5 * grid.dt;

  /* Perform an integration step */
  /* .--------------------------------------------------------------.
   * | Step 1: Apply inner boundary conditions to Phi, a, and alpha |
   * .--------------------------------------------------------------.
   */
  Phi.level_np1[0]   = 0.0;
  a.level_np1[0]     = 1.0;
  alpha.level_np1[0] = 1.0;
  
  /* .-----------------------------------.
   * | Step 2: Integrate phi, Phi and Pi |
   * .-----------------------------------.
   */
  evolution::time_step_scalarfield_gridfunctions( 1, grid, 
						  phi.level_n, Phi.level_n  , Pi.level_n  , a.level_n  , alpha.level_n  , 
						  Phi.level_nm1, Pi.level_nm1, a.level_nm1, alpha.level_nm1,
						  Phi.level_np1, Pi.level_np1, phi.level_np1 );
  /* Step 2.d: Apply boundary conditions */
  evolution::apply_outgoing_radiation_bdry_cond( 1, grid,
						 phi.level_nm1, Pi.level_nm1,
						 phi.level_n  , Phi.level_n  , a.level_n, alpha.level_n,
						 phi.level_np1, Phi.level_np1, Pi.level_np1);
  
  Pi.level_np1[0] = - Pi.level_n[0] + Pi.level_np1[1] + Pi.level_n[1];
  /* .-------------------------------.
   * | Step 3: Integrate a and alpha |
   * .-------------------------------.
   */
  // TODO: can this be parallelized in some way? Maybe SIMD?
  LOOP(1,grid.Nx0Total) {
    /* Step 3.a: Compute a */
    a.level_np1[j] = evolution::pointwise_solution_of_the_Hamiltonian_constraint(j,grid,Phi.level_np1,Pi.level_np1,a.level_np1);
    
    /* Step 3.b: Compute alpha */
    alpha.level_np1[j] = evolution::pointwise_solution_of_the_polar_slicing_condition( j, grid, a.level_np1, alpha.level_np1 );
  }
  /* Step 3.d: Now rescale alpha */
  evolution::rescaling_of_the_lapse(grid,a.level_np1,alpha.level_np1);

  /* .-------------------------.
   * | Step 4: Update the time |
   * .-------------------------.
   */
  grid.t += 0.5 * grid.dt;

  cout << scientific << setprecision(15) << "After: " << alpha.level_np1[0] << endl;

  /* Shift time levels appropriately */
  phi.level_n   = phi.level_nm1;
  Phi.level_n   = Phi.level_nm1;
  Pi.level_n    = Pi.level_nm1;
  a.level_n     = a.level_nm1;
  alpha.level_n = alpha.level_nm1;

  /* Finally, update all grid parameters */
  grid.Nx0      = Nx0_new;
  grid.x0_max   = x0_max_new;
  grid.x[0]     = x_new;
  grid.r_ito_x0 = r_star;
  grid.ds_min   = grid.r_ito_x0[1] - grid.r_ito_x0[0];
  grid.dt       = CFL_FACTOR * grid.ds_min;
  grid.Nt       = (int)grid.t_final/grid.dt;
#if( COORD_SYSTEM == SINH_SPHERICAL )
  grid.sinhA             = sinhA_new;
  grid.sinhW             = sinhW_new;
  grid.inv_sinhW         = 1.0/grid.sinhW;
  grid.sinh_inv_W        = sinh(grid.inv_sinhW);
  grid.A_over_sinh_inv_W = grid.sinhA / grid.sinh_inv_W;
#endif
  grid.current_regrid_level++;

}

/* Lagrange interpolator */
void utilities::Lagrange_interpolator( const int interp_index, const int interp_stencil_size, const realvec x,
				       const realvec phi, const realvec Phi, const realvec Pi, const realvec a, const realvec alpha,
				       realvec &phi_star, realvec &Phi_star, realvec &Pi_star, realvec &a_star, realvec &alpha_star,
				       const real x_star ) {

  /* This function implements a simple Lagrange polynomial interpolator.
   * It follows closely the discussion in: https://en.wikipedia.org/wiki/Lagrange_polynomial
   *
   * We get as inputs the values { x_i } and the function values { y_i = y(x_i) }. We then
   * want to compute the value of the function y_star = y(x_star), assuming that x_star is
   * inside the interval [x_{0},x_{max}].
   */

  /* .------------------------------------------.
   * | Step 1: Find the appropriate array index |
   * .------------------------------------------.
   *
   * Step 1.a: Set the size of x and bisect the closest value of x_star from x */
  int idx       = utilities::bisection_index_finder( x, x_star );
  int idx_min   = max(0,idx-interp_stencil_size/2-1);

  /* Step 2.b: At this point idx_min holds the integer for which | x[idx_min] - x_star |
   *           is minimal. We now want to perform the interpolation. Ideally, we will
   *           choose points on both sides of idx_min so that we can bracket the value
   *           for interpolation. To this end, we need to make sure we are within the 
   *           bounds of the array x.
   */
  // while( idx_min < 0 ) idx_min++;
  // while( idx_min + interp_stencil_size > size_of_x ) idx_min--;
  // int idx_max = idx_min + interp_stencil_size;

  /* .------------------------------------------------.
   * | Step 2: Compute the Lagrange basis polynomials |
   * .------------------------------------------------.
   */
  realvec l_i_of_x_star(interp_stencil_size);
  for(int i=0;i<interp_stencil_size;i++) {
    real numer = 1.0;
    real denom = 1.0;
    for(int j=0;j<i;j++) {
      numer *=     x_star     - x[idx_min + j];
      denom *= x[idx_min + i] - x[idx_min + j];
    }
    for(int j=i+1;j<interp_stencil_size;j++) {
      numer *=     x_star     - x[idx_min + j];
      denom *= x[idx_min + i] - x[idx_min + j];
    }
    l_i_of_x_star[i] = numer/denom;
  }

  /* .-----------------------------------.
   * | Step 3: Perform the interpolation |
   * .-----------------------------------.
   */
  phi_star[interp_index]   = 0.0;
  Phi_star[interp_index]   = 0.0;
  Pi_star[interp_index]    = 0.0;
  a_star[interp_index]     = 0.0;
  alpha_star[interp_index] = 0.0;
  LOOP(0,interp_stencil_size) {
    phi_star[interp_index]   += phi[   idx_min + j] * l_i_of_x_star[j];
    Phi_star[interp_index]   += Phi[   idx_min + j] * l_i_of_x_star[j];
    Pi_star[interp_index]    += Pi[    idx_min + j] * l_i_of_x_star[j];
    a_star[interp_index]     += a[     idx_min + j] * l_i_of_x_star[j];
    alpha_star[interp_index] += alpha[ idx_min + j] * l_i_of_x_star[j];
  }
  
}
/* Bisection index finder */
int utilities::bisection_index_finder( const realvec x, const real x_star ) {

  /* Set the size of x */
  int size_of_x = x.size();
  
  /* Find the initial minimum and maximum indices */
  int j1 = 0;
  int j2 = size_of_x-1;

  /* Find x1 and x2 */
  real x1 = x_star - x[j1];
  real x2 = x_star - x[j2];

  /* Check if x_star is indeed inside the interval [x1,x2] */
  if( x1*x2 > 0 ) {
    cerr << setprecision(5) << "\n(bisection_index_finder ERROR) j1 = " << j1 << " | j2 = " << j2 << endl;
    cerr << setprecision(5) << "(bisection_index_finder ERROR) x1 = " << x1 << " | x2 = " << x2 << " | x_star = " << x_star << endl;
    utilities::SFcollapse1D_error(BISECTION_INTERVAL_ERROR);
  }

  /* Perform the bisection */
  for(int j=0; j<size_of_x; j++) {

    /* Compute the midpoint of j1 and j2 */
    const int  j_midpoint = 0.5*(j1 + j2);

    /* Compute x_star - x_{j_midpoint} */
    const real x_midpoint = x_star - x[j_midpoint];
    if( x_midpoint*x1 < 0 ) {
      j2 = j_midpoint;
      x2 = x_midpoint;
    }
    else {
      j1 = j_midpoint;
      x1 = x_midpoint;
    }
    /* Since we are bisecting the elements of an array, the search will end
     * when | j2-j1 | = 1, i.e. we have two consecutive elements */
    if( abs( j1 - j2 ) == 1 ) {
      /* Now check whether x[j1] or x[j2] is closer to x_star */
      if( fabs(x1) < fabs(x2) )
	return j1;
      else
	return j2;
    }
  }

  /* If we reach this point, then the algorithm failed. Return an error */
  utilities::SFcollapse1D_error(BISECTION_CONVERGENCE_ERROR);

  /* This is added to get rid of a compiler warning, but the
   * function is never able to actually reach this point
   */
  return BISECTION_CONVERGENCE_ERROR;

}

/* Print parameter information to the user and to file */
void utilities::parameter_information( const real phi0, grid::parameters grid ) {

  DECLARE_GRID_PARAMETERS;

  /* Compute information about the run to share with the user */
  const real rmax = r_ito_x0[Nx0-1];
  int Nr_bet_01 = 0, Nr_bet_05 = 0;
  LOOP(0,Nx0Total) {
    const real rlocal = r_ito_x0[j];
    if( rlocal < 1.0 ) Nr_bet_01++; // Counts points for which 0<r<1
    if( rlocal < 5.0 ) Nr_bet_05++; // Counts points for which 0<r<5
  }

  /* Trick based on Nawaz's answer in Stackoverflow. This avoids code duplication.
   * Source: https://stackoverflow.com/questions/10150468/how-to-redirect-cin-and-cout-to-files
   */
  ofstream fileout; // Set fileout with the name of the output file
  if( t==0 ) {
    fileout.open("out_parameters.txt");
  }
  else {
    fileout.open("out_parameters.txt",ios_base::app);
  }

  auto *coutbuf = cout.rdbuf();           // Save current cout buf
  int iter = 0;
  
  while( iter < 2 ) {
    
    if( iter == 1 ) cout.rdbuf(fileout.rdbuf()); // Redirect cout to the file opened by fileout
    
    cout << "\n";
    if( t == 0 ) {
      cout << ".------------------------."       << endl;
      cout << "| Parameters of this run |"       << endl;
      cout << ".------------------------."       << endl;
    }
    else {
      cout << ".----------------."             << endl;
      cout << "| New parameters |"             << endl;
      cout << ".----------------."             << endl;
      cout << "Regrid level: " << current_regrid_level << endl;
    }
#if( COORD_SYSTEM == SPHERICAL )
    cout << "Coordinate system: Spherical"     << endl;
#elif( COORD_SYSTEM == SINH_SPHERICAL )
    cout << "Coordinate system: SinhSpherical" << endl;
#endif
    cout << "Initial condition information:"   << endl;
    cout << fixed      << setprecision(13) << "phi_{0}    = " << phi0            << endl;
    cout << fixed      << setprecision(2)  << "r_{0}      = " << R0              << endl;
    cout << fixed      << setprecision(2)  << "delta      = " << DELTA           << endl;
    cout << "\nRadial grid information:"       << endl;
    cout << fixed      << setprecision(0)  << "N_{r}      = " << Nx0             << endl;
    cout << fixed      << setprecision(0)  << "r_{max}    = " << rmax            << endl;
    cout << scientific << setprecision(5)  << "dr_{min}   = " << ds_min          << endl;
#if( COORD_SYSTEM == SINH_SPHERICAL )
    cout << fixed      << setprecision(2)  << "sinhA      = " << sinhA           << endl;
    cout << fixed      << setprecision(5)  << "sinhW      = " << sinhW           << endl;
#endif
    cout << fixed      << setprecision(0)  << "Points in 0<r<1: " << Nr_bet_01   << endl;
    cout << fixed      << setprecision(0)  << "Points in 0<r<5: " << Nr_bet_05   << endl;
    cout << "\nTime evolution information:"    << endl;
    cout << fixed      << setprecision(0)  << "N_{t}      = " << Nt              << endl;
    cout << fixed      << setprecision(2)  << "t_{final}  = " << t_final         << endl;
    cout << scientific << setprecision(3)  << "dt         = " << dt              << endl;
    cout << fixed      << setprecision(2)  << "CFL factor = " << CFL_FACTOR      << endl;
    cout << "\n";

    if( iter == 1 ) cout.rdbuf(coutbuf); // Reset to standard output again

    iter++;

  }

}

/* Output the mass-aspect function */
void utilities::compute_and_output_mass_aspect_function( const int which_level, const int n, const grid::parameters grid, const gridfunction a ) {

  DECLARE_GRID_PARAMETERS;

  /* Compute the mass aspect function based on which_level */
  ofstream outfile;
  const int number_of_digits = 8;
  outfile.open("out/mass_"+string(number_of_digits - to_string(n).length(),'0')+to_string(n)+".dat");
  outfile.precision(15);
  LOOP(0,Nx0Total) {
    const real a_local = (which_level == -1) * a.level_nm1[j] + (which_level == 0) * a.level_n[j] + (which_level == 1) * a.level_np1[j];
    const real mass    = 0.5 * r_ito_x0[j] * ( 1.0 - 1.0/SQR(a_local) );
    outfile << scientific << x[0][j] << " " << r_ito_x0[j]  << " " << mass << endl;
  }
  outfile.close();

}

/* Output central values (r=0) of the gridfunctions */
void utilities::output_gridfunctions_central_values( const int n, const grid::parameters grid,
						     const realvec phi, const realvec Phi, const realvec Pi, const realvec a, const realvec alpha ) {

  DECLARE_GRID_PARAMETERS;

  ofstream out_central;
  if( t > 0.0 ) {
    out_central.open("out_central_values.dat",ios_base::app);
  }
  else {
    out_central.open("out_central_values.dat");
  }

  out_central << scientific << setprecision(15)
	      << t        << " "
	      << alpha[0] << " "
	      << phi[0]   << " "
	      << Phi[0]   << " "
	      << Pi[0]    << " "
	      << a[0]     << endl;

  out_central.close();

}

/* Check whether or not a black hole has formed by checking whether or not the lapse has collapsed */
bool utilities::check_for_collapse_of_the_lapse( const grid::parameters grid, const gridfunction alpha ) {

  DECLARE_GRID_PARAMETERS;
  
  LOOP(0,Nx0Total) {
    if( alpha.level_np1[j] < LAPSE_COLLAPSE_CRITERION ) return true;
  }

  return false;

}

/* NaN checker: verifies the gridfunctions still have valid values */
void utilities::NaN_checker( const int n, grid::parameters grid, gridfunction phi, gridfunction Phi, gridfunction Pi, gridfunction a, gridfunction alpha ) {

  DECLARE_GRID_PARAMETERS;

  /* Loop over gridfunctions and check for NaN's */
  LOOP(0,Nx0Total) {
    if( isnan( phi.level_np1[j] ) || isnan( Phi.level_np1[j] ) || isnan( Pi.level_np1[j] ) || isnan( a.level_np1[j] ) || isnan( alpha.level_np1[j] ) ) {
      cerr << "\n(NaN_checker) Iteration: " << n
	   << " | time = " << t
	   << " | j = "    << j
	   << endl;
      utilities::SFcollapse1D_error( NAN_ERROR );
    }
  }

}

/* Various errors that may occur during execution */
void utilities::SFcollapse1D_error( const int error ) {

  switch (error) {

    case SPHERICAL_USAGE_ERROR:
      cerr << "(SFcollapse1D ERROR) Incorrect usage of the program!\n";
      cerr << "(SFcollapse1D INFO) Spherical coordinates selected.\n";
      cerr << "(SFcollapse1D INFO) Correct usage is: ./SFcollapse1D Nx Domain_size t_final phi0\n";
      cerr << "(SFcollapse1D INFO) Terminating the program...\n";
      exit(SPHERICAL_USAGE_ERROR);
      break;

    case SINH_SPHERICAL_USAGE_ERROR:
      cerr << "(SFcollapse1D ERROR) Incorrect usage of the program!\n";
      cerr << "(SFcollapse1D INFO) SinhSpherical coordinates selected.\n";
      cerr << "(SFcollapse1D INFO) Correct usage is: ./SFcollapse1D Nx Domain_size t_final sinhW phi0\n";
      cerr << "(SFcollapse1D INFO) Terminating the program...\n";
      exit(SINH_SPHERICAL_USAGE_ERROR);
      break;

    case GRID_STRUCTURE_ERROR:
      cerr << "(initialize_parameters ERROR) Unknown grid structure! Please check the macros.hpp file.\n";
      cerr << "(SFcollapse1D INFO) Terminating the program...\n";
      exit(GRID_STRUCTURE_ERROR);
      break;

    case COORD_SYSTEM_ERROR:
      cerr << "(initialize_parameters ERROR) Unknown coordinate system! Please check the macros.hpp file.\n";
      cerr << "(SFcollapse1D INFO) Terminating the program...\n";
      exit(COORD_SYSTEM_ERROR);
      break;

    case REGRIDDING_OPTION_ERROR:
      cerr << "(regrid ERROR) Unknown regridding option! Please check the macros.hpp file.\n";
      cerr << "(SFcollapse1D INFO) Terminating the program...\n";
      exit(NAN_ERROR);
      break;

    case GRIDFUNCTION_TIMELEVEL_ERROR:
      cerr << "(shift_timelevels ERROR) Unknown value of which_level inside shift_timelevels!\n";
      cerr << "(SFcollapse1D INFO) Terminating the program...\n";
      exit(GRIDFUNCTION_TIMELEVEL_ERROR);
      break;

    case BISECTION_INTERVAL_ERROR:
      cerr << "(Bisection_index_finder ERROR) Bad interval!\n";
      cerr << "(SFcollapse1D INFO) Terminating the program...\n";
      exit(BISECTION_INTERVAL_ERROR);
      break;

    case BISECTION_CONVERGENCE_ERROR:
      cerr << "(Bisection_index_finder ERROR) Did not converge!\n";
      cerr << "(SFcollapse1D INFO) Terminating the program...\n";
      exit(BISECTION_CONVERGENCE_ERROR);
      break;

    case NAN_ERROR:
      cerr << "(NaN_checker ERROR) NaN found!\n";
      cerr << "(SFcollapse1D INFO) Terminating the program...\n";
      exit(NAN_ERROR);
      break;

  }

}
