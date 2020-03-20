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
#include <vector>
#include "macros.hpp"
#include "grid.hpp"
#include "utilities.hpp"
#include "gridfunction.hpp"

using namespace std;

/* Check whether or not a regrid is necessary */
bool utilities::check_regrid_criterion( const grid::parameters grid, const vector<REAL> phi ) {

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
  LOOP(0,Nx0Total) {
    /* First, compute phi_{x} and phi_{xx} */
    REAL phi_x, phi_xx;
    if( (j != 0) && (j != Nx0Total-1) ) {
      phi_x  = inv_dx0*0.5  * ( phi[j+1] - phi[j-1] );
      phi_xx = inv_dx0_sqrd * ( phi[j+1] - 2.0*phi[j] + phi[j-1] );
    }
    else if( j==0 ) {
      phi_x  = inv_dx0*0.5  * ( -3.0*phi[j] + 4.0*phi[j+1] - phi[j+2] );
      phi_xx = inv_dx0_sqrd * ( 2.0*phi[j] - 5.0*phi[j+1] + 4.0*phi[j+2] - phi[j+3] );
    }
    else {
      phi_x  = inv_dx0*0.5  * ( 3.0*phi[j] - 4.0*phi[j-1] + phi[j-2] );
      phi_xx = inv_dx0_sqrd * ( 2.0*phi[j] - 5.0*phi[j-1] + 4.0*phi[j-2] - phi[j-3] );
    }
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
    const REAL tmp0   = cosh( x[0][j] * inv_sinhW );
    const REAL tmp1   = tanh( x[0][j] * inv_sinhW ) * inv_sinhW;
    const REAL tmp2   = 1.0 / SQR(tmp0);
    const REAL tmp3   = SQR(sinhW) / SQR(A_over_sinh_inv_W);
    const REAL phi_rr = tmp3 * tmp2 * ( phi_xx - tmp1 * phi_x );
    const REAL l      = 1.0 / sqrt( abs(phi_rr) );

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
  const int  max_idx    = Nx0_new;
  const REAL x0_max_new = x0_max;
#if( COORD_SYSTEM == SINH_SPHERICAL )
  const REAL sinhA_new  = sinhA;
  const REAL sinhW_new  = sinhW;
#endif

#elif( REGRID_OPTION == REGRID_OUTER_BOUNDARY )

  /* .-------------------------------------------------------.
   * | Option 2: Decrease the position of the outer boundary |
   * .-------------------------------------------------------.
   *
   * This option sets x0_max -> x0_max_new, with x0_max_new < x0_max.
   */
  const int  Nx0_new    = Nx0Total;
  const int  max_idx    = Nx0_new - 1;
  const REAL x0_max_new = x0_max * REGRID_FACTOR;
#if( COORD_SYSTEM == SINH_SPHERICAL )
  const REAL sinhA_new  = sinhA * REGRID_FACTOR;
  const REAL sinhW_new  = sinhW;
#endif

#elif( (REGRID_OPTION == REGRID_POINT_DENSITY) && (COORD_SYSTEM == SINH_SPHERICAL) )

  /* .------------------------------------------------------.
   * | Option 3: Increase the point density near the origin |
   * .------------------------------------------------------.
   *
   * This option sets sinhW -> sinhW_new, with sinhW_new < sinhW.
   */
  const int  Nx0_new    = Nx0Total;
  const int  max_idx    = Nx0_new - 1;
  const REAL x0_max_new = x0_max;
  const REAL sinhA_new  = sinhA;
  const REAL sinhW_new  = sinhW * REGRID_FACTOR;

#else
  utilities::SFcollapse1D_error( REGRIDDING_OPTION_ERROR );
#endif

  /* Set arrays for the values of r_star and x_star */
  vector<REAL> x_new(Nx0_new);
  vector<REAL> r_star(Nx0_new);
  vector<REAL> x_star(Nx0_new);

  /* Compute the values of r_star from the new parameters */
  const REAL dx0_new = (x0_max_new - x0_min)/((REAL)Nx0_new-1.0);
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

  /* Now we know which points we need to find our gridfunctions at, so we
   * define new gridfunctions to be defined at the points x_star.
   */
  gridfunction phi_star(Nx0_new);
  gridfunction Phi_star(Nx0_new);
  gridfunction Pi_star(Nx0_new);
  gridfunction a_star(Nx0_new);
  gridfunction alpha_star(Nx0_new);

  /* Set the gridfunctions at the origin and outer boundary */
#if( REGRID_OPTION == REGRID_RADIAL_POINTS )
  const vector<int> indices = {0};
#else
  const vector<int> indices = {0,max_idx};
#endif

  for( int which_idx=0;which_idx<(int)indices.size();which_idx++) {
    
    const int idx = indices[which_idx];
    
    phi_star.level_nm1[idx] = phi.level_nm1[idx];
    phi_star.level_n[idx]   = phi.level_n[idx];
    phi_star.level_np1[idx] = phi.level_np1[idx];

    Phi_star.level_nm1[idx] = Phi.level_nm1[idx];
    Phi_star.level_n[idx]   = Phi.level_n[idx];
    Phi_star.level_np1[idx] = Phi.level_np1[idx];

    Pi_star.level_nm1[idx] = Pi.level_nm1[idx];
    Pi_star.level_n[idx]   = Pi.level_n[idx];
    Pi_star.level_np1[idx] = Pi.level_np1[idx];

    a_star.level_nm1[idx] = a.level_nm1[idx];
    a_star.level_n[idx]   = a.level_n[idx];
    a_star.level_np1[idx] = a.level_np1[idx];

    alpha_star.level_nm1[idx] = alpha.level_nm1[idx];
    alpha_star.level_n[idx]   = alpha.level_n[idx];
    alpha_star.level_np1[idx] = alpha.level_np1[idx];

  }

  /* Now perform the interpolations, one point at a time */
  LOOP(1,max_idx) {
    utilities::Lagrange_interpolator( j, REGRID_INTERP_STENCIL_SIZE, grid.x[0],
				      phi     , Phi     , Pi     , a     , alpha     ,
				      phi_star, Phi_star, Pi_star, a_star, alpha_star,
				      x_star[j]);
  }

  /* Update all gridfunctions */
  /* phi */
  phi.level_nm1 = phi_star.level_nm1;
  phi.level_n   = phi_star.level_n;
  phi.level_np1 = phi_star.level_np1;
  /* Phi */
  Phi.level_nm1 = Phi_star.level_nm1;
  Phi.level_n   = Phi_star.level_n;
  Phi.level_np1 = Phi_star.level_np1;
  /* Pi */
  Pi.level_nm1 = Pi_star.level_nm1;
  Pi.level_n   = Pi_star.level_n;
  Pi.level_np1 = Pi_star.level_np1;
  /* a */
  a.level_nm1 = a_star.level_nm1;
  a.level_n   = a_star.level_n;
  a.level_np1 = a_star.level_np1;
  /* alpha */
  alpha.level_nm1 = alpha_star.level_nm1;
  alpha.level_n   = alpha_star.level_n;
  alpha.level_np1 = alpha_star.level_np1;

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
void utilities::Lagrange_interpolator( const int interp_index, const int interp_stencil_size, const vector<REAL> x,
				       gridfunction phi      , gridfunction Phi      , gridfunction Pi      , gridfunction a      , gridfunction alpha,
				       gridfunction &phi_star, gridfunction &Phi_star, gridfunction &Pi_star, gridfunction &a_star, gridfunction &alpha_star,
				       const REAL x_star ) {

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
  int size_of_x = x.size();
  int idx_min   = utilities::bisection_index_finder( x, x_star ) - interp_stencil_size/2;

  /* Step 2.b: At this point idx_min holds the integer for which | x[idx_min] - x_star |
   *           is minimal. We now want to perform the interpolation. Ideally, we will
   *           choose points on both sides of idx_min so that we can bracket the value
   *           for interpolation. To this end, we need to make sure we are within the 
   *           bounds of the array x.
   */
  while( idx_min < 0 ) idx_min++;
  while( idx_min + interp_stencil_size > size_of_x ) idx_min--;
  int idx_max = idx_min + interp_stencil_size;

  /* .------------------------------------------------.
   * | Step 2: Compute the Lagrange basis polynomials |
   * .------------------------------------------------.
   */
  vector<REAL> l_j_of_x_star(interp_stencil_size);
  for(int j=idx_min;j<idx_max;j++) {
    REAL numer = 1.0;
    REAL denom = 1.0;
    for(int m=idx_min;m<idx_max;m++) {
      numer *= ( (m == j) ? 1.0 : x_star - x[m] ); // If m=j, multiply by 1
      denom *= ( (m == j) ? 1.0 : x[j]   - x[m] ); // If m=j, multiply by 1
    }
    l_j_of_x_star[j-idx_min] = numer/denom;
  }

  /* .-----------------------------------.
   * | Step 3: Perform the interpolation |
   * .-----------------------------------.
   */
  LOOP(0,interp_stencil_size) {
    /* phi */
    phi_star.level_nm1[interp_index]   += phi.level_nm1[idx_min + j]   * l_j_of_x_star[j];
    phi_star.level_n[interp_index]     += phi.level_n[idx_min + j]     * l_j_of_x_star[j];
    phi_star.level_np1[interp_index]   += phi.level_np1[idx_min + j]   * l_j_of_x_star[j];
    /* Phi */
    Phi_star.level_nm1[interp_index]   += Phi.level_nm1[idx_min + j]   * l_j_of_x_star[j];
    Phi_star.level_n[interp_index]     += Phi.level_n[idx_min + j]     * l_j_of_x_star[j];
    Phi_star.level_np1[interp_index]   += Phi.level_np1[idx_min + j]   * l_j_of_x_star[j];
    /* Pi */
    Pi_star.level_nm1[interp_index]    += Pi.level_nm1[idx_min + j]    * l_j_of_x_star[j];
    Pi_star.level_n[interp_index]      += Pi.level_n[idx_min + j]      * l_j_of_x_star[j];
    Pi_star.level_np1[interp_index]    += Pi.level_np1[idx_min + j]    * l_j_of_x_star[j];
    /* a */
    a_star.level_nm1[interp_index]     += a.level_nm1[idx_min + j]     * l_j_of_x_star[j];
    a_star.level_n[interp_index]       += a.level_n[idx_min + j]       * l_j_of_x_star[j];
    a_star.level_np1[interp_index]     += a.level_np1[idx_min + j]     * l_j_of_x_star[j];
    /* alpha */
    alpha_star.level_nm1[interp_index] += alpha.level_nm1[idx_min + j] * l_j_of_x_star[j];
    alpha_star.level_n[interp_index]   += alpha.level_n[idx_min + j]   * l_j_of_x_star[j];
    alpha_star.level_np1[interp_index] += alpha.level_np1[idx_min + j] * l_j_of_x_star[j];
  }
  
}
/* Bisection index finder */
int utilities::bisection_index_finder( const vector<REAL> x, const REAL x_star ) {

  /* Set the size of x */
  int size_of_x = x.size();
  
  /* Find the initial minimum and maximum indices */
  int j1 = 0;
  int j2 = size_of_x-1;

  /* Find x1 and x2 */
  REAL x1 = x_star - x[j1];
  REAL x2 = x_star - x[j2];

  /* Check if x_star is indeed inside the interval [x1,x2] */
  if( x1*x2 >= 0 ) {
    cerr << "\n(bisection_index_finder ERROR) j1 = " << j1 << " | j2 = " << j2 << endl;
    cerr << "(bisection_index_finder ERROR) x1 = " << x1 << " | x2 = " << x2 << " | x_star = " << x_star << endl;
    utilities::SFcollapse1D_error(BISECTION_INTERVAL_ERROR);
  }

  /* Perform the bisection */
  for(int j=0; j<size_of_x; j++) {

    /* Compute the midpoint of j1 and j2 */
    const int  j_midpoint = 0.5*(j1 + j2);

    /* Compute x_star - x_{j_midpoint} */
    const REAL x_midpoint = x_star - x[j_midpoint];
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
void utilities::parameter_information( grid::parameters grid ) {

  DECLARE_GRID_PARAMETERS;

  /* Compute information about the run to share with the user */
  const REAL rmax = r_ito_x0[Nx0-1];
  int Nr_bet_01 = 0, Nr_bet_05 = 0;
  LOOP(0,Nx0Total) {
    const REAL rlocal = r_ito_x0[j];
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
    cout << fixed      << setprecision(0) << "phi_{0}    = " << PHI0            << endl;
    cout << fixed      << setprecision(0) << "r_{0}      = " << R0              << endl;
    cout << fixed      << setprecision(0) << "delta      = " << DELTA           << endl;
    cout << "\nRadial grid information:"       << endl;
    cout << fixed      << setprecision(0) << "N_{r}      = " << Nx0             << endl;
    cout << fixed      << setprecision(0) << "r_{max}    = " << rmax            << endl;
#if( COORD_SYSTEM == SINH_SPHERICAL )
    cout << fixed      << setprecision(2) << "sinhA      = " << sinhA           << endl;
    cout << fixed      << setprecision(5) << "sinhW      = " << sinhW           << endl;
#endif
    cout << fixed      << setprecision(0) << "Points in 0<r<1: " << Nr_bet_01   << endl;
    cout << fixed      << setprecision(0) << "Points in 0<r<5: " << Nr_bet_05   << endl;
    cout << "\nTime evolution information:"    << endl;
    cout << fixed      << setprecision(0) << "N_{t}      = " << Nt              << endl;
    cout << fixed      << setprecision(2) << "t_{final}  = " << t_final         << endl;
    cout << scientific << setprecision(3) << "dt         = " << dt              << endl;
    cout << fixed      << setprecision(2) << "CFL factor = " << CFL_FACTOR      << endl;
    cout << "\n";

    if( iter == 1 ) cout.rdbuf(coutbuf); // Reset to standard output again

    iter++;

  }

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
      cerr << "(SFcollapse1D INFO) Correct usage is: ./SFcollapse1D Nx Domain_size t_final\n";
      cerr << "(SFcollapse1D INFO) Terminating the program...\n";
      exit(SPHERICAL_USAGE_ERROR);
      break;

    case SINH_SPHERICAL_USAGE_ERROR:
      cerr << "(SFcollapse1D ERROR) Incorrect usage of the program!\n";
      cerr << "(SFcollapse1D INFO) SinhSpherical coordinates selected.\n";
      cerr << "(SFcollapse1D INFO) Correct usage is: ./SFcollapse1D Nx Domain_size t_final sinhW\n";
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
