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
#include "utilities.hpp"

using namespace std;

/* Lagrange interpolator */
REAL utilities::Lagrange_interpolator( const int interp_stencil_size, const vector<REAL> x, const vector<REAL> y, const REAL x_star ) {

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
    l_j_of_x_star[j] = numer/denom;
  }

  /* .-----------------------------------.
   * | Step 3: Perform the interpolation |
   * .-----------------------------------.
   */
  REAL y_star = 0.0;
  LOOP(0,interp_stencil_size) {
    y_star += y[idx_min + j] * l_j_of_x_star[j];
  }

  return y_star;
  
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
  if( x1*x2 >= 0 ) utilities::SFcollapse1D_error(BISECTION_INTERVAL_ERROR);

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
	return x1;
      else
	return x2;
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
  ofstream fileout("out_parameters.txt"); // Set fileout with the name of the output file
  auto *coutbuf = cout.rdbuf();           // Save current cout buf
  int iter = 0;
  
  while( iter < 2 ) {
    
    if( iter == 1 ) cout.rdbuf(fileout.rdbuf()); // Redirect cout to the file opened by fileout
    
    cout << "\n";
    cout << ".------------------------."       << endl;
    cout << "| Parameters of this run |"       << endl;
    cout << ".------------------------."       << endl;
#if( COORD_SYSTEM == SPHERICAL )
    cout << "Coordinate system: Spherical"     << endl;
#elif( COORD_SYSTEM == SINH_SPHERICAL )
    cout << "Coordinate system: SinhSpherical" << endl;
#endif
    cout << "Initial condition information:"   << endl;
    cout << "phi_{0}    = " << PHI0            << endl;
    cout << "r_{0}      = " << R0              << endl;
    cout << "delta      = " << DELTA           << endl;
    cout << "\nRadial grid information:"       << endl;
    cout << "N_{r}      = " << Nx0             << endl;
    cout << "r_{max}    = " << rmax            << endl;
#if( COORD_SYSTEM == SINH_SPHERICAL )
    cout << "sinhA      = " << sinhA           << endl;
    cout << "sinhW      = " << sinhW           << endl;
#endif
    cout << "Points in 0<r<1: " << Nr_bet_01   << endl;
    cout << "Points in 0<r<5: " << Nr_bet_05   << endl;
    cout << "\nTime evolution information:"    << endl;
    cout << "N_{t}      = " << Nt              << endl;
    cout << "t_{final}  = " << t_final         << endl;
    cout << "dt         = " << dt              << endl;
    cout << "CFL factor = " << CFL_FACTOR      << endl;
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
      cerr << "(NaN_checker) Iteration: " << n
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
      cerr << "(initialize_parameters ERROR) Unknown coordinate system! Please check the macros.hpp file!\n";
      cerr << "(SFcollapse1D INFO) Terminating the program...\n";
      exit(COORD_SYSTEM_ERROR);
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
