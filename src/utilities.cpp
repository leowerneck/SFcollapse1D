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
#include <cmath>
#include <vector>
#include "macros.hpp"
#include "utilities.hpp"

using namespace std;

/* Lagrange interpolator */
REAL utilities::Lagrange_interpolator( const int interp_stencil_size, const vector<REAL> x, const vector<REAL> y, const REAL x_star ) {

  /* This function implements a simple Lagrange polynomial interpolator.
   * It follows closely the discussion in: https://en.wikipedia.org/wiki/Lagrange_polynomial
   *
   * We get as inputs the values { x_i } and the function values { y_i = y(x_i) }. We then
   * want to compute the value of the function y_star = y(x_star).
   *
   * First find the appropriate indices to use in the interpolation */
  int size_of_x    = x.size();
  int idx_min = utilities::Bisection_index_finder( x, x_star ) - interp_stencil_size/2;

  /* At this point idx_min holds the integer for which | x[idx_min] - x_star | is minimal.
   * We now want to perform the interpolation. Ideally, we will choose points on both sides
   * of idx_min so that we can bracket the value for interpolation. 
   *
   * Make sure we are within the bounds of the array x */
  while( idx_min < 0 ) idx_min++;
  while( idx_min + interp_stencil_size > size_of_x ) idx_min--;
  int idx_max = idx_min + interp_stencil_size;

  /* Then compute the Lagrange basis polynomials */
  vector<REAL> l_j_of_x_star(interp_stencil_size);
  REAL numer = 1.0;
  REAL denom = 1.0;
  for(int j=idx_min;j<idx_max;j++) {
    for(int m=idx_min;m<idx_max;m++) {
      numer *= ( (m == j) ? 1.0 : x_star - x[m] ); // If m=j, multiply by 1
      denom *= ( (m == j) ? 1.0 : x[j]   - x[m] ); // If m=j, multiply by 1
      cout << x_star << " " << x[m] << " | " << x[j] << " " << x[m] << endl;
    }
    cout << "numer = " << numer << endl;
    cout << "denom = " << denom << endl;
    l_j_of_x_star[j] = numer/denom;
  }

  /* Perform the interpolation */
  REAL y_star = 0.0;
  LOOP(0,interp_stencil_size) {
    y_star += y[idx_min + j] * l_j_of_x_star[j];
  }

  return y_star;
  
}
/* Bisection index finder */
int utilities::Bisection_index_finder( const vector<REAL> x, const REAL x_star ) {

  /* Set the size of x */
  int size_of_x = x.size();
  
  /* Find the initial minimum and maximum indices */
  int j1 = 0;
  int j2 = size_of_x-1;

  /* Find x1 and x2 */
  REAL x1 = x_star - x[j1];
  REAL x2 = x_star - x[j2];

  cout << j1 << " " << j2 << " | " << x1 << " " << x2 << endl;

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

}

/* Various errors that may occur during execution */
void utilities::SFcollapse1D_error( const int error ) {

  switch (error) {

    case SPHERICAL_USAGE_ERROR:
      cerr << "(SFcollapse1D) ERROR: Incorrect usage of the program!\n";
      cerr << "(SFcollapse1D) Spherical coordinates selected.\n";
      cerr << "(SFcollapse1D) Correct usage is: ./SFcollapse1D Nx Domain_size t_final\n";
      exit(SPHERICAL_USAGE_ERROR);
      break;

    case SINH_SPHERICAL_USAGE_ERROR:
      cerr << "(SFcollapse1D) ERROR: Incorrect usage of the program!\n";
      cerr << "(SFcollapse1D) SinhSpherical coordinates selected.\n";
      cerr << "(SFcollapse1D) Correct usage is: ./SFcollapse1D Nx Domain_size t_final sinhW\n";
      exit(SINH_SPHERICAL_USAGE_ERROR);
      break;

    case GRID_STRUCTURE_ERROR:
      cerr << "(initialize_parameters) ERROR: Unknown grid structure! Please check the macros.hpp file.\n";
      exit(GRID_STRUCTURE_ERROR);
      break;

    case COORD_SYSTEM_ERROR:
      cerr << "(initialize_parameters) ERROR: Unknown coordinate system! Please check the macros.hpp file!\n";
      exit(COORD_SYSTEM_ERROR);
      break;

    case GRIDFUNCTION_TIMELEVEL_ERROR:
      cerr << "(shift_timelevels) ERROR: Unknown value of which_level inside shift_timelevels!\n";
      exit(GRIDFUNCTION_TIMELEVEL_ERROR);
      break;

    case BISECTION_INTERVAL_ERROR:
      cerr << "(Bisection_index_finder) ERROR: bad interval!\n";
      exit(BISECTION_INTERVAL_ERROR);
      break;

    case BISECTION_CONVERGENCE_ERROR:
      cerr << "(Bisection_index_finder) ERROR: did not converge!\n";
      exit(BISECTION_CONVERGENCE_ERROR);
      break;

  }

}
