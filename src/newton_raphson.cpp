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
#include "grid.hpp"
#include "gridfunction.hpp"

using namespace std;

/* Function prototypes needed by this file */
REAL pointwise_Newton_method_Spherical( const int j, grid::parameters grid, const std::vector<REAL> Phi, const std::vector<REAL> Pi, const std::vector<REAL> a ) {

  DECLARE_GRID_PARAMETERS;

  /* Set auxiliary variables */
  const REAL A         = log(a[j-1]);               // A^{n+1}_{j}
  const REAL avgPhi    = 0.5*( Phi[j] + Phi[j-1] ); // 0.5*( Phi^{n+1}_{j+1} + Phi^{n+1}_{j+1} )
  const REAL avgPi     = 0.5*( Pi[j]  + Pi[j-1]  ); // 0.5*(  Pi^{n+1}_{j+1} +  Pi^{n+1}_{j+1} )
  const REAL PhiSqr    = SQR(avgPhi);
  const REAL PiSqr     = SQR(avgPi);
  const REAL PhiPiTerm = M_PI * ( x[0][j] + x[0][j-1] ) * ( PhiSqr + PiSqr );
  const REAL half_invr = 0.5*inv_dx0/(j-0.5);

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
    const REAL f  = inv_dx0*(A_old - A) + tmp0 - half_invr - PhiPiTerm;
    const REAL df = inv_dx0 + tmp0;

    /* Update A_new */
    A_new = A_old - f/df;

    /* Increment the iterator */
    iter++;

  } while( ( fabs(A_new - A_old) > NEWTON_TOL ) && ( iter <= NEWTON_MAX_ITER ) );

  /* Check for convergence */
  if( iter > NEWTON_MAX_ITER ) cerr << "WARNING: Newton's method did not converge to a root! j = " << j << " | iter = " << iter << endl;

  /* Return the value of a */
  return( exp(A_new) );

}

REAL pointwise_Newton_method_SinhSpherical( const int j, grid::parameters grid, const std::vector<REAL> Phi, const std::vector<REAL> Pi, const std::vector<REAL> a ) {

  DECLARE_GRID_PARAMETERS;

  /* Set auxiliary variables */
  const REAL A         = log(a[j-1]);               // A^{n+1}_{j}
  const REAL avgPhi    = 0.5*( Phi[j] + Phi[j-1] ); // 0.5*( Phi^{n+1}_{j+1} + Phi^{n+1}_{j+1} )
  const REAL avgPi     = 0.5*( Pi[j]  + Pi[j-1]  ); // 0.5*(  Pi^{n+1}_{j+1} +  Pi^{n+1}_{j+1} )
  const REAL PhiSqr    = SQR(avgPhi);
  const REAL PiSqr     = SQR(avgPi);
  const REAL midx0     = ( x[0][j] + x[0][j-1] ) / 2.0;
  const REAL PhiPiTerm = 2.0 * M_PI * (SQR(sinhA) / sinhW) * sinh(midx0*inv_sinhW) * cosh(midx0*inv_sinhW) / SQR(sinh_inv_W) * ( PhiSqr + PiSqr );
  const REAL half_invr = 0.5/( sinhW * tanh(x[0][j]) );

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
    const REAL f    = inv_dx0*(A_old - A) + tmp0 - half_invr - PhiPiTerm;
    const REAL df   = inv_dx0 + tmp0;

    /* Update A_new */
    A_new = A_old - f/df;

    /* Increment the iterator */
    iter++;

  } while( ( fabs(A_new - A_old) > NEWTON_TOL ) && ( iter <= NEWTON_MAX_ITER ) );

  /* Check for convergence */
  if( iter > NEWTON_MAX_ITER ) cerr << "WARNING: Newton's method did not converge to a root! j = " << j << " | iter = " << iter << endl;

  /* Return the value of a */
  return( exp(A_new) );

}
