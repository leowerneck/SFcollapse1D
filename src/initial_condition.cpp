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
REAL pointwise_Newton_method( const int j, grid::parameters grid, const std::vector<REAL> Phi, const std::vector<REAL> Pi, const std::vector<REAL> a );
void rescaling_of_the_lapse( grid::parameters grid, const std::vector<REAL> a, std::vector<REAL> &alpha );

void initial_condition( grid::parameters grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha ) {

  DECLARE_GRID_PARAMETERS;

  LOOP(0,Nx0Total) {
    /* Set some useful auxiliary variables */
    const REAL r = r_ito_x0[j];
    const REAL factor = (r-R0)/SQR(DELTA);
    const REAL expfactor = (r-R0)*factor;
    const REAL exp_rmr0_over_deltasqrd = exp(-expfactor);

    /* Set the initial condition for Phi and Pi */
    phi.level_nm1[j] = PHI0*exp_rmr0_over_deltasqrd;
    Phi.level_nm1[j] = -2.0*factor*PHI0*exp_rmr0_over_deltasqrd;
    Pi.level_nm1[j]  = 0.0;

    if( j>0 ) {
      /* Compute a */
      a.level_nm1[j] = pointwise_Newton_method(j,grid,Phi.level_nm1,Pi.level_nm1,a.level_nm1);

      /* Compute auxiliary quantities */
      const REAL b = a.level_nm1[j] + a.level_nm1[j-1];
      const REAL c = a.level_nm1[j] - a.level_nm1[j-1];
#if( COORD_SYSTEM == SPHERICAL )      
      const REAL midway_r = 0.5 * (x[0][j] + x[0][j-1]);
#elif( COORD_SYSTEM == SINH_SPHERICAL )
      const REAL midway_r = sinhW * tanh( inv_sinhW * 0.5 * (x[0][j] + x[0][j-1]) );
#else
      cerr << "ERROR! Unknown coordinate system! Check the macros.hpp file!\n";
      exit(5);
#endif
      const REAL d = ( 1.0 - 0.25 * SQR(b) )/( 2.0 * midway_r ) - inv_dx0 * c / b;

      /* Compute alpha */
      alpha.level_nm1[j] = alpha.level_nm1[j-1]*( 1.0 - d*dx0 )/( 1.0 + d*dx0 );
    }
    else {
      /* If at inner boundary, impose a = 1 */
      a.level_nm1[j] = 1.0;
      /* If at inner boundary, impose alpha = 1 */
      alpha.level_nm1[j] = 1.0;
    }
  }

  /* Now rescale alpha */
  rescaling_of_the_lapse(grid,a.level_nm1,alpha.level_nm1);

}
