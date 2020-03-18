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

/* Necessary includes */
#include <iostream>
#include <cmath>
#include "macros.hpp"
#include "grid.hpp"

using namespace std;

/* As a reminder, the parameters class contains the following members (all public):
 * .------------------------------------------------------------------------------------------.
 * | Nx0                -> Number of interior grid points in the x0 direction (r)             |
 * .------------------------------------------------------------------------------------------.
 * | Ngz0               -> Number of exterior grid points in the x0 direction (r)             |
 * .------------------------------------------------------------------------------------------.
 * | Nx0Total           -> Nx0+2*Ngz0: Total number of grid points in the x0 direction        |
 * .------------------------------------------------------------------------------------------.
 * | x0_min, x0_max     -> Range of x0, i.e. x0 in [x0_min,x0_max]                            |
 * .------------------------------------------------------------------------------------------.
 * | dx0                -> Step sizes in the x0 direction                                     |
 * .------------------------------------------------------------------------------------------.
 * | t_initial, t_final -> Initial and final values of t. By default, t_initial = 0           |
 * .------------------------------------------------------------------------------------------.
 */
void grid::parameters::initialize_parameters(char *argv[]) {

  /* .-----------------------------------.
   * | Initialize all spatial parameters |
   * .-----------------------------------.
   *
   * Set DIM */
  DIM = 1;
  /* Read in Nx0 */
  Nx0 = atoi(argv[1])+1;
  /* Set Ngz0 */
  Ngz0 = NGHOSTS0;
  /* Set Nx0Total */
  Nx0Total = Nx0 + 2*Ngz0;

  /* For Spherical/SinhSpherical coordinates, we set the following values:
   * .----------.--------------.---------------------.
   * | Variable |    Value     | Physical Coordinate |
   * .----------.--------------.---------------------.
   * |  x0_min  |      0       |        r_min        |
   * .----------.--------------.---------------------.
   * |  x0_max  | +domain_size |        r_max        |
   * .----------.--------------.---------------------.
   */
#if( COORD_SYSTEM == SPHERICAL )
  x0_max = fabs(atof(argv[2]));
#elif( COORD_SYSTEM == SINH_SPHERICAL )
  x0_max            = 1.0;
  sinhA             = fabs(atof(argv[2]));
  sinhW             = fabs(atof(argv[4]));
  inv_sinhW         = 1.0/sinhW;
  sinh_inv_W        = sinh(inv_sinhW);
  A_over_sinh_inv_W = sinhA / sinh_inv_W;
#endif
  x0_min = 0.0;

  /* Set dx0,dx1,dx2 */
  dx0 = (x0_max - x0_min)/((REAL)Nx0-1.0);
  /* Set inv_dx0,inv_dx1,inv_dx2 */
  inv_dx0 = 1.0/dx0;
  /* Set inv_dx0_sqrd,inv_dx1_sqrd,inv_dx2_sqrd */
  inv_dx0_sqrd = inv_dx0*inv_dx0;
  /* Allocate memory for x,y,z */
  x.resize(1);
  x[0].resize(Nx0Total,0.0);

  /* Populate the x,y,z */
#if( GRID_STRUCTURE == CELL_CENTERED )
  /* Cell centered grid structure:
   *
   *      x          x           x           x        ...
   * 0  0.5dr  dr  1.5dr  2dr  2.5dr  3dr  3.5dr  4dr ...
   */
  LOOP(0,Nx0Total) x[0][j] = x0_min + (j-0.5-NGHOSTS0+1)*dx0;

#elif( GRID_STRUCTURE == VERT_CENTERED )
  /* Vertex centered grid structure:
   *
   * x         x           x           x          x   ...
   * 0  0.5dr  dr  1.5dr  2dr  2.5dr  3dr  3.5dr  4dr ...
   */
  LOOP(0,Nx0Total) x[0][j] = x0_min + (j-Ngz0)*dx0;

#else
  utilities::SFcollapse1D_error(GRID_STRUCTURE_ERROR);
#endif

  /* Populate r in terms of x0 */
  r_ito_x0.resize(Nx0Total,0.0);
#if( COORD_SYSTEM == SPHERICAL )
  r_ito_x0 = x[0];
#elif( COORD_SYSTEM == SINH_SPHERICAL )
  LOOP(0,Nx0Total) r_ito_x0[j] = A_over_sinh_inv_W * sinh( x[0][j] * inv_sinhW );
#endif

  /* .--------------------------------.
   * | Initialize all time parameters |
   * .--------------------------------.
   *
   * Set t_initial; read in t_final */
  t_initial = 0.0;
  t_final   = atof(argv[3]);
  t         = t_initial;

  /* Set dt based on CFL factor and min(dr_{j}) */
  ds_min = r_ito_x0[1] - r_ito_x0[0];
  dt = CFL_FACTOR * ds_min;

  /* Set Nt; add 0.5 because C++ rounds down */
  Nt = (int)(t_final/dt + 0.5);

  /* Initialize regrid variables */
  current_regrid_level = 0;
  max_regrid_levels    = MAX_REGRID_LEVELS;

}

/* Spherical coordinates: x(r) and r(x) */
REAL grid::compute_x_of_r( const REAL r ) {
  return r;
}
REAL grid::compute_r_of_x( const REAL x ) {
  return x;
}

/* SinhSpherical coordinates: x(r) and r(x) */
REAL grid::compute_x_of_r( const REAL r, const REAL A, const REAL w ) {
  return( w * asinh( sinh( 1.0/w ) / A * r ) );
}
REAL grid::compute_r_of_x( const REAL x, const REAL A, const REAL w ) {
  return( A * sinh(x/w) / sinh(1.0/w) );
}

/* .--------------------------------------.
 * | The compute_xyz_Cartesian() function |
 * .--------------------------------------.
 */
void grid::compute_xyz_Cartesian(const REAL x0, const REAL x1, const REAL x2, REAL xyz[3]) {

  /* Standard conversion from Spherical
   * to Cartesian coordinates:
   * .-------------------------------.
   * | x = r * sin(theta) * cos(phi) |
   * .-------------------------------.
   * | y = r * sin(theta) * sin(phi) |
   * .-------------------------------.
   * | z = r * cos(theta)            |
   * .-------------------------------.
   */
  xyz[0] = x0 * sin(x1) * cos(x2);
  xyz[1] = x0 * sin(x1) * sin(x2);
  xyz[2] = x0 * cos(x1);
  
}
