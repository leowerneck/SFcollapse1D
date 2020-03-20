/* .-----------------------------------------------------------------------.
 * | SFcollapse1D                                                          |
 * | Gravitational collapse of scalar fields in spherical symmetry         |
 * |                                                                       |
 * | Copyright (c) 2020 Leonardo Werneck                                   |
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

#ifndef __MACROS_HPP__
#define __MACROS_HPP__

/* Set the REAL variable type */
#define REAL double

/* Set the coordinate system */
#define SPHERICAL      (0)
#define SINH_SPHERICAL (1)
#define COORD_SYSTEM   SINH_SPHERICAL

/* Set cell or vertex centered grid macros */
#define CELL_CENTERED (0)
#define VERT_CENTERED (1)
#define GRID_STRUCTURE VERT_CENTERED

/* Set CFL factor */
#define CFL_FACTOR (0.5)

/* Various error codes */
#define SPHERICAL_USAGE_ERROR        (1)
#define SINH_SPHERICAL_USAGE_ERROR   (2)
#define GRID_STRUCTURE_ERROR         (3)
#define COORD_SYSTEM_ERROR           (4)
#define REGRIDDING_OPTION_ERROR      (5)
#define GRIDFUNCTION_TIMELEVEL_ERROR (6)
#define BISECTION_INTERVAL_ERROR     (7)
#define BISECTION_CONVERGENCE_ERROR  (8)
#define NAN_ERROR                    (9)

/* Set ghostzones. The default is the same
 * number of ghostzones in every direction
 */
#define NGHOSTS  (0)
#define NGHOSTS0 NGHOSTS

/* Scalar field collapse parameters */
#define PHI0_EXAMPLE_WEAK_FIELD (0.100)
#define PHI0_EXAMPLE_INTR_FIELD (0.300)
#define PHI0_EXAMPLE_STRG_FIELD (0.400)
#define PHI0  PHI0_EXAMPLE_STRG_FIELD
#define R0    (0.0)
#define DELTA (1.0)

/* Square macro */
#define SQR(x) ( (x)*(x) )

/* Cube macro */
#define CBD(x) ( (x)*(x)*(x) )

/* Newton's method parameters */
#define NEWTON_TOL      (1e-7)
#define NEWTON_MAX_ITER (100)

/* Regridding parameters */
#define REGRID_RADIAL_POINTS       (1)
#define REGRID_OUTER_BOUNDARY      (2)
#define REGRID_POINT_DENSITY       (3)
#define REGRID_OPTION              REGRID_POINT_DENSITY
#define REGRID_FACTOR              (0.9) // Should be < 1 for options 2 and 3, and > 1 for option 1
#define MAX_REGRID_LEVELS          (10)
#define REGRID_INTERP_STENCIL_SIZE (4)

/* Checkpoints */
#define OUTPUT_CHECKPOINT         (200)
#define REGRID_CHECKER_CHECKPOINT (20)
#define NAN_CHECKER_CHECKPOINT    (20)
#define INFORMATION_CHECKPOINT    (50)

/* Set the declare grid parameters macro. This macro
 * is useful so that we don't need to keep appending
 * "grid." to the variables inside the grid_parameters
 * class. This also allows us to have more control over
 * when the parameters in the grid_parameters class are
 * modified, since we must then append "grid." to the
 * variable name, otherwise it won't be updated properly.
 */
#define DECLARE_GRID_PARAMETERS                                                                 \
  const int DIM                            __attribute__((unused)) = grid.DIM;                  \
  const int Nx0                            __attribute__((unused)) = grid.Nx0;                  \
  const int Ngz0                           __attribute__((unused)) = grid.Ngz0;                 \
  const int Nx0Total                       __attribute__((unused)) = grid.Nx0Total;             \
  const int Nt                             __attribute__((unused)) = grid.Nt;                   \
  const int current_regrid_level           __attribute__((unused)) = grid.current_regrid_level; \
  const int max_regrid_levels              __attribute__((unused)) = grid.max_regrid_levels;    \
  const REAL x0_min                        __attribute__((unused)) = grid.x0_min;               \
  const REAL x0_max                        __attribute__((unused)) = grid.x0_max;               \
  const REAL dx0                           __attribute__((unused)) = grid.dx0;                  \
  const REAL ds_min                        __attribute__((unused)) = grid.ds_min;               \
  const REAL dt                            __attribute__((unused)) = grid.dt;                   \
  const REAL inv_dx0                       __attribute__((unused)) = grid.inv_dx0;              \
  const REAL inv_dx0_sqrd                  __attribute__((unused)) = grid.inv_dx0_sqrd;         \
  const REAL t_initial                     __attribute__((unused)) = grid.t_initial;            \
  const REAL t_final                       __attribute__((unused)) = grid.t_final;              \
  const REAL t                             __attribute__((unused)) = grid.t;                    \
  const REAL sinhA                         __attribute__((unused)) = grid.sinhA;                \
  const REAL sinhW                         __attribute__((unused)) = grid.sinhW;                \
  const REAL inv_sinhW                     __attribute__((unused)) = grid.inv_sinhW;            \
  const REAL sinh_inv_W                    __attribute__((unused)) = grid.sinh_inv_W;           \
  const REAL A_over_sinh_inv_W             __attribute__((unused)) = grid.A_over_sinh_inv_W;    \
  const std::vector<REAL> r_ito_x0         __attribute__((unused)) = grid.r_ito_x0;             \
  const std::vector< std::vector<REAL> > x __attribute__((unused)) = grid.x;

/* Set the LOOP macro */
#define LOOP(jmin,jmax) for(int j=jmin;j<jmax;j++)

/* Information macro */
#define INTEGRATION_INFO					                                                                 \
  {                                                                                                                              \
    auto current_time = std::chrono::high_resolution_clock::now();                                                               \
    auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>( current_time - start_time ).count();                   \
    REAL time_elapsed = (REAL)elapsed_time;                                                                                      \
    REAL time_left    = ((REAL)grid.Nt/(REAL)n - 1.0) * time_elapsed;		                                                 \
    std::cout << fixed                                                                                                           \
	      << "\r(SFcollapse1D INFO) " << "Iter " << setfill('0') << setw(to_string(grid.Nt).length()) << n << "/" << grid.Nt \
	      << " | t = " << setprecision(3) << grid.t                                                                          \
	      << " |  Runtime: " << setprecision(0) << time_elapsed                                                              \
	      << " seconds | ETA: " << setprecision(0) << time_left << " seconds" << std::flush;                                 \
  }
  
#endif // __MACROS_HPP__
