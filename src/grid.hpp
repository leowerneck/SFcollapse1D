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

#ifndef __GRID__
#define __GRID__

#include <vector>
#include "macros.hpp"

/* Define parameters class */
namespace grid {

  class parameters {
  public:
    /* .--------------------.
     * | Spatial Parameters |
     * .--------------------.
     *
     * Number of spatial dimensions */
    int DIM;
    /* Number of points in each spatial direction */
    int Nx0;
    /* Number of ghostzones in each spatial direction */
    int Ngz0;
    /* Total number of points in each spatial direction */
    int Nx0Total;
    /* Domain size */
    REAL x0_min,x0_max;
    /* Spatial step sizes */
    REAL dx0;
    /* Vector of vectors for x0,x1,x2 */
    std::vector< std::vector<REAL> > x;
    /* Inverse spatial step sizes */
    REAL inv_dx0;
    /* Inverse spatial step sizes squared */
    REAL inv_dx0_sqrd;
    /* Sinh Spherical coordinate system variables */
    REAL sinhA, sinhW, inv_sinhW, sinh_inv_W, A_over_sinh_inv_W;
    std::vector< REAL > r_ito_x0;

    /* .-----------------.
     * | Time Parameters |
     * .-----------------.
     *
     * Number of time steps */
    int Nt;
    /* Initial and final times */
    REAL t_initial,t_final;
    /* Time and time step */
    REAL t,dt;

    /* .-----------------.
     * | The constructor |
     * .-----------------.
     */
    parameters(char *argv[]) {
      initialize_parameters(argv);
    }

  private:
    /* .------------------------------------.
     * | The initialize_parameters function |
     * .------------------------------------.
     */
    void initialize_parameters(char *argv[]);
  };

  /* .-----------------------------.
   * | The compute_xyz_Cartesian() |
   * .-----------------------------.
   */
  void compute_xyz_Cartesian(const REAL x0, const REAL x1, const REAL x2, REAL xyz[3]);

  /* .-----------------------------.
   * | The compute_xyz_Cartesian() |
   * .-----------------------------.
   */
  REAL compute_r_from_x0(const int j);

}

#endif // __GRID__
