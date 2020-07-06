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

#ifndef __GRID_HPP__
#define __GRID_HPP__

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
    real x0_min,x0_max;
    /* Vector of vectors for x0,x1,x2 */
    realmat x;
    /* Step sizes */
    real dx0, ds_min;
    /* Inverse spatial step sizes */
    real inv_dx0;
    /* Inverse spatial step sizes squared */
    real inv_dx0_sqrd;
    /* Sinh Spherical coordinate system variables */
    real sinhA, sinhW, inv_sinhW, sinh_inv_W, A_over_sinh_inv_W;
    realvec r_ito_x0;
    /* Regridding parameters */
    int current_regrid_level, max_regrid_levels;
    /* Amplitude of initial condition */
    real phi0;

    /* .-----------------.
     * | Time Parameters |
     * .-----------------.
     *
     * Number of time steps */
    int Nt;
    /* Initial and final times */
    real t_initial,t_final;
    /* Time and time step */
    real t,dt;

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

  /* Set functions to compute x(r) and r(x). Note that
   * we choose to overload the functions that compute
   * these quantities with multiple definitions, one
   * that works for Spherical coordinates and one that
   * works for SinhSpherical coordinates. There is no
   * confusion because the function call is clearly
   * different.
   */
  real compute_x_of_r( const real r );
  real compute_r_of_x( const real x );
  real compute_x_of_r( const real r, const real A, const real w );
  real compute_r_of_x( const real x, const real A, const real w );
  
  /* .-----------------------------.
   * | The compute_xyz_Cartesian() |
   * .-----------------------------.
   */
  void compute_xyz_Cartesian(const real x0, const real x1, const real x2, real xyz[3]);

  /* .-----------------------------.
   * | The compute_xyz_Cartesian() |
   * .-----------------------------.
   */
  real compute_r_from_x0(const int j);

}

#endif // __GRID_HPP__
