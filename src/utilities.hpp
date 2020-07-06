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

#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__

/* Basic includes */
#include <vector>
#include "macros.hpp"
#include "grid.hpp"
#include "gridfunction.hpp"

namespace utilities {

  /* Compute mass-aspect function */
  void compute_and_output_mass_aspect_function( const int, const int, const grid::parameters, const gridfunction );

  /* NaN checker function ) */
  void NaN_checker( const int, grid::parameters, gridfunction, gridfunction, gridfunction, gridfunction, gridfunction );

  /* Function to check whether regridding is necessary or not */
  bool check_regrid_criterion( const grid::parameters, const realvec );

  /* Regridding function */
  void regrid( grid::parameters &, gridfunction &, gridfunction &, gridfunction &, gridfunction &, gridfunction & );

  /* Lagrange interpolator */
  void Lagrange_interpolator( const int interp_index, const int interp_stencil_size, const realvec x,
			      const realvec phi, const realvec Phi, const realvec Pi, const realvec a, const realvec alpha,
			      realvec &phi_star, realvec &Phi_star, realvec &Pi_star, realvec &a_star, realvec &alpha_star,
			      const real x_star );

  /* Bisection index finder */
  int bisection_index_finder( const realvec, const real );

  /* Printing useful grid information to the user */
  void parameter_information( grid::parameters );

  /* Various error messages for SFcollapse1D */
  void SFcollapse1D_error( const int );

  /* Print central values */
  void output_gridfunctions_central_values( const int, const grid::parameters,
					    const realvec, const realvec, const realvec, const realvec, const realvec );

  /* Output energy density to file */
  void output_energy_density_to_file( const grid::parameters, const realvec, const realvec, const realvec, const int );

  /* Lapse collapse checker */
  bool check_for_collapse_of_the_lapse( const grid::parameters, const gridfunction );

}

#endif // __UTILITIES_HPP__
