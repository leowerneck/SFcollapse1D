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

#ifndef __EVOLUTION_HPP__
#define __EVOLUTION_HPP__

/* Basic includes */
#include <cmath>
#include "macros.hpp"
#include "grid.hpp"
#include "gridfunction.hpp"

namespace evolution {

  /* Function to set the initial condition for all gridfunctions: phi, Phi, Pi, a, and alpha */
  void initial_condition( grid::parameters, gridfunction &, gridfunction &, gridfunction &, gridfunction &, gridfunction & );

  /* Function to step phi, Phi, and Pi forward in time */
  void time_step_scalarfield_gridfunctions( const int, const grid::parameters,
					    const realvec, const realvec, const realvec, const realvec, const realvec,
					    const realvec, const realvec, const realvec, const realvec, 
					    realvec &, realvec &, realvec & );

  /* Function to apply outgoing radiation boundary conditions to phi, Phi, and Pi */
  void apply_outgoing_radiation_bdry_cond( const int, grid::parameters,
					   const realvec, const realvec,
					   const realvec, const realvec, const realvec, const realvec,
					   realvec &, realvec &, realvec & );

  /* Function to solve the Hamiltonian constraint */
  real pointwise_solution_of_the_Hamiltonian_constraint( const int, grid::parameters, const realvec, const realvec, const realvec );

  /* Function to solve the polar slicing condition */
  real pointwise_solution_of_the_polar_slicing_condition( const int, grid::parameters, const realvec, const realvec );

  /* Function to perform the rescaling of the lapse function */
  void rescaling_of_the_lapse( grid::parameters, const realvec, realvec & );

}

#endif // __EVOLUTION_HPP__
