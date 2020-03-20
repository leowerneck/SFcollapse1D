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
#include <vector>
#include "macros.hpp"
#include "grid.hpp"
#include "gridfunction.hpp"

namespace evolution {

  /* Function to set the initial condition for all gridfunctions: phi, Phi, Pi, a, and alpha */
  void initial_condition( grid::parameters, gridfunction &, gridfunction &, gridfunction &, gridfunction &, gridfunction & );

  /* Function to step phi, Phi, and Pi forward in time */
  void time_step_scalarfield_gridfunctions( const int, const grid::parameters,
					    const std::vector<REAL>, const std::vector<REAL>, const std::vector<REAL>, const std::vector<REAL>, const std::vector<REAL>,
					    const std::vector<REAL>, const std::vector<REAL>, const std::vector<REAL>, const std::vector<REAL>, 
					    std::vector<REAL> &, std::vector<REAL> &, std::vector<REAL> & );

  /* Function to apply outgoing radiation boundary conditions to phi, Phi, and Pi */
  void apply_outgoing_radiation_bdry_cond( const int, grid::parameters,
					   const std::vector<REAL>, const std::vector<REAL>,
					   const std::vector<REAL>, const std::vector<REAL>, const std::vector<REAL>, const std::vector<REAL>,
					   std::vector<REAL> &, std::vector<REAL> &, std::vector<REAL> & );

  /* Function to solve the Hamiltonian constraint */
  REAL pointwise_solution_of_the_Hamiltonian_constraint( const int, grid::parameters, const std::vector<REAL>, const std::vector<REAL>, const std::vector<REAL> );

  /* Function to solve the polar slicing condition */
  REAL pointwise_solution_of_the_polar_slicing_condition( const int, grid::parameters, const std::vector<REAL>, const std::vector<REAL> );

  /* Function to perform the rescaling of the lapse function */
  void rescaling_of_the_lapse( grid::parameters, const std::vector<REAL>, std::vector<REAL> & );

  /* Function to compute and output the mass aspect ratio */
  void compute_and_output_mass_aspect_ratio( const int, const int, const grid::parameters, const gridfunction );

}

#endif // __EVOLUTION_HPP__
