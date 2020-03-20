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

/* Function prototypes needed by this file */

namespace evolution {

  /* Function to set the initial condition for all gridfunctions: phi, Phi, Pi, a, and alpha */
  void initial_condition( grid::parameters grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha );

  /* Function to step phi, Phi, and Pi forward in time */
  void time_step_scalarfield_gridfunctions( const int n, const grid::parameters grid,
					    const std::vector<REAL> phi_n, const std::vector<REAL> Phi_n, const std::vector<REAL> Pi_n, const std::vector<REAL> a_n, const std::vector<REAL> alpha_n,
					    const std::vector<REAL> Phi_nm1, const std::vector<REAL> Pi_nm1, const std::vector<REAL> a_nm1, const std::vector<REAL> alpha_nm1, 
					    std::vector<REAL> &Phi_np1, std::vector<REAL> &Pi_np1, std::vector<REAL> &phi_np1 );

  /* Function to apply outgoing radiation boundary conditions to phi, Phi, and Pi */
  void apply_outgoing_radiation_bdry_cond( const int n, grid::parameters grid,
					   const std::vector<REAL> phi_nm1, const std::vector<REAL> Pi_nm1,
					   const std::vector<REAL> phi_n  , const std::vector<REAL> Phi_n, const std::vector<REAL> a_n, const std::vector<REAL> alpha_n,
					   std::vector<REAL> &phi_np1, std::vector<REAL> &Phi_np1, std::vector<REAL> &Pi_np1 );

  /* Function to solve the Hamiltonian constraint */
  REAL pointwise_solution_of_the_Hamiltonian_constraint( const int j, grid::parameters grid, const std::vector<REAL> Phi, const std::vector<REAL> Pi, const std::vector<REAL> a );

  /* Function to solve the polar slicing condition */
  REAL pointwise_solution_of_the_polar_slicing_condition( const int j, grid::parameters grid, const std::vector<REAL> a, const std::vector<REAL> alpha );

  /* Function to compute and output the mass aspect ratio */
  void compute_and_output_mass_aspect_ratio( const int, const int, const grid::parameters, const gridfunction );

  /* Function to perform the rescaling of the lapse function */
  void rescaling_of_the_lapse( grid::parameters, const std::vector<REAL>, std::vector<REAL> & );

}

#endif // __EVOLUTION_HPP__
