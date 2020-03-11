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
#include <cmath>
#include <vector>
#include "macros.hpp"
#include "grid.hpp"
#include "gridfunction.hpp"

/* Function prototypes needed by this file */

namespace evolution {

  /* First time step driver function prototype */
  void first_time_step( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha );

  /* Generic time step driver function prototype */
  void time_step( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha );

  /* First time step in Spherical coordinates prototype */
  void first_time_step_Spherical( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha );

  /* Generic time step in Spherical coordinates prototype */
  void time_step_Spherical( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha );

  /* First time step in SinhSpherical coordinates prototype */
  void first_time_step_SinhSpherical( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha );

  /* Generic time step in SinhSpherical coordinates prototype */
  void time_step_SinhSpherical( grid::parameters &grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha );

  /* .---------------------------.
   * | Compute mass-aspect ratio |
   * .---------------------------.
   */
  inline REAL mass_aspect_ratio( const REAL, const REAL );

  /* Print mass_aspect_ratio */
  void output_mass_aspect_ratio( const int, const int, const grid::parameters, const gridfunction );

}
