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

void rescaling_of_the_lapse( grid::parameters grid, const std::vector<REAL> a, std::vector<REAL> &alpha ) {

  DECLARE_GRID_PARAMETERS;

  /* Set the initial value of kappa */
  REAL kappa = a[0]/alpha[0];

  /* Loop over the grid, updating kappa if needed */
  LOOP(1,Nx0Total) {
    REAL kappa_new = a[j]/alpha[j];
    if( kappa_new < kappa )
      kappa = kappa_new;
  }

  /* Rescale the lapse */
  LOOP(0,Nx0Total) alpha[j] *= kappa;

}
