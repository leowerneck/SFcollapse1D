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

namespace utilities {

  /* Lagrange interpolator */
  REAL Lagrange_interpolator( const int, const std::vector<REAL>, const std::vector<REAL>, const REAL );

  /* Lagrange basis polynomials*/
  void Lagrange_basis_polynomials( const int, const REAL, const std::vector<REAL>, std::vector<REAL> & );

  /* Bisection index finder */
  int Bisection_index_finder( const std::vector<REAL>, const REAL );

  /* Various error messages for SFcollapse1D */
  void SFcollapse1D_error( const int );

}

#endif // __UTILITIES_HPP__
