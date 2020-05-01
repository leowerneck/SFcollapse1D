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

#ifndef __GRIDFUNCTION_HPP__
#define __GRIDFUNCTION_HPP__

/* Basic includes */
#include <string>
#include "macros.hpp"
#include "grid.hpp"

/* Define the gridfunction class */
class gridfunction {
public:
  /* .--------------------.
   * | Gridfunction array |
   * .--------------------.
   *
   * We need three time levels
   */
  realvec level_np1;
  realvec level_n;
  realvec level_nm1;
    
  /* .-----------------.
   * | The constructor |
   * .-----------------.
   */
  gridfunction( const int N ) {
    initialize_gridfunction(N);
  }

  /* .--------------------------------.
   * | The shift time levels function |
   * .--------------------------------.
   */
  void shift_timelevels( const int );

  /* .----------------------.
   * | Output level to file |
   * .----------------------.
   */
  void output_to_file( const grid::parameters, const std::string, const int, const int );

private:
  /* .---------------------------------------.
   * | The initialize_gridfunctions function |
   * .---------------------------------------.
   */
  void initialize_gridfunction( const int );
};

#endif // __GRIDFUNCTION_HPP__
