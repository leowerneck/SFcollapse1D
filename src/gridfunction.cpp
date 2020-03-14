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
#include <iostream>
#include <fstream>
#include <string>
#include "grid.hpp"
#include "gridfunction.hpp"

using namespace std;

/* .---------------------------------------.
 * | The initialize_gridfunctions function |
 * .---------------------------------------.
 *
 * This function is used by the constructor to initialize
 * the gridfunction class. Quite simply, it receives a
 * vector size and resizes the three vectors in the
 * gridfunctions to be of that size, initializing the
 * entries of the vectors to zero.
 *
 * The value of size in our program is typically
 * grid.NxTotal (See simulation_parameters.hpp).
 */
void gridfunction::initialize_gridfunction( const grid::parameters grid ) {

  /* Resize the gridfunction to size and initialize them to zero */
  level_np1.resize(grid.Nx0Total,0.0);
  level_n.resize(  grid.Nx0Total,0.0);
  level_nm1.resize(grid.Nx0Total,0.0);

}

void gridfunction::shift_timelevels( const int which_levels ) {

  if( which_levels == 2 ) {
    level_nm1 = level_n;
    level_n   = level_np1;
  }
  else if( which_levels == 1 ) {
    level_n = level_np1;
  }
  else{
    utilities::SFcollapse1D_error( GRIDFUNCTION_TIMELEVEL_ERROR );
  }

}

void gridfunction::output_to_file( const grid::parameters grid, const string basename, const int which_level, const int n ) {

  DECLARE_GRID_PARAMETERS;

  ofstream outfile;
  outfile.open("out/"+basename+"_"+string(to_string(Nt).length() - to_string(n).length(),'0')+to_string(n)+".dat");
  outfile.precision(15);
  LOOP(0,Nx0Total) {
    outfile << scientific << x[0][j] << " " << r_ito_x0[j]  << " " << (which_level == -1) * level_nm1[j] + (which_level == 0) * level_n[j] + (which_level == 1) * level_np1[j] << endl;
  }
  outfile.close();

}
