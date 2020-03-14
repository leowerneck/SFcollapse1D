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
#include <cmath>
#include <chrono>
#include "macros.hpp"
#include "grid.hpp"
#include "gridfunction.hpp"
#include "evolution.hpp"

using namespace std;

/* Function prototypes needed by this file */
void initial_condition( grid::parameters grid, gridfunction &phi, gridfunction &Phi, gridfunction &Pi, gridfunction &a, gridfunction &alpha );

int main( int argc, char *argv[] ) {

  /* Print logo to the user */
#include "logo.hpp"

  /* Check correct usage */
#if( COORD_SYSTEM == SPHERICAL )
  if( argc != 4 ) {
    cerr << "ERROR: Incorrect usage of the program!\n";
    cerr << "Spherical coordinates selected.\n";
    cerr << "Correct usage is: ./SFcollapse1D Nx Domain_size t_final" << endl;
    exit(1);
  }
#elif( COORD_SYSTEM == SINH_SPHERICAL )
  if( argc != 5 ) {
    cerr << "ERROR: Incorrect usage of the program!\n";
    cerr << "SinhSpherical coordinates selected.\n";
    cerr << "Correct usage is: ./SFcollapse1D Nx Domain_size t_final sinhW" << endl;
    exit(1);
  }
#endif

  /* Start the timer */
  auto start_time = std::chrono::high_resolution_clock::now();

  /* Set the grid of the program */
  grid::parameters grid(argv);

  /* Declare all needed gridfunctions */
  gridfunction phi(grid), Phi(grid), Pi(grid), a(grid), alpha(grid);

  /* Set the initial condition */
  initial_condition( grid, phi, Phi, Pi, a, alpha );

  /* Print information to the user */
  phi.output_to_file(grid,"scalarfield",-1,0);
  Phi.output_to_file(grid,"Phi",-1,0);
  Pi.output_to_file(grid,"Pi",-1,0);
  a.output_to_file(grid,"a",-1,0);
  alpha.output_to_file(grid,"alpha",-1,0);
  evolution::output_mass_aspect_ratio(-1,0,grid,a);

  /* The first time step is special, perform it seperately */
  evolution::first_time_step( grid, phi, Phi, Pi, a, alpha );

  /* Print information to the user */
  phi.output_to_file(grid,"scalarfield",1,1);
  Phi.output_to_file(grid,"Phi",1,1);
  Pi.output_to_file(grid,"Pi",1,1);
  a.output_to_file(grid,"a",1,1);
  alpha.output_to_file(grid,"alpha",1,1);
  evolution::output_mass_aspect_ratio(1,1,grid,a);

  /* Now shift the time levels */
  phi.shift_timelevels(1);
  Phi.shift_timelevels(1);
  Pi.shift_timelevels(1);
  a.shift_timelevels(1);
  alpha.shift_timelevels(1);

  /* Begin evolving the system */
  for(int n=2;n<=grid.Nt;n++) {

    /* Perform an integration step */
    evolution::time_step( grid, phi, Phi, Pi, a, alpha );

    /* Print information to the user */
    if( n%200 == 0 ) {
      phi.output_to_file(grid,"scalarfield",1,n);
      Phi.output_to_file(grid,"Phi",1,n);
      Pi.output_to_file(grid,"Pi",1,n);
      a.output_to_file(grid,"a",1,n);
      alpha.output_to_file(grid,"alpha",1,n);
      evolution::output_mass_aspect_ratio(1,n,grid,a);
    }
    
    /* Shift time levels appropriately */
    phi.shift_timelevels(2);
    Phi.shift_timelevels(2);
    Pi.shift_timelevels(2);
    a.shift_timelevels(2);
    alpha.shift_timelevels(2);

    if( n%20 == 0 ) {
      INTEGRATION_INFO;
    }
    
  }
  
  return 0;

}
