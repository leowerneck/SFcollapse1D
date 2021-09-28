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
#include <iomanip>
#include <fstream>
#include <cmath>
#include <chrono>
#include "macros.hpp"
#include "utilities.hpp"
#include "grid.hpp"
#include "gridfunction.hpp"
#include "evolution.hpp"

using namespace std;

int main( int argc, char *argv[] ) {

  /* Print logo to the user */
#include "logo.hpp"

  /* Check correct usage */
#if( COORD_SYSTEM == SPHERICAL )
  if( argc != 5 ) utilities::SFcollapse1D_error(SPHERICAL_USAGE_ERROR);
#elif( COORD_SYSTEM == SINH_SPHERICAL )
  if( argc != 6 ) utilities::SFcollapse1D_error(SINH_SPHERICAL_USAGE_ERROR);
#endif

  /* Start the timer */
  auto start_time = std::chrono::high_resolution_clock::now();

  /* Construct the base grid */
  grid::parameters grid(argv);
  
  /* Print information about the run to the user */
  utilities::parameter_information(grid);

  /* Declare all needed gridfunctions */
  gridfunction phi(grid.Nx0Total), Phi(grid.Nx0Total), Pi(grid.Nx0Total), a(grid.Nx0Total), alpha(grid.Nx0Total);

  /* Set the initial condition */
  evolution::initial_condition( grid, phi, Phi, Pi, a, alpha );

  // Estimate ADM mass for BSSN computation
  // real M_ADM = 0.0;
  // cout << "M_ADM = " << scientific << setprecision(15) << M_ADM << endl;
  // LOOP(0,grid.Nx0Total-1) {
  //   const real r_left    = grid.r_ito_x0[j];
  //   const real r_right   = grid.r_ito_x0[j+1];
  //   const real dr        = r_right - r_left;
  //   const real tmp_left  = sqrt(2.0*M_PI) * r_left / a.level_nm1[j];
  //   const real tmp_right = sqrt(2.0*M_PI) * r_right / a.level_nm1[j+1];
  //   const real X_left    = tmp_left * Phi.level_nm1[j];
  //   const real X_right   = tmp_right * Phi.level_nm1[j+1];
  //   const real Y_left    = tmp_left * Pi.level_nm1[j];
  //   const real Y_right   = tmp_right * Pi.level_nm1[j+1];
  //   const real I_left    = SQR(X_left) + SQR(Y_left);
  //   const real I_right   = SQR(X_right) + SQR(Y_right);

  //   M_ADM += 0.5 * dr * ( I_left + I_right );

  //   cout << scientific << setprecision(15) << "r = " << 0.5*(r_left+r_right) << " | M_ADM = " << M_ADM << endl;
  // }

  // utilities::compute_and_output_mass_aspect_function(-1,0,grid,a);
  
  // return 0;

  /* Print information to the user */
  phi.output_to_file(grid,"scalarfield",-1,0);
  Phi.output_to_file(grid,"Phi",-1,0);
  Pi.output_to_file(grid,"Pi",-1,0);
  a.output_to_file(grid,"a",-1,0);
  alpha.output_to_file(grid,"alpha",-1,0);
  utilities::compute_and_output_mass_aspect_function(-1,0,grid,a);
  // utilities::output_energy_density_to_file( grid, Phi.level_nm1, Pi.level_nm1, a.level_nm1, 0 );

  utilities::output_gridfunctions_central_values( 0, grid, phi.level_nm1, Phi.level_nm1, Pi.level_nm1, a.level_nm1, alpha.level_nm1 );

  /* The first time step is special, requires two half-integrations */

  /* Perform an integration step */
  /* .--------------------------------------------------------------.
   * | Step 1: Apply inner boundary conditions to Phi, a, and alpha |
   * .--------------------------------------------------------------.
   */
  Phi.level_n[0]   = 0.0;
  a.level_n[0]     = 1.0;
  alpha.level_n[0] = 1.0;

  /* .-----------------------------------.
   * | Step 2: Integrate phi, Phi and Pi |
   * .-----------------------------------.
   */
  evolution::time_step_scalarfield_gridfunctions( 0, grid, 
						  phi.level_nm1, Phi.level_nm1, Pi.level_nm1, a.level_nm1, alpha.level_nm1, 
						  Phi.level_nm1, Pi.level_nm1, a.level_nm1, alpha.level_nm1,
						  Phi.level_n  , Pi.level_n  , phi.level_n );
  /* Step 2.d: Apply boundary conditions */
  evolution::apply_outgoing_radiation_bdry_cond( 0, grid,
						 phi.level_nm1, Pi.level_nm1,
						 phi.level_nm1, Phi.level_nm1, a.level_nm1, alpha.level_nm1,
						 phi.level_n  , Phi.level_n  , Pi.level_n );

  Pi.level_n[0] = - Pi.level_nm1[0] + Pi.level_n[1] + Pi.level_nm1[1];
  /* .-------------------------------.
   * | Step 3: Integrate a and alpha |
   * .-------------------------------.
   */
  // TODO: can this be parallelized in some way? Maybe SIMD?
  LOOP(1,grid.Nx0Total) {
    /* Step 3.a: Compute a */
    a.level_n[j] = evolution::pointwise_solution_of_the_Hamiltonian_constraint(j,grid,Phi.level_n,Pi.level_n,a.level_n);
    
    /* Step 3.b: Compute alpha */
    alpha.level_n[j] = evolution::pointwise_solution_of_the_polar_slicing_condition( j, grid, a.level_n, alpha.level_n );
  }
  /* Step 3.d: Now rescale alpha */
  evolution::rescaling_of_the_lapse(grid,a.level_n,alpha.level_n);
  
  /* .-------------------------.
   * | Step 4: Update the time |
   * .-------------------------.
   */
  grid.t += 0.5 * grid.dt;
  
  /* Perform an integration step */
  /* .--------------------------------------------------------------.
   * | Step 1: Apply inner boundary conditions to Phi, a, and alpha |
   * .--------------------------------------------------------------.
   */
  Phi.level_np1[0]   = 0.0;
  a.level_np1[0]     = 1.0;
  alpha.level_np1[0] = 1.0;
  
  /* .-----------------------------------.
   * | Step 2: Integrate phi, Phi and Pi |
   * .-----------------------------------.
   */
  evolution::time_step_scalarfield_gridfunctions( 1, grid, 
						  phi.level_n, Phi.level_n  , Pi.level_n  , a.level_n  , alpha.level_n  , 
						  Phi.level_nm1, Pi.level_nm1, a.level_nm1, alpha.level_nm1,
						  Phi.level_np1, Pi.level_np1, phi.level_np1 );
  /* Step 2.d: Apply boundary conditions */
  evolution::apply_outgoing_radiation_bdry_cond( 1, grid,
						 phi.level_nm1, Pi.level_nm1,
						 phi.level_n  , Phi.level_n  , a.level_n, alpha.level_n,
						 phi.level_np1, Phi.level_np1, Pi.level_np1);
  
  Pi.level_np1[0] = - Pi.level_n[0] + Pi.level_np1[1] + Pi.level_n[1];
  /* .-------------------------------.
   * | Step 3: Integrate a and alpha |
   * .-------------------------------.
   */
  // TODO: can this be parallelized in some way? Maybe SIMD?
  LOOP(1,grid.Nx0Total) {
    /* Step 3.a: Compute a */
    a.level_np1[j] = evolution::pointwise_solution_of_the_Hamiltonian_constraint(j,grid,Phi.level_np1,Pi.level_np1,a.level_np1);
    
    /* Step 3.b: Compute alpha */
    alpha.level_np1[j] = evolution::pointwise_solution_of_the_polar_slicing_condition( j, grid, a.level_np1, alpha.level_np1 );
  }
  /* Step 3.d: Now rescale alpha */
  evolution::rescaling_of_the_lapse(grid,a.level_np1,alpha.level_np1);

  /* .-------------------------.
   * | Step 4: Update the time |
   * .-------------------------.
   */
  grid.t += 0.5 * grid.dt;

  /* Print information to the user */
  // phi.output_to_file(grid,"scalarfield",1,1);
  // Phi.output_to_file(grid,"Phi",1,1);
  // Pi.output_to_file(grid,"Pi",1,1);
  // a.output_to_file(grid,"a",1,1);
  // alpha.output_to_file(grid,"alpha",1,1);
  // utilities::compute_and_output_mass_aspect_function(1,1,grid,a);
  // utilities::output_energy_density_to_file( grid, Phi.level_n, Pi.level_n, a.level_n, 1 );

  utilities::output_gridfunctions_central_values( 1, grid, phi.level_n, Phi.level_n, Pi.level_n, a.level_n, alpha.level_n );

  /* Define the central density */
  real max_central_density = 0.0;
  {
    const real central_Phi_sqrd = SQR( Phi.level_np1[0] );
    const real central_Pi_sqrd  = SQR( Pi.level_np1[0]  );
    const real central_a_sqrd   = SQR( a.level_np1[0]   );
    const real central_density  = 0.5 * ( central_Phi_sqrd + central_Pi_sqrd ) / central_a_sqrd;
    if( central_density > max_central_density) {
      max_central_density = central_density;
    }
  }

  /* Now shift the time levels */
  phi.shift_timelevels(1);
  Phi.shift_timelevels(1);
  Pi.shift_timelevels(1);
  a.shift_timelevels(1);
  alpha.shift_timelevels(1);

  int n = 2;

  bool lapse_collapsed = false;

  /* Begin evolving the ADM+EKG equations */
  while( grid.t < grid.t_final ) {

    /* Perform an integration step */
    /* .--------------------------------------------------------------.
     * | Step 1: Apply inner boundary conditions to Phi, a, and alpha |
     * .--------------------------------------------------------------.
     */
    Phi.level_np1[0]   = 0.0;
    a.level_np1[0]     = 1.0;
    alpha.level_np1[0] = 1.0;

    /* .-----------------------------------.
     * | Step 2: Integrate phi, Phi and Pi |
     * .-----------------------------------.
     */
    evolution::time_step_scalarfield_gridfunctions( n, grid, 
						    phi.level_n, Phi.level_n  , Pi.level_n  , a.level_n  , alpha.level_n  , 
						    Phi.level_nm1, Pi.level_nm1, a.level_nm1, alpha.level_nm1,
						    Phi.level_np1, Pi.level_np1, phi.level_np1 );
    /* Step 2.a: Apply boundary conditions */
    evolution::apply_outgoing_radiation_bdry_cond( n, grid,
						   phi.level_nm1, Pi.level_nm1,
						   phi.level_n  , Phi.level_n  , a.level_n, alpha.level_n,
						   phi.level_np1, Phi.level_np1, Pi.level_np1);

    Pi.level_np1[0] = - Pi.level_n[0] + Pi.level_np1[1] + Pi.level_n[1];
    /* .-------------------------------.
     * | Step 3: Integrate a and alpha |
     * .-------------------------------.
     */
    // TODO: can this be parallelized in some way? Maybe SIMD?
    LOOP(1,grid.Nx0Total) {
      /* Step 3.a: Compute a */
      a.level_np1[j] = evolution::pointwise_solution_of_the_Hamiltonian_constraint(j,grid,Phi.level_np1,Pi.level_np1,a.level_np1);
    
      /* Step 3.b: Compute alpha */
      alpha.level_np1[j] = evolution::pointwise_solution_of_the_polar_slicing_condition( j, grid, a.level_np1, alpha.level_np1 );
    }
    /* Step 3.d: Now rescale alpha */
    evolution::rescaling_of_the_lapse(grid,a.level_np1,alpha.level_np1);

    /* .-------------------------.
     * | Step 4: Update the time |
     * .-------------------------.
     */
    grid.t += grid.dt;

    /* Check whether or not a regrid is necessary */
    // if( n%REGRID_CHECKER_CHECKPOINT == 0 && grid.current_regrid_level < grid.max_regrid_levels ) {
    // if( grid.current_regrid_level < grid.max_regrid_levels ) {
    //   const bool need_to_regrid = utilities::check_regrid_criterion( grid, phi.level_np1 );
    //   if( need_to_regrid == true ) {
    // 	cout << "\n(SFcollapse1D INFO) Regridding at iteration " << n << endl;
    // 	utilities::regrid( grid, phi, Phi, Pi, a, alpha );
    // 	utilities::parameter_information(phi0,grid);
    // 	n++;
    //   }
    // }

    /* Check for NaNs */
    if( n%NAN_CHECKER_CHECKPOINT == 0 ) utilities::NaN_checker( n, grid, phi, Phi, Pi, a, alpha );

    // /* Compute the central density */
    const real central_Phi_sqrd = SQR( Phi.level_np1[0] );
    const real central_Pi_sqrd  = SQR( Pi.level_np1[0]  );
    const real central_a_sqrd   = SQR( a.level_np1[0]   );
    const real central_density  = 0.5 * ( central_Phi_sqrd + central_Pi_sqrd ) / central_a_sqrd;
    if( central_density > max_central_density ) {
      max_central_density = central_density;
    }

    /* Check whether or not the lapse function has collapsed */
    if( utilities::check_for_collapse_of_the_lapse( grid, alpha ) ) {
      cout << "\n(SFcollapse1D INFO) The lapse function has collapsed!\n";
      cout <<   "(SFcollapse1D INFO) Iteration = " << n << " | t = " << setprecision(4) << grid.t << endl;
      cout <<   "(SFcollapse1D INFO) Stopping time integration...\n";
      lapse_collapsed = true;
      break;
    }

    /* Print information to the user */
    if( n%OUTPUT_CHECKPOINT == 0 ) {
      // phi.output_to_file(grid,"scalarfield",1,n);
      // Phi.output_to_file(grid,"Phi",1,n);
      // Pi.output_to_file(grid,"Pi",1,n);
      // a.output_to_file(grid,"a",1,n);
      // alpha.output_to_file(grid,"alpha",1,n);
      // utilities::compute_and_output_mass_aspect_function(1,n,grid,a);
      // utilities::output_energy_density_to_file( grid, Phi.level_np1, Pi.level_np1, a.level_np1, n );
    }

    // Output central values
    if( grid.t < 4.0 ) {
      if( n%500 == 0 ) {
    	utilities::output_gridfunctions_central_values( n, grid, phi.level_np1, Phi.level_np1, Pi.level_np1, a.level_np1, alpha.level_np1 );
      }
    }
    else if( grid.t > 4.0 && grid.t < 5.0 ) {
      if( n%200 == 0 ) {
    	utilities::output_gridfunctions_central_values( n, grid, phi.level_np1, Phi.level_np1, Pi.level_np1, a.level_np1, alpha.level_np1 );
      }
    }
    else {
      utilities::output_gridfunctions_central_values( n, grid, phi.level_np1, Phi.level_np1, Pi.level_np1, a.level_np1, alpha.level_np1 );
    }
    
    /* Shift time levels appropriately */
    phi.shift_timelevels(2);
    Phi.shift_timelevels(2);
    Pi.shift_timelevels(2);
    a.shift_timelevels(2);
    alpha.shift_timelevels(2);

    /* Print integration information to the user */
    INTEGRATION_INFO;

    /* Update the iteration number */
    n++;
    
  }

  /* Output max central density */
  // Gaussian shell
  // const real eta_weak   = 0.3364266156435;
  // const real eta_strong = 0.3364266156436;
  // Gaussian shell v2
  // const real eta_weak   = 1.591203885327e-2;
  // const real eta_strong = 1.591203885328e-2;
  // Tanh shell
  // const real eta_weak   = 2.914370451806e-1;
  // const real eta_strong = 2.914370451807e-1;
  // // Set eta_{*}
  // const real eta_c      = 0.5*( eta_weak + eta_strong );
  // ofstream out_central_density;
  // out_central_density.open("max_central_density_values_tanh.dat",ios_base::app);
  // out_central_density << scientific << setprecision(15)
  // 		      << grid.phi0                << " "
  // 		      << max_central_density      << " "
  // 		      << eta_c - grid.phi0        << " "
  // 		      << log(eta_c - grid.phi0)   << " "
  // 		      << log(max_central_density) << endl;
  // out_central_density.close();

  if( lapse_collapsed ) {
    cout << "(SFcollapse1D INFO) Run result: Gravitational collapse! Strong field regime.\n";
    //   cout << setprecision(0) << "1\n";
  }
  else {
    cout << "\n(SFcollapse1D INFO) Run result: Full dispersion! Weak field regime.\n";
    //   cout << setprecision(0) << "0\n";
  }

  cout << "\n(SFcollapse1D INFO) Program terminated without errors!\n";

  return 0;

}
