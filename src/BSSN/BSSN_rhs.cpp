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

/* We will now implement the RHSs of the BSSN evolution equations */

#include <vector>
#include "macros.hpp"
#include "grid.hpp"

void compute_BSSN_rhs( grid::parameters grid, vector<REAL> in_gfs, vector<REAL> &rhs ) {

  DECLARE_GRID_PARAMETERS;

  /* .-----------------------------.
   * | Step 1: The 1+log condition |
   * .-----------------------------.
   *
   * The evolution equation for the lapse function
   * will be the 1+log condition:
   * .---------------------------------------------------------------.
   * | \partial_{t}\alpha = \alpha\partial_{i}\beta^{i} - 2 \alpha K |
   * .---------------------------------------------------------------.
   */
#pragma omp parallel for
  LOOP3D(Ngz0,Nx0Total-Ngz0, Ngz1,Nx1Total-Ngz1, Ngz2,Nx2Total-Ngz2) {

    /* \bar\gamma_{00} auxiliary variables */
    const REAL gammabarDD00_i0m2  = in_gfs[INDEX4D(IDX_HDD00,i0-2,i1  ,i2  )] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00_i0m1  = in_gfs[INDEX4D(IDX_HDD00,i0-1,i1  ,i2  )] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00       = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1  ,i2  )] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00_i0p1  = in_gfs[INDEX4D(IDX_HDD00,i0+1,i1  ,i2  )] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00_i0p2  = in_gfs[INDEX4D(IDX_HDD00,i0+2,i1  ,i2  )] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00_i0p3  = in_gfs[INDEX4D(IDX_HDD00,i0+3,i1  ,i2  )] * reference_metric::scaleDD[0][0];

    const REAL gammabarDD00_i1m2  = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1-2,i2  )] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00_i1m1  = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1-1,i2  )] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00       = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1  ,i2  )] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00_i1p1  = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1+1,i2  )] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00_i1p2  = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1+2,i2  )] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00_i1p3  = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1+3,i2  )] * reference_metric::scaleDD[0][0];

    const REAL gammabarDD00_i2m2  = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1  ,i2-2)] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00_i2m1  = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1  ,i2-1)] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00       = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1  ,i2  )] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00_i2p1  = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1  ,i2+1)] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00_i2p2  = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1  ,i2+2)] * reference_metric::scaleDD[0][0];
    const REAL gammabarDD00_i2p3  = in_gfs[INDEX4D(IDX_HDD00,i0  ,i1  ,i2+3)] * reference_metric::scaleDD[0][0];

    /* \bar\gamma_{01} auxiliary variables */
    const REAL gammabarDD01_i0m2  = in_gfs[INDEX4D(IDX_HDD01,i0-2,i1  ,i2  )] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01_i0m1  = in_gfs[INDEX4D(IDX_HDD01,i0-1,i1  ,i2  )] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01       = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1  ,i2  )] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01_i0p1  = in_gfs[INDEX4D(IDX_HDD01,i0+1,i1  ,i2  )] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01_i0p2  = in_gfs[INDEX4D(IDX_HDD01,i0+2,i1  ,i2  )] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01_i0p3  = in_gfs[INDEX4D(IDX_HDD01,i0+3,i1  ,i2  )] * reference_metric::scaleDD[0][1];

    const REAL gammabarDD01_i1m2  = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1-2,i2  )] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01_i1m1  = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1-1,i2  )] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01       = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1  ,i2  )] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01_i1p1  = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1+1,i2  )] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01_i1p2  = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1+2,i2  )] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01_i1p3  = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1+3,i2  )] * reference_metric::scaleDD[0][1];

    const REAL gammabarDD01_i2m2  = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1  ,i2-2)] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01_i2m1  = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1  ,i2-1)] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01       = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1  ,i2  )] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01_i2p1  = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1  ,i2+1)] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01_i2p2  = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1  ,i2+2)] * reference_metric::scaleDD[0][1];
    const REAL gammabarDD01_i2p3  = in_gfs[INDEX4D(IDX_HDD01,i0  ,i1  ,i2+3)] * reference_metric::scaleDD[0][1];

    /* \bar\gamma_{02} auxiliary variables */
    const REAL gammabarDD02_i0m2  = in_gfs[INDEX4D(IDX_HDD02,i0-2,i1  ,i2  )] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02_i0m1  = in_gfs[INDEX4D(IDX_HDD02,i0-1,i1  ,i2  )] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02       = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1  ,i2  )] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02_i0p1  = in_gfs[INDEX4D(IDX_HDD02,i0+1,i1  ,i2  )] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02_i0p2  = in_gfs[INDEX4D(IDX_HDD02,i0+2,i1  ,i2  )] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02_i0p3  = in_gfs[INDEX4D(IDX_HDD02,i0+3,i1  ,i2  )] * reference_metric::scaleDD[0][2];

    const REAL gammabarDD02_i1m2  = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1-2,i2  )] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02_i1m1  = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1-1,i2  )] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02       = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1  ,i2  )] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02_i1p1  = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1+1,i2  )] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02_i1p2  = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1+2,i2  )] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02_i1p3  = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1+3,i2  )] * reference_metric::scaleDD[0][2];

    const REAL gammabarDD02_i2m2  = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1  ,i2-2)] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02_i2m1  = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1  ,i2-1)] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02       = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1  ,i2  )] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02_i2p1  = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1  ,i2+1)] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02_i2p2  = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1  ,i2+2)] * reference_metric::scaleDD[0][2];
    const REAL gammabarDD02_i2p3  = in_gfs[INDEX4D(IDX_HDD02,i0  ,i1  ,i2+3)] * reference_metric::scaleDD[0][2];

    /* \bar\gamma_{11} auxiliary variables */
    const REAL gammabarDD11_i0m2  = in_gfs[INDEX4D(IDX_HDD11,i0-2,i1  ,i2  )] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11_i0m1  = in_gfs[INDEX4D(IDX_HDD11,i0-1,i1  ,i2  )] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11       = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1  ,i2  )] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11_i0p1  = in_gfs[INDEX4D(IDX_HDD11,i0+1,i1  ,i2  )] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11_i0p2  = in_gfs[INDEX4D(IDX_HDD11,i0+2,i1  ,i2  )] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11_i0p3  = in_gfs[INDEX4D(IDX_HDD11,i0+3,i1  ,i2  )] * reference_metric::scaleDD[1][1];

    const REAL gammabarDD11_i1m2  = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1-2,i2  )] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11_i1m1  = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1-1,i2  )] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11       = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1  ,i2  )] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11_i1p1  = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1+1,i2  )] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11_i1p2  = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1+2,i2  )] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11_i1p3  = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1+3,i2  )] * reference_metric::scaleDD[1][1];

    const REAL gammabarDD11_i2m2  = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1  ,i2-2)] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11_i2m1  = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1  ,i2-1)] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11       = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1  ,i2  )] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11_i2p1  = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1  ,i2+1)] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11_i2p2  = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1  ,i2+2)] * reference_metric::scaleDD[1][1];
    const REAL gammabarDD11_i2p3  = in_gfs[INDEX4D(IDX_HDD11,i0  ,i1  ,i2+3)] * reference_metric::scaleDD[1][1];

    /* \bar\gamma_{12} auxiliary variables */
    const REAL gammabarDD12_i0m2  = in_gfs[INDEX4D(IDX_HDD12,i0-2,i1  ,i2  )] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12_i0m1  = in_gfs[INDEX4D(IDX_HDD12,i0-1,i1  ,i2  )] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12       = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1  ,i2  )] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12_i0p1  = in_gfs[INDEX4D(IDX_HDD12,i0+1,i1  ,i2  )] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12_i0p2  = in_gfs[INDEX4D(IDX_HDD12,i0+2,i1  ,i2  )] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12_i0p3  = in_gfs[INDEX4D(IDX_HDD12,i0+3,i1  ,i2  )] * reference_metric::scaleDD[1][2];

    const REAL gammabarDD12_i1m2  = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1-2,i2  )] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12_i1m1  = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1-1,i2  )] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12       = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1  ,i2  )] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12_i1p1  = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1+1,i2  )] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12_i1p2  = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1+2,i2  )] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12_i1p3  = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1+3,i2  )] * reference_metric::scaleDD[1][2];

    const REAL gammabarDD12_i2m2  = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1  ,i2-2)] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12_i2m1  = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1  ,i2-1)] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12       = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1  ,i2  )] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12_i2p1  = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1  ,i2+1)] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12_i2p2  = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1  ,i2+2)] * reference_metric::scaleDD[1][2];
    const REAL gammabarDD12_i2p3  = in_gfs[INDEX4D(IDX_HDD12,i0  ,i1  ,i2+3)] * reference_metric::scaleDD[1][2];

    /* \bar\gamma_{22} auxiliary variables */
    const REAL gammabarDD22_i0m2  = in_gfs[INDEX4D(IDX_HDD22,i0-2,i1  ,i2  )] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22_i0m1  = in_gfs[INDEX4D(IDX_HDD22,i0-1,i1  ,i2  )] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22       = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1  ,i2  )] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22_i0p1  = in_gfs[INDEX4D(IDX_HDD22,i0+1,i1  ,i2  )] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22_i0p2  = in_gfs[INDEX4D(IDX_HDD22,i0+2,i1  ,i2  )] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22_i0p3  = in_gfs[INDEX4D(IDX_HDD22,i0+3,i1  ,i2  )] * reference_metric::scaleDD[2][2];

    const REAL gammabarDD22_i1m2  = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1-2,i2  )] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22_i1m1  = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1-1,i2  )] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22       = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1  ,i2  )] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22_i1p1  = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1+1,i2  )] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22_i1p2  = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1+2,i2  )] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22_i1p3  = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1+3,i2  )] * reference_metric::scaleDD[2][2];

    const REAL gammabarDD22_i2m2  = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1  ,i2-2)] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22_i2m1  = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1  ,i2-1)] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22       = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1  ,i2  )] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22_i2p1  = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1  ,i2+1)] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22_i2p2  = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1  ,i2+2)] * reference_metric::scaleDD[2][2];
    const REAL gammabarDD22_i2p3  = in_gfs[INDEX4D(IDX_HDD22,i0  ,i1  ,i2+3)] * reference_metric::scaleDD[2][2];
    
    /* alpha auxiliary variables */
    const REAL alpha_i0m2  = in_gfs[INDEX4D(IDX_ALPHA,i0-2,i1  ,i2  )];
    const REAL alpha_i0m1  = in_gfs[INDEX4D(IDX_ALPHA,i0-1,i1  ,i2  )];
    const REAL alpha       = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1  ,i2  )];
    const REAL alpha_i0p1  = in_gfs[INDEX4D(IDX_ALPHA,i0+1,i1  ,i2  )];
    const REAL alpha_i0p2  = in_gfs[INDEX4D(IDX_ALPHA,i0+2,i1  ,i2  )];
    const REAL alpha_i0p3  = in_gfs[INDEX4D(IDX_ALPHA,i0+3,i1  ,i2  )];

    const REAL alpha_i1m2  = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1-2,i2  )];
    const REAL alpha_i1m1  = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1-1,i2  )];
    const REAL alpha       = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1  ,i2  )];
    const REAL alpha_i1p1  = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1+1,i2  )];
    const REAL alpha_i1p2  = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1+2,i2  )];
    const REAL alpha_i1p3  = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1+3,i2  )];

    const REAL alpha_i2m2  = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1  ,i2-2)];
    const REAL alpha_i2m1  = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1  ,i2-1)];
    const REAL alpha       = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1  ,i2  )];
    const REAL alpha_i2p1  = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1  ,i2+1)];
    const REAL alpha_i2p2  = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1  ,i2+2)];
    const REAL alpha_i2p3  = in_gfs[INDEX4D(IDX_ALPHA,i0  ,i1  ,i2+3)];

    /* beta^{0} auxiliary variables */
    const REAL betaU0_i0m1 = in_gfs[INDEX4D(IDX_VETU0,i0-1,i1  ,i2  )] * reference_metric::scaleU[0];
    const REAL betaU0      = in_gfs[INDEX4D(IDX_VETU0,i0  ,i1  ,i2  )] * reference_metric::scaleU[0];
    const REAL betaU0_i0p1 = in_gfs[INDEX4D(IDX_VETU0,i0+1,i1  ,i2  )] * reference_metric::scaleU[0];
    const REAL betaU0_i0p2 = in_gfs[INDEX4D(IDX_VETU0,i0+2,i1  ,i2  )] * reference_metric::scaleU[0];
    const REAL betaU0_i0p3 = in_gfs[INDEX4D(IDX_VETU0,i0+3,i1  ,i2  )] * reference_metric::scaleU[0];

    const REAL betaU0_i1m1 = in_gfs[INDEX4D(IDX_VETU0,i0  ,i1-1,i2  )] * reference_metric::scaleU[0];
    const REAL betaU0      = in_gfs[INDEX4D(IDX_VETU0,i0  ,i1  ,i2  )] * reference_metric::scaleU[0];
    const REAL betaU0_i1p1 = in_gfs[INDEX4D(IDX_VETU0,i0  ,i1+1,i2  )] * reference_metric::scaleU[0];
    const REAL betaU0_i1p2 = in_gfs[INDEX4D(IDX_VETU0,i0  ,i1+2,i2  )] * reference_metric::scaleU[0];
    const REAL betaU0_i1p2 = in_gfs[INDEX4D(IDX_VETU0,i0  ,i1+3,i2  )] * reference_metric::scaleU[0];

    const REAL betaU0_i2m1 = in_gfs[INDEX4D(IDX_VETU0,i0  ,i1  ,i2-1)] * reference_metric::scaleU[0];
    const REAL betaU0      = in_gfs[INDEX4D(IDX_VETU0,i0  ,i1  ,i2  )] * reference_metric::scaleU[0];
    const REAL betaU0_i2p1 = in_gfs[INDEX4D(IDX_VETU0,i0  ,i1  ,i2+1)] * reference_metric::scaleU[0];
    const REAL betaU0_i2p2 = in_gfs[INDEX4D(IDX_VETU0,i0  ,i1  ,i2+2)] * reference_metric::scaleU[0];
    const REAL betaU0_i2p2 = in_gfs[INDEX4D(IDX_VETU0,i0  ,i1  ,i2+3)] * reference_metric::scaleU[0];

    /* beta^{1} auxiliary variables */
    const REAL betaU1_i0m1 = in_gfs[INDEX4D(IDX_VETU1,i0-1,i1  ,i2  )] * reference_metric::scaleU[1];
    const REAL betaU1      = in_gfs[INDEX4D(IDX_VETU1,i0  ,i1  ,i2  )] * reference_metric::scaleU[1];
    const REAL betaU1_i0p1 = in_gfs[INDEX4D(IDX_VETU1,i0+1,i1  ,i2  )] * reference_metric::scaleU[1];
    const REAL betaU1_i0p2 = in_gfs[INDEX4D(IDX_VETU1,i0+2,i1  ,i2  )] * reference_metric::scaleU[1];
    const REAL betaU1_i0p3 = in_gfs[INDEX4D(IDX_VETU1,i0+3,i1  ,i2  )] * reference_metric::scaleU[1];

    const REAL betaU1_i1m1 = in_gfs[INDEX4D(IDX_VETU1,i0  ,i1-1,i2  )] * reference_metric::scaleU[1];
    const REAL betaU1      = in_gfs[INDEX4D(IDX_VETU1,i0  ,i1  ,i2  )] * reference_metric::scaleU[1];
    const REAL betaU1_i1p1 = in_gfs[INDEX4D(IDX_VETU1,i0  ,i1+1,i2  )] * reference_metric::scaleU[1];
    const REAL betaU1_i1p2 = in_gfs[INDEX4D(IDX_VETU1,i0  ,i1+2,i2  )] * reference_metric::scaleU[1];
    const REAL betaU1_i1p3 = in_gfs[INDEX4D(IDX_VETU1,i0  ,i1+3,i2  )] * reference_metric::scaleU[1];

    const REAL betaU1_i2m1 = in_gfs[INDEX4D(IDX_VETU1,i0  ,i1  ,i2-1)] * reference_metric::scaleU[1];
    const REAL betaU1      = in_gfs[INDEX4D(IDX_VETU1,i0  ,i1  ,i2  )] * reference_metric::scaleU[1];
    const REAL betaU1_i2p1 = in_gfs[INDEX4D(IDX_VETU1,i0  ,i1  ,i2+1)] * reference_metric::scaleU[1];
    const REAL betaU1_i2p2 = in_gfs[INDEX4D(IDX_VETU1,i0  ,i1  ,i2+2)] * reference_metric::scaleU[1];
    const REAL betaU1_i2p3 = in_gfs[INDEX4D(IDX_VETU1,i0  ,i1  ,i2+3)] * reference_metric::scaleU[1];

    /* beta^{2} auxiliary variables */
    const REAL betaU2_i0m1 = in_gfs[INDEX4D(IDX_VETU2,i0-1,i1  ,i2  )] * reference_metric::scaleU[2];
    const REAL betaU2      = in_gfs[INDEX4D(IDX_VETU2,i0  ,i1  ,i2  )] * reference_metric::scaleU[2];
    const REAL betaU2_i0p1 = in_gfs[INDEX4D(IDX_VETU2,i0+1,i1  ,i2  )] * reference_metric::scaleU[2];
    const REAL betaU2_i0p2 = in_gfs[INDEX4D(IDX_VETU2,i0+2,i1  ,i2  )] * reference_metric::scaleU[2];
    const REAL betaU2_i0p3 = in_gfs[INDEX4D(IDX_VETU2,i0+3,i1  ,i2  )] * reference_metric::scaleU[2];
    
    const REAL betaU2_i1m1 = in_gfs[INDEX4D(IDX_VETU2,i0  ,i1-1,i2  )] * reference_metric::scaleU[2];
    const REAL betaU2      = in_gfs[INDEX4D(IDX_VETU2,i0  ,i1  ,i2  )] * reference_metric::scaleU[2];
    const REAL betaU2_i1p1 = in_gfs[INDEX4D(IDX_VETU2,i0  ,i1+1,i2  )] * reference_metric::scaleU[2];
    const REAL betaU2_i1p2 = in_gfs[INDEX4D(IDX_VETU2,i0  ,i1+2,i2  )] * reference_metric::scaleU[2];
    const REAL betaU2_i1p3 = in_gfs[INDEX4D(IDX_VETU2,i0  ,i1+3,i2  )] * reference_metric::scaleU[2];

    const REAL betaU2_i2m1 = in_gfs[INDEX4D(IDX_VETU2,i0  ,i1  ,i2-1)] * reference_metric::scaleU[2];
    const REAL betaU2      = in_gfs[INDEX4D(IDX_VETU2,i0  ,i1  ,i2  )] * reference_metric::scaleU[2];
    const REAL betaU2_i2p1 = in_gfs[INDEX4D(IDX_VETU2,i0  ,i1  ,i2+1)] * reference_metric::scaleU[2];
    const REAL betaU2_i2p2 = in_gfs[INDEX4D(IDX_VETU2,i0  ,i1  ,i2+2)] * reference_metric::scaleU[2];
    const REAL betaU2_i2p3 = in_gfs[INDEX4D(IDX_VETU2,i0  ,i1  ,i2+3)] * reference_metric::scaleU[2];

    /* B^{0} auxiliary variables */
    const REAL BU0_i0m1    = in_gfs[INDEX4D(IDX_BETU0,i0-1,i1  ,i2  )] * reference_metric::scaleU[0];
    const REAL BU0         = in_gfs[INDEX4D(IDX_BETU0,i0  ,i1  ,i2  )] * reference_metric::scaleU[0];
    const REAL BU0_i0p1    = in_gfs[INDEX4D(IDX_BETU0,i0+1,i1  ,i2  )] * reference_metric::scaleU[0];
    const REAL BU0_i0p2    = in_gfs[INDEX4D(IDX_BETU0,i0+2,i1  ,i2  )] * reference_metric::scaleU[0];
    const REAL BU0_i0p3    = in_gfs[INDEX4D(IDX_BETU0,i0+3,i1  ,i2  )] * reference_metric::scaleU[0];
   
    const REAL BU0_i1m1    = in_gfs[INDEX4D(IDX_BETU0,i0  ,i1-1,i2  )] * reference_metric::scaleU[0];
    const REAL BU0         = in_gfs[INDEX4D(IDX_BETU0,i0  ,i1  ,i2  )] * reference_metric::scaleU[0];
    const REAL BU0_i1p1    = in_gfs[INDEX4D(IDX_BETU0,i0  ,i1+1,i2  )] * reference_metric::scaleU[0];
    const REAL BU0_i1p2    = in_gfs[INDEX4D(IDX_BETU0,i0  ,i1+2,i2  )] * reference_metric::scaleU[0];
    const REAL BU0_i1p2    = in_gfs[INDEX4D(IDX_BETU0,i0  ,i1+3,i2  )] * reference_metric::scaleU[0];

    const REAL BU0_i2m1    = in_gfs[INDEX4D(IDX_BETU0,i0  ,i1  ,i2-1)] * reference_metric::scaleU[0];
    const REAL BU0         = in_gfs[INDEX4D(IDX_BETU0,i0  ,i1  ,i2  )] * reference_metric::scaleU[0];
    const REAL BU0_i2p1    = in_gfs[INDEX4D(IDX_BETU0,i0  ,i1  ,i2+1)] * reference_metric::scaleU[0];
    const REAL BU0_i2p2    = in_gfs[INDEX4D(IDX_BETU0,i0  ,i1  ,i2+2)] * reference_metric::scaleU[0];
    const REAL BU0_i2p2    = in_gfs[INDEX4D(IDX_BETU0,i0  ,i1  ,i2+3)] * reference_metric::scaleU[0];

    /* B^{1} auxiliary variables */
    const REAL BU1_i0m1    = in_gfs[INDEX4D(IDX_BETU1,i0-1,i1  ,i2  )] * reference_metric::scaleU[1];
    const REAL BU1         = in_gfs[INDEX4D(IDX_BETU1,i0  ,i1  ,i2  )] * reference_metric::scaleU[1];
    const REAL BU1_i0p1    = in_gfs[INDEX4D(IDX_BETU1,i0+1,i1  ,i2  )] * reference_metric::scaleU[1];
    const REAL BU1_i0p2    = in_gfs[INDEX4D(IDX_BETU1,i0+2,i1  ,i2  )] * reference_metric::scaleU[1];
    const REAL BU1_i0p3    = in_gfs[INDEX4D(IDX_BETU1,i0+3,i1  ,i2  )] * reference_metric::scaleU[1];

    const REAL BU1_i1m1    = in_gfs[INDEX4D(IDX_BETU1,i0  ,i1-1,i2  )] * reference_metric::scaleU[1];
    const REAL BU1         = in_gfs[INDEX4D(IDX_BETU1,i0  ,i1  ,i2  )] * reference_metric::scaleU[1];
    const REAL BU1_i1p1    = in_gfs[INDEX4D(IDX_BETU1,i0  ,i1+1,i2  )] * reference_metric::scaleU[1];
    const REAL BU1_i1p2    = in_gfs[INDEX4D(IDX_BETU1,i0  ,i1+2,i2  )] * reference_metric::scaleU[1];
    const REAL BU1_i1p3    = in_gfs[INDEX4D(IDX_BETU1,i0  ,i1+3,i2  )] * reference_metric::scaleU[1];

    const REAL BU1_i2m1    = in_gfs[INDEX4D(IDX_BETU1,i0  ,i1  ,i2-1)] * reference_metric::scaleU[1];
    const REAL BU1         = in_gfs[INDEX4D(IDX_BETU1,i0  ,i1  ,i2  )] * reference_metric::scaleU[1];
    const REAL BU1_i2p1    = in_gfs[INDEX4D(IDX_BETU1,i0  ,i1  ,i2+1)] * reference_metric::scaleU[1];
    const REAL BU1_i2p2    = in_gfs[INDEX4D(IDX_BETU1,i0  ,i1  ,i2+2)] * reference_metric::scaleU[1];
    const REAL BU1_i2p3    = in_gfs[INDEX4D(IDX_BETU1,i0  ,i1  ,i2+3)] * reference_metric::scaleU[1];

    /* B^{2} auxiliary variables */
    const REAL BU2_i0m1    = in_gfs[INDEX4D(IDX_BETU2,i0-1,i1  ,i2  )] * reference_metric::scaleU[2];
    const REAL BU2         = in_gfs[INDEX4D(IDX_BETU2,i0  ,i1  ,i2  )] * reference_metric::scaleU[2];
    const REAL BU2_i0p1    = in_gfs[INDEX4D(IDX_BETU2,i0+1,i1  ,i2  )] * reference_metric::scaleU[2];
    const REAL BU2_i0p2    = in_gfs[INDEX4D(IDX_BETU2,i0+2,i1  ,i2  )] * reference_metric::scaleU[2];
    const REAL BU2_i0p3    = in_gfs[INDEX4D(IDX_BETU2,i0+3,i1  ,i2  )] * reference_metric::scaleU[2];
    
    const REAL BU2_i1m1    = in_gfs[INDEX4D(IDX_BETU2,i0  ,i1-1,i2  )] * reference_metric::scaleU[2];
    const REAL BU2         = in_gfs[INDEX4D(IDX_BETU2,i0  ,i1  ,i2  )] * reference_metric::scaleU[2];
    const REAL BU2_i1p1    = in_gfs[INDEX4D(IDX_BETU2,i0  ,i1+1,i2  )] * reference_metric::scaleU[2];
    const REAL BU2_i1p2    = in_gfs[INDEX4D(IDX_BETU2,i0  ,i1+2,i2  )] * reference_metric::scaleU[2];
    const REAL BU2_i1p3    = in_gfs[INDEX4D(IDX_BETU2,i0  ,i1+3,i2  )] * reference_metric::scaleU[2];

    const REAL BU2_i2m1    = in_gfs[INDEX4D(IDX_BETU2,i0  ,i1  ,i2-1)] * reference_metric::scaleU[2];
    const REAL BU2         = in_gfs[INDEX4D(IDX_BETU2,i0  ,i1  ,i2  )] * reference_metric::scaleU[2];
    const REAL BU2_i2p1    = in_gfs[INDEX4D(IDX_BETU2,i0  ,i1  ,i2+1)] * reference_metric::scaleU[2];
    const REAL BU2_i2p2    = in_gfs[INDEX4D(IDX_BETU2,i0  ,i1  ,i2+2)] * reference_metric::scaleU[2];
    const REAL BU2_i2p3    = in_gfs[INDEX4D(IDX_BETU2,i0  ,i1  ,i2+3)] * reference_metric::scaleU[2];

    /* TrK auxiliary variables */
    const REAL trK         = in_gfs[INDEX4D(IDX_TRK  ,i0  ,i1  ,i2  )];

    /* Compute gammabarUU */
    REAL gammabarDD[3][3];
    gammabarDD[0][0] = gammabarDD00;
    gammabarDD[0][1] = gammabarDD01;
    gammabarDD[0][2] = gammabarDD02;

    gammabarDD[1][0] = gammabarDD01;
    gammabarDD[1][1] = gammabarDD11;
    gammabarDD[1][2] = gammabarDD12;

    gammabarDD[2][0] = gammabarDD02;
    gammabarDD[2][1] = gammabarDD12;
    gammabarDD[2][2] = gammabarDD22;
    
    REAL gammabarUU[3][3];
    utilities::invert_3x3_symmetric_matrix( gammabarDD, gammabarUU );

    /* Compute the derivatives of \bar\gamma_{ij} */
    REAL gammaDD_dD[3][3][3];
    gammaDD_dD[0][0][0] = inv_dx0 * ( (1.0/12.0)*gammabarDD00_i0m2 - (3.0/4.0)*gammabarDD00_i0m1 + (3.0/4.0)*gammabarDD00_i0p1 - (1.0/12.0)*gammabarDD00_i0p2 );
    gammaDD_dD[0][0][1] = inv_dx1 * ( (1.0/12.0)*gammabarDD00_i1m2 - (3.0/4.0)*gammabarDD00_i1m1 + (3.0/4.0)*gammabarDD00_i1p1 - (1.0/12.0)*gammabarDD00_i1p2 );
    gammaDD_dD[0][0][2] = inv_dx2 * ( (1.0/12.0)*gammabarDD00_i2m2 - (3.0/4.0)*gammabarDD00_i2m1 + (3.0/4.0)*gammabarDD00_i2p1 - (1.0/12.0)*gammabarDD00_i2p2 );

    gammaDD_dD[0][1][0] = inv_dx0 * ( (1.0/12.0)*gammabarDD01_i0m2 - (3.0/4.0)*gammabarDD01_i0m1 + (3.0/4.0)*gammabarDD01_i0p1 - (1.0/12.0)*gammabarDD01_i0p2 );
    gammaDD_dD[0][1][1] = inv_dx1 * ( (1.0/12.0)*gammabarDD01_i1m2 - (3.0/4.0)*gammabarDD01_i1m1 + (3.0/4.0)*gammabarDD01_i1p1 - (1.0/12.0)*gammabarDD01_i1p2 );
    gammaDD_dD[0][1][2] = inv_dx2 * ( (1.0/12.0)*gammabarDD01_i2m2 - (3.0/4.0)*gammabarDD01_i2m1 + (3.0/4.0)*gammabarDD01_i2p1 - (1.0/12.0)*gammabarDD01_i2p2 );

    gammaDD_dD[0][2][0] = inv_dx0 * ( (1.0/12.0)*gammabarDD02_i0m2 - (3.0/4.0)*gammabarDD02_i0m1 + (3.0/4.0)*gammabarDD02_i0p1 - (1.0/12.0)*gammabarDD02_i0p2 );
    gammaDD_dD[0][2][1] = inv_dx1 * ( (1.0/12.0)*gammabarDD02_i1m2 - (3.0/4.0)*gammabarDD02_i1m1 + (3.0/4.0)*gammabarDD02_i1p1 - (1.0/12.0)*gammabarDD02_i1p2 );
    gammaDD_dD[0][2][2] = inv_dx2 * ( (1.0/12.0)*gammabarDD02_i2m2 - (3.0/4.0)*gammabarDD02_i2m1 + (3.0/4.0)*gammabarDD02_i2p1 - (1.0/12.0)*gammabarDD02_i2p2 );

    gammaDD_dD[1][1][0] = inv_dx0 * ( (1.0/12.0)*gammabarDD11_i0m2 - (3.0/4.0)*gammabarDD11_i0m1 + (3.0/4.0)*gammabarDD11_i0p1 - (1.0/12.0)*gammabarDD11_i0p2 );
    gammaDD_dD[1][1][1] = inv_dx1 * ( (1.0/12.0)*gammabarDD11_i1m2 - (3.0/4.0)*gammabarDD11_i1m1 + (3.0/4.0)*gammabarDD11_i1p1 - (1.0/12.0)*gammabarDD11_i1p2 );
    gammaDD_dD[1][1][2] = inv_dx2 * ( (1.0/12.0)*gammabarDD11_i2m2 - (3.0/4.0)*gammabarDD11_i2m1 + (3.0/4.0)*gammabarDD11_i2p1 - (1.0/12.0)*gammabarDD11_i2p2 );

    gammaDD_dD[1][2][0] = inv_dx0 * ( (1.0/12.0)*gammabarDD12_i0m2 - (3.0/4.0)*gammabarDD12_i0m1 + (3.0/4.0)*gammabarDD12_i0p1 - (1.0/12.0)*gammabarDD12_i0p2 );
    gammaDD_dD[1][2][1] = inv_dx1 * ( (1.0/12.0)*gammabarDD12_i1m2 - (3.0/4.0)*gammabarDD12_i1m1 + (3.0/4.0)*gammabarDD12_i1p1 - (1.0/12.0)*gammabarDD12_i1p2 );
    gammaDD_dD[1][2][2] = inv_dx2 * ( (1.0/12.0)*gammabarDD12_i2m2 - (3.0/4.0)*gammabarDD12_i2m1 + (3.0/4.0)*gammabarDD12_i2p1 - (1.0/12.0)*gammabarDD12_i2p2 );

    gammaDD_dD[2][2][0] = inv_dx0 * ( (1.0/12.0)*gammabarDD22_i0m2 - (3.0/4.0)*gammabarDD22_i0m1 + (3.0/4.0)*gammabarDD22_i0p1 - (1.0/12.0)*gammabarDD22_i0p2 );
    gammaDD_dD[2][2][1] = inv_dx1 * ( (1.0/12.0)*gammabarDD22_i1m2 - (3.0/4.0)*gammabarDD22_i1m1 + (3.0/4.0)*gammabarDD22_i1p1 - (1.0/12.0)*gammabarDD22_i1p2 );
    gammaDD_dD[2][2][2] = inv_dx2 * ( (1.0/12.0)*gammabarDD22_i2m2 - (3.0/4.0)*gammabarDD22_i2m1 + (3.0/4.0)*gammabarDD22_i2p1 - (1.0/12.0)*gammabarDD22_i2p2 );

    gammaDD_dD[1][0][0] = gammaDD_dD[0][1][0];
    gammaDD_dD[1][0][1] = gammaDD_dD[0][1][1];
    gammaDD_dD[1][0][2] = gammaDD_dD[0][1][2];

    gammaDD_dD[2][0][0] = gammaDD_dD[0][2][0];
    gammaDD_dD[2][0][1] = gammaDD_dD[0][2][1];
    gammaDD_dD[2][0][2] = gammaDD_dD[0][2][2];

    gammaDD_dD[2][1][0] = gammaDD_dD[1][2][0];
    gammaDD_dD[2][1][1] = gammaDD_dD[1][2][1];
    gammaDD_dD[2][1][2] = gammaDD_dD[1][2][2];

    /* Then compute GammabarUDD */
    REAL GammabarUDD[3][3][3];
    for(int k=0;k<3;k++)
      for(int j=0;j<3;j++)
	for(int i=0;i<3;i++)
	  GammabarUDD[k][i][j] = 0;

    for(int k=0;k<3;k++)
      for(int j=0;j<3;j++)
	for(int i=0;i<3;i++)
	  for(int l=0;l<3;l++)
	    GammabarUDD[k][i][j] += 0.5*gammabarUU[k][l] * (gammabarDD_dD[j][l][i] + gammabarDD_dD[i][l][j] - gammabarDD_dD[i][j][l]);
    

    /* Compute the upwinded derivatives of alpha */
    const REAL alpha_dupD0 = inv_dx0 * ( -(1.0/4.0)*alpha_i0m1 - (5.0/6.0)*alpha + (3.0/2.0)*alpha_i0p1 - (1.0/2.0)*alpha_i0p2 + (1.0/12.0)*alpha_i0p3 );
    const REAL alpha_dupD1 = inv_dx1 * ( -(1.0/4.0)*alpha_i1m1 - (5.0/6.0)*alpha + (3.0/2.0)*alpha_i1p1 - (1.0/2.0)*alpha_i1p2 + (1.0/12.0)*alpha_i1p3 );
    const REAL alpha_dupD2 = inv_dx2 * ( -(1.0/4.0)*alpha_i2m1 - (5.0/6.0)*alpha + (3.0/2.0)*alpha_i2p1 - (1.0/2.0)*alpha_i2p2 + (1.0/12.0)*alpha_i2p3 );

    /* Compute the upwinded derivatives of beta^{i} */
    const REAL betaU0_dupD0 = inv_dx0 * ( -(1.0/4.0)*betaU0_i0m1 - (5.0/6.0)*betaU0 + (3.0/2.0)*betaU0_i0p1 - (1.0/2.0)*betaU0_i0p2 + (1.0/12.0)*betaU0_i0p3 );
    const REAL betaU0_dupD1 = inv_dx1 * ( -(1.0/4.0)*betaU0_i1m1 - (5.0/6.0)*betaU0 + (3.0/2.0)*betaU0_i1p1 - (1.0/2.0)*betaU0_i1p2 + (1.0/12.0)*betaU0_i1p3 );
    const REAL betaU0_dupD2 = inv_dx2 * ( -(1.0/4.0)*betaU0_i2m1 - (5.0/6.0)*betaU0 + (3.0/2.0)*betaU0_i2p1 - (1.0/2.0)*betaU0_i2p2 + (1.0/12.0)*betaU0_i2p3 );
    
    const REAL betaU1_dupD0 = inv_dx0 * ( -(1.0/4.0)*betaU1_i0m1 - (5.0/6.0)*betaU1 + (3.0/2.0)*betaU1_i0p1 - (1.0/2.0)*betaU1_i0p2 + (1.0/12.0)*betaU1_i0p3 );
    const REAL betaU1_dupD1 = inv_dx1 * ( -(1.0/4.0)*betaU1_i1m1 - (5.0/6.0)*betaU1 + (3.0/2.0)*betaU1_i1p1 - (1.0/2.0)*betaU1_i1p2 + (1.0/12.0)*betaU1_i1p3 );
    const REAL betaU1_dupD2 = inv_dx2 * ( -(1.0/4.0)*betaU1_i2m1 - (5.0/6.0)*betaU1 + (3.0/2.0)*betaU1_i2p1 - (1.0/2.0)*betaU1_i2p2 + (1.0/12.0)*betaU1_i2p3 );
    
    const REAL betaU2_dupD0 = inv_dx0 * ( -(1.0/4.0)*betaU2_i0m1 - (5.0/6.0)*betaU2 + (3.0/2.0)*betaU2_i0p1 - (1.0/2.0)*betaU2_i0p2 + (1.0/12.0)*betaU2_i0p3 );
    const REAL betaU2_dupD1 = inv_dx1 * ( -(1.0/4.0)*betaU2_i1m1 - (5.0/6.0)*betaU2 + (3.0/2.0)*betaU2_i1p1 - (1.0/2.0)*betaU2_i1p2 + (1.0/12.0)*betaU2_i1p3 );
    const REAL betaU2_dupD2 = inv_dx2 * ( -(1.0/4.0)*betaU2_i2m1 - (5.0/6.0)*betaU2 + (3.0/2.0)*betaU2_i2p1 - (1.0/2.0)*betaU2_i2p2 + (1.0/12.0)*betaU2_i2p3 );

    /* Compute the upwinded derivatives of B^{i} */
    const REAL BU0_dupD0    = inv_dx0 * ( -(1.0/4.0)*BU0_i0m1 - (5.0/6.0)*BU0 + (3.0/2.0)*BU0_i0p1 - (1.0/2.0)*BU0_i0p2 + (1.0/12.0)*BU0_i0p3 );
    const REAL BU0_dupD1    = inv_dx1 * ( -(1.0/4.0)*BU0_i1m1 - (5.0/6.0)*BU0 + (3.0/2.0)*BU0_i1p1 - (1.0/2.0)*BU0_i1p2 + (1.0/12.0)*BU0_i1p3 );
    const REAL BU0_dupD2    = inv_dx2 * ( -(1.0/4.0)*BU0_i2m1 - (5.0/6.0)*BU0 + (3.0/2.0)*BU0_i2p1 - (1.0/2.0)*BU0_i2p2 + (1.0/12.0)*BU0_i2p3 );
    
    const REAL BU1_dupD0    = inv_dx0 * ( -(1.0/4.0)*BU1_i0m1 - (5.0/6.0)*BU1 + (3.0/2.0)*BU1_i0p1 - (1.0/2.0)*BU1_i0p2 + (1.0/12.0)*BU1_i0p3 );
    const REAL BU1_dupD1    = inv_dx1 * ( -(1.0/4.0)*BU1_i1m1 - (5.0/6.0)*BU1 + (3.0/2.0)*BU1_i1p1 - (1.0/2.0)*BU1_i1p2 + (1.0/12.0)*BU1_i1p3 );
    const REAL BU1_dupD2    = inv_dx2 * ( -(1.0/4.0)*BU1_i2m1 - (5.0/6.0)*BU1 + (3.0/2.0)*BU1_i2p1 - (1.0/2.0)*BU1_i2p2 + (1.0/12.0)*BU1_i2p3 );
    
    const REAL BU2_dupD0    = inv_dx0 * ( -(1.0/4.0)*BU2_i0m1 - (5.0/6.0)*BU2 + (3.0/2.0)*BU2_i0p1 - (1.0/2.0)*BU2_i0p2 + (1.0/12.0)*BU2_i0p3 );
    const REAL BU2_dupD1    = inv_dx1 * ( -(1.0/4.0)*BU2_i1m1 - (5.0/6.0)*BU2 + (3.0/2.0)*BU2_i1p1 - (1.0/2.0)*BU2_i1p2 + (1.0/12.0)*BU2_i1p3 );
    const REAL BU2_dupD2    = inv_dx2 * ( -(1.0/4.0)*BU2_i2m1 - (5.0/6.0)*BU2 + (3.0/2.0)*BU2_i2p1 - (1.0/2.0)*BU2_i2p2 + (1.0/12.0)*BU2_i2p3 );

    /* Compute the contraction beta^{i}\partial_{i}alpha */
    const REAL contraction_betaU__alpha_dupD = betaU0 * alpha_dupD0 + betaU1 * alpha_dupD1 + betaU2 * alpha_dupD2;

    /* Compute the contractions beta^{j}\partial_{j}\beta^{i} */
    const REAL contraction_betaU__betaU_dupD0 = betaU0 * betaU0_dupD0 + betaU1 * betaU0_dupD1 + betaU2 * betaU0_dupD2;
    const REAL contraction_betaU__betaU_dupD1 = betaU0 * betaU1_dupD0 + betaU1 * betaU1_dupD1 + betaU2 * betaU1_dupD2;
    const REAL contraction_betaU__betaU_dupD2 = betaU0 * betaU2_dupD0 + betaU1 * betaU2_dupD1 + betaU2 * betaU2_dupD2;

    /* Compute the contractions beta^{j}\partial_{j}B^{i} */
    const REAL contraction_betaU__BU_dupD0 = betaU0 * BU0_dupD0 + betaU1 * BU0_dupD1 + betaU2 * BU0_dupD2;
    const REAL contraction_betaU__BU_dupD1 = betaU0 * BU1_dupD0 + betaU1 * BU1_dupD1 + betaU2 * BU1_dupD2;
    const REAL contraction_betaU__BU_dupD2 = betaU0 * BU2_dupD0 + betaU1 * BU2_dupD1 + betaU2 * BU2_dupD2;
    
    /* Compute the RHS of alpha */
    const REAL rhs_alpha = contraction_betaU__alpha_dupD - 2.0 * alpha * trK;

    /* Compute the RHS of beta^{i} */
    const REAL rhs_betaU0 = contraction_betaU__betaU_dupD0 + BU0;
    const REAL rhs_betaU1 = contraction_betaU__betaU_dupD1 + BU1;
    const REAL rhs_betaU2 = contraction_betaU__betaU_dupD2 + BU2;

    /* Compute the RHS of B^{i} */
    const REAL rhs_BU0    = contraction_betaU__BU_dupD0 + (3.0/4.0) * ( rhs_LambdabarU0 - contraction_betaU__LambdabarU_dupD0 ) - BSSN_ETA * B0;
    const REAL rhs_BU1    = contraction_betaU__BU_dupD0 + (3.0/4.0) * ( rhs_LambdabarU1 - contraction_betaU__LambdabarU_dupD1 ) - BSSN_ETA * B1;
    const REAL rhs_BU2    = contraction_betaU__BU_dupD0 + (3.0/4.0) * ( rhs_LambdabarU2 - contraction_betaU__LambdabarU_dupD2 ) - BSSN_ETA * B2;

    /* .-------------------------------------------------.
     * | Step 2: Covariant Gamma-driving shift condition |
     * .-------------------------------------------------.
     *
     * For the shift vector we use Brown's covariant Gamma-driving shift
     * condition, which reads
     * .------------------------------------------------------------------------------.
     * | \partial_{t}\beta^{i} = beta^{j}\partial_{j}beta^{i}                         |
     * |                       + beta^{j}barGamma^{i}_{mj}\alpha\partial_{i}\beta^{i} |
     * |                       + B^{i}                                                |
     * .------------------------------------------------------------------------------.
     */
    /* beta^{0} auxiliary variables */
    

    /* beta^{1} auxiliary variables */
    

    
  }
  
}
