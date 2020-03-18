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

#include "BSSN.hpp"
#include "macros.hpp"
#include "utilities.hpp"
#include "gridfunctions.hpp"

void BSSN::tensors::compute_BSSN_tensors_from_gridfunctions( const int i0, const int i1, const int i2, const gridfunctions in_gfs ) {

  /* Start by initializing the variables inside the class that are already defined */
  alpha            = in_gfs[INDEX3DGF(IDX_ALPHA   ,i0,i1,i2)];
  cf               = in_gfs[INDEX3DGF(IDX_CF      ,i0,i1,i2)];
  gammabarDD[0][0] = in_gfs[INDEX3DGF(IDX_HDD00   ,i0,i1,i2)] * BSSN::reference_metric::scaleDD[0][0] + BSSN::reference_metric::gammahatDD[0][0];
  gammabarDD[0][1] = in_gfs[INDEX3DGF(IDX_HDD00   ,i0,i1,i2)] * BSSN::reference_metric::scaleDD[0][1] + BSSN::reference_metric::gammahatDD[0][1];
  gammabarDD[0][2] = in_gfs[INDEX3DGF(IDX_HDD00   ,i0,i1,i2)] * BSSN::reference_metric::scaleDD[0][2] + BSSN::reference_metric::gammahatDD[0][2];
  gammabarDD[1][1] = in_gfs[INDEX3DGF(IDX_HDD00   ,i0,i1,i2)] * BSSN::reference_metric::scaleDD[1][1] + BSSN::reference_metric::gammahatDD[1][1];
  gammabarDD[1][2] = in_gfs[INDEX3DGF(IDX_HDD00   ,i0,i1,i2)] * BSSN::reference_metric::scaleDD[1][2] + BSSN::reference_metric::gammahatDD[1][2];
  gammabarDD[2][2] = in_gfs[INDEX3DGF(IDX_HDD00   ,i0,i1,i2)] * BSSN::reference_metric::scaleDD[2][2] + BSSN::reference_metric::gammahatDD[2][2];
  gammabarDD[1][0] = gammabarDD[0][1];
  gammabarDD[2][0] = gammabarDD[0][2];
  gammabarDD[2][1] = gammabarDD[1][2];
  AbarDD[0][0]     = in_gfs[INDEX3DGF(IDX_ADD00   ,i0,i1,i2)] * BSSN::reference_metric::scaleDD[0][0];
  AbarDD[0][1]     = in_gfs[INDEX3DGF(IDX_ADD00   ,i0,i1,i2)] * BSSN::reference_metric::scaleDD[0][1];
  AbarDD[0][2]     = in_gfs[INDEX3DGF(IDX_ADD00   ,i0,i1,i2)] * BSSN::reference_metric::scaleDD[0][2];
  AbarDD[1][1]     = in_gfs[INDEX3DGF(IDX_ADD00   ,i0,i1,i2)] * BSSN::reference_metric::scaleDD[1][1];
  AbarDD[1][2]     = in_gfs[INDEX3DGF(IDX_ADD00   ,i0,i1,i2)] * BSSN::reference_metric::scaleDD[1][2];
  AbarDD[2][2]     = in_gfs[INDEX3DGF(IDX_ADD00   ,i0,i1,i2)] * BSSN::reference_metric::scaleDD[2][2];
  AbarDD[1][0]     = AbarDD[0][1];
  AbarDD[2][0]     = AbarDD[0][2];
  AbarDD[2][1]     = AbarDD[1][2];
  trK              = in_gfs[INDEX3DGF(IDX_TRK     ,i0,i1,i2)];
  betaU[0]         = in_gfs[INDEX3DGF(IDX_VETU0   ,i0,i1,i2)] * BSSN::reference_metric::scaleU[0];
  betaU[1]         = in_gfs[INDEX3DGF(IDX_VETU1   ,i0,i1,i2)] * BSSN::reference_metric::scaleU[1];
  betaU[2]         = in_gfs[INDEX3DGF(IDX_VETU2   ,i0,i1,i2)] * BSSN::reference_metric::scaleU[2];
  BU[0]            = in_gfs[INDEX3DGF(IDX_BETU0   ,i0,i1,i2)] * BSSN::reference_metric::scaleU[0];
  BU[1]            = in_gfs[INDEX3DGF(IDX_BETU1   ,i0,i1,i2)] * BSSN::reference_metric::scaleU[1];
  BU[2]            = in_gfs[INDEX3DGF(IDX_BETU2   ,i0,i1,i2)] * BSSN::reference_metric::scaleU[2];
  LambdabarU[0]    = in_gfs[INDEX3DGF(IDX_LAMBDAU0,i0,i1,i2)] * BSSN::reference_metric::scaleU[0];
  LambdabarU[1]    = in_gfs[INDEX3DGF(IDX_LAMBDAU1,i0,i1,i2)] * BSSN::reference_metric::scaleU[1];
  LambdabarU[2]    = in_gfs[INDEX3DGF(IDX_LAMBDAU2,i0,i1,i2)] * BSSN::reference_metric::scaleU[2];

  /* Now set \bar\gamma^{ij} and \bar\gamma := det( \bar\gamma_{ij} ) */
  utilities::invert_3x3_symmetric_matrix__compute_det( gammabarDD, gammabarDET, gammabarUU );

  
  
}

namespace BSSN {

  class tensors {
  public:
    /* alpha and beta^{i} */
    REAL alpha, betaU[3];

    /* Conformal factor */
    REAL cf;
    
    /* \bar\gamma_{ij} , det(gamma_{ij}), and \bar\gamma^{ij} */
    REAL gammabarDD[3][3], gammabarUU[3][3], gammabarDET;

    /* Derivatives of the metric */
    REAL gammabarDD_dD[3][3][3], gammabarDD_dupD[3][3][3], gammabarDD_dDD[3][3][3][3];

    /* Christoffel symbols, Gammabar^{k}_{ij} */
    REAL GammabadUDD[3][3][3];

    /* A_{ij} and A^{ij} */
    REAL AbarDD[3][3], AbarUU[3][3];

    /* K_{ij} and K^{ij} */
    REAL KDD[3][3], KUU[3][3], trK;

    /* \bar\Lambda^{i} */
    REAL LambdabarU[3];

    /* Constructor */
    tensors( const int i0, const int i1, const int i2, const gridfunctions in_gfs ) {
      compute_BSSN_tensors_from_gridfunctions( i0, i1, i2, in_gfs );
    }

  private:
    /* Initializer function */
    void compute_BSSN_tensors_from_gridfunctions( const int i0, const int i1, const int i2, const gridfunctions in_gfs );
  };

}

#endif // __BSSN__
