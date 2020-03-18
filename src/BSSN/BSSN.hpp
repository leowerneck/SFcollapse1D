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

#ifndef __BSSN__
#define __BSSN__

#include "macros.hpp"
#include "gridfunction.hpp"

namespace BSSN {

  class tensors {
  public:
    /* alpha and beta^{i} */
    REAL alpha, betaU[3];
    
    /* \bar\gamma_{ij} , det(gamma_{ij}), and \bar\gamma^{ij} */
    REAL gammabarDD[3][3], gammabarUU[3][3], gammaDET;

    /* Derivative of the metric */
    REAL gammabarDD_dD[3][3][3];

    /* Christoffel symbols, Gammabar^{k}_{ij} */
    REAL GammabadUDD[3][3][3];

    /* A_{ij} and A^{ij} */
    REAL AbarDD[3][3], AbarUU[3][3];

    /* K_{ij} and K^{ij} */
    REAL KDD[3][3], KUU[3][3], trK;

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
