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

#ifndef __BSSN_VARS__
#define __BSSN_VARS__

/* We use the following BSSN variables in the evolution:
 *
 * - alpha, the lapse function
 * - cf, the conformal factor
 * - vet^{i}, the rescaled version of \beta^{i}, the shift vector
 * - bet^{i}, the rescaled version of the B^{i} vector
 * - lambda^{i}, the rescaled version of \bar\Lambda^{i}
 * - h_{ij}, the rescaled version of eps_{ij} = \bar\gamma_{ij} - \hat\gamma_{ij}
 * - a_{ij}, the rescaled version of A_{ij} = e^{-4phi}( K_{ij} - (1/3)*K*gamma_{ij} )
 */

/* Here we define the BSSN variable macros */
#define IDX_ALPHA    (0)
#define IDX_CF       (1)
#define IDX_TRK      (2)
#define IDX_BETU0    (3)
#define IDX_BETU1    (4)
#define IDX_BETU2    (5)
#define IDX_VETU0    (6)
#define IDX_VETU1    (7)
#define IDX_VETU2    (8)
#define IDX_LAMBDAU0 (9)
#define IDX_LAMBDAU1 (10)
#define IDX_LAMBDAU2 (11)
#define IDX_HDD00    (12)
#define IDX_HDD01    (13)
#define IDX_HDD02    (14)
#define IDX_HDD11    (15)
#define IDX_HDD12    (16)
#define IDX_HDD22    (17)
#define IDX_ADD00    (18)
#define IDX_ADD01    (19)
#define IDX_ADD02    (20)
#define IDX_ADD11    (21)
#define IDX_ADD12    (22)
#define IDX_ADD22    (23)

#endif // __BSSN_VARS__
