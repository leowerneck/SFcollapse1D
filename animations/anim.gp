# .-----------------------------------------------------------------------.
# | SFcollapse1D                                                          |
# | Gravitational collapse of scalar fields in spherical symmetry         |
# |                                                                       |
# | Copyright (c) 2020, Leonardo Werneck                                  |
# |                                                                       |
# | This program is free software: you can redistribute it and/or modify  |
# | it under the terms of the GNU General Public License as published by  |
# | the Free Software Foundation, either version 3 of the License, or     |
# | (at your option) any later version.                                   |
# |                                                                       |
# | This program is distributed in the hope that it will be useful,       |
# | but WITHOUT ANY WARRANTY; without even the implied warranty of        |
# | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
# | GNU General Public License for more details.                          |
# |                                                                       |
# | You should have received a copy of the GNU General Public License     |
# | along with this program.  If not, see <https://www.gnu.org/licenses/>.|
# .-----------------------------------------------------------------------.

set term gif enhanced animate delay 10 size 900,600
set output sprintf("%s.gif",which_var)
set encoding utf8

if ( which_var eq 'scalarfield' ) {
  ymin  = -0.04
  ymax  = +0.04
  label = "phi(t,r)"
}
if( which_var eq 'Phi' ) {
  ymin = -0.01
  ymax = +0.01
  label = "Phi(t,r)"
}
if( which_var eq 'Pi' ) {
  ymin = -0.02
  ymax = +0.02
  label = "Pi(t,r)"
}
if( which_var eq 'a' ) {
  ymin = 0.9999
  ymax = 1.002
  set ytics (1.0, 1.0005,1.001,1.0015,1.002)
  label = "a(t,r)"
}
if( which_var eq 'alpha' ) {
  ymin = 0.980
  ymax = 1.010
  label = "alpha(t,r)"
}
if( which_var eq 'mass' ) {
  ymin = 0.0
  ymax = 0.006
  label = "Mass(t,r)"
}

set xrange [0:50]
set yrange [ymin:ymax]
set grid
unset key
set xlabel "r"
set ylabel label

do for[i=0:120] {
  filename = sprintf("../out/%s_%04d.dat",which_var,10*i)
  plot filename using 1:3 w l lw 2 lc rgb 'blue'
}
