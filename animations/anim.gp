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
  ymin  = -0.4
  ymax  = +0.4
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
  ymax = 1.005
  label = "a(t,r)"
}
if( which_var eq 'alpha' ) {
  ymin = 0.0
  ymax = 1.2
  label = "alpha(t,r)"
}
if( which_var eq 'mass' ) {
  ymin = 0.0
  ymax = 0.5
  label = "Mass(t,r)"
}

set xrange [0.0:5.0]
set yrange [ymin:ymax]
set grid
unset key
set xlabel "r"
set ylabel label

last_file_number  = 27600/2
output_multiplier = 200
number_of_files   = last_file_number/output_multiplier

do for[i=0:number_of_files] {
  filename = sprintf("../out/cfl_05/%s_%05d.dat",which_var,i*output_multiplier)
  plot filename using 2:3 w l lw 2 lc rgb 'blue'
}
