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
  ymin  = -0.5
  ymax  = +0.5
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
  ymax = 1.0
  label = "alpha(t,r)"
}
if( which_var eq 'mass' ) {
  ymin = 0.0
  ymax = 0.5
  label = "Mass(t,r)"
}
  if( which_var eq 'rho' ) {
  ymin = 0.0
  ymax = 0.25
  label = "rho(t,r)"
}

set xrange [0.0:5.0]
set yrange [ymin:ymax]
set grid
unset key
set xlabel "r"
set ylabel label

last_file_number  = 12900
output_multiplier = 100
number_of_files   = last_file_number/output_multiplier
dt                = 8.484e-04
do for[i=0:number_of_files] {
  n = i*output_multiplier
  filename = sprintf("../out/%s_%08d.dat",which_var,n)
  time = n*dt
  set title sprintf("Time = %.2f",time)
  plot filename using 2:3 w l lw 2 lc rgb 'blue'
}
