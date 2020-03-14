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

set term epslatex
set output sprintf("%s_%s.tex",which_var,which_regime)

if ( which_var eq 'scalarfield' ) {
  if ( which_regime eq 'weak' ) {
    ymin_tic   = -0.04
    ymax_tic   = +0.04
    tic_format = '%.2f'
  }
  if ( which_regime eq 'inter' ) {
    ymin_tic   = -0.4
    ymax_tic   = +0.2
    tic_format = '%.2f'
    }
  if ( which_regime eq 'strong' ) {
    ymin_tic   = -0.4
    ymax_tic   = +0.4
    tic_format = '%.2f'
  }
  label      = "$\\phi(t,r)$"
}
if( which_var eq 'Phi' ) {
  ymin_tic = -0.01
  ymax_tic = +0.01
  label = "$\\Phi(t,r)$"
}
if( which_var eq 'Pi' ) {
  ymin_tic = -0.02
  ymax_tic = +0.02
  label = "$\\Pi(t,r)$"
}
if( which_var eq 'a' ) {
  ymin_tic = 0.9999
  ymax_tic = 1.005
  label = "$a(t,r)$"
}
if( which_var eq 'alpha' ) {
  if ( which_regime eq 'weak' ) {
    ymin_tic   = 0.980
    ymax_tic   = 1.002
    tic_format = '%.3f'
  }
  if ( which_regime eq 'inter' ) {
    ymin_tic   = 0.3
    ymax_tic   = 1.1
    tic_format = '%.2f'
  }
  if ( which_regime eq 'strong' ) {
    ymin_tic   = 0.0
    ymax_tic   = 1.2
    tic_format = '%.1f'
  }
  label      = "$\\alpha(t,r)$"
}
if( which_var eq 'mass' ) {
  if ( which_regime eq 'inter' ) {
    ymin_tic   = 0.0
    ymax_tic   = 0.07
    tic_format = '%.3f'
  }
  if ( which_regime eq 'inter' ) {
    ymin_tic   = 0.0
    ymax_tic   = 2.0
    tic_format = '%.1f'
  }
  if ( which_regime eq 'strong' ) {
    ymin_tic   = 0.0
    ymax_tic   = 0.5
    tic_format = '%.3f'
  }
  label      = "$M(t,r)$"
}

increment = (ymax_tic - ymin_tic)/4.0
ymin = ymin_tic - increment
ymax = ymax_tic + increment

set multiplot layout 3,3 margins 0.1,0.98,0.1,0.98 spacing 0,0

set xtics font ",6pt"
set ytics font ",6pt"

# Plot number 1:
# .---.---.---.
# | x |   |   |
# .---.---.---.
# |   |   |   |
# .---.---.---.
# |   |   |   |
# .---.---.---.
# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-1:6]
unset xlabel
set xtics (0,1,2,3,4,5)
set xtics format ''
# Configure the y axis
set yrange [ymin:ymax]
set ylabel label
set ytics ymin_tic,increment,ymax_tic
set ytics format tic_format
# Generate the plot
set label 1 '$t = 0$' at 2.5,ymax_tic+increment/3.0 center font ",6pt"
filename = sprintf("../out/%s_000000.dat",which_var)
plot filename using 2:3 w l lw 2 lc rgb "blue"

# Plot number 2:
# .---.---.---.
# |   | x |   |
# .---.---.---.
# |   |   |   |
# .---.---.---.
# |   |   |   |
# .---.---.---.
# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-1:6]
unset xlabel
set xtics (0,1,2,3,4,5)
set xtics format ''
# Configure the y axis
set yrange [ymin:ymax]
unset ylabel
set ytics ymin_tic,increment,ymax_tic
set ytics format ''
# Generate the plot
set label 1 '$t = 1.1\times10^{5}\Delta t$' at 2.5,ymax_tic+increment/3.0 center font ",6pt"
filename = sprintf("../out/%s_011000.dat",which_var)
plot filename using 2:3 w l lw 2 lc rgb "blue"

# Plot number 3:
# .---.---.---.
# |   |   | x |
# .---.---.---.
# |   |   |   |
# .---.---.---.
# |   |   |   |
# .---.---.---.
# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-1:6]
unset xlabel
set xtics (0,1,2,3,4,5)
set xtics format ''
# Configure the y axis
set yrange [ymin:ymax]
unset ylabel
set ytics ymin_tic,increment,ymax_tic
set ytics format ''
# Generate the plot
set label 1 '$t = 2.2\times10^{5}\Delta t$' at 2.5,ymax_tic+increment/3.0 center font ",6pt"
filename = sprintf("../out/%s_022000.dat",which_var)
plot filename using 2:3 w l lw 2 lc rgb "blue"

# Plot number 4:
# .---.---.---.
# |   |   |   |
# .---.---.---.
# | x |   |   |
# .---.---.---.
# |   |   |   |
# .---.---.---.
# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-1:6]
unset xlabel
set xtics (0,1,2,3,4,5)
set xtics format ''
# Configure the y axis
set yrange [ymin:ymax]
set ylabel label
set ytics ymin_tic,increment,ymax_tic
set ytics format tic_format
# Generate the plot
set label 1 '$t = 3.3\times10^{5}\Delta t$' at 2.5,ymax_tic+increment/3.0 center font ",6pt"
filename = sprintf("../out/%s_033000.dat",which_var)
plot filename using 2:3 w l lw 2 lc rgb "blue"

# Plot number 5:
# .---.---.---.
# |   |   |   |
# .---.---.---.
# |   | x |   |
# .---.---.---.
# |   |   |   |
# .---.---.---.
# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-1:6]
unset xlabel
set xtics (0,1,2,3,4,5)
set xtics format ''
# Configure the y axis
set yrange [ymin:ymax]
unset ylabel
set ytics ymin_tic,increment,ymax_tic
set ytics format ''
# Generate the plot
set label 1 '$t = 4.4\times10^{5}\Delta t$' at 2.5,ymax_tic+increment/3.0 center font ",6pt"
filename = sprintf("../out/%s_044000.dat",which_var)
plot filename using 2:3 w l lw 2 lc rgb "blue"

# Plot number 6:
# .---.---.---.
# |   |   |   |
# .---.---.---.
# |   |   | x |
# .---.---.---.
# |   |   |   |
# .---.---.---.
# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-1:6]
unset xlabel
set xtics (0,1,2,3,4,5)
set xtics format ''
# Configure the y axis
set yrange [ymin:ymax]
unset ylabel
set ytics ymin_tic,increment,ymax_tic
set ytics format ''
# Generate the plot
set label 1 '$t = 5.5\times10^{5}\Delta t$' at 2.5,ymax_tic+increment/3.0 center font ",6pt"
filename = sprintf("../out/%s_055000.dat",which_var)
plot filename using 2:3 w l lw 2 lc rgb "blue"

# Plot number 7:
# .---.---.---.
# |   |   |   |
# .---.---.---.
# |   |   |   |
# .---.---.---.
# | x |   |   |
# .---.---.---.
# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-1:6]
set xlabel "$r$"
set xtics (0,1,2,3,4,5)
set xtics format
# Configure the y axis
set yrange [ymin:ymax]
set ylabel label
set ytics ymin_tic,increment,ymax_tic
set ytics format tic_format
# Generate the plot
set label 1 '$t = 6.6\times10^{5}\Delta t$' at 2.5,ymax_tic+increment/3.0 center font ",6pt"
filename = sprintf("../out/%s_066000.dat",which_var)
plot filename using 2:3 w l lw 2 lc rgb "blue"

# Plot number 8:
# .---.---.---.
# |   |   |   |
# .---.---.---.
# |   |   |   |
# .---.---.---.
# |   | x |   |
# .---.---.---.
# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-1:6]
set xlabel "$r$"
set xtics (0,1,2,3,4,5)
set xtics format
# Configure the y axis
set yrange [ymin:ymax]
unset ylabel
set ytics ymin_tic,increment,ymax_tic
set ytics format ''
# Generate the plot
set label 1 '$t = 7.7\times10^{5}\Delta t$' at 2.5,ymax_tic+increment/3.0 center font "0,6pt"
filename = sprintf("../out/%s_077000.dat",which_var)
plot filename using 2:3 w l lw 2 lc rgb "blue"

# Plot number 9:
# .---.---.---.
# |   |   |   |
# .---.---.---.
# |   |   |   |
# .---.---.---.
# |   |   | x |
# .---.---.---.
# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-1:6]
set xlabel "$r$"
set xtics (0,1,2,3,4,5)
set xtics format
# Configure the y axis
set yrange [ymin:ymax]
unset ylabel
set ytics ymin_tic,increment,ymax_tic
set ytics format ''
# Generate the plot
set label 1 '$t = 8.8\times10^{5}\Delta t$' at 2.5,ymax_tic+increment/3.0 center font "0,6pt"
filename = sprintf("../out/%s_088000.dat",which_var)
plot filename using 2:3 w l lw 2 lc rgb "blue"

unset multiplot
