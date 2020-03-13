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

set term epslatex size 4,5
set output "origin_problem.tex"

set multiplot layout 3,1 margins 0.1,0.98,0.1,0.98 spacing 0,0

set xtics font ",6pt"
set ytics font ",6pt"

# Plot number 1:
# .---.
# | x |
# .---.
# |   |
# .---.
# |   |
# .---.
# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-0.2:1.2]
unset xlabel
set xtics (0,0.2,0.4,0.6,0.8,1.0)
set xtics format ''
# Configure the y axis
ymin = 0.0
ymax = 1.2e-8
increment = (ymax - ymin)/4.0
set yr [ymin-increment:ymax+increment]
set ylabel "$\\Phi_{\\rm weak}(0,r)$"
set ytics ymin,increment,ymax
set ytics format '%.2e'
# Generate the plot
filename = "Phi_weak.dat"
plot filename using 2:3 w l lw 2 lc rgb "blue"

# Plot number 2:
# .---.
# |   |
# .---.
# | x |
# .---.
# |   |
# .---.
# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-0.2:1.2]
unset xlabel
set xtics (0,0.2,0.4,0.6,0.8,1.0)
set xtics format ''
# Configure the y axis
ymin = 0.0
ymax = 6e-8
increment = (ymax - ymin)/4.0
set yr [ymin-increment:ymax+increment]
set ylabel "$\\Phi_{\\rm inter}(0,r)$"
set ytics ymin,increment,ymax
set ytics format '%.2e'
# Generate the plot
filename = "Phi_inter.dat"
plot filename using 2:3 w l lw 2 lc rgb "blue"

# Plot number 3:
# .---.
# |   |
# .---.
# |   |
# .---.
# | x |
# .---.
# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-0.2:1.2]
set xlabel "$r$"
set xtics (0,0.2,0.4,0.6,0.8,1.0)
set xtics format '%.1f'
# Configure the y axis
ymin = 0.0
ymax = 1.2e-7
increment = (ymax - ymin)/4.0
set yr [ymin-increment:ymax+increment]
set ylabel "$\\Phi_{\\rm strong}(0,r)$"
set ytics ymin,increment,ymax
set ytics format '%.2e'
# Generate the plot
filename = "Phi_strong.dat"
plot filename using 2:3 w l lw 2 lc rgb "blue"

unset multiplot
