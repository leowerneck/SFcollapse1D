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

set term epslatex size 5,3
set output "resources/origin_problem.tex"

set xtics font ",6pt"
set ytics font ",6pt"

# Basic setup
set grid
unset key
# Configure the x axis
set xrange [-0.2:1.2]
set xlabel "$r$"
set xtics (0,0.2,0.4,0.6,0.8,1.0)
set xtics format "%.1f"
# Configure the y axis
ymin = 0.0
ymax = 8e-8
increment = (ymax - ymin)/4.0
set yr [ymin-increment:ymax+increment]
set ylabel "$\\Phi(0,r)$"
set ytics ymin,increment,ymax
# Generate the plot
filename = "../out/Phi_00000000.dat"
plot filename using 2:3 w l lw 2 lc rgb "blue"
