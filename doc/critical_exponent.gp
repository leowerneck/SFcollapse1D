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

# Set terminal and output
set terminal epslatex
set output "resources/critical_exponent_BSSN.tex"

# Set axis ranges
set xrange [-26:-7]
set yrange [2:18]

set xtics (-25,-20,-15,-10)
set ytics (4,8,12,16)

# Set axis labels
set xlabel "$\\ln\\left(\\eta_{*}-\\eta\\right)$"
set ylabel "$\\ln\\left(\\rho^{\\max}_{\\rm central}\\right)$"

# Set the key
key_data = "Data points"
key_fit  = "Fitting curve"

# Set the data file
data_file = "../rho_BSSN.dat"
fit_file  = "../fitting_curve_critical_exponent_BSSN.dat"

set grid
set key width box lc rgb "#BDBDBD"

# Generate the plot
plot fit_file u 1:2 w l lw 3 lc rgb "red" notitle, data_file u 4:5 w p ls 4 lw 3 ps 2 lc rgb "blue" notitle, data_file u 4:(1000*$5) w p ls 4 lw 3 ps 2 lc rgb "blue" t key_data, fit_file u 1:(1000*$2) w l lw 3 lc rgb "red" t key_fit
