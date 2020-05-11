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
set terminal epslatex size 5,3
set output "resources/critical_exponent_oscillation.tex"

# Set axis ranges
set xrange [-34:-6]
set yrange [-1:+1]

# Set xtics and ytics positions
set xtics (-30,-25,-20,-15,-10)
set ytics (-0.75,-0.50,-0.25,0.00,+0.25,+0.50,+0.75)

# Set the grid
set grid

# Set the data files
data_file = "../numerical_data_critical_exponent.dat"
fit_file  = "../fitting_curve_critical_exponent.dat"

# Set axis labels
set xlabel "$\\ln\\left(\\eta_{*}-\\eta\\right)$" offset 12,0
set ylabel "$\\ln\\left(\\rho^{\\max}_{\\rm central}\\right) - C + 2\\gamma\\ln\\left|\\eta_{*}-\\eta\\right|$"

# Set multiplot environment
set multiplot layout 1,2 margins 0.1,0.98,0.1,0.98 spacing 0,0

# ---------- START OF PLOT 1 ----------

# Set x position of vertical lines
x1 = -28.0
x2 = x1 + 2.3
x3 = x2 + 2.3
x4 = x3 + 2.3
x5 = x4 + 2.3
x6 = x5 + 2.3
x7 = x6 + 2.3
x8 = x7 + 2.3
x9 = x8 + 2.3

y_i = -0.4
y_f = +0.4

y_label_bot = -0.5
y_label_top = +0.5

# Draw vertical lines to estimate the period
set arrow 1 from x1,y_i to x1,y_f nohead lw 2 dt 2
set arrow 2 from x2,y_i to x2,y_f nohead lw 2 dt 2
set arrow 3 from x3,y_i to x3,y_f nohead lw 2 dt 2
set arrow 4 from x4,y_i to x4,y_f nohead lw 2 dt 2
set arrow 5 from x5,y_i to x5,y_f nohead lw 2 dt 2
set arrow 6 from x6,y_i to x6,y_f nohead lw 2 dt 2
set arrow 7 from x7,y_i to x7,y_f nohead lw 2 dt 2
set arrow 8 from x8,y_i to x8,y_f nohead lw 2 dt 2
set arrow 9 from x9,y_i to x9,y_f nohead lw 2 dt 2

# Set ytic format
set ytics format '$%.2f$'

# Remove the key for the first plot
unset key

# Generate the plot
plot data_file u 1:2 w p ls 4 lw 2 ps 1.5 lc rgb "blue" notitle

# ----------  END OF PLOT 1  ----------
# -------------------------------------
# ---------- START OF PLOT 2 ----------

# Set vertical lines for second plot
unset arrow 1
unset arrow 2
unset arrow 3
unset arrow 4
unset arrow 5
unset arrow 6
unset arrow 7
unset arrow 8
unset arrow 9

# Remove x and y axis label for second plot
unset xlabel
unset ylabel
unset label

# Remove ytics labels for second plot
set ytics format ''

# Set key and draw a box around it
set key width -2.5 box lc rgb "#BDBDBD"

# Set the key
key_data = "Data points"
key_fit  = "Fitting curve"

# Generate the plot
plot data_file u 1:2 w p ls 4 lw 2 ps 1.5 lc rgb "blue" t key_data, fit_file u 1:3 w l lw 3 lc rgb "red" t key_fit

# ----------  END OF PLOT 2  ----------

unset multiplot
