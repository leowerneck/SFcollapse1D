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
set output "resources/accumulation_tau.tex"

# Set the accumulation time
tau_star = 1.22958674

# Set axis ranges
set xrange [-0.1:+1.3]
set yrange [-0.7:+0.7]

set xtics (0,0.3,0.6,0.9,1.2) format '$%.1f$'
set ytics (-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6) format '$%.1f$'

# Set axis labels
set xlabel "$\\tau$"
set ylabel "$\\phi_{\\rm central}$"

# Set the key
key_data = "Data points"
key_fit  = "Fitting curve"

# Set the data file
data_file = "../out_combined.dat"

set grid
unset key

# Generate the plot
plot data_file u 2:4 w l lw 3 lc rgb "blue" notitle

set output "resources/accumulation_Lambda.tex"

# Set axis ranges
set xrange [-2:12]
set yrange [-0.7:+0.7]

set xtics (0,2.5,5.0,7.5,10.0) format '$%.1f$'
set ytics (-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6) format '$%.1f$'

# Set axis labels
set xlabel "$\\Lambda$"
set ylabel "$\\phi_{\\rm central}$"

# Set the data file
data_file = "../out_combined.dat"

set grid
unset key

# tau = 0.625914 | diff = 0.483540 | Lambda = 0.504722
# tau = 1.109454 | diff = 0.098750 | Lambda = 2.119155
# tau = 1.208204 | diff = 0.017499 | Lambda = 3.845173
# tau = 1.225703 | diff = 0.003190 | Lambda = 5.550879
# tau = 1.228893 | diff = 0.000565 | Lambda = 7.273271
# tau = 1.229458 | diff = 0.000106 | Lambda = 8.958907
# tau = 1.229564 | diff = 0.000054 | Lambda = 10.689805
x1 = 0.504722
x2 = 2.119155
x3 = 3.845173
x4 = 5.550879
x5 = 7.273271
x6 = 8.958907
x7 = 10.689805

y_end     = +0.0
y_startb  = -0.2
y_startt  = +0.2
y_labelb  = -0.25
y_labelt  = +0.25
x_offset = 0.5
x_offset_label = 0.6

# Set a few useful arrows
set arrow 1 from (x1-x_offset),y_startb to x1,y_end filled front lw 3
set arrow 2 from (x2-x_offset),y_startt to x2,y_end filled front lw 3
set arrow 3 from (x3-x_offset),y_startb to x3,y_end filled front lw 3
set arrow 4 from (x4-x_offset),y_startt to x4,y_end filled front lw 3
set arrow 5 from (x5-x_offset),y_startb to x5,y_end filled front lw 3
set arrow 6 from (x6-x_offset),y_startt to x6,y_end filled front lw 3
set arrow 7 from (x7-x_offset),y_startb to x7,y_end filled front lw 3

# Set a few useful labels
set label 1 "$n=0$" at (x1-x_offset_label),y_labelb center
set label 2 "$n=1$" at (x2-x_offset_label),y_labelt center
set label 3 "$n=2$" at (x3-x_offset_label),y_labelb center
set label 4 "$n=3$" at (x4-x_offset_label),y_labelt center
set label 5 "$n=4$" at (x5-x_offset_label),y_labelb center
set label 6 "$n=5$" at (x6-x_offset_label),y_labelt center
set label 7 "$n=6$" at (x7-x_offset_label),y_labelb center

# Generate the plot
plot data_file u (-log(tau_star-$2)):4 w l lw 3 lc rgb "blue" notitle
