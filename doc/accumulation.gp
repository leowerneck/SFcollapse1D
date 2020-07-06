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
set output "resources/accumulation_tau_BSSN.tex"

set border #@WHITE

# Set axis ranges
set xrange [-0.1:+1.7]
set yrange [-0.7:+0.7]

# set xtics (0,0.3,0.6,0.9,1.2) format '$%.1f$'
# set xtics (0,0.5,1.0,1.5,2.0,2.5) format '$%.1f$'
# set xtics (0,1,2,3,4,5) format '$%.1f$'
set xtics (0,0.4,0.8,1.2,1.6) format '$%.1f$'
set ytics (-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6) format '$%.1f$'

# Set axis labels
set xlabel "$\\tau$"
set ylabel "$\\psi_{\\rm central}$"

# Set the data file
data_file = "../proper_time_central_BSSN.dat"

set grid #@WHITE
unset key

# Generate the plot
plot data_file u 1:3 w l lw 3 @BLUE notitle

##########################################################################################
##########################################################################################
##########################################################################################

set output "resources/accumulation_Lambda_BSSN.tex"

# Set axis ranges
set xrange [-2:10]
set yrange [-0.7:+0.7]

set border #@WHITE

set xtics (0,2,4,6,8) format '$%.0f$'
set ytics (-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6) format '$%.1f$'

# Set axis labels
set xlabel "$\\Lambda$"
set ylabel "$\\psi_{\\rm central}$"

# Set the data file
data_file = "../proper_time_central_BSSN.dat"

set grid #@WHITE
unset key

# Gaussian shell
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

# Gaussian shell v2

# tau = 3.372778 | diff = 1.149750 | Lambda = -0.370178
# tau = 4.522528 | diff = 0.245291 | Lambda = 1.209850
# tau = 4.767819 | diff = 0.043332 | Lambda = 2.938391
# tau = 4.811151 | diff = 0.007898 | Lambda = 4.644055
# tau = 4.819049 | diff = 0.001398 | Lambda = 6.364845
# tau = 4.820447 | diff = 0.000266 | Lambda = 8.037696
# tau = 4.820712 | diff = 0.000086 | Lambda = 9.762958

# x1 = -0.370178
# x2 = 1.209850
# x3 = 2.938391
# x4 = 4.644055
# x5 = 6.364845
# x6 = 8.037696
# x7 = 9.762958

# Tanh shell

# tau = 1.407498 | diff = 0.962253 | Lambda = -0.170959
# tau = 2.369752 | diff = 0.184436 | Lambda = 1.495268
# tau = 2.554187 | diff = 0.032531 | Lambda = 3.225081
# tau = 2.586718 | diff = 0.005930 | Lambda = 4.930668
# tau = 2.592648 | diff = 0.001050 | Lambda = 6.651730
# tau = 2.593698 | diff = 0.000200 | Lambda = 8.325791
# tau = 2.593898 | diff = 0.000047 | Lambda = 10.083830

# x1 = -0.170959
# x2 = 1.495268
# x3 = 3.225081
# x4 = 4.930668
# x5 = 6.651730
# x6 = 8.325791
# x7 = 10.083830

y_end     = +0.0
y_startb  = -0.2
y_startt  = +0.2
y_labelb  = -0.25
y_labelt  = +0.25
x_offset = 0.5
x_offset_label = 0.6

# Set a few useful arrows
# set arrow 1 from (x1-x_offset),y_startb to x1,y_end filled front lw 3
# set arrow 2 from (x2-x_offset),y_startt to x2,y_end filled front lw 3
# set arrow 3 from (x3-x_offset),y_startb to x3,y_end filled front lw 3
# set arrow 4 from (x4-x_offset),y_startt to x4,y_end filled front lw 3
# set arrow 5 from (x5-x_offset),y_startb to x5,y_end filled front lw 3
# set arrow 6 from (x6-x_offset),y_startt to x6,y_end filled front lw 3
# set arrow 7 from (x7-x_offset),y_startb to x7,y_end filled front lw 3

# Set a few useful labels
# set label 1 "$n=0$" at (x1-x_offset_label),y_labelb center
# set label 2 "$n=1$" at (x2-x_offset_label),y_labelt center
# set label 3 "$n=2$" at (x3-x_offset_label),y_labelb center
# set label 4 "$n=3$" at (x4-x_offset_label),y_labelt center
# set label 5 "$n=4$" at (x5-x_offset_label),y_labelb center
# set label 6 "$n=5$" at (x6-x_offset_label),y_labelt center
# set label 7 "$n=6$" at (x7-x_offset_label),y_labelb center

# Generate the plot
plot data_file u 2:3 w l lw 3 @BLUE notitle
