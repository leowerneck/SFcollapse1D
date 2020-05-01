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
set output "resources/critical_lapse_best.tex"

# Set grid
set grid

# Set plot keys
title_weak = "$\\eta = 0.3364266156435$"
title_strg = "$\\eta = 0.3364266156436$"
# title_weak = "$\\eta = 0.33634459617$"
# title_strg = "$\\eta = 0.33634459618$"

# Multiplot layout
# .-----------------.
# |     Plot 1      |
# .--------.--------.
# | Plot 3 | Plot 2 |
# .--------.--------.

BORDER_L = 0.10
BORDER_R = 0.07
BORDER_B = 0.20
BORDER_T = 0.07

SPACING_X = 0.10
SPACING_Y = 0.10

NUM_ROWS = 2
NUM_COLS = 2

PLT_W = (1 - BORDER_L - BORDER_R - (NUM_COLS-1)*SPACING_X)/NUM_COLS
PLT_H = (1 - BORDER_B - BORDER_T - (NUM_ROWS-1)*SPACING_Y)/NUM_ROWS

set multiplot

# .--------.
# | Plot 1 |
# .--------.

# Configure key
set key width -2.5 height 0.5 box lc rgb "#BDBDBD"

set xlabel "$t$" offset 0,-9
set ylabel "$\\alpha_{\\rm central}$"

file_w = "../out_weak_best.dat"
file_s = "../out_strong_best.dat"

set tmargin at screen 1 - BORDER_T
set bmargin at screen 1 - BORDER_T - PLT_H
set lmargin at screen BORDER_L
set rmargin at screen 1 - BORDER_R

# Set x and y axis ranges for plot 1
set xrange [-0.5:6.5]
set yrange [-0.05:0.65]

# Set x and y tics for plot 1
set xtics (0,1,2,3,4,5,6)
set ytics ("0.0" 0,0.2,0.4,0.6)

# Draw "zoom" objects that connect plot 1 to plot 2
# First the rectangle
set object 1 rectangle from 5.05,-0.02 to 5.35,0.25 lw 2 dt 2 fs empty border lc rgb "#BCBCBC"

# Set auxiliary variables
end_x1 = 1 - BORDER_R - PLT_W
end_y1 = 1 - BORDER_T - PLT_H - SPACING_Y
end_x2 = 1 - BORDER_R
end_y2 = 1 - BORDER_T - PLT_H - SPACING_Y

# Draw lines that link to plot 2
# set arrow 2 from 5.4,-0.02 to screen end_x2,end_y2 lw 2 dt 2 lc rgb "#BCBCBC" nohead
set arrow 1 from 5.05,-0.02 to screen end_x1,end_y1 nohead filled back lw 2 lc rgb "#BABABA" dt 2
set arrow 2 from 5.35,-0.02 to screen end_x2,end_y2 nohead filled back lw 2 lc rgb "#BBBBBB" dt 2

# Generate plot 1
plot file_w using 1:2 with lines linewidth 3 @BLUE title title_weak, file_s using 1:2 with lines linewidth 3 dashtype 2 @ORANGE title title_strg

# .--------.
# | Plot 2 |
# .--------.
unset xlabel
unset ylabel
unset key
unset arrow 1
unset arrow 2
unset object 1
set tmargin at screen 1 - BORDER_T - PLT_H - SPACING_Y
set bmargin at screen BORDER_B
set lmargin at screen 1 - BORDER_R - PLT_W
set rmargin at screen 1 - BORDER_R

# Draw "zoom" objects that connect plot 2 to plot 3
# First the rectangle
set object 1 rectangle from 5.3,-0.01 to 5.32,0.14 lw 2 dt 2 fs empty border lc rgb "#BCBCBC"

# Set auxiliary variables
end_x1 = BORDER_L + PLT_W
end_y1 = 1 - BORDER_T - PLT_H - SPACING_Y
end_x2 = BORDER_L + PLT_W
end_y2 = BORDER_B

# Draw lines that link to plot 3
set arrow 1 from 5.3,+0.14 to screen end_x1,end_y1 nohead filled back lw 2 lc rgb "#BABABA" dt 2
set arrow 2 from 5.3,-0.01 to screen end_x2,end_y2 nohead filled back lw 2 lc rgb "#BBBBBB" dt 2

# Set x and y axis ranges for plot 2
set xrange [5.05:5.35]
set yrange [-0.02:0.25]

# Set x and y tics for plot 2
set xtics ("5.08" 5.08,"5.20" 5.2,"5.32" 5.32)
set ytics ("0.0" 0,0.1,0.2)

plot file_w using 1:2 with lines linewidth 3 @BLUE title title_weak, file_s using 1:2 with lines linewidth 3 dashtype 2 @ORANGE title title_strg

# .--------.
# | Plot 3 |
# .--------.
unset xlabel
unset ylabel
unset key
unset arrow 1
unset arrow 2
unset object 1
set tmargin at screen 1 - BORDER_T - PLT_H - SPACING_Y
set bmargin at screen BORDER_B
set lmargin at screen BORDER_L
set rmargin at screen BORDER_L + PLT_W

# Set x and y axis ranges for plot 3
set xrange [5.300:5.320]
set yrange [-0.02:0.14]

# Set x and y tics for plot 3
set xtics ("5.302" 5.302,"5.310" 5.310,"5.318" 5.318)
set ytics ("0.00" 0,0.06,"0.12" 0.12)

plot file_w using 1:2 with lines linewidth 3 @BLUE title title_weak, file_s using 1:2 with lines linewidth 3 dashtype 2 @ORANGE title title_strg
