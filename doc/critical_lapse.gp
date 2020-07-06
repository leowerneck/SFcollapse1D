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
set output "resources/critical_lapse_BSSN_black.tex"

# Set grid
set grid @WHITE
set border @WHITE

# Set plot keys
# Gaussian shell
# title_weak = "$\\eta = 0.3364266156435$"
# title_strg = "$\\eta = 0.3364266156436$"

# Gaussian shell v2
# title_weak = "$\\eta = 0.01591203885327$"
# title_strg = "$\\eta = 0.01591203885328$"

# Tanh shell
# title_weak = "$\\eta = 0.2914370451806$"
# title_strg = "$\\eta = 0.2914370451807$"

# Gaussian shell (BSSN)
title_weak = "$\\eta = 0.30332394090$"
title_strg = "$\\eta = 0.30332394095$"

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
# set key width -1.5 height 0.5 at 4,1.2 box lc rgb "#BDBDBD"
set key width -1.5 height 0.5 at 4,1.2 box @WHITE

set xlabel "$t$" offset 0,-9
set ylabel "$\\alpha_{\\rm central}$"

# Gaussian shell
# file_w = "../out_weak_best.dat"
# file_s = "../out_strong_best.dat"

# Gaussian shell v2
# file_w = "../out_weak_v2.dat"
# file_s = "../out_strong_v2.dat"

# Tanh shell
# file_w = "../out_weak_tanh.dat"
# file_s = "../out_strong_tanh.dat"

# Gaussian shell (BSSN)
file_w = "../weak_BSSN.dat"
file_s = "../strong_BSSN.dat"

set tmargin at screen 1 - BORDER_T
set bmargin at screen 1 - BORDER_T - PLT_H
set lmargin at screen BORDER_L
set rmargin at screen 1 - BORDER_R

# Set x and y axis ranges for plot 1
# Gaussian shell
set xrange [-0.4:7.0]
set yrange [-0.2:1.2]

# Gaussian shell v2
# set xrange [-0.5:14.5]
# set yrange [-0.05:0.75]

# Tanh shell
# set xrange [-0.5:10.5]
# set yrange [-0.05:0.75]

# Set x and y tics for plot 1
# Gaussian shell
set xtics (0,1,2,3,4,5,6) format "$%.0f$"
set ytics (0,0.2,0.4,0.6,0.8,1) format "$%.1f$"

# Gaussian shell v2
# set xtics (0,2,4,6,8,10,12,14) format "$%.0f$"
# set ytics (0,0.35,0.7) format "$%.2f$"

# Tanh shell
# set xtics (0,2,4,6,8,10)
# set ytics ("0.00" 0, "0.35" 0.35, "0.70" 0.70)

# Draw "zoom" objects that connect plot 1 to plot 2
# First the rectangle
# Gaussian shell
# rectangle_left_border  = +5.06
# rectangle_right_border = +5.34
# rectangle_top_border   = +0.23
# rectangle_bot_border   = -0.03

# Gaussian shell v2
# rectangle_left_border  = +12.95
# rectangle_right_border = +13.55
# rectangle_top_border   = +0.23
# rectangle_bot_border   = -0.03

# Tanh shell
# rectangle_left_border  = +9.3
# rectangle_right_border = +9.9
# rectangle_top_border   = +0.23
# rectangle_bot_border   = -0.03

# Gaussian shell (BSSN)
rectangle_left_border  = +5.6
rectangle_right_border = +6.6
rectangle_top_border   = +0.9
rectangle_bot_border   = -0.1

set object 1 rectangle from rectangle_left_border,rectangle_bot_border to rectangle_right_border,rectangle_top_border lw 2 dt 2 fs empty border lc rgb "#BCBCBC"

# Set auxiliary variables
end_x1 = 1 - BORDER_R - PLT_W
end_y1 = 1 - BORDER_T - PLT_H - SPACING_Y
end_x2 = 1 - BORDER_R
end_y2 = 1 - BORDER_T - PLT_H - SPACING_Y

# Draw lines that link to plot 2
set arrow 1 from rectangle_left_border,rectangle_bot_border to screen end_x1,end_y1 nohead filled back lw 2 lc rgb "#BABABA" dt 2
set arrow 2 from rectangle_right_border,rectangle_bot_border to screen end_x2,end_y2 nohead filled back lw 2 lc rgb "#BBBBBB" dt 2

# Generate plot 1
plot file_w using 1:2 with lines linewidth 3 @CYAN title title_weak, file_s using 1:2 with lines linewidth 3 dashtype 2 @ORANGE title title_strg

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

# Set x and y axis ranges for plot 2
set xrange [rectangle_left_border:rectangle_right_border]
set yrange [rectangle_bot_border:rectangle_top_border]

# Draw "zoom" objects that connect plot 2 to plot 3
# First the rectangle
# Gaussian shell
# rectangle_left_border  = +5.3
# rectangle_right_border = +5.32
# rectangle_top_border   = +0.14
# rectangle_bot_border   = -0.02

# Gaussian shell
# rectangle_left_border  = +13.465
# rectangle_right_border = +13.505
# rectangle_top_border   = +0.12
# rectangle_bot_border   = -0.02

# Tanh shell
# rectangle_left_border  = +9.850
# rectangle_right_border = +9.882
# rectangle_top_border   = +0.12
# rectangle_bot_border   = -0.02

# Gaussian shell (BSSN)
rectangle_left_border  = +6.53
rectangle_right_border = +6.59
rectangle_top_border   = +0.55
rectangle_bot_border   = -0.05

set object 1 rectangle from rectangle_left_border,rectangle_bot_border to rectangle_right_border,rectangle_top_border lw 2 dt 2 fs empty border lc rgb "#BCBCBC"

# Set auxiliary variables
end_x1 = BORDER_L + PLT_W
end_y1 = 1 - BORDER_T - PLT_H - SPACING_Y
end_x2 = BORDER_L + PLT_W
end_y2 = BORDER_B

# Draw lines that link to plot 3
set arrow 1 from rectangle_left_border,rectangle_top_border to screen end_x1,end_y1 nohead filled back lw 2 lc rgb "#BABABA" dt 2
set arrow 2 from rectangle_left_border,rectangle_bot_border to screen end_x2,end_y2 nohead filled back lw 2 lc rgb "#BBBBBB" dt 2

# Set x and y tics for plot 2
# Gaussian shell
# set xtics (5.08,5.2,5.32) format '$%.2f$'
# set ytics (0,0.1,0.2) format '$%.1f$'

# Gaussian shell v2
# set xtics (13,13.25,13.5) format '$%.2f$'
# set ytics (0,0.1,0.2) format '$%.1f$'

# Tanh shell
# set xtics (9.4,9.6,9.8) format '$%.1f$'
# set ytics ("0.0" 0,0.1,0.2)

# Gaussian shell (BSSN)
set xtics (5.9,6.2,6.5) format '$%.1f$'
set ytics (0,0.4,0.8) format '$%.1f$'

plot file_w using 1:2 with lines linewidth 3 @CYAN title title_weak, file_s using 1:2 with lines linewidth 3 dashtype 2 @ORANGE title title_strg

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
set xrange [rectangle_left_border:rectangle_right_border]
set yrange [rectangle_bot_border:rectangle_top_border]

# Set x and y tics for plot 3
# Gaussian shell
# set xtics (5.302,5.31,5.318) format '$%.3f$'
# set ytics (0,0.06,0.12) format '$%.2f$'

# Gaussian shell
# set xtics (13.47,13.485,13.5) format '$%.3f$'
# set ytics (0,0.05,0.10) format '$%.2f$'

# Tanh shell
# set xtics (9.854,9.866,9.878) format '$%.3f$'
# set ytics ("0.00" 0,0.05,"0.10" 0.10)

# Gaussian shell (BSSN)
set xtics (6.54,6.56,6.58) format '$%.2f$'
set ytics (0,0.25,0.5) format '$%.2f$'

plot file_w using 1:2 with lines linewidth 3 @CYAN title title_weak, file_s using 1:2 with lines linewidth 3 dashtype 2 @ORANGE title title_strg
