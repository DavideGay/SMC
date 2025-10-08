#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.2 patchlevel 4    last modified 2018-06-01 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2018
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal qt 0 font "Sans,9"
# set output
unset clip points
set clip one
unset clip two
set errorbars front 1.000000 
set border 31 front lt black linewidth 1.000 dashtype solid
set zdata 
set ydata 
set xdata 
set y2data 
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc  bgnd fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02 
set style ellipse size graph 0.05, 0.03 angle 0 units xy
set dummy x, y
set format x "% h" 
set format y "% h" 
set format x2 "% h" 
set format y2 "% h" 
set format z "% h" 
set format cb "% h" 
set format r "% h" 
set ttics format "% h"
set timefmt "%d/%m/%y,%H:%M"
set angles radians
set tics back
unset grid
unset raxis
set theta counterclockwise right
set style parallel front  lt black linewidth 2.000 dashtype solid
set key title "" center
set key fixed left top vertical Right noreverse enhanced autotitle nobox
set key noinvert samplen 2 spacing 1.2 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
set style line 1  linecolor 1 linewidth 2.000 dashtype solid pointtype 1 pointsize default pointinterval 0 pointnumber 0
set style line 2  linecolor 2 linewidth 2.000 dashtype solid pointtype 2 pointsize default pointinterval 0 pointnumber 0
set style line 3  linecolor 3 linewidth 2.000 dashtype solid pointtype 3 pointsize default pointinterval 0 pointnumber 0
set style line 4  linecolor 4 linewidth 2.000 dashtype solid pointtype 4 pointsize default pointinterval 0 pointnumber 0
set style line 5  linecolor 5 linewidth 2.000 dashtype solid pointtype 5 pointsize default pointinterval 0 pointnumber 0
set style line 6  linecolor 6 linewidth 2.000 dashtype solid pointtype 6 pointsize default pointinterval 0 pointnumber 0
set style line 7  linecolor 7 linewidth 2.000 dashtype solid pointtype 7 pointsize default pointinterval 0 pointnumber 0
set style line 8  linecolor 8 linewidth 2.000 dashtype solid pointtype 8 pointsize default pointinterval 0 pointnumber 0
set style line 9  linecolor 9 linewidth 2.000 dashtype solid pointtype 9 pointsize default pointinterval 0 pointnumber 0
set style line 10  linecolor 10 linewidth 2.000 dashtype solid pointtype 10 pointsize default pointinterval 0 pointnumber 0
set style line 11  linecolor 11 linewidth 2.000 dashtype solid pointtype 11 pointsize default pointinterval 0 pointnumber 0
set style line 12  linecolor 12 linewidth 2.000 dashtype solid pointtype 12 pointsize default pointinterval 0 pointnumber 0
set style line 13  linecolor 13 linewidth 2.000 dashtype solid pointtype 13 pointsize default pointinterval 0 pointnumber 0
set style line 14  linecolor 14 linewidth 2.000 dashtype solid pointtype 14 pointsize default pointinterval 0 pointnumber 0
set style line 15  linecolor 15 linewidth 2.000 dashtype solid pointtype 15 pointsize default pointinterval 0 pointnumber 0
set style line 16  linecolor 16 linewidth 2.000 dashtype solid pointtype 16 pointsize default pointinterval 0 pointnumber 0
set style line 17  linecolor 17 linewidth 2.000 dashtype solid pointtype 17 pointsize default pointinterval 0 pointnumber 0
set style line 18  linecolor 18 linewidth 2.000 dashtype solid pointtype 18 pointsize default pointinterval 0 pointnumber 0
set style line 19  linecolor 19 linewidth 2.000 dashtype solid pointtype 19 pointsize default pointinterval 0 pointnumber 0
set style line 20  linecolor 20 linewidth 2.000 dashtype solid pointtype 20 pointsize default pointinterval 0 pointnumber 0
unset style arrow
set style histogram clustered gap 2 title textcolor lt -1
unset object
set style textbox transparent margins  1.0,  1.0 border  lt -1 linewidth  1.0
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding utf8
unset polar
unset parametric
unset decimalsign
unset micro
unset minussign
set view 60, 30, 1, 1
set view azimuth 0
set rgbmax 255
set samples 100, 100
set isosamples 10, 10
set surface 
unset contour
set cntrlabel  format '%8.3g' font '' start 5 interval 20
set mapping cartesian
set datafile separator whitespace
set datafile commentschars '#!%'
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5 unsorted
set cntrparam firstlinetype 0
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data lines
set style function lines
unset xzeroaxis
unset yzeroaxis
unset zzeroaxis
unset x2zeroaxis
unset y2zeroaxis
set xyplane relative 0.5
set tics scale  1, 0.5, 1, 1, 1
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set mrtics default
set nomttics
set xtics border in scale 1,0.5 mirror norotate  autojustify
set xtics  norangelimit autofreq 
set ytics border in scale 1,0.5 mirror norotate  autojustify
set ytics  norangelimit autofreq 
set ztics border in scale 1,0.5 nomirror norotate  autojustify
set ztics  norangelimit autofreq 
unset x2tics
unset y2tics
set cbtics border in scale 1,0.5 mirror norotate  autojustify
set cbtics  norangelimit autofreq 
set rtics axis in scale 1,0.5 nomirror norotate  autojustify
set rtics  norangelimit autofreq 
unset ttics
set title "" 
set title  font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  font "" norotate
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "" 
set xlabel  font "" textcolor lt -1 norotate
set x2label "" 
set x2label  font "" textcolor lt -1 norotate
set xrange [ * : * ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set ylabel "" 
set ylabel  font "" textcolor lt -1 rotate
set y2label "" 
set y2label  font "" textcolor lt -1 rotate
set yrange [ * : * ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zlabel "" 
set zlabel  font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse writeback
set cblabel "" 
set cblabel  font "" textcolor lt -1 rotate
set cbrange [ * : * ] noreverse writeback
set rlabel "" 
set rlabel  font "" textcolor lt -1 norotate
set rrange [ * : * ] noreverse writeback
unset logscale
unset jitter
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "en_US.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles noborder corners2color mean
set pm3d nolighting
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2 size screen 0.05, 0.6 front  noinvert bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit brief errorvariables nocovariancevariables errorscaling prescale nowrap v5
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x)
round2(x, n) = round(x*10**n)*10.0**(-n)
bin(x,s) = s*floor(x/s)
echo(x) = x
GNUTERM = "qt"
pdfplot = "set term push; set term pdfcairo enhanced color rounded size 8,6 lw 1.5 font 'Arial,28'; set border lw 1.5 front; set out 'temp.pdf' ; replot; set out; set term pop"
epsplot = "set term push; set term epscairo enhanced color rounded size 8,6 lw 1.5 font 'Arial,28'; set border lw 1.5 front; set out 'temp.eps' ; replot; set out; set term pop"
pngplot = "set term push; set term pngcairo enhanced color rounded size 1600,1200 lw 1.5 font 'Arial,28'; set border lw 1.5 front; set out 'temp.png' ; replot; set out; set term pop"
matlab = "set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB;           set palette defined  (0  0.0 0.0 0.5,                                 1  0.0 0.0 1.0,                                 2  0.0 0.5 1.0,                                 3  0.0 1.0 1.0,                                 4  0.5 1.0 0.5,                                 5  1.0 1.0 0.0,                                 6  1.0 0.5 0.0,                                 7  1.0 0.0 0.0,                                 8  0.5 0.0 0.0 )"
rainbow = "set palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV;            set palette defined ( 0 0 1 1, 1 0.85 1 1 )"
mandelbrot = "set sample 500; set isosample 500; set pm3d map;               complex(x,y) = x*{1,0}+y*{0,1};               mandel(x,y,z,n) = (abs(z)>2.0 || n>=100) ? n : mandel(x,y,z*z+complex(x,y),n+1);               splot mandel(x,y,{0,0},0) w pm3d not"
b = 1334.6226171233
FIT_CONVERGED = 1
FIT_NDF = 109
FIT_STDFIT = 4813.27205232354
FIT_WSSR = 2525267075.61499
FIT_P = 0.0
FIT_NITER = 3
b_err = 10.5283161526558
a = 1100.92254070091
a_err = 7.27586086402037
## Last datafile plotted: "recursion2.dat"
p "recursion1.dat" u ($1/1000) w p t "recursion events", a/1000*x w lines lw 2 t "ax, a=1.10092", "recursion2.dat" u ($1/1000) w p t "recursion events, 2", b/1000*x w l lw 2 t "bx, b=1.33462"
## fit a*x "recursion1.dat" via a
#    EOF
