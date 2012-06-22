reset

f(x) = a * x + b
fit f(x) 'last-regression-slope.dat' using 1:3 via a, b

set terminal push
set terminal epslatex standalone monochrome solid size 14cm, 8cm

set output 'results1.tex'

set grid
unset key
set xlabel 'Time'
set ylabel 'Slope of linear regression'
plot 'time-vs-slopes.dat' with lines

set output 'results2.tex'

unset grid
set key below
set title 'Final error'
set xlabel '$\Delta t$ (Log. base 10)'
set ylabel 'Error (Log. base 10)'
plot 'last-regression-slope.dat' using 1:3 with lines linetype 0 linewidth 2 \
     title sprintf('Least-squares line of slope %02.02f', a), \
     '' using 1:2 with points pointtype 7 pointsize 2 title 'Actual data'

unset output
set terminal pop

system('latex results1.tex && dvips results1.dvi && ps2eps -f results1.ps && dvipdf results1.dvi')
system('latex results2.tex && dvips results2.dvi && ps2eps -f results2.ps && dvipdf results2.dvi')
