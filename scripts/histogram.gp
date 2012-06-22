reset

set terminal push
set terminal epslatex standalone monochrome solid size 12cm, 8cm
set output 'histogram.tex'

unset key

set xrange [-pi:pi]

set xlabel 'Angle'
set ylabel 'Frequency'
plot 'histogram.dat' with lines linewidth 2

unset output
set terminal pop

system("latex histogram.tex && dvips histogram.dvi && ps2eps -f histogram.ps && dvipdf histogram.dvi")
