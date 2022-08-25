#!/bin/gnuplot

set term epscairo
set output "out.eps"

set title("H_2O Dyson calculations")
set xlabel("Ionisation energy (a.u.)")
set ylabel("Dyson Norm^2 (arb. units)")
set grid

p "svp.csv" w l t "SVP", "tzvpp.csv" w l t "TZVPP", "qzvpp.csv" w l t "QZVPP"
