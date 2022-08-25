#!/bin/gnuplot

set term epscairo
set output "out.eps"

set title("H_2O Dyson calculations")
set xlabel("Ionisation energy (a.u.)")
set ylabel("Dyson Norm^2 (arb. units)")
set grid

p "nstate3.csv" w l t "N_{state} = 3", "nstate6.csv" w l t "N_{state} = 6"
