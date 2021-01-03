set terminal postscript eps color enhanced "Times" 20 size 4, 3

set xlabel "q/k_F"
set ylabel "{/Symbol c}_{Lindhard}(q,{/Symbol w}=0)/{/Symbol r}({/Symbol e}_F)"
set grid

set output "result/q_chi_w=0.eps"

plot "result/q_chi_w=0.txt" using 1:2 with lines notitle
