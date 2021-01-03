set terminal postscript eps color enhanced "Times" 20 size 4, 3

set title "Im {/Symbol c}_{Lindhard}(q,{/Symbol w})/{/Symbol r}({/Symbol e}_F)"
set xlabel "q/k_F"
set ylabel "{/Helvetica-Oblique @^-h}{/Symbol w}/{/Symbol e}_F"

set pm3d 
set pm3d map


set output "result/q_w_Im_chi.eps"

splot "result/q_w_chi.txt" using 1:2:4 with pm3d notitle
