
qz = 10.
rf = 5e2
rn = 1e1
al = 1.1
bt = 100
#f(x) = (rn-rf)/(1 + bt*exp(-al*(x-qz)))+rf
f(x) = 1/(1 + bt*exp(-al*(x-qz)))

plot [:] 'Mem1.txt' u 1:(f($6)) t 'memrstor status w', 'Mem1.txt' u 1:6 t 'total charge Q'
pause -1
plot [:] 'Mem1.txt' u 1:(f($6)) t 'memrstor status w'
pause -1
set terminal png
set output 'Mem1.png'
rep
q

