set title "Approximate solution"
set key off
splot 'data/solution.data'using 1:2:3 with lines palette lw 2
