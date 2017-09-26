#er homedir and remove previous file
rm -f ‘$2′
gnuplot << EOP
 ### set data source file
datafile = '$1'

### set graph type and size
#set terminal jpeg size 640,480
set terminal pngcairo size 640,480 
set autoscale
set logscale y
set xtic auto
set ytic auto
#set yrange [1:1000000]
### set titles
set grid x y
set xlabel "Vertices"
set ylabel "Degree"
unset key
### set output filename
set output '$2'

### build graph
# plot datafile with lines
plot datafile using 1:2 with lines axes x1y1  
#plot datafile using 1:2 title "user" with lines axes x1y1

EOP
