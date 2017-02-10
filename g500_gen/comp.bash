for((k=0;k<4;k++))
do
for((i=0;i<4;i++))
do
	for((j=0;j<4;j++))
	do
	echo "thread $k, row $i, col $j"
	diff scale-22-rank-$k-of-4-par-rowcol-$i-"$j".txt scale-22-rank-$k-of-4-par-rowcol-$i-"$j".dat
	done
done

done
