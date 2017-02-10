for((i=0;i<=7;i++))
do
	cat level-"$i"* | sort -g | uniq | wc -l
done
