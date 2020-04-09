awk 'BEGIN {OFS="\t";}
    { if (NR > 3) {
	sum = 0;
	for (i=1;i<=NF;i++){
		sum += $i;
	}
	avg = sum/NF;
	sum = 0;
	for (i=0;i<=NF;i++){
		sum += ($i - avg)^2;
	}
	sd = sqrt(sum/NF);
	print $0, avg, sd; }
}' learnAWK.txt 