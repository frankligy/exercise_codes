# This script can automatically run DoubletDecon main function n times and automatically take intersection
# of the doublets coming back from each run, end up with a single file. ./output/intersection.txt. In addtion
# you can also directly change the key arguments for running DoubletDecon, including filename, removeCC, rhop,
# centroid, num_doubs, only50, etc via passing corresponding arguments to this .sh script.

# example usage:

# ./multiRun.sh 20 0 1 1 100 0 4

# 20 means run 20 times with subsequent arguments
# 0 means removeCC=FALSE   (rule of thumb: 0 means FALSE, 1 means TRUE)
# 1 means rhop=1
# 1 means centroid=TRUE    (rule of thumb: 0 means FALSE, 1 means TRUE)
# 100 means num_doubs=100
# 0 means only50=FALSE     (rule of thumb: 0 means FALSE, 1 means TRUE)
# 4 means min_uniq=4



cd /Users/ligk2e/Downloads/DoubletDecon-master/R  # change to the folder with .R files


i=1
while [ ${i} -le $1 ]   
do 
filename="round${i}"
../unitRun.R ${filename} $2 $3 $4 $5 $6 $7   # run DoubletDecon
cat ../output/Final_doublets_groups_${filename}.txt | cut -f 1 | sort > ../output/${filename}doublets.txt   # extract doublets from each run
if [ -e ../output/doubletsResult_round$[i-1].txt ]
then 
    comm -12 ../output/${filename}doublets.txt ../output/doubletsResult_round$[i-1].txt | sort > ../output/doubletsResult_${filename}.txt  # take intersection
else
    mv ../output/${filename}doublets.txt ../output/doubletsResult_${filename}.txt
fi
i=$[$i+1]
done
mv ../output/doubletsResult_round$1.txt ../output/intersection.txt  # rename, done!



