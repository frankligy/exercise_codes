
# usage: ./copy.sh /Users/ligk2e/Desktop/temp/folder1 /Users/ligk2e/Desktop/temp/folder2     
# please don't include final forward slash /
# recommened to use full path
# the command means: I want to only copy unique file in folder1 to folder2

# comment out one version, make sure only run one version.


# stable version
for i in ${1}/*; do file1=$(basename $i); echo ${file1} >> folder1filename.txt; done
for j in ${2}/*; do file2=$(basename $j); echo ${file2} >> folder2filename.txt; done 
comm -23 folder1filename.txt folder2filename.txt > unique1.txt
cat unique1.txt | while read line; do cp ${1}/${line} ${2}/; done
rm -i folder1filename.txt folder2filename.txt unique1.txt


# easier version
ls ${1} > folder1filename.txt
ls ${2} > folder2filename.txt
comm -23 folder1filename.txt folder2filename.txt > unique1.txt
cat unique1.txt | while read line; do cp ${1}/${line} ${2}/; done
rm -i folder1filename.txt folder2filename.txt unique1.txt
