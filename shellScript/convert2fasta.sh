#! /bin/bash


task=0702

cat $task.txt | while read line; do echo ${line//[:space:]}; done > $task.new.txt
rm $task.txt

touch $task.fasta

counter=0
cat $task.new.txt | while read line
do 
counter=$[$counter+1]
echo ">peptide$counter" >> $task.fasta
echo $line >> $task.fasta
done

task1=${task:0:2}
task2=${task:2:4}

./netCTLpan -a HLA-C$task1:$task2 $task.fasta > $task.result.txt



cat $task.result.txt | tail -n +35 | grep peptide | grep -v identi | sort -u -k2,2 | awk '{print $8}'

cat $task.result.txt | tail -n +35 | grep peptide | grep -v identi | sort -u -k2,2 | awk '{print $8}' | wc -l

# cat 5701.result.txt | tail -n +35 | grep peptide | grep -v identi | sort -u -k2,2 | awk '{print $2}' > get.txt
# for i in {1..45}; do echo "peptide$i"; done > all.txt
# comm -23 all.txt get.txt 
# cat 5101.result.txt | tail -n +35 | grep peptide | grep -v identi | sort -u -k2,2 > whole.txt



