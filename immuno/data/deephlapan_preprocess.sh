#! /bin/bash

task='new_testing'

cat $task.txt | cut -f 1 | tail -n +2 > ${task}_peptide.txt
cat $task.txt | cut -f 2 | tail -n +2 | while read line; do f=${line:0:7};b=${line:7:9};c=$f:$b; echo $c; done > ${task}_hla.txt
cat $task.txt | cut -f 3 | tail -n +2 > ${task}_immuno.txt



touch ${task}_annotation.txt
length=`cat ${task}_peptide.txt | wc -l`
echo $length
for i in {1..1399}; do echo "new_testing" >> ${task}_annotation.txt; done


paste -d "," ${task}_annotation.txt ${task}_hla.txt ${task}_peptide.txt > ${task}_combine.txt

{ echo "Annotation,HLA,peptide"; cat ${task}_combine.txt; } > ${task}_final.txt