#BSUB -W 5:00
#BSUB -M 250000
#BSUB -n 1
#BSUB -J OptiType
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

cd $(pwd)
module load singularity/3.1.0
mkdir TCGA-A8-A08I
singularity run -B /data/salomonis2/LabFiles/Frank-Li/OptiType/TCGA-A8-A08I:/mnt my_software.sif \
    -i TCGA-A8-A08I-01A-11R-A00Z-07.qsort.read1.fastq TCGA-A8-A08I-01A-11R-A00Z-07.qsort.read2.fastq --rna -v -o /mnt

function program() {

sample=$1

mkdir $(pwd)/${sample}

echo 'copying files'
for i in /data/salomonis2/FASTQs/TCGA-Breast-Cancer/*.fq.gz 
do 
    base1=`basename ${i}`
    case ${i} in
        *"${sample}"*)
        cp ${i} $(pwd) 
        gzip -d $(pwd)/${base1}       
    esac
done

echo 'rename files for Optitype'
declare -a int_files
for j in $(pwd)/*.fq
do
    if [[ "${j}" == *"${sample}"* ]]
    then
        base2=`basename -s .fq ${j}`
        mv ${j} $(pwd)/${base2}.fastq
        int_files+=(${base2}.fastq)
    fi
done

echo 'run Optitype in Singularity'
singularity run -B /data/salomonis2/LabFiles/Frank-Li/OptiType/${sample}:/mnt my_software.sif -i ${int_files[0]} ${int_files[1]} --rna -v -o /mnt

echo 'remove the fastq files'
for m in $(pwd)/*.fastq
do
    if [[ "${m}" == *"${sample}"* ]]
    then
        rm ${m}
    fi
done

}


cat remain.txt | while read line; do echo "starting to process ${line} ......"; program ${line}; echo "finished processing ${line} ......"; done










