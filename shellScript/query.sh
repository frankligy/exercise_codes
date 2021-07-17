#BSUB -W 5:00
#BSUB -M 80000
#BSUB -n 1
#BSUB -J HLAtype
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

cd /data/salomonis2/LabFiles/Frank-Li/OptiType



function program() {

echo ">$1" >> result.txt
find $1 | grep .tsv | xargs cat | awk 'NR==2 {print $2,$3,$4,$5,$6,$7}'

}


touch result.txt
cat already.txt | while read line; do program ${line}; done 

# export PATH=$PATH:/usr/bin/find

​# find $1 | grep .tsv | xargs cat | awk 'NR==2 {print $2,$3,$4,$5,$6,$7}' >> result.txt