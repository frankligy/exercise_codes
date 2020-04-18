######## shortcut
# using VS code as default editor for shell script, command+K+C for block comment, command+K+U can undo that
# In unix system: Ctrl+A goes back to the first letter, Ctrl+E goes to the last letter
# Cont'd: Ctrl+ D terminate the current command, Ctrl+L clear the screen

###### resources
# personal onedrive, linux_courseware/geeksforgeeks/


#### configuration
# change the theme: open ~/.bash_profile(get executed when you are logged in) or ~/.bash_rc(everytime you open a new terminmal window, you should 
# modify /.bash_rc on cluster, mac should modify /.bash_profile), add the following to said files.

#export PS1="\e[0;31m[\u@\h \W]\$ \e[m"​
export CLICOLOR=1​
export LSCOLORS=GxFxCxDxBxegedabagaced (each two pairs – text color and bc color; 11 pairs In total representing 11 types of files)​
alias l="ls -al"   no spaces are allowed

# above PS1 would result longish command line wrap around the line, alterntive configuration
export PS1='\[\033[1;34m\]\$\[\033[0m\] '​
export CLICOLOR=1​
export LSCOLORS=GxFxCxDxBxegedabagaced





# variable, the scope of a varibale is the shell in which the variable get defined
var="test"    # no space is allowd !!!!!
echo $var   
echo ${var}   # line3 and line4 are the same, in case of confusion, line4 is recommended, they all mean the value of this variable
echo $(echo $var)   # $() means execute the command in the parenthesis and will pass the internal output to the exteral command


# bash arithmetic operation
echo $((5+10))
x=1
y=2
ans=$((x+y))
ans=$[$a/$b]   # bash doesn't support float-point operation, awk does.


# array
my_array=(4 5 6)
echo ${my_array}  # only return my_array[0], which is 4
echo ${my_array[2]} # will return 6


# aforementioned variables are only alive in this shell
# solution to make them permanent is to use export 
# then it will be saved in .bash_profile through which it will load automatically everytime you are logged into
export PATH=/data/home/mjchen/app/package:$PATH   # this is also how to add path to environment 




# while loop, in-line format and formal format
cat list.txt | while read line; do echo $line; done
cat list.txt | while read line
do
echo $line
done

# for loop, in-line format and formal format is the same as above
for i in ./srafiles/*; do echo $i; done    
for i in {1..500}


# different way to traverse files
for i in ./srafiles/*           # will return the relative path of every file in this folder
for i in ./srafiles              # will return the file name of every file in this folder

# backend
nohup echo "go to the backend"    # will throw the command to backend, stout will be appending into nohup.out
ps -es  # check all the processes

# sub-shell 
cat list.txt | while read line; do echo $line &;done  # amperand will throw each line into a sub-shell

# pipe
ls -al | cut -d "\t" -f 1    # the output of ls, ps, etc will also be used as input for other command
ps -ef | cut -d "\t" -f 1
ls | xargs touch            # for command that is not able to directly accept the last output, xargs could teach them process the last output line by line

# awk could handle all simple operations / arithematic calculation in each line
awk 'BEGIN {OFS="\t";F=[:-"\t"];}
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

    else ($2 ~!/unix/){
        $4 = ""    # delete fourth column
        printf "var%s%s" $2 $3
    }
}' learnAWK.txt 


# sed could perform all kinds of substitution in each line
sed 's/unix/linux/g;s/1/2/g' file       # consecutive substitution


# interact with input
grep ${1} mRNA-ExonIDs.txt > mannual.txt    # ${1} menans first argument

echo -n "Enter matched position:"   # -n option means don't print the trailing newline \n
read pos        # it generate a variable pos, the value is equal to your input
echo $pos

awk -v awk_pos=$pos 'BEGIN{print "The matched one is"}    # pass the shell variable into awk namespace
{if (NR==awk_pos) print NR,$1,$2,$3,$4}
END{print "See aboved"}' mannual.txt

# if condition
if [ $another != 0 ]
then
    grep $another mRNA-ExonIDs.txt > mannual1.txt
    grep $another Hs_Ensembl_exon.txt > exonlist1.txt
else
    echo "OK,that's it"
fi


# conda
module load anaconda3   # only on cluster
proxy_on   # only on cluster, before creating new environment or install packages
conda create -n env python = 3   # = means fuzzy match, == means exact match
conda env list  # all the env
conda init --all bash  # only on cluster
conda/source activate env  # on cluster only use source
conda list -n env # show all the packages
conda list -n env scipy  # show specific package
conda remove -n env scipy # remove specific package
conda remove -n env --all # remove all env
conda deactivate # deactivate current env, no argument needed
proxy_on
pip install pyscenic -user  # -user only on cluster, install depends on the channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge   # these are settings for anaocnda, not specific environment
conda install -c bioconda scanpy
conda install -c  anaconda git

