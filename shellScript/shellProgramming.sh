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


# shebang to specify the which shell to use
cat /etc/shells  # will output all available shells

# if one line is too long, use \ to newline, there is a space before \

# variable, the scope of a varibale is the shell in which the variable get defined
# the data type: file, line, numeric, string, array
var="test"    # no space is allowd !!!!!
echo $var   
echo ${var}   # line3 and line4 are the same, in case of confusion, line4 is recommended, they all mean the value of this variable
echo "$var"   # this is a better practice to use double quotation
echo $(echo $var)   # $() means execute the command in the parenthesis and will pass the internal output to the exteral command
echo `echo $var`    # equivalent to above command

# in-line arguments, 
$1, $2 .... ${10}
$*  # all the arguments were viewed as a whole, for i in "$*" , i will be the whole
$@ # all the arguments were viewed as individual, for i in "$@", i will be each argument
$#   # the number of arguments

# special character
$$ # current PID
$! # last process's PID in backend
$? # the result of last process, will be 0 if executed correctly, will be non-zero if not

# all the varibales
set  # display all the variables
unset var # delete variable

# bash arithmetic operation
echo $((5+10))
x=1
y=2
ans=$((x+y))
ans=$[$a/$b]   # bash doesn't support float-point operation, awk does.
temp=`expr 2 \* 3`   # \* means multiply, / means divide, % means take remainder


# array
my_array=(4 5 6)
echo ${my_array}  # only return the my_array[0], just 4
echo ${my_array[2]} # will return 6
my_array[2]=7   # change the value
my_array[3]=8   # add new value
unset my_array[2]  # then the array will be (4 5 NULL 8)

# string
x="pattern"
echo ${#x}  # return the length of string
echo ${x:2:5}  # slicing the string


# aforementioned variables are only alive in this shell
# solution to make them permanent is to use export 
# then it will be saved in .bash_profile through which it will load automatically everytime you are logged into
export PATH=/data/home/mjchen/app/package:$PATH   # this is also how to add path to environment 
source /etc/profile   # activate the configuration file that we just made change to




# while loop, in-line format and formal format
# example1
cat list.txt | while read line; do echo $line; done
cat list.txt | while read line
do
echo $line
done

# example2
SUM=0
i=0
while [ $i -le $1]
do
	SUM=$[$SUM+$i]
	i=$[$i+1]
done
echo $SUM

# example3: until
ans="yes"
until [ $ans != "yes"]; do echo "Hi"; done





# for loop, in-line format and formal format is the same as above
# example 1
for i in ./srafiles/*; do echo $i; done    
for i in {1..500}

# example 2
SUM=0
for ((i=1;i<=100;i++))
do 
	SUM=$[$SUM+$i]
done
echo $SUM

# example 3
for i in {5..15..3}; do echo $i; done    # seems doesn't work


# different way to traverse files
for i in ./srafiles/*           # will return the relative path of every file in this folder
for i in ./srafiles              # will return the file name of every file in this folder

# backend
nohup echo "go to the backend"    # will throw the command to backend, stout will be appending into nohup.out
ps -es  # check all the processes
# then using PID
cd /proc/PID/fd 
ls -l    # will display the process's stdin, stdout, stderr, and their designated destination.

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

# if condition   :    a space is required after [ and before ]
# example 1:
if [ $another != 0 ]
then
    grep $another mRNA-ExonIDs.txt > mannual1.txt
    grep $another Hs_Ensembl_exon.txt > exonlist1.txt
elif []
then
    echo "OK,that's it"
else
	echo "Hi"
fi

# example2:
if [ 5 -gt 10 -o 10 -gt 4 ]   # -ge, -lt, -le      # -o means or, -a means and, ! means negation
then
	echo "Hi"
fi

# example3
# [ -f file ] means if it is a file ; [ -d file ] means if the file is a folder; [ -s file ] will be true when file exist and it is not empty; [! -s file]
if [ -e /root/shell/a.txt ] || [10 -gt 5]   # if exist this file   # && means and, || means or
then
	echo "Hi"
fi

# example 4
case $1 in
	"1")
	echo "Monday"
	;;
	"2")
	echo "Tuesday"
	;;
esac

# example 5
if [ -z "$str1" ]     # whether it is empty string
if [ -n "$str2" ]     # whether it is not empty string
if [ "$str1" = "$str2" ]
if [ "$str1" != "$str2" ]
if [[ "$str1" =~ "str2" ]]      


# function

# built-in
basename include/stdio.h .h   # will return stdio, since here we specify the suffix, which will be trimmed off
dirname /usr/bin/    # will delete all content after last /, include / itself, if no /, then return '.' current path

# customized 
# in shell, every variable will be viewed as global variable unless you use local to declare in your function scope
function getSum(){
		local var="hi"
		SUM=$[$n1+$n2]
		echo "sum=$SUM"
		return 0     # return just used to denote whether this function is successfully executed or not
}

read -p "please type in first parameter n1: " n1
read -p "please type in second parameter n2 " n2

getSum $n1 $n2

# log file
./test.sh > log.dat 2>&1    # 0 means stdin, 1 means stdout, 2 means stderr, it means direct stderr to stdout, which is log.dat as well.


# parse the arguments
usage="Usage: $0 -S sample_name -I input_dir -O output_dir"
while getopts S:I:O: flag; do
    case $flag in
    S)
    sample=$OPTARG
    ;;
    I)
    input_dir=$OPTARG
    ;;
    O)
    output_dir=$OPTARG
    ;;
    \?)
    echo $usage >& 2    # redirect to stderr
    exit -1
    ;;
    esac
done
shift $(( OPTIND - 1));

# calculate running time
master_time_start=`date +%s`
# your program
master_time_end=`date +%s`
(master_time_exec=`expr $(( $master_time_end - $master_time_start ))`


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

