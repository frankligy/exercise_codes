# using VS code as default editor for shell script, command+K+C for block comment, command+K+U can undo that
# In unix system: Ctrl+A goes back to the first letter, Ctrl+E goes to the last letter
# Cont'd: Ctrl+ D terminate the current command, Ctrl+L clear the screen

# In order to run this file, ./check.sh

# If using VIM, type i to enter editing mode, then ESC to control mode, colon to input command, :q means exit
# :q! don't save and exit, :wq save and exit

# If using nano, directly into edting mode, using Ctrl(~) and Alt(M) key to operate. More info please refer to any tutorials

# change the theme: open ~/.bash_profile(get executed when you are logged in) or ~/.bash_rc(everytime you open a new terminmal window, you should 
# modify /.bash_rc on cluster, mac should modify /.bash_profile), add the following to said files.

#export PS1="\e[0;31m[\u@\h \W]\$ \e[m"​
#export CLICOLOR=1​
#export LSCOLORS=GxFxCxDxBxegedabagaced (each two pairs – text color and bc color; 11 pairs In total representing 11 types of files)​
#alias l="ls -al"   no spaces are allowed

# chmod 755 file   7:user itself, read(4)+write(2)+execute(1)   5: groups   5: groups

# echo ${1}
grep ${1} mRNA-ExonIDs.txt > mannual.txt
grep ${1} Hs_Ensembl_exon.txt > exonlist.txt


# echo "$(awk 'NR==2' mannual.txt)"
# EnsID=${1}
# echo $EnsID




echo -n "Enter matched position:"   # -n option means don't print the trailing newline \n
read pos

# awk -v awk_pos=$pos 'NR==awk_pos' mannual.txt
awk -v awk_pos=$pos 'BEGIN{print "The matched one is"}
{if (NR==awk_pos) print NR,$1,$2,$3,$4}
END{print "See aboved"}' mannual.txt

echo -n "Enter the training ENST digits for query:"
read trailing

awk -v awk_trailing=$trailing$ 'BEGIN{print "The matched one is"}
{if ($2 ~ awk_trailing) print NR,$1,$2,$3,$4}          
END{print "see aboved"}' mannual.txt

# pre-defining the pattern, then ~ pattern.