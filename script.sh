# Another way to execute .sh file, bash -f script.sh
# tr -d '\n' < output.txt > temp.txt && mv temp.txt output.txt     how to use tr and mv and < > symbol
# tr usages:
# 1. lower case to upper case
# 2. translate white-space to tab
# 3. translate braces to parenthesis
# 4. squeeze continuous repetition of characters using -s
# 5. (most) delete specified charactoers using -d  
# 6. remove all digits/non-digits from string, negation can be achived using -c
# An easy version of sed which could fulfill a lot of daily task, refer to https://www.geeksforgeeks.org/tr-command-in-unix-linux-with-examples/


echo -n "Enter the EnsID;"
read EnsID
grep $EnsID mRNA-ExonIDs.txt > mannual.txt
grep $EnsID Hs_Ensembl_exon.txt > exonlist.txt
echo -n "Enter another EnsID:"
read another
if [ $another != 0 ]
then
    grep $another mRNA-ExonIDs.txt > mannual1.txt
    grep $another Hs_Ensembl_exon.txt > exonlist1.txt
else
    echo "OK,that's it"
fi
