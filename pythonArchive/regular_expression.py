#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 11:58:43 2020

@author: ligk2e
"""

# This script and comments are for understanding python3 re package and regex package, also brush up regular expression

# detailed regular expression please refer to two html I downloaded recently, cureently sit on Desktop:
# 1. Python 3 - Regular Expressions - Tutorialspoint
# 2. Python RegEx
# 3. A snapshot in my onedrive (learning linux from now)

# we are using raw string, r'', meaning, there is no escape character \, it is quite confusing. So in a standard string,
# \ is an escape character, a leading \ will make the following letter has special meaning, like \n, \t. So in raw string, 
# \ is just \. However, in regular expression, raw string '\' is used to escape. Try to understand this logic.

# Summary of regular expression:
# 1. ordinary pattern
# 2. metacharacter, use \ to escape if you wanna their liternal meaning
# 3. wild card special character
# 4. Anchor character
# 5. Quantify character
# 6. character set (metacharacters are not special in character set, caret ^ will have "negation" meaning in character set),[az-] here - means literal, [a-z] here - means range
# 7. Alternatives(TAA|TGG|TGA), [] character set can only handle alternatives in single leagth. Also if you wanna capture them (point9)
# 8. Grouping in conjuction with quantifier (\w+\b)*
# 9. Grouping in conjunction with back reference ()()  ->> \1 refer to the first paranthesis, \2 refer to the second, 
#    It can be used in re.sub and also even used in one pattern, (P|p)ython%\1erl will be either Python%Perl or python%perl
# 10. greedy and non-greedy: default of .* are greedy, .*? will be confined to non-greedy
# 11. (?=exp): match portion before exp; (?<=exp): match portion after exp
# 12. metacharacters: . ^ $ * + ? {} [] \ | ()

import re
import regex # regex has some additional function, I will use its fuzzy matching

text = 'I\'m dancing'
re.search(r'\b\w+(ing)',text).group(1)   # the match is dancing, return 'ing' since we specifying subgroup1
re.search(r'\b\w+(?=ing)',text)          # the match will be danc

# re:
# 1. re.match or re.search, match can only happen if at the beginning of the string. re.match().group(), groups() will return all the 
#    subgroups, not all matches, match function and search function can only return the first occurence. group(1) return first subgroup
#    span will return the (start,end) in a tuple. Matchobj and Searchobj cound be used as condition in IF command
# 2. re.findall() could return all matches in a list
# 3. re.finditer() could return all index information, it is a iterator
# 4. split, sub function
# 5. re.fullmatch()
# 6. re.complile()
# 7. in each method, there is an optional argument called flags, flags=re.I(ignore case), flags = re.M(multiple line matching)
# 8. re.sub()


# regex:

a = 'ATTTCRRR'
b = 'YYYYYYYYYYYYYYYTTTRRR'
c = 'TTTRRRP'

pattern = regex.compile('(%s){d<=2}' % c)    # if you delet two of a, will get a match, {i <= 4, s <= 4, d <= 4, e <= 5}
pattern.search(b)

# the result, fuzzy_counts = (substitution, deletion, insertion)