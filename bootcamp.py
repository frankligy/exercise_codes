#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 11:00:13 2020

@author: ligk2e
"""


#### copy module
# Assignmnent means let the left point to the right
# Let's discuss when they will interfere and when they are seemingly not interfering with each other
g = 5
h = g
h = 6   # let h points to a new object 6

a = [1,2,3,4]
b = a       #################### avoid directly assigning of mutable object  #############################
b[0] = 'test'  # list,set, dict and customized classes are mutable, allowing in-place change

import copy
a = [1,2,[3,4],4]
c = copy.copy(a)
c[0] = 'test'  # shallow copy copied the whole original object' reference?
c[2][0] = 'test'  # compound object

a = [1,2,[3,4],4]
d = copy.deepcopy(a)  # deep copy copied the whole original object
d[0] = 'test'
d[2][0] = 'test'


#### itertools module
'''
0. Class is a code template for creating objects, objects' type is the same as the name of the Class

1. iterables: list,tuple,dict,set,string,range object, generator object, enumerate object, iterator object ...
    They can be iterated in for loop, and last four can only be visible using for loop and print, or print(next())

2. itertools provide a plethora of methods related to iterables

    2.1 infinite iterators:(remember to use break in you if condition)
        2.1.1 itertools.count(5,5)  return count object: 5,10,15,20,25.......   (range function is sufficent?)
        2.1.2 itertools.cycle('ABC') return cycle object: 'A','B','C','A','B','C'....  
        2.1.3 itertools.repeat(4,5) return repeat object: 4,4,4,4,4
    2.2 Combinatoric iterators:
        2.2.1 itertools.product('AB','CD') return product object: 'AC','AD','BC','BD', cartesian product, here product might be nested inner class under itertoos class
        2.2.2 itertools.permutations()
        2.2.3 itertools.Combinations()
        2.2.4 itertools.Combinations_with_replacement()
    2.3 Terminating iterators:
        2.3.1 accumulate(iter, func)
        2.3.2 chain(iter1, iter2..)
        2.3.3 chain.from_iterable()
        2.3.4 compress(iter, selector)
        2.3.5 dropwhile(func, seq)
        2.3.6 filterfalse(func, seq)
        2.3.7 islice(iterable, start, stop, step)
        2.3.8 starmap(func., tuple list)
        2.3.9 takewhile(func, iterable)
        2.3.10 tee(iterator, count)
        2.3.11 zip_longest( iterable1, iterable2, fillval)
    2.4 Itertools.groupby()
                
'''
#range object:
e = range(10)  # e is a range object

# enumerate object:
strings = ['import','is','with','if','file','exception']
D = {key: index for index,key in enumerate(strings)} 
t = enumerate(strings)   # t is a enumerate object

#generator object
def fibonacci(max):
    n, a, b = 0, 0, 1
    while n < max:
        yield b  # same as print, but lead the function to a special function
        a, b = b, a + b  # assigned simutaneously, using original a and b value
        n = n + 1
    return 'done'
f = fibonacci(6)   # f is a generator object
print(next(f))
    
# list_itertor object
a = [1,2,3]
a1 = iter(a)



#### Collections Module
'''
four additional containers: counter, DefaultDict, OrderedDict,NamedTuple, they are all subclass of dict, NamedTuple is a function 
'''
import collections
from collections import Counter
c=Counter('Hello')
c.update('bfg')
c['f']  # access certain element
c.elements() # return iterators for each element
c.most_common(2) # n most frequent elements
# arithmatic: c1+c2, c1&c2 take the intersection

'''
collection.deque is preferred when you need quicker pop and append over list object, because list is a array but deque is 
a linked table. In same sense, set is a hash table, it is faster to query one element, time complexity is only O(1)

'''










