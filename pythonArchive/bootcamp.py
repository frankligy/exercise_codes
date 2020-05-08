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

# python, function call by reference, so original copy will be changed 
# but, for immutable type, it doesn't matter, it won't be changed but raise an error
# for mutable type, it will be changed, but we almost never really need to change the original list,
# we mostly read the original list and construct new stuff. 
# when we need to made change on original list, it is dangerous


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

# advanced usage
items1 = list(map(lambda x: x ** 2, filter(lambda x: x % 2, range(1, 10))))
colors = list(map(lambda x:x%2,[random.randint(0,200) for i in range(200)]))
occurence = [k for k in range(len(hlaAllele)) if hlaAllele[i] == hlaQuery]  # find all occurence of a item in a list

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


#### file IO
# modes: r -- read; w -- write(file no need to exist); a-- append(file no need to exist); b -- binary file 
f = open('file3.txt','r',encoding = 'utf-8')  # f will be a TextIOWrapper type
content = f.read()   # it will be a string with newline in between   (read all file in one go)
f.close()

with open('file3.txt','r',encoding = 'utf-8') as f1:   # f is still TextIOWrapper type
    line1 = f1.readline()     # only read one line with newline   
    line2 = f1.readline()
    
with open('file3.txt','r',encoding = 'utf-8') as f2:    
    lines = f2.readlines()    # all the lines will be stored in a list, each item with trailing newline

with open('file3.txt','r',encoding = 'utf-8') as f3:
    for line in f3: print(line)
    
with open('file1.txt','r') as f4, open('file2.txt','w') as f5:
    
    
    
# GUI programming ****** Don't run the following code, please   **********
import tkinter   # come with python base installer, no need to pip again
                 # other professional packages: wxPython、PyQt、PyGTK


def main():
    flag = True

    # change the word in label
    def change_label_text():
        nonlocal flag   # indicate this variable's scope is belonging to an outer layer but not the global 
        flag = not flag
        color, msg = ('red', 'Hello, world!')\
            if flag else ('blue', 'Goodbye, world!')
        label.config(text=msg, fg=color)

    # confirm exit
    def confirm_to_quit():
        if tkinter.messagebox.askokcancel('prompt', 'sure to quit'):
            top.quit()

    # generate top bar
    top = tkinter.Tk()
    # configure the size of the top bar
    top.geometry('240x160')
    # configure the title of top bar
    top.title('example')
    # the whole body
    label = tkinter.Label(top, text='Hello, world!', font='Arial -32', fg='red')
    label.pack(expand=1)
    # a container for buttons
    panel = tkinter.Frame(top)
    # generate button object, assign to a container, bind with functions using command argument
    button1 = tkinter.Button(panel, text='modify', command=change_label_text)
    button1.pack(side='left')
    button2 = tkinter.Button(panel, text='exit', command=confirm_to_quit)
    button2.pack(side='right')
    panel.pack(side='bottom')
    # main loop
    tkinter.mainloop()


if __name__ == '__main__':
    main()


#### game development: pygame, panda3D
# below snippet is not a complete code, it is representative of the mindset, how to capture the event, how to communicate 
while running:
        # obtain event from event queue and process that 
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        # obtain mousebuttondown event
        if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
            # obtain mouse position
            x, y = event.pos
            radius = randint(10, 100)
            sx, sy = randint(-10, 10), randint(-10, 10)
            color = Color.random_color()
            # instantiate a ball on the mouse position
            ball = Ball(x, y, radius, sx, sy, color)
  

#### List, tuple, set
# List
a = [1,2,3,4,5,6]
a.insert(1,8)  # insert 8, 8 will be the 1st element in the new list, which means insert 8 before second element
a.pop(3)   # pop the third element
a.remove(6)   # remove element 6
a.clear()    # become a empty list

# Tuple: useful in multi-thread programming since it is immutable, secure

# set: 
a = {1,2,3,4}
a.add(5)   # add 5 to the set
b = {3,4,5,6,7}
a.update(b)   # add b to a and remove overlapped elements

#### OOP
#example1: leading double underscore
class A:
    __tt = 'd'    # when create the classA, name mangling will happen

    def __init__(self, foo):
        self.__foo = foo   # after instantiate the object, name mangling will happen

    def __bar(self):    # when create the class A, name mangling will happen
        print(self.__foo)
        print('__bar')


test = A('hi')
test.__foo     # AttributeError
test._A__foo   # for __foo, __bar, if they are attribute or method in class, interpreter will do name mangling, will change
               # __foo to _A__foo, __bar to _A__bar.

#example2: @property
class Person(object):

    __slots__ = ('_name', '_age', '_gender')   # this class can only have these three attributes, can not add other attributes
    
    
    def __init__(self, name, age):
        self._name = name
        self._age = age


    @property      # it means, every time when we are trying to access the value of object.name, we are not gonna look up in __dict__, instead, we will trigger the below function
    def name(self):
        return self._name

    @property
    def age(self):
        return self._age

    @age.setter     # it means, every time when we are trying to assign a new velue to object.age, we are going to trigger the following function 
    def age(self, age):
        self._age = age

    @staticmethod    # this method is belong to the class, not any instance of this class
    def is_valid(a, b, c):
        return a + b > c and b + c > a and a + c > b

    @classmethod
    def now(cls):
        ctime = localtime(time())
        return cls(ctime.tm_hour, ctime.tm_min, ctime.tm_sec)
    
# example3: abstract class, abstract can not be instantiated, it can only be inherited
from abc import ABCMeta, abstractmethod


class Pet(object, metaclass=ABCMeta):


    def __init__(self, nickname):
        self._nickname = nickname

    @abstractmethod
    def make_voice(self):
        pass


class Dog(Pet):

    def make_voice(self):
        print('%s: wong...' % self._nickname)


class Cat(Pet):


    def make_voice(self):
        print('%s: miao...' % self._nickname)


# example4: how to add magic method to a class
# meaning, we all know that a class'method could be implemented by A.method1(), A.method2()
# then, is it possible to define the method like slicing A[4,5], with A as f:, those non-standard behavior are handley in class definitin
# for instance, wanna accomplist slicing, you should define __getitem__(self,key) method, it will be called when you try to do A[3,4]
# a full list could be found in my onedrive: learning python from now
        
# example5: singleton, __new__ method
class Singleton(object):
    instance = None

    def __new__(cls,age,name):
        if not cls.instance:
            cls.instance = object.__new__(cls)
            return cls.instance

a = Singleton(18,'John')
b = Singleton(48,'Tom')
        

#### process and thread, asyschronous I/O
# example1: multi-process in python
from multiprocessing import Process
from os import getpid
from random import randint
from time import time, sleep


def download_task(filename):
    print('Start downloading, pid is {0}.'.format(getpid()))
    print('Start downloading {}...'.format(filename))
    time_to_download = randint(5, 10)
    sleep(time_to_download)
    print('Finished the {0} downloading, consume {1} second'.format(filename, time_to_download))


def main():
    start = time()
    p1 = Process(target=download_task, args=('file1.pdf',))
    p1.start()
    p2 = Process(target=download_task, args=('file2.pdf',))
    p2.start()
    p1.join()   # wait until p1 is finished
    p2.join()
    end = time()
    print('Consume {} seconds in total'.format(end - start))


if __name__ == '__main__':
    main()


# example2: multi-thread in python
from random import randint
from threading import Thread
from time import time, sleep


def download(filename):
    print('Start to download {}'.format(filename))
    time_to_download = randint(5, 10)
    sleep(time_to_download)
    print('{} has been downloaded, consume {} seconds'.format(filename, time_to_download))


def main():
    start = time()
    t1 = Thread(target=download, args=('file1.pdf',))
    t1.start()
    t2 = Thread(target=download, args=('file2.pdf',))
    t2.start()
    t1.join()
    t2.join()
    end = time()
    print('Consume {} seconds in total'.format(end - start))


if __name__ == '__main__':
    main()

# example3: multi-thread with lock
from time import sleep
from threading import Thread, Lock


class Account(object):

    def __init__(self):
        self._balance = 0
        self._lock = Lock()  # there is gonna be a lock

    def deposit(self, money):
        self._lock.acquire()   # if the resource is available, current thread will proceed, acquire the key and re-lock the resource, so other thread can not access at this moment
        try:
            new_balance = self._balance + money
            sleep(0.01)
            self._balance = new_balance
        finally:
            self._lock.release()   # open the lock, so the resource is avilable now, next thread could enter

    @property
    def balance(self):
        return self._balance


class AddMoneyThread(Thread):

    def __init__(self, account, money):
        super().__init__()
        self._account = account
        self._money = money

    def run(self):
        self._account.deposit(self._money)


def main():
    account = Account()
    threads = []
    for _ in range(100):
        t = AddMoneyThread(account, 1)
        threads.append(t)
        t.start()
    for t in threads:
        t.join()
    print('balance is: $%d' % account.balance)


if __name__ == '__main__':
    main()

#### Heap sort, stack and queue
# heap sort: an efficient sorting algorithm using heap data strcture
import heapq

list1 = [34, 25, 12, 99, 87, 63, 58, 78, 88, 92]

print(heapq.nlargest(3, list1))
print(heapq.nsmallest(3, list1))


# Stack and queue
# they are abstarct data structure, describing a concept
# Stack: Last-in-First-Out
# Queue: Fisrt-in-First-Out

# A simple class stack that only allows pop and push operations
class Stack:

    def __init__(self):
        self.stack = []

    def pop(self):
        if len(self.stack) < 1:
            return None
        return self.stack.pop()

    def push(self, item):
        self.stack.append(item)

    def size(self):
        return len(self.stack)

# And a queue that only has enqueue and dequeue operations
class Queue:

    def __init__(self):
        self.queue = []

    def enqueue(self, item):
        self.queue.append(item)

    def dequeue(self):
        if len(self.queue) < 1:
            return None
        return self.queue.pop(0)

    def size(self):
        return len(self.queue) 


#### sorting
#example1: operate on list, sorted has return value, sort doesn't have
s = ['ab','abc','a','djkj']
b = sorted(s,key=lambda x: len(x),reverse = True)
s.sort(key=len)


# example2: sort based on key of a dict
dic = {'name':'zs','sex':'man','city':'bj'}
foo = zip(dic.keys(),dic.values())   # zip object, iterator, each is a tuple ('name','zs'), ('sex','man')
foo = [i for i in foo]
b = sorted(foo,key=lambda x: x[0])
new_dic = {i[0]:i[1] for i in b}


# example3: lambda function with condition
foo = [-5,8,0,4,9,-4,-20,-2,8,2,-4]
a = sorted(foo,key=lambda x: (x<0,abs(x)))   # first sort by boolean value, then abs value for negative number
                                             # lambda function return a function object, this example, it returns(True,value)


# example4: max function also has key argument, and see how max could return the key which has largest value
count = {'a':100,'b':3000,'d':3000,'e':67}
maxKey = max(count,key=lambda x:count[x])
maxTogether = max(count.items(),key=lambda x: count[x])



#### *args **kwargs
def myFun(*args):  
    for arg in args:  
        print (arg) 
    
myFun('Hello', 'Welcome', 'to', 'GeeksforGeeks') 


def myFun(**kwargs):  
    for key, value in kwargs.items(): 
        print ("{0} == {1}".format(key,value)) 
  
myFun(first ='Geeks', mid ='for', last='Geeks') 


#### subprocess module
from subprocess import run,Popen,PIPE

# subprocess.run is a high-level wrapper and also the latest one, compared to subprecess.call etc
# return a CompleteProcess object, text means the stdout, stderr will be string instead bytes b''
# capture_output can not be used with stdout, stderr, only one way
c = subprocess.run([path,'-p',intFile,'-BA','-a','HLA-A01:01,HLA-A03:01,HLA-B07:02'],capture_output=True,text=True)
# stdout is a file
with open('log2.txt','w') as f2:   # give stdout a TextIOWrapper object
    subprocess.run([path,'-p',intFile,'-BA','-a','HLA-A01:01,HLA-A03:01,HLA-B07:02'],stdout=f2)
    
# low-level, conferring a lot of flexibility,return a Popen object
p = subprocess.Popen([path,'-p',intFile,'-BA','-a','HLA-A01:01,HLA-A03:01,HLA-B07:02'],stdout=subprocess.PIPE,text=True)
(stdout,stderr) = p.communicate

# stdout is a file
with open('log2.txt','w') as f2:   # give stdout a TextIOWrapper object
    subprocess.Popen([path,'-p',intFile,'-BA','-a','HLA-A01:01,HLA-A03:01,HLA-B07:02'],stdout=f2)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    