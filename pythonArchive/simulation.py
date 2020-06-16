#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 22:07:18 2020

@author: ligk2e
"""

'''
Problem setting:

Students are submitting jobs to printer, probability of having job submitted is 1/180 per secs.    

printer has two mode:
    low quality mode, pagePerSecond = 10
    high quality mode, pagePersecond = 5

Want to decide the mose time-saving mode as default mode such that the total waiting time of jobs sitting in the queue
can be minimized. 

We Simulate the whole process second by second.

'''


from Linear import Queue

import random

class Printer():
    def __init__(self,ppm):
        self.pagerate = ppm
        self.currentTask = None
        self.timeRemaining = 0
        
    def tick(self):
        if self.currentTask != None:
            self.timeRemaining = self.timeRemaining - 1
            if self.timeRemaining <= 0:
                self.currentTask = None
                
    def busy(self):
        if self.currentTask != None:
            return True
        else:
            return False
        
    def startNext(self,newtask):
        self.currentTask = newtask
        self.timeRemaining = newtask.getPages() + 60/self.pagerate
        
        
class Task:
    def __init__(self,time):
        self.timestamp = time
        self.pages = random.randrange(1,21) 
        
    def getStamp(self):
        return self.timestamp
    
    def getPages(self):
        return self.pages

    def waitTime(self,currenttime):
        return currenttime - self.timestamp
    
def newPrintTask():
    num = random.randrange(1,181)
    if num == 180:
        return True
    else:
        return False
    
def simulation(numSeconds,pagesPerMinute):
    
    labprinter = Printer(pagesPerMinute)
    printQueue = Queue()
    waitingtimes = []
    
    for currentSecond in range(numSeconds):
        if newPrintTask():
            task = Task(currentSecond)
            printQueue.enqueue(task)
            
    if (not labprinter.busy()) and (not printQueue.isEmpty()):
        nexttask = printQueue.dequeue()
        waitingtimes.append(nexttask.waittime(currentSecond))
        labprinter.startNext(nexttask)
        
    labprinter.tick()
    
    averageWait = sum(waitingtimes)/len(waitingtimes)
    print('Average Wait %6.2 secs %3d tasks remaining.' %(averageWait,printQueue.size()))  # %6.2 means output
    # will reserve 6 character in the console, rounding off to 2 decimal places
    # %3d means output will reserve 3 character in the console
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    