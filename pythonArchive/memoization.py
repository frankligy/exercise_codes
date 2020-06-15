#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 18:19:33 2020

@author: ligk2e
"""


# Method 1: primitive function decorator

def memoize(func):
    cache = dict()

    def memoized_func(*args):
        if args in cache:   # I first time know that this will work.
            return cache[args]
        result = func(*args)
        cache[args] = result
        return result

    return memoized_func   # return a function

def fibonacci(n):
    if n == 0:
        return 0
    elif n == 1:
        return 1
    return fibonacci(n - 1) + fibonacci(n - 2)

memoized_fibonacci = memoize(fibonacci)   # memoized_fibonacci is a function

memoized_fibonacci(35)    # *args will capture all non-specified arguments in a list named args, so args=[35,]
memoized_fibonacci(35)    # second time, the speed up could manifest.


# Method 2: using built-in functools
import functools

@functools.lru_cache(maxsize=128)
def fibonacci(n):
    if n == 0:
        return 0
    elif n == 1:
        return 1
    return fibonacci(n - 1) + fibonacci(n - 2)

fibonacci(35)  # even faster than method 1 because we actually turn it into a dynamic programming problem.