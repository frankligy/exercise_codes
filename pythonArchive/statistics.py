#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 13:13:36 2020

@author: ligk2e
"""

# example1: 

# calculate mean, standard deviation, standard error of mean, confidence Interval, errorbar
# for CI: a list of data, we assume they are drawn from a normal distribution, in that sense, this sample
# will follow t-distribution with df=len(list)-1, loc = mean, scale = sem(a).

def confidenceInterval(lis):
    import numpy as np, scipy.stats as st
    a = np.array(lis)
    me = np.mean(a)
    confInt = st.t.interval(0.95, len(a)-1, loc=me, scale=st.sem(a))  # will return a tuple (lower,upper)
    errorBar1 = np.std(a)   # using standard deviation as errorbar, yerr=errorBar1, lower error = upper error = errorBar1
    errorBar2 = st.sem(a)   # using standard error of mean as errorbar, same as above
    errorBar3 = [me-confInt[0],confInt[1]-me]  # using confidence interval as errorbar, yerr=errorBar3, lower error = errorBar3[0], upper error = errorBar3[1]
    return me, confInt, errorBar3   # these three will be combined as tuple automatically, if function return multiple values



# example2: 
    
# calculate frequently-used metrics of a distribution
# scipy.stats.rv_continous is a object, norm, etc are inherited from that and inherit all the method
# below is a normal distribution with loc=10,scale=3

from scipy.stats import norm

norm.pdf(60,10,3)   # the point distribution when x=60
norm.ppf(0.99,10,3)  # the quantile when cdf(quantile)=0.99, meaning, P(x<quantile)=0.99
norm.cdf(12,10,3)    # p(x<12)
norm.rvs(size=100,loc=10,scale=3,random_state=1)  # sample 100 points from this distribution
norm.sf(12,10,3)   # complementary of cdf, p(x>12)
norm.isf(0.99,10,3)  # complentatary of ppf, p(x>quantile)=0.99


# example3:

# normarlize the data, self-reasoning if input is just a list

def normalize(matrix):
    from sklearn.preprocessing import normalize
    matrixNew = normalize(matrix,norm='l1',axis=0)   # normalize by column using L1 normalization
    return matrixNew