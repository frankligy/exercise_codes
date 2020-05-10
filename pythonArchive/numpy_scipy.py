#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 15:31:19 2020

@author: ligk2e
"""

import numpy as np
from scipy.sparse import isspmatrix,csr_matrix,csc_matrix,dok_matrix

# ndarray and matrix object, but behavior-wise, matrix doesn't seem to have any differences

a = np.array([[3,4,5],[5,6,7]])    # numpy.ndarray
a1 = np.matrix([[3,6,5],[9,6,7]])   # numpy.matrix

# illustrate reshape, 1D array, 2D row and column vector
d = np.array([1,2,3,4])
d.shape   # (4,) this is a column vector, 1-d array
d1 = d.reshape(1,-1)
d1.shape   # (1,4) this is a row vector, 2-d array
d0 = d.reshape(-1,1)
d0.shape  # (4,1)
d2 = d1.reshape(-1,1)
d2.shape   # (4,1) this is a column vector again, 2-d array
e = np.array([1])
e.shape    # 1-d array as well

#python adopt row-major as default when filling in, R is column major as default.
a = np.array([[2,4],[5,6],[2,7],[7,8]])
a.shape
a.reshape(2,4)

# in julia, we specify column vector and row vector, [1 2 3 4] means row vector, [1,2,3,4] means column vector
# [1 2 3 4; 5,6,7,8] means a matrix 2*4


# initialization

b = np.zeros([3,4],dtype=np.int64)
b1 = np.ones([3,4],dtype=np.int64)

# raise power, multiplication
a2  = np.multiply(a,a1)   # element-wise multiplication
a1_t = np.transpose(a1)   # transpose
a3 = np.matmul(a,a1_t)    # matrix multiplication
a3 = np.dot(a,a1_t)
a4 = np.power(a3,3)       # element-wise power

np.add(A,B,out=B)     # B value will change, B will store the output of this command
np.divide(A,2,out=A)
np.negative(A,out=A)

# some useful and frequently-used function:
x = np.random.random((3,3))   # [0.0,1.0)
np.diag(x)
np.diag(x,k=1)   # the oblique line in parallel with main diagonal but above the diagonal
np.diag(x,k=-1)  # the oblique line in parallel with main diagonal but below the diagonal
y0 = np.diag([1,1,1],k=1)

y1 = np.diag([1,1,1])
y2 = np.eye(3)    # these two are all identity matrix

a = np.nonzero([1,2,0,1])   # index of all non-zero

x1 = np.pad(x,pad_width=1,mode='constant',constant_values=0)

x = np.unravel_index(99,(6,7,8))   # the order of unravelling, first matrix (row-wise), then second matrix...

Z = np.tile( np.array([[0,1],[1,0]]), (4,4))  # repeat a pattern 4*4 times

a = np.random.random(5)
np.trunc(a)    # the fractional part will be discarded
a.astype(int)   # change to int element-wise


# convert a cartesian coordinates to polar coordinates
Z = np.random.random((10,2))
X,Y = Z[:,0], Z[:,1]
R = np.sqrt(X**2+Y**2)
T = np.arctan2(Y,X)
print(R)
print(T)

# meshgrid
nx, ny = (3, 2)
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
xv, yv = np.meshgrid(x, y)

# iterate all value and their index
Z = np.arange(9).reshape(3,3)
for index, value in np.ndenumerate(Z):
    print(index, value)
    
X = np.random.randn(100)  # sample length 100 vector from normal distribution
Y = np.random.randint(0,100,(1000,100))  # sample (1000*100) matrix from [0,100)
confint = np.percentile(X, [2.5, 97.5])   # 95% confidence interval

x = np.array([1, 2, 4, 7, 0])
np.diff(x)


# scipy could convert to .mat format, special function, linear algebra, integration, etc....

aSparseC = csc_matrix(a)     # convert ndarray to csc object
aSparseR = aSparseC.tocsr()   # convert to csr object
aSparseD = aSparseC.todok()   # convert to dok object, it is less expensive to operate compare to csc, csr object
a_ori = aSparseC.toarray()   # convert back
