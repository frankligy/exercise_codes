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
a4 = np.power(a3,3)       # element-wise power


# scipy could convert to .mat format, special function, linear algebra, integration, etc....

aSparseC = csc_matrix(a)     # convert ndarray to csc object
aSparseR = aSparseC.tocsr()   # convert to csr object
aSparseD = aSparseC.todok()   # convert to dok object, it is less expensive to operate compare to csc, csr object
a_ori = aSparseC.toarray()   # convert back
