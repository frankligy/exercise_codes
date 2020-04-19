#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 15:31:19 2020

@author: ligk2e
"""

import numpy as np
from scipy.sparse import isspmatrix,csr_matrix,csc_matrix,dok_matrix

a = np.array([[3,4,5],[5,6,7]])    # numpy.ndarray
a1 = np.matrix([[3,6,5],[9,6,7]])   # numpy.matrix

# in julia, we specify column vector and row vector, [1 2 3 4] means row vector, [1,2,3,4] means column vector
# [1 2 3 4; 5,6,7,8] means a matrix 2*4

b = np.zeros([3,4],dtype=np.int64)
b1 = np.ones([3,4],dtype=np.int64)

a2  = np.multiply(a,a1)   # element-wise multiplication
a1_t = np.transpose(a1)   # transpose
a3 = np.matmul(a,a1_t)    # matrix multiplication
a4 = np.power(a3,3)       # element-wise power


# scipy could convert to .mat format, special function, linear algebra, integration, etc....

aSparseC = csc_matrix(a)     # convert ndarray to csc object
aSparseR = aSparseC.tocsr()   # convert to csr object
aSparseD = aSparseC.todok()   # convert to dok object, it is less expensive to operate compare to csc, csr object
a_ori = aSparseC.toarray()   # convert back
