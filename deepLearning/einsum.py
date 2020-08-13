#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 22:02:53 2020

@author: ligk2e
"""

'''
rules:
1. if a dimension disappear in output, it means sum along this dimension




'''




import torch

x = torch.rand([2,3])

# permutaion, the general case for transpose
torch.einsum('ij -> ji',x)

# Summation
torch.einsum('ij -> ', x)  

# column sum
torch.einsum('ij -> j', x)

# row sum
torch.einsum('ij -> i', x)

# matrix-matrix
v = torch.rand([1,3])
torch.einsum('ij,kj -> ij',x,v)

# dot product first row with first row of matrix
torch.einsum('i,i -> ',x[0],x[0])

# hadamard product
torch.einsum('ij,ij -> ij', x,x)

# outer product
torch.einsum('i,j -> ij', x[0],x[0])

# batch matrix multiplication
a = torch.rand((3,2,5))
b = torch.rand((3,5,3))
torch.einsum('ijk,ikl -> ijl',a,b)

# matrix diagonal
x = torch.rand((3,3))
torch.einsum('ii -> i',x)

# matrix trace
torch.einsum('ii -> ', x)