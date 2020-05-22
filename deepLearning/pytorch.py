#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 16:11:25 2020

@author: ligk2e
"""

# conda create -n pytorch python=3.6
# conda activate pytorch
# conda install pytorch torchvision -c pytorch



################ Chapter1: tensor object ############################
import torch


'''
float: single precision value, 32bit, up to 7 significant point
double: double precision value, 64bit, up to 14 significant point
short: 16bit
int: 32bit
long: 64bit

'''
 
x = torch.empty(5, 3)   # uninitialized matrix, temperary memory value will be assigned to it
x = torch.zeros(5, 3, dtype=torch.long)   # random initialized matrix
x = torch.randn(4, 4)
x = torch.tensor([2,3])

y = x.view(16)
z = x.view(-1, 8)  # the size -1 is inferred from other dimensions
print(x.size(), y.size(), z.size())


# convert between numpy.darray and torch.tensor, they are connected, change one will also change another one.
b = torch.from_numpy(x)   # convert ndarray x to torch.tensor b
b = x.numpy()              # convert tensor to x

device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')  # device object



################# Chapter2: Autograd() ###############################
x = torch.ones(2, 2, requires_grad=True)
y = x + 2
z = y * y * 3
out = z.mean()

out.backward()
x.grad

x = torch.randn(3, requires_grad=True)
y = x * 2

a = torch.randn(2, 2)
a.requires_grad_(True)

with torch.no_grad():
    print((x ** 2).requires_grad)   # requires_grad will void
    
y = x.detach()       # turn requires_grad from on to off
print(y.requires_grad)


'''
tensor object:
     .grad_fn   # how this tensor is generated from last tensor
     .require_grad # boolean
     .grad    # when loss.backward(),all net.parameters will has a grad
     .dtype
     .data
'''


################ Chapter3: Neural Network #############################
import torch
import torch.nn as nn
import torch.nn.functional as F


class Net(nn.Module):

    def __init__(self):
        super(Net, self).__init__()
        # 1 input image channel, 6 output channels, 3x3 square convolution
        # kernel
        self.conv1 = nn.Conv2d(1, 6, 3)
        self.conv2 = nn.Conv2d(6, 16, 3)
        # an affine operation: y = Wx + b
        self.fc1 = nn.Linear(16 * 6 * 6, 120)  # 6*6 from image dimension
        self.fc2 = nn.Linear(120, 84)
        self.fc3 = nn.Linear(84, 10)

    def forward(self, x):
        # Max pooling over a (2, 2) window
        x = F.max_pool2d(F.relu(self.conv1(x)), (2, 2))
        # If the size is a square you can only specify a single number
        x = F.max_pool2d(F.relu(self.conv2(x)), 2)
        x = x.view(-1, self.num_flat_features(x))
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x

    def num_flat_features(self, x):
        size = x.size()[1:]  # all dimensions except the batch dimension
        num_features = 1
        for s in size:
            num_features *= s
        return num_features


net = Net()
print(net)

params = list(net.parameters())  # net.parameters() is a generator object
print(len(params))
print(params[0].size())  # conv1's .weight

input = torch.randn(1, 1, 32, 32)  # one image greyscale for above-built network, (batch_size, channel, height, width)
output = net(input)                   # get the output when feed input to the net
net.zero_grad()                    # zero all gradient buffer
output.backward(torch.randn(1, 10))   # autograd

target = torch.randn(1,10)
target = target.view(1, -1)  # make it the same shape as output
criterion = nn.MSELoss()   # MSTLoss object
loss = criterion(output, target)   # loss is a tensor - scalar
loss.grad_fn.next_functions[0][0].next_functions[0][0]

net.zero_grad() 
print(net.conv1.bias.grad)
loss.backward()
print(net.conv1.bias.grad)

import torch.optim as optim

# create your optimizer
optimizer = optim.SGD(net.parameters(), lr=0.01)
optimizer.step()    # Does the update





















