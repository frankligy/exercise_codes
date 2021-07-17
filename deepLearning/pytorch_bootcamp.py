#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 10:35:34 2020

@author: ligk2e
"""




import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

input_1d = torch.tensor([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype = torch.float)
input_1d.size()

input_2d = torch.tensor([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]], dtype = torch.float)
input_2d.size()

input_2d_img = torch.tensor([[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]], [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]], [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]], dtype = torch.float)
input_2d_img.size()

cnn1d_4 = nn.Conv1d(in_channels=1, out_channels=5, kernel_size=3, stride=1)    # cnn1d_4 is an object
# it requires the input to be [batch_size, in_channel, length]
input_1d = torch.unsqueeze(input_1d, 0)
input_1d = torch.unsqueeze(input_1d, 0)   # now it is [1,1,10]
output_1d = cnn1d_4(input_1d)    # it is [1,5,8]



cnn2d_4 = nn.Conv2d(in_channels=3, out_channels=5, kernel_size=3, stride=1)
# it requires input to be [batch_size,in_channel,in_height,in_width]
input_2d_img = torch.unsqueeze(input_2d_img,0)   # [1,3,3,10]
output_2d_img = cnn2d_4(input_2d_img)    # it is [1,5,1,8]









rnn_input = torch.from_numpy(np.arange(20)).view(4,5).unsqueeze(2).float()  # [batch_size,seq_length,input_size]
rnn = nn.RNN(input_size=1,hidden_size = 2, batch_first=True, num_layers=3,bidirectional=True)
out,h_n = rnn(rnn_input)

# hidden size means number of output for each cell
# out: [4,5,4]   output of all time-steps from last layer [batch,seq_len,direction*hidden_size]
# h_n: [6,4,2]   last time-step output from all layers [layer*directions,batch,hidden_size]







torch.arange(4) + 1
torch.randn(4)



# seq2seq model, Encoder/Deconder, Attention
import torch
import torch.nn as nn
from torch import optim
import torch.nn.functional as F

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")   # device object


embedding = nn.Embedding(10, 3)




import torch
import torch.nn as nn
from torch import optim
import torch.nn.functional as F

class Encoder(nn.Module):
  def __init__(self, input_size, hidden_size, bidirectional = True):
    super(Encoder, self).__init__()
    self.hidden_size = hidden_size
    self.input_size = input_size
    self.bidirectional = bidirectional
    
    self.lstm = nn.LSTM(input_size, hidden_size, bidirectional = bidirectional)
  
  def forward(self, inputs, hidden):
    
    output, hidden = self.lstm(inputs.view(1, 1, self.input_size), hidden)  # ([seq_len,batch_size,input_size],(h0,c0))
    return output, hidden   # here hidden is a list of length 2, hidden[0], hidden[1] has the same dimension as before
    
  def init_hidden(self):
    return (torch.zeros(1 + int(self.bidirectional), 1, self.hidden_size),  # the same dimension as h_n, [layer*directions,batch,hidden_size]
      torch.zeros(1 + int(self.bidirectional), 1, self.hidden_size))



class AttentionDecoder(nn.Module):
  
  def __init__(self, hidden_size, output_size, vocab_size):
    super(AttentionDecoder, self).__init__()
    self.hidden_size = hidden_size
    self.output_size = output_size
    
    self.attn = nn.Linear(hidden_size + output_size, 1)
    self.lstm = nn.LSTM(hidden_size + vocab_size, output_size) #if we are using embedding hidden_size should be added with embedding of vocab size
    self.final = nn.Linear(output_size, vocab_size)
  
  def init_hidden(self):
    return (torch.zeros(1, 1, self.output_size),
      torch.zeros(1, 1, self.output_size))
  
  def forward(self, decoder_hidden, encoder_outputs, input):
    
    weights = []
    for i in range(len(encoder_outputs)):
      print(decoder_hidden[0][0].shape)
      print(encoder_outputs[0].shape)
      weights.append(self.attn(torch.cat((decoder_hidden[0][0], 
                                          encoder_outputs[i]), dim = 1)))
    normalized_weights = F.softmax(torch.cat(weights, 1), 1)
    
    attn_applied = torch.bmm(normalized_weights.unsqueeze(1),
                             encoder_outputs.view(1, -1, self.hidden_size))
    
    input_lstm = torch.cat((attn_applied[0], input[0]), dim = 1) #if we are using embedding, use embedding of input here instead
    
    output, hidden = self.lstm(input_lstm.unsqueeze(0), decoder_hidden)
    
    output = self.final(output[0])
    
    return output, hidden, normalized_weights


bidirectional = True
c = Encoder(10, 20, bidirectional)
a, b = c.forward(torch.randn(10), c.init_hidden())
print(a.shape)
print(b[0].shape)
print(b[1].shape)

x = AttentionDecoder(20 * (1 + bidirectional), 25, 30)
y, z, w = x.forward(x.init_hidden(), torch.cat((a,a)), torch.zeros(1,1, 30)) #Assuming <SOS> to be all zeros
print(y.shape)
print(z[0].shape)
print(z[1].shape)
print(w)
































