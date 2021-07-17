#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copy from Aladdin, not my own coding!!

Just for use of learning and reference!

"""


"""
let's talk about character-level LSTM

1. identiy the tensor size, tensor.Size([]) = tensor(22), 0-D tensor, you can conver it to tensor([22]) = tensor.Size([1]) by view(1)
2. view(-1) is used for unrolling
3. idea: we have a corpus, let's see a very long string, you first get_random_sample (length of 250) by randomly pick a start_index, then end_index is determined,
    afte that, the tensor input would be this 250-long substring[:-1], tensor target would be [1:], each character's LSTM output is the predicted next character.
    Let's say batch=5, character_total =100, so the input would be tensor.Size([5,250]), each entry has been converted to a index, the nth character in 100 corpus.
    we feed each column into model one by one, for each column, it is tensor.Size([5]), embedding layer (100, 256), each entry will be a index, so it will return a 
    256 of length vector, result in tensor.Size([5,256]). before fed into LSTM, unsqueeze(1) to tensor.Size([5,1,256]). the hidden size is 256, but we connect a FC
    layer, so the output size become 100. Each character(column) will have a loss, we use crossEntropyLoss, so loss(y_pred,y), y_pred is tensor.Size([5,1,100]) 
    [$$$$$$$$$ I think it should be [5,100]$$$$$$$$$$$$$$$$], y is tensor.Size([5]), which is exactly the tensor target's corresponding column. We will accumulate 
    all individual character's loss together, backprop to finish the training. Remember every character (LSTM cell) use the (hidden,cell) from last cell.

    In the prediction phase, we have initial several character, we only iterate say we have 5 characters, we only iterate first4, we only care the (hidden,cell),
    the prediction starts with the last character, the output of the last character will be served as next input. we can define predict_Length, here, we will temparature
    and torch.multinomial to get the most possible character.

    All the time, there's only one LSTM cell. (hidden, cell) is contagious and get passed on.
"""

import torch
import torch.nn as nn
import string
import random
import sys
import unidecode

# Device configuration
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Get characters from string.printable
all_characters = string.printable
n_characters = len(all_characters)   # 100
file = unidecode.unidecode(open("/Users/ligk2e/Desktop/Machine-Learning-Collection/ML/Projects/text_generation_babynames/data/names.txt").read())


class RNN(nn.Module):
    def __init__(self, input_size, hidden_size, num_layers, output_size):
        super(RNN, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers

        self.embed = nn.Embedding(input_size, hidden_size)
        self.lstm = nn.LSTM(hidden_size, hidden_size, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_size, output_size)

    def forward(self, x, hidden, cell):
        out = self.embed(x)
        out, (hidden, cell) = self.lstm(out.unsqueeze(1), (hidden, cell))
        out = self.fc(out.reshape(out.shape[0], -1))
        return out, (hidden, cell)

    def init_hidden(self, batch_size):
        hidden = torch.zeros(self.num_layers, batch_size, self.hidden_size).to(device)
        cell = torch.zeros(self.num_layers, batch_size, self.hidden_size).to(device)
        return hidden, cell
    

class Generator:
    def __init__(self):
        self.chunk_len = 250
        self.num_epochs = 5000
        self.batch_size = 1
        self.print_every = 50
        self.hidden_size = 256
        self.num_layers = 2
        self.lr = 0.003

    def char_tensor(self, string):
        tensor = torch.zeros(len(string)).long()
        for c in range(len(string)):
            tensor[c] = all_characters.index(string[c])
        return tensor

    def get_random_batch(self):
        start_idx = random.randint(0, len(file) - self.chunk_len)
        end_idx = start_idx + self.chunk_len + 1
        text_str = file[start_idx:end_idx]
        text_input = torch.zeros(self.batch_size, self.chunk_len)
        text_target = torch.zeros(self.batch_size, self.chunk_len)

        for i in range(self.batch_size):
            text_input[i, :] = self.char_tensor(text_str[:-1])
            text_target[i, :] = self.char_tensor(text_str[1:])

        return text_input.long(), text_target.long()

    def generate(self, initial_str="A", predict_len=100, temperature=0.85):
        hidden, cell = self.rnn.init_hidden(batch_size=self.batch_size)
        initial_input = self.char_tensor(initial_str)
        predicted = initial_str

        for p in range(len(initial_str) - 1):
            _, (hidden, cell) = self.rnn(
                initial_input[p].view(1).to(device), hidden, cell
            )

        last_char = initial_input[-1]

        for p in range(predict_len):
            output, (hidden, cell) = self.rnn(
                last_char.view(1).to(device), hidden, cell
            )
            output_dist = output.data.view(-1).div(temperature).exp()
            top_char = torch.multinomial(output_dist, 1)[0]
            predicted_char = all_characters[top_char]
            predicted += predicted_char
            last_char = self.char_tensor(predicted_char)

        return predicted

    # input_size, hidden_size, num_layers, output_size
    def train(self):
        self.rnn = RNN(
            n_characters, self.hidden_size, self.num_layers, n_characters
        ).to(device)

        optimizer = torch.optim.Adam(self.rnn.parameters(), lr=self.lr)
        criterion = nn.CrossEntropyLoss()
        #writer = SummaryWriter(f"runs/names0")  # for tensorboard

        print("=> Starting training")

        for epoch in range(1, self.num_epochs + 1):
            inp, target = self.get_random_batch()
            hidden, cell = self.rnn.init_hidden(batch_size=self.batch_size)

            self.rnn.zero_grad()
            loss = 0
            inp = inp.to(device)
            target = target.to(device)

            for c in range(self.chunk_len):
                output, (hidden, cell) = self.rnn(inp[:, c], hidden, cell)
                loss += criterion(output, target[:, c])

            loss.backward()
            optimizer.step()
            loss = loss.item() / self.chunk_len

            if epoch % self.print_every == 0:
                print(f"Loss: {loss}")
                print(self.generate())

            #writer.add_scalar("Training loss", loss, global_step=epoch)
            
            
            
gennames = Generator()
gennames.train()           
            
            
            
            
            
            
            
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        