#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 18:14:10 2020

@author: ligk2e
"""
import numpy as np

class ProteinFamily():
    
    def __init__(self,familyID,seq,ss):
        self.ID = familyID
        self.seq = seq
        self.ss = ss
    
    def ambientOneHotEncoding(self,length):   # well, assume sliding window length will be an odd number
        result = []
        dic = {'H':0,'C':1,'E':2}
        for i in range(0,len(self.seq)-length+1):
            windowSeq = self.seq[i:i+length]
            windowSS = self.ss[i:i+length]
            middle = (length-1)//2
            label = dic[windowSS[middle]]
            oneHot = ProteinFamily.oneHotEncoding(windowSeq)
            result.append((oneHot,label))   # [0,0,0.....1,0,0],'0'
        return result
        
        
        
    @staticmethod
    def oneHotEncoding(seq):
        '''
        A-R-N-D-C-Q-E-G-H-I-L-K-M-F-P-S-T-W-Y-V
        
        '''
        # Cysteine(C) will result in [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        result = []
        template = np.zeros(20)
        dic = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,
               'T':16,'W':17,'Y':18,'V':19}
        for letter in seq:
            template[dic[letter]] = 1   # assign the certain position to 1
            result.extend(template)     # extend it to result array
            template = np.zeros(20)     # reset the template to all zeros
        return result
            
        
        
        
    @staticmethod
    def parseFile(path):
        with open(path,'r') as file1:
            content = file1.readlines()
        familyID, seq, ss = [],[],[]
        for index,value in enumerate(content):
            if index % 4 == 0: familyID.append(value.lstrip('>').rstrip('\n'))  # >1cix
            elif index % 4 == 1: seq.append(value.rstrip('\n'))   # consensus protein sequence
            elif index % 4 == 2: ss.append(value.rstrip('\n'))   # consensus secondary structure notation
        familyInfo = [(familyID[i],seq[i],ss[i]) for i in range(len(familyID))]
        # filter items with 'X' in the sequence
        familyInfoNew = list(filter(lambda x: 'X' not in x[1],familyInfo))
        return familyInfoNew
                
            

    @staticmethod
    def kFoldSplit(lis,i):    # well, only for 5 fold
        import random
        random.Random(4).shuffle(lis)  # specify random state as 4
        training, testing = [],[]
        for index, value in enumerate(lis):
            if index % 5 == i: testing.append(value)
            else: training.append(value)
        return training,testing




class KNNmachine():

    
    def __init__(self,X1,Y1,testingData):
        self.X1 = X1
        self.Y1 = Y1
        self.testingData = testingData
    
    def constructModel(self,k,mode='distance'):  # mode = 'uniform' so don't account for distance
        from sklearn.neighbors import KNeighborsClassifier
        clf = KNeighborsClassifier(k,weights=mode,algorithm='kd_tree')
        clf.fit(self.X1,self.Y1)
        self.clf = clf
       
    
    def predict(self):  #[   [   ([],0)  ,  ([],1)   ],[],[]]
        testing = len(testingData)   # how many protein in testing set
        finalResult = []
        for eachProtein in testingData:   # [  (  [],0     ),(),()   ]
            window = len(eachProtein)   # how many windows in a protein
            stat = []
            for eachWindow in eachProtein:   # (   [],0  )
                temp = np.array(eachWindow[0]).reshape(1,-1)   # from column vector to row vector
                #print(temp.shape)
                prediction = list(self.clf.predict(temp))[0]
                #print(prediction,type(prediction))

                stat.append(1) if prediction == eachWindow[1] else stat.append(0)
            from functools import reduce
            percentage = reduce(lambda a,b:a+b,stat)/window   # 80% of the residue' secondary structure are correctly predicted
            #if percentage >= 0.6: print(eachProtein,stat)
            finalResult.append(round(percentage,2))
        self.prediction = finalResult   #[80%,45%,34%,98%...]
        from statistics import mean
        print(finalResult)
        print('This round yield {0} average accuracy'.format(mean(finalResult)))
        

    @staticmethod
    def decoder(oneHotEncodingProtein,stat):
        dic1= {0:'A',1:'R',2:'N',3:'D',4:'C',5:'Q',6:'E',7:'G',8:'H',9:'I',10:'L',11:'K',12:'M',13:'F',14:'P',15:'S',
                16:'T',17:'W',18:'Y',19:'V'}
        dic2= {0:'H',1:'C',2:'E'}
        protein,secondStructure = [],[]
        for eachWindow in oneHotEncodingProtein:
            ss = dic2[eachWindow[1]]
            secondStructure.append(ss)
            seqOri = eachWindow[0]
            seqDeco = ''
            for i in range(0,len(seq),20):
                aa = dic[seq[i,i+20].index(1.0)]
                seqDeco += aa
            protein.append(seqDeco)
        # protein: [RTYDY,TYDYD,...,FETGD]
        # secondStructure: [H,C,C,E...E]
        reconstruct = ''
        for each in protein:
            if protein.index(each) == len(protein) - 1: remain = each
            else: remain = each[0]
            reconstruct += remain
        print(reconstruct,secondStructure,stat)
            
        
            
            
                

        
            
        
        





if __name__ == '__main__':
    familyInfoNew = ProteinFamily.parseFile('/Users/ligk2e/Desktop/ssprotein/sec_stru_benchmark_sable135.txt')
    for i in range(5):
        trainingData,testingData = [],[]
        training,testing = ProteinFamily.kFoldSplit(familyInfoNew,i)
        for item in training:
            member = ProteinFamily(item[0],item[1],item[2])
            result = member.ambientOneHotEncoding(11)
            trainingData.extend(result)    #[([],0),([],1)...]
        for item in testing:
            member = ProteinFamily(item[0],item[1],item[2])
            result = member.ambientOneHotEncoding(11)
            testingData.append(result)    # for testing set, I preserve where each sliding window comes from
        X1,Y1= [],[]
        for j in trainingData:
            X1.append(j[0])
            Y1.append(j[1])
        X1,Y1 = np.array(X1),np.array(Y1)  # [[0,0,0,1,0,0,1,....],[0,1,0,1,0,0,0,...]], [0,1,2,1,2,1,...]
        #print(X1,Y1,type(X1),type(Y1),X1.shape,Y1.shape)
        machine = KNNmachine(X1,Y1,testingData)
        machine.constructModel(15)
        machine.predict()
        
            
    
    
    
    












    