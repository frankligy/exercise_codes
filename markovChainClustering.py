#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 15:24:53 2020

@author: ligk2e
"""

'''
Disclaimer: I modified the code from https://github.com/GuyAllard/markov_clustering

'''

import numpy as np
from scipy.sparse import csr_matrix,csc_matrix,dok_matrix,isspmatrix



def constructAdjacency(edges,adjacency,dic):
    for edge in edges:
        row,column,weight = edge[0],edge[1],edge[2]
        rowInd,colInd = dic[row],dic[column]
        adjacency[rowInd,colInd] = weight
    return adjacency
    
    




def addSelfLoop(matrix,loopValue): # adjacency matrix has to be a square
    # change to dok_matrix for subsequent calculation because it is expensive to change the structure of sparse matrix
    matrixParsimony = matrix.todok()   # dok object, (0,0):1.0   (3,4):2.0
    shape = matrixParsimony.shape
    for i in range(shape[0]):
        matrixParsimony[i,i] = loopValue
    return matrixParsimony.tocsc()




def normalize(matrix):
    from sklearn.preprocessing import normalize
    matrixNew = normalize(matrix,norm='l1',axis=0)   # normalize by column using L1 normalization
    return matrixNew



def iteration(matrix,times,powerExpand,powerInflate):
    for i in range(times):   # time of iterations
        lastMatrix = matrix.copy()
        # expansion
        matrix = matrix ** powerExpand
        # inflation
        matrix = normalize(matrix.power(powerInflate))
        condition = matrixConverge(matrix,lastMatrix)    # True means has converged
        if condition: break
    return matrix




def matrixConverge(matrix1,matrix2,rtol=1e-5,atol=1e-8):
    c = np.abs(matrix1-matrix2) - rtol*np.abs(matrix2)
    return c.max() <= atol
    


def getCluster(matrix):
    attractors = matrix.diagonal().nonzero()[0]  # diagonal return an array, nonzero return an array with the index where its element is non-zero
    # above return value will be like array([3,4,7,8],), that's way we specify the first element[0]
    clusters = set()
    for attractor in attractors:  
        cluster = tuple(matrix.getrow(attractor).nonzero()[1].tolist())# getrow will return a 1*colnum csc object
                                                       # nonzero will return 2-D array, first for row, second for column, that's way we specify [1]
                                                       # set is a hashtable, key has to be immutable element, so convert list to tuple
        clusters.add(cluster)
    return list(clusters) 

                                                    




def visualization(matrix,clusters,inflation):
    import networkx as nx
    from matplotlib import cm
    import matplotlib.pyplot as plt
    graph = nx.Graph(matrix)   # Graph object
    dicNode = {node:cluster for cluster,nodes in enumerate(clusters) for node in nodes}
    nodeAssign = [dicNode[i] for i in range(len(graph.nodes()))]  # [0,1,0,0,1,0]
    print(nodeAssign)
    nx.draw_networkx(graph,node_color = nodeAssign, node_size=50, with_labels=True, edge_color="silver",cmap=cm.tab20)
    plt.title('inflation parameter = {}'.format(inflation))
    plt.savefig('clustering with inflation as {}.pdf'.format(inflation),bbox_inches = 'tight')
    plt.show()


def roundMatrix(matrix):     # 2d ndarray
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            matrix[i,j] = round(matrix[i,j],1)
    return matrix
    
    
def get_clusters(matrix):
    """
    Retrieve the clusters from the matrix
    
    :param matrix: The matrix produced by the MCL algorithm
    :returns: A list of tuples where each tuple represents a cluster and
              contains the indices of the nodes belonging to the cluster
    """
    if not isspmatrix(matrix):
        # cast to sparse so that we don't need to handle different 
        # matrix types
        matrix = csc_matrix(matrix)

    # get the attractors - non-zero elements of the matrix diagonal
    attractors = matrix.diagonal().nonzero()[0]

    # somewhere to put the clusters
    clusters = set()

    # the nodes in the same row as each attractor form a cluster
    for attractor in attractors:
        cluster = tuple(matrix.getrow(attractor).nonzero()[1].tolist())
        clusters.add(cluster)

    return sorted(list(clusters))



if __name__ == '__main__':
    edges = [('A','B',2),('A','C',2),('B','C',3),('B','D',1),('D','E',2),('D','G',4),('E','G',3),('E','F',2),('F','G',4),('E','H',1),('F','M',1),('H','M',2),('H','J',3),
         ('H','K',4),('M','J',4),('J','K',3),('J','S',1),('M','Q',1),('M','L',3),('K','L',4),('K','N',1),('Q','N',2),('Q','P',3),('N','P',3),('P','R',4)]

    adjacencyOri = np.zeros([17,17],dtype=np.float64)

    dic = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'J':8,'K':9,'L':10,'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16}

    adjacency = constructAdjacency(edges,adjacencyOri,dic)  # numpy.ndarray
    adjacency = csc_matrix(adjacency)    # csc object
    adjacencySelfLoop = addSelfLoop(adjacency,1)
    adjLoopNorm = normalize(adjacencySelfLoop)
    adjFinal = iteration(adjLoopNorm,10,2,2.1)
    
    
    
    
    
    adjArray = roundMatrix(adjFinal.toarray())
    clusters = get_clusters(adjFinal) 
    visualization(adjacency,clusters,1.1) 




    import markov_clustering as mc
    import networkx as nx
    import random
    result = mc.run_mcl(adjacency,inflation=2.1,iterations=10,expansion=2)           # run MCL with default parameters
    clusters = mc.get_clusters(result)    # get clusters








    