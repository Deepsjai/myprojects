# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import numpy as np
from scipy.signal import convolve2d

def gridGame(Grid,k, Rules):
    for i in range (0,k):
        num_rows, num_cols = Grid.shape
       # Grid = [[0 for x in range(num_cols)] for y in range(num_rows)]
        Rule_indication = {}
        Rule_indication.update({idx:0 if val =='dead' else 1 for idx, val in enumerate(Rules)})
        Grid_Neighbors = convolve2d(Grid,np.ones((3,3),dtype=int),'same') - Grid
        Grid_Transformed = np.vectorize(Rule_indication.get)(Grid_Neighbors)
        Grid = Grid_Transformed
    return Grid

#g = [[1,0,1,0],[0,1,0,1]]

#a = ['dead', 'alive', 'dead', 'alive', 'dead']
    
n = int(input("Enter the number of rows in a matrix for Grid: "))
c = int(input("Enter the number of columns in a matrix for Grid : "))
g = [[int(input("type elements: ")) for x in range (c)] for y in range(n)]
g_n = np.array(g)
r = int(input("Please enter no. of rounds: "))
a = list()
num = input("no. of elements for Rules :")
print ('Enter dead or alive in array: ')
for i in range(int(num)):
    no = input("enter element: ")
    a.append(no)
result = gridGame(g_n,r, a)
print(result)