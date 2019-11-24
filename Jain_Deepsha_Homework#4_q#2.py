#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:51:37 2019

@author: deepsha
"""

import pandas as pd
import numpy as np
import heapq 

def findIndex(array, Rank):
    #x2 = np.random.randint(10, size=(5, 4)) 
    Sum_Array =np.sum(array, axis = 1)
    l = len(Sum_Array)
    indexes = np.arange(l)
    find_i = heapq.nlargest(Rank, indexes, key=lambda i: Sum_Array[i])
    find_index = find_i[-1]
    return find_index

n = int(input("Enter the number of rows in a matrix: "))
c = int(input("Enter the number of columns in a matrix: "))
array = [[int(input("type elements: ")) for x in range (c)] for y in range(n)]
#array = [[2,2],[21,11],[5,5],[6,7]]
Rank = int(input("please put rank: "))
result = findIndex(array, Rank)
print(result)