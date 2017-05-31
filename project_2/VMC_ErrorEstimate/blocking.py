#!/usr/bin/env python3

#Use first argument with Importance Sampling

import sys
import numpy as np
#import matplotlib.pyplot as plt

dataset = np.fromfile(sys.argv[1])


for i in range(2, len(sys.argv)):
    np.concatenate((dataset, np.fromfile(sys.argv[i])), axis=0)

def sumEverySecond(dataset):
   dataset = [sum(dataset[current: current+2]) for current in range(0, len(dataset), 2)]
   return np.divide(dataset,2)

def blockingAnalisys(dataset):
   std = []
   block = []
   i = 1
   sizeOfInitialDataset = len(dataset)
   sizeOfDataset =  sizeOfInitialDataset
   std.append(np.var(dataset)/sizeOfDataset)
   block.append(i)
   while (sizeOfDataset > 20):
      tempArray = sumEverySecond(dataset)
      dataset = tempArray
      sizeOfDataset = len(dataset)
      std.append(np.var(dataset)/sizeOfInitialDataset)
      i += 1
      block.append(i)
   return std, block

std_one,block_one = blockingAnalisys(dataset)

np.savetxt('blocking.csv', (block_one,std_one))

