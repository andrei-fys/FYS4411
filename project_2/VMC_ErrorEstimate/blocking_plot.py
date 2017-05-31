#!/usr/bin/env python3

#Use first argument with Importance Sampling

import sys
import numpy as np
import matplotlib.pyplot as plt

input_one = sys.argv[1]
input_two = sys.argv[2]

dataset_one = np.loadtxt(input_one)
dataset_two = np.loadtxt(input_two)

def sumEverySecond(dataset):
   dataset = [sum(dataset[current: current+2]) for current in range(0, len(dataset), 2)]
   return np.divide(dataset,2)

def blockingAnalisys(dataset):
   std = []
   block = []
   i = 1
   sizeOfDataset = len(dataset)
   std.append(np.var(dataset)/len(dataset))
   block.append(i)
   while (sizeOfDataset > 256):
      tempArray = sumEverySecond(dataset)
      dataset = tempArray
      sizeOfDataset = len(dataset)
      std.append(np.var(dataset)/len(dataset))
      i += 1
      block.append(i)
   return std, block

std_one,block_one = blockingAnalisys(dataset_one)
std_two,block_two = blockingAnalisys(dataset_two)

#print(np.var(dataset_one)/len(dataset_one))
#print(np.var(dataset_two)/len(dataset_two))

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(block_one, std_one, color="black",label=r'$Imp S$')
ax1.plot(block_two, std_two, color="red",label=r'$St$')
plt.grid()

plt.legend(loc="lower right", fontsize=18)
#
plt.xlabel('block', fontsize=20)
plt.ylabel(r'$\sigma^2$', fontsize=20)
#
plt.draw()
plt.show()
##plt.savefig(filename3)
