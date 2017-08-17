#
# Python file to receive an .ac file and output a feed-forward .ac file
# to calculate marginal probabilities
# 
# By: Andrew Jemin Choi 
#

import sys
import math

# Node Class
class node:
    # Node type can be 'n' or 'v' for leaf nodes
    # Node type can be '+' or '*' for non-leaf
    def __init__(self, ntype, index, value, derivative):
        self.ntype = ntype
        self.index = index # denotes the n'th variable
        self.vr = value
        self.dr = derivative
        self.numChild = 0 # number of child nodes of non-leaf
        self.childList = [] # indices of children
        self.prL = [] # product cache register
        self.prR = []

    # Helper functions
    def update_vr(self, newVal):
        self.vr = newVal

    def update_dr(self, newDr):
        self.dr = newDr

    def add_child(self, child):
        self.childList.append(child)

    def update_numChild(self, num):
        self.numChild = num

    def update_prL(self, cache):
        self.prL = cache

    def update_prR(self, cache):
        self.prR = cache


# Check if file was passed by command line
FILE_NAME = ''
FILE_SIZE = 0
if (len(sys.argv) < 3):
    FILE_NAME = input("Please input an AC file: ")
    FILE_SIZE = input("Please input the AC size: ")
else:    
    FILE_NAME = sys.argv[1]
    FILE_SIZE = sys.argv[2]

print("File name:", FILE_NAME)
print("size:", FILE_SIZE)

# List to add nodes
circuit = []
index = 0

# Read the AC file
AC_FILE = open(FILE_NAME, 'r')
for line in AC_FILE:
    #print(line)
    if line[0] == '(':
        print("\t ... reading file ...")
    elif line[0] == 'E':
        print("\t ... done reading file ...")
        index = index - 1
        circuit[index].dr = 1
        break
    else:
        parsed = line.split()
        if parsed[0] == 'n':
            # index and dr equals 0
            val = float(parsed[1])
            circuit.append(node(parsed[0], 0, val, 0))
            index = index + 1
            
        elif parsed[0] == 'v':
            ind = int(parsed[1])
            val = float(parsed[2])
            circuit.append(node(parsed[0], ind, val, 0))
            index = index + 1

        elif parsed[0] == '+':
            value = 0
            addNode = node(parsed[0], 0, 0, 0)
            for child in parsed:
                if child == '+':
                    #nothing
                    value = 0
                else:
                    child = int(child)
                    addNode.add_child(child)
                    value = value + circuit[child].vr
            addNode.update_vr(value)
            circuit.append(addNode)
            index = index + 1

        elif parsed[0] == '*':
            value = 1.0
            numChild = 0
            multNode = node(parsed[0], 0, 0, 0)
            childValues = [1.0]
            prL = []
            prR = []
            for child in parsed:
                if child == '*':
                    value = 1.0
                else:
                    child = int(child)
                    multNode.add_child(child)
                    childValues.append(circuit[child].vr)
                    #value = value * circuit[child].vr
                    numChild = numChild + 1

            multNode.update_numChild(numChild)
            prL.append(1.0)
            prR.append(1.0)
            
            k = numChild
            for j in range(1, (numChild+1), 1):
                productL = childValues[j] * prL[j-1]
                prL.append(productL)
                productR = childValues[k] * prR[j-1]
                prR.append(productR)
                k = k - 1
            print("L", prL)
            print("R", prR)
            value = prL[numChild]
            
            multNode.update_prL(prL)
            multNode.update_prR(prR)
            multNode.update_vr(value)
            circuit.append(multNode)
            index = index + 1

    #print(index)
for i in circuit:
    #print(i.vr)
    break

outputVal = circuit[index].vr

print("output:", outputVal)
print("output log", math.log10(outputVal))

# Write to a test file

# Close AC file
AC_FILE.close()
