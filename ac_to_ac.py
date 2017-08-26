#
# Python file to receive an .ac file and output a feed-forward .ac file
# to calculate marginal probabilities
# 
# By: Andrew Jemin Choi 
#

import sys
import math

# Global Variables
# List to add nodes
circuit = []
index = 0

# Node Class
class node:
    # Node type can be 'n' or 'v' for leaf nodes
    # Node type can be '+' or '*' for non-leaf
    def __init__(self, ntype, index, var, value, derivative):
        self.ntype = ntype
        self.nIndex = index
        self.var = var # denotes the n'th variable
        self.vr = value
        self.dr = derivative
        self.numChild = 0 # number of child nodes of non-leaf
        self.childList = [] # indices of children
        self.prL = [] # product cache register
        self.prR = []
        self.childValues = [] # values of children
        self.line = 0 # Line in output AC

    # Helper functions
    def update_vr(self, newVal):
        self.vr = newVal

    def update_dr(self, newDr):
        self.dr = newDr

    def update_childList(self, child):
        self.childList = child

    def add_child(self, child):
        self.childList.append(child)

    def update_numChild(self, num):
        self.numChild = num

    def update_prL(self, cache):
        self.prL = cache

    def update_prR(self, cache):
        self.prR = cache

    def update_childValues(self, cVal):
        self.childValues = cVal

# Stores constant pairs to be inserted into AC
class ACNode:
     def __init__(self, n, i):
         self.n = n # node
         self.index = i # line number in AC

# Functions

        
# Check if file was passed by command line
FILE_NAME = ''
FILE_SIZE = 0
if (len(sys.argv) < 2):
    FILE_NAME = input("Please input an AC file: ")
    # FILE_SIZE = input("Please input the AC size: ")
else:    
    FILE_NAME = sys.argv[1]
    #FILE_SIZE = sys.argv[2]

print("File name:", FILE_NAME)
#print("size:", FILE_SIZE)



# Read the AC file
AC_FILE = open(FILE_NAME, 'r')
for line in AC_FILE:
    #print(line)
    if line[0] == '(':
        print("\t ... reading file ...")
    elif line[0] == 'E':
        print("\t ... done reading file ...")
        index = index - 1
        #circuit[index].update_dr(1)
        break
    else:
        parsed = line.split()
        
        if parsed[0] == 'n':
            # index and dr equals 0
            val = float(parsed[1])
            circuit.append(node(parsed[0], index, 0, val, 0))
            index = index + 1

        elif parsed[0] == 'v':
            var = int(parsed[1])
            val = float(parsed[2])
            circuit.append(node(parsed[0], index, var, val, 0))
            index = index + 1

        elif parsed[0] == '+':
            value = 0
            addNode = node(parsed[0], index, 0, 0, 0)
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
            multNode = node(parsed[0], index, 0, 0, 0)
            childValues = [1.0] # first value is a placeholder
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

            value = prL[numChild]
            
            multNode.update_prL(prL)
            multNode.update_prR(prR)
            multNode.update_vr(value)
            multNode.update_childValues((childValues.pop(0)))
            circuit.append(multNode)
            index = index + 1

    #print(index)

# Close AC file
AC_FILE.close()

outputVal = circuit[index].vr
circuit[index].update_dr(1)

print("output:", outputVal)
print("output log", math.log10(outputVal))

''' 
Start cache-backpropagation
The leaf index is the index of the node in the circuit
We will find the marginal for this index
'''
leafIndex = int(input("Find marginal of leaf node index: "))
counter = 0 # counts the line/node number in the ac
isDone = False # breaks out of loop if goal node is reached
outputAC = '' # Output string

# Reverse the original circuit for "backpropagation"
reversedCircuit = list(circuit)
reversedCircuit.reverse()

outputAC += "(TBA)\n" # Placeholder

for parent in reversedCircuit:
    if parent.ntype == '+':
        parent.line = counter
        outputAC += 'n ' + str(parent.dr) + '\n'
        counter += 1
        for childIndex in parent.childList:
            childNode = circuit[childIndex]

            # End program if we reached goal node
            if leafIndex in parent.childList:
                # Create a node for child dr
                if childNode.line == 0: # child was not assigned a line
                    childNode.line = counter
                    outputAC += 'n ' + str(childNode.dr) + '\n'
                    counter += 1
                childNode.dr += parent.dr
                parentLine = str(parent.line)
                childLine = str(childNode.line)
                outputAC += '+ ' + parentLine + ' ' + childLine + '\n'
                counter += 1
                isDone = True
                break

            # Continue - we haven't reached goal node
            else: 
                childNode.dr += parent.dr
                
                # Check if the child node is in output AC
                if childNode.line == 0:
                    # Saves one edge by redirecting the child dr to parent
                    childNode.line = parent.line
                else:
                    # The child node is already in the AC
                    # Insert the add node in output AC
                    parentLine = str(parent.line)
                    childLine = str(childNode.line)
                    outputAC += '+ ' + parentLine + ' ' + childLine + '\n'
                    childNode.line = counter # Change the line to the add node
                    counter += 1

    elif parent.ntype == '*':
        pos = 1 # Keep track of node position in cache
        
        parent.line = counter
        outputAC += 'n ' + str(parent.dr) + '\n'
        counter += 1
        

        for childIndex in parent.childList:
            childNode = circuit[childIndex]
            w = parent.numChild

            # Calculate dr
            rightVal = parent.prR[(w-pos)]
            leftVal = parent.prL[(pos-1)]
            childNode.dr += parent.dr * rightVal * leftVal
            pos += 1

            # Insert constant nodes for cache value
            prLLine = str(counter)
            counter += 1
            outputAC += 'n ' + str(rightVal) + '\n'
            prRLine = str(counter)
            counter += 1
            outputAC += 'n ' + str(leftVal) + '\n'

            parentLine = str(parent.line)

            # Insert multiplication node
            outputAC += '* ' + prLLine + ' ' + prRLine + ' ' + parentLine + '\n'
            counter += 1
        
    if isDone:
        break

    #index -= 1 # decrement the index (node position in circuit)
    
outputAC += 'EOF'

#print(outputAC)

# Write to a test file
f = open('out.ac', 'w')
f.write(outputAC)
f.close()
