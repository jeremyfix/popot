
import sys
import json
import glob
import matplotlib.pyplot as plt


# We list all the .json files in the current directory
# and collect for all the files
# the problem name, algorithm name, dimension and the associated json file
# 

# We can now propose a Text interface with 
# 1) a selection of the problem to work on (
#    -> display a list of dimensions
#    -> display a list of algorithm


fp = open("results.json")
results = json.load(fp)
print(results.keys())

plt.figure()
#plt.show()
