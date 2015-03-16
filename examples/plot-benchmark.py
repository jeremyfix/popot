
import sys
import json
import glob
import matplotlib.pyplot as plt
import numpy as np


json_files = []

if(len(sys.argv) == 1):
    # We list all the .json files in the current directory
    json_files = glob.glob('./*.json')
    use_selector = True
else:
    json_files = sys.argv[1:]
    use_selector = False


if(use_selector):
    print("Using an interface to select the files to plot")

# We can now propose a Text interface with 
# 1) a selection of the problem to work on (
#    -> display a list of dimensions
#    -> display a list of algorithm

for f in json_files:
    fp = open(f)
    results = json.load(fp)

    error = results['mean']
    std = results['std']
    FE = results['FE']
    algo = results['algorithm']
    problem = results['problem']
    title = "%s, d = %i" % (problem['name'], problem['dimension'])

    plt.plot(FE, np.log(error), label=algo)
    #plt.errorbar(FE, np.log(error), np.log(std), linestyle='None', marker='^', label=algo)

#plt.title(title)
plt.ylabel('Mean function value (log)')
plt.xlabel('Function evaluations')

plt.legend()


plt.show()
