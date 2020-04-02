import matplotlib.pyplot as plt
import numpy as np
import pandas
import pylab

def plot_profile(directory, number):
    filename = directory + "/model_" + str(number) + ".profile"
    f = open(filename)
    vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10:
            spl = line.split()
            vals.append(float(spl[-1]))


    pylab.figure(1, figsize=(10,6))
    pylab.xlabel('Residue position')
    pylab.ylabel('DOPE per-residue score')
    pylab.plot(vals, color='red', linewidth=2, label='Model ' + str(number))
    pylab.legend()
    pylab.savefig(directory + "/model_" + str(number) + "_profile.png", dpi=65)
