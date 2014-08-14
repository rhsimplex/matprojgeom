import pymatgen as pm
import matplotlib.pyplot as plt
from numpy import empty
def periodicTablePlot(elements, values, title):
    """
    Builds a plot over the periodic table
    """
    table = empty((9,18))
    table[:] = None
    entries = zip(elements, values)
    for entry in entries:
        A = pm.Element(entry[0])
        table[A.row - 1, A.group - 1] = entry[1]
        plt.text(A.group - 1, A.row - 1, entry[0], horizontalalignment = 'center', verticalalignment = 'bottom')
        if(entry[1] is not None):
            plt.text(A.group - 1, A.row - 1, entry[1], horizontalalignment = 'center', verticalalignment = 'top')
        else:
            plt.text(A.group - 1, A.row - 1, '-', horizontalalignment = 'center', verticalalignment = 'top')
    plt.axis('off')
    plt.title(title)
    plt.imshow(table, interpolation='None') 
