import numpy
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import spline
from sklearn.linear_model import Lasso


def main():
    dist = [1, 2, 3]
    U = [3, 3, 2]
    V = [4, 4, 4]
    W = [5, 5, 5]

    sp1 = plt.figure(1).add_subplot(111)
    sp2 = plt.figure(2).add_subplot(111)

    xnew = np.linspace(min(dist), max(dist), 300)
    sp1.plot(dist, U, 'ko')

    sp1.plot(xnew, spline(dist, U, xnew), c=plt.cm.get_cmap('hsv', 15)(13), label='U')
    sp2.plot(dist, V, 'ko')
    sp2.plot(dist, V, 'b', label='V')
    sp2.plot(dist, W, 'ko')
    sp2.plot(dist, W, 'r', label='W')
    sp1.legend()
    sp2.legend()
    plt.show()

if __name__ == '__main__':
    l = Lasso()
    print(l.random_state)