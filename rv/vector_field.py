import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats

from GC2019.gc2019 import read_gc2019


def main():
    layer = (-100, 100)

    xs = []
    ys = []
    zs = []
    vxs = []
    vys = []
    vzs = []

    first_line = True
    counter = 0

    r = 5000

    # f = open(f'xyz{r}.tsv', 'w')
    cat = read_gc2019()

    # for line in open('xyz.tsv', 'r').readlines():
    #     if first_line:
    #         first_line = False
    #         continue
    #     parts = line.strip().split('\t')
    #     x = float(parts[0])
    #     y = float(parts[1])
    #     z = float(parts[2])
    #     vx = float(parts[3])
    #     vy = float(parts[4])
    #     vz = float(parts[5])
    for index, row in cat.iterrows():
        x = row['x']
        y = row['y']
        z = row['z']
        vx = row['vx']
        vy = row['vy']
        vz = row['vz']

        counter += 1
        if counter % 10000 == 0:
            print(counter)

        # if z < layer[0] or z > layer[1]:
        #     continue
        if -r < x < r and -r < y < r and -r < z < r:
            xs.append(x)
            ys.append(y)
            zs.append(z)
            # vxs.append(vx - -10.34848874410621)
            # vys.append(vy - -17.554795954741007)
            # vzs.append(vz - -7.233114343092995)
            vxs.append(vx * 3)
            vys.append((vy)  * 3)
            vzs.append(vz * 3)
        # if -r < x < r and -r < y < r and -r < z < r:
        #     f.write(f'{x}\t{y}\t{z}\t{vx}\t{vy}\t{vz}\n')

    # COUNT = plt.hist2d(xs,
    #            ys,
    #            bins=(100, 100),
    #            range=((-r, r), (-r, r)),
    #            cmap='Blues')[0]
    #
    #
    # # Add a color bar to the histogram
    #
    # # Add labels, title, and display the plot
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.title(f'{COUNT.sum()} stars (min {COUNT.min()}, max {COUNT.max()})')
    # plt.show()
    # return

    # f.close()

    bins = 8
    # r = 200.0

    EDGES = np.linspace(-r, r, bins + 1)
    x = np.linspace(-r + r / bins, r - r / bins, bins)
    X, Y, Z = np.meshgrid(x, x, x)
    M = stats.binned_statistic_dd([xs, ys, zs], [vxs, vys, vzs], statistic='median', bins=[EDGES, EDGES, EDGES]).statistic
    c = COUNT = stats.binned_statistic_dd([xs, ys, zs], values=None, statistic='count', bins=[EDGES, EDGES, EDGES]).statistic

    c = c.ravel() / c.max()
    # Repeat for each body line and two head lines
    c = np.concatenate((c, np.repeat(c, 2)))
    # Colormap
    c = plt.cm.hsv(c)

    fig = plt.figure()
    ax = Axes3D(fig)
    # q = ax.quiver(X, Y, Z, *M, length=bins / 4, colors=c)
    q = ax.quiver(xs, ys, zs, vxs, vys, vzs, colors=c)
    # plt.colorbar()

    # avex = np.median(M[0])
    # avey = np.median(M[1])
    # avez = np.median(M[2])
    # print(avex)
    # print(avey)
    # print(avez)
    # print(np.median(vxs))
    # print(np.median(vys))
    # print(np.median(vzs))
    plt.title(f'{int(COUNT.sum())} (min {int(COUNT.min())}, max {int(COUNT.max())}) starts in cube with edge {2 * r}')
    plt.draw()
    plt.show()



if __name__ == '__main__':
    main()