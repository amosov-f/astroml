from common.gaia.with_rv import read_raw_gaia_with_rv, read_raw_gaia_with_rv_no_errors
from matplotlib import pyplot as plt

def main():
    df = read_raw_gaia_with_rv_no_errors()

    px = 1 / df['parallax']
    plt.hist(px, bins=100, range=(-0.1, 10))
    # plt.show()
    plt.xlabel('Расстояние, кпк')
    plt.ylabel('Число звезд')
    plt.title('Распределение расстояний звезд GAIA DR3 with RV')
    plt.savefig('distance_dr3.png')

if __name__ == '__main__':
    main()