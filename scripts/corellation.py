import numpy as np

def main():
    xs = np.random.randn(100, 5)
    noise = np.random.randn(100)  # white noise
    xs[:, 4] = 2 * xs[:, 0] + 3 * xs[:, 2] + .5 * noise  # collinearity


    corr = np.corrcoef(xs, rowvar=0)

    print(corr)

    w, v = np.linalg.eig(corr)
    print(w)
    print(v[:,0])

if __name__ == '__main__':
    main()