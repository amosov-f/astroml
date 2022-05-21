def main():
    t211 = 14.8
    s310 = -12.1
    v310 = -8.5

    r2 = 0.933

    dM12dr1_s = s310 / (-0.26) / (r2)

    dM12dr1_v = v310 / (-0.22) / (r2)

    print(dM12dr1_v)
    print(dM12dr1_s)

    dM12dr1 = (dM12dr1_v + dM12dr1_s) / 2
    print(dM12dr1)

    dw3dr1 = (t211 / (0.37 * (r2)) - dM12dr1) / 3

    print(dw3dr1)

    print(dM12dr1 * 0.37)

    A4 = 20.6
    A3 = 17.5
    A2 = 15
    A1 = 11.1

    r4 = 1.447
    r3 = 1.168
    r2 = 0.933
    r1 = 0.757

    print((A2 - A1) / (r2 - r1))
    print((A3 - A2) / (r3 - r2))

    dM12dr1_2 = ((A2 - A1) / (r2 - r1))
    dM12dr1_4 = ((A4 - A3) / (r4 - r3))

    s310_2 = -0.26 * (dM12dr1_2) * r2
    s310_4 = -0.26 * (dM12dr1_4) * r4

    print(s310_2)
    print(s310_4)

    # dM12dr1_t = t211 / (-0.22) / (0.933)



if __name__ == '__main__':
    main()