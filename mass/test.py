from math import sin, cos, acos, asin, atan2, radians, degrees, pi

L0 = 0.5747703990741704  # 32.931918056


def angle(l, b):
    if not l or not b:
        return None
    l = radians(l)
    b = radians(b)
    ro = acos(cos(b) * cos(l - L0))
    return degrees(ro)


def main():
    print(angle(32.9, 22.5) < 1.0)


if __name__ == '__main__':
    main()
