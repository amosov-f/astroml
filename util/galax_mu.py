import numpy


def galax_mu(mua, mud, l, b, d):
    L0 = 32.931918056
    ci = 0.45598379779  # sin 62.871748611
    si = 0.88998807641  # cos 62.871748611

    cd = numpy.cos(d)
    sfi = si * numpy.cos(l - L0) / cd
    cfi = (numpy.cos(b) * ci - numpy.sin(b) * si * numpy.sin(l - L0)) / cd
    mul = cfi * mua + sfi * mud
    mub = -sfi * mua + cfi * mud
    return (mul, mub)

if __name__ == '__main__':
    print galax_mu(0, 0, 0, 0, 0)