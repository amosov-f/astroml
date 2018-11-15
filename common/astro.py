from numpy import sin, cos, arcsin, arctan2, pi
import numpy as np

Leo = 4.936829260965283  # 282.85948083°
L0 = 0.5747703990741704  # 32.931918056°
si = 0.88998807641_8  # sin62.871748611°
ci = 0.45598379779_8  # cos62.871748611°


def galaxy(a, d):
    al = a - Leo
    sa = sin(al)
    ca = cos(al)
    sd = sin(d)
    cd = cos(d)
    b = arcsin(sd * ci - cd * si * sa)
    l = arctan2(sd * si + cd * ci * sa, cd * ca) + L0
    l = np.where(l < 0, l + 2 * pi, l)
    # if l < 0:
    #     l += 2 * pi
    return l, b


def galaxy_mu(mua, mud, l, b, d):
    cd = cos(d)
    sfi = si * cos(l - L0) / cd
    cfi = (cos(b) * ci - sin(b) * si * sin(l - L0)) / cd
    mul = cfi * mua + sfi * mud
    mub = - sfi * mua + cfi * mud
    return mul, mub