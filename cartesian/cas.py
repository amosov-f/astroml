import numbers
from dataclasses import dataclass, field
from typing import List, Dict, OrderedDict


@dataclass
class Summand:
    coeff: float = 1
    vars: Dict[str, int] = field(default_factory=OrderedDict)

    def __mul__(self, other):
        if isinstance(other, numbers.Number):
            return Summand(coeff=self.coeff * other, vars=self.vars)
        result = {}
        for k in set(self.vars.keys()).union(other.vars.keys()):
            new_value = self.vars.get(k, 0) + other.vars.get(k, 0)
            if new_value != 0:
                result[k] = new_value
        return Summand(coeff=self.coeff * other.coeff, vars=result)

    def __truediv__(self, other):
        return self * Summand(coeff=1/other.coeff, vars={k: -v for k, v in other.vars.items()})

    def __str__(self):
        res = ''
        if self.coeff != 1:
            res += str(self.coeff)
        for k, v in self.vars.items():
            res += f' {k}'
            if v != 1:
                res += f'**{v}'
        return res

    def __pow__(self, power, modulo=None):
        result = Summand()
        for i in range(power):
            result *= self
        return result


def plus(**kwargs):
    return Summand(vars=kwargs)

def minus(**kwargs):
    return Summand(coeff=-1, vars=kwargs)

def sum(*args):
    return Summands(list(args))

@dataclass
class Summands:
    summands: List[Summand]

    def __mul__(self, other):
        result = []
        for s1 in self.summands:
            for s2 in other.summands:
                result.append(s1 * s2)
        return Summands(result)

    def __pow__(self, power, modulo=None):
        result = self
        for i in range(power - 1):
            result *= self
        return result

    def __truediv__(self, other: Summand):
        result = []
        for s in self.summands:
            result.append(s / other)
        return Summands(result)

    def __add__(self, other):
        return Summands(self.summands + other.summands)

    def __sub__(self, other):
        return Summands(self.summands + [s * -1 for s in other.summands])

    def __str__(self):
        res = ''
        for s in self.summands:
            res += f' + {s}'
        return res

    def replace(self, var: str, expr):
        if isinstance(expr, Summand):
            expr = Summands([expr])
        result = []
        for s in self.summands:
            if var in s.vars:
                pwd = s.vars[var]
                lt = Summands([s / Summand(vars={var: pwd})]) * (expr ** pwd)
                result.extend(lt.summands)
            else:
                result.append(s)
        return Summands(result)

def replace_all(expr, replace2d=False):
    expr = expr.replace('x', plus(r=1, cosb=1, cosl=1))
    expr = expr.replace('y', plus(r=1, cosb=1, sinl=1))
    expr = expr.replace('z', plus(r=1, sinb=1))

    ax = plus(M11=1)
    ay = sum(plus(M12=1), minus(wz=1))
    az = sum(plus(M13=1), plus(wy=1))
    bx = sum(plus(M12=1), plus(wz=1))
    by = plus(M22=1)
    bz = sum(plus(M23=1), minus(wx=1))
    cx = sum(plus(M13=1), minus(wy=1))
    cy = sum(plus(M23=1), plus(wx=1))
    cz = plus(M33=1)
    
    expr = (
        expr.replace('ax', ax)
            .replace('ay', ay)
            .replace('az', az)
            .replace('bx', bx)
            .replace('by', by)
            .replace('bz', bz)
            .replace('cx', cx)
            .replace('cy', cy)
            .replace('cz', cz)
    )

    if replace2d:
        axx = sum(plus(dM11dx=1))
        ayy = sum(plus(dM12dy=1), minus(dwzdy=1))
        azz = sum(plus(dM13dz=1), plus(dwydz=1))
        axy = sum(plus(dM11dy=1), plus(dM12dx=1), minus(dwzdx=1))
        axz = sum(plus(dM11dz=1), plus(dM13dx=1), plus(dwydx=1))
        ayz = sum(plus(dM13dy=1), plus(dwydy=1), plus(dM12dz=1), minus(dwzdz=1))

        expr = (
            expr.replace('axx', axx)
            .replace('ayy', ayy)
            .replace('azz', azz)
            .replace('axy', axy)
            .replace('axz', axz)
            .replace('ayz', ayz)
        )

        bxx = sum(plus(dM12dx=1), plus(dwzdx=1))
        byy = sum(plus(dM22dy=1))
        bzz = sum(plus(dM23dz=1), minus(dwxdz=1))
        bxy = sum(plus(dM22dx=1), plus(dM12dy=1), minus(dwzdy=1))
        bxz = sum(plus(dM23dx=1), minus(dwxdx=1), plus(dM12dz=1), plus(dwzdz=1))
        byz = sum(plus(dM23dy=1), minus(dwxdy=1), plus(dM22dz=1))

        expr = (
            expr.replace('bxx', bxx)
            .replace('byy', byy)
            .replace('bzz', bzz)
            .replace('bxy', bxy)
            .replace('bxz', bxz)
            .replace('byz', byz)
        )

        cxx = sum(plus(dM13dx=1), minus(dwydx=1))
        cyy = sum(plus(dM23dy=1), plus(dwxdy=1))
        czz = plus(dM33dz=1)
        cxy = sum(plus(dM23dx=1), plus(dwxdx=1), plus(dM13dy=1), minus(dwydy=1))
        cxz = sum(plus(dM33dx=1), plus(dM13dz=1), minus(dwydz=1))
        cyz = sum(plus(dM33dy=1), plus(dM23dz=1), plus(dwxdz=1))

        expr = (
            expr.replace('cxx', cxx)
            .replace('cyy', cyy)
            .replace('czz', czz)
            .replace('cxy', cxy)
            .replace('cxz', cxz)
            .replace('cyz', cyz)
        )

    return expr

def main():
    x = Summands([plus(x=1)])
    # vy = Summands([Summand(coeff=-1, vars=dict(V=1)), Summand(vars=dict(b1=1, x=1)), Summand(vars=dict(b2=1, y=1)), Summand(vars=dict(b3=1, z=1))])
    vy = Summands([Summand(coeff=-1, vars=dict(V=1)),
                   plus(bx=1, x=1),
                   plus(by=1, y=1),
                   plus(bz=1, z=1),
                   plus(bxx=1, x=2),
                   plus(byy=1, y=2),
                   plus(bzz=1, z=2),
                   plus(bxy=1, x=1, y=1),
                   plus(bxz=1, x=1, z=1),
                   plus(byz=1, y=1, z=1),
                   ])
    y = Summands([plus(y=1)])
    # vx = Summands([Summand(coeff=-1, vars=dict(U=1)), Summand(vars=dict(a1=1, x=1)), Summand(vars=dict(a2=1, y=1)), Summand(vars=dict(a3=1, z=1))])
    vx = Summands([Summand(coeff=-1, vars=dict(U=1)),
                   plus(ax=1, x=1),
                   plus(ay=1, y=1),
                   plus(az=1, z=1),
                   plus(axx=1, x=2),
                   plus(ayy=1, y=2),
                   plus(azz=1, z=2),
                   plus(axy=1, x=1, y=1),
                   plus(axz=1, x=1, z=1),
                   plus(ayz=1, y=1, z=1),
                   ])
    z = Summands([plus(z=1)])
    vz = Summands([minus(W=1),
                   plus(cx=1, x=1),
                   plus(cy=1, y=1),
                   plus(cz=1, z=1),
                   plus(cxx=1, x=2),
                   plus(cyy=1, y=2),
                   plus(czz=1, z=2),
                   plus(cxy=1, x=1, y=1),
                   plus(cxz=1, x=1, z=1),
                   plus(cyz=1, y=1, z=1),
                   ])
    k_mul_cosb = (x * vy - y * vx) / plus(r=2, cosb=1)

    # print(k_mul_cosb)

    k_mul_cosb = replace_all(k_mul_cosb)

    # print(k_mul_cosb)

    kmub = ((x ** 2 + y ** 2) * vz - z * (x * vx + y * vy)) / plus(r=3, cosb=1)
    print(kmub)
    return

    kmub = replace_all(kmub)

    vrr = (x * vx + y * vy + z * vz) / plus(r=2)
    vrr = replace_all(vrr)

    print('kmub')
    print(kmub)

    print('kmulcosb')
    print(k_mul_cosb)

    print('vrr')
    print(vrr)

    # print('kmulcosb')
    # print(k_mul_cosb)


if __name__ == '__main__':
    main()