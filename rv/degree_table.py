from dataclasses import dataclass

from numpy import sin, cos
from setuptools._vendor.ordered_set import OrderedSet

full_print = False


@dataclass
class Deg:
    column: str
    coeff: float
    sinb: int
    cosb: int
    sinl: int
    cosl: int
    row_type: str = ""
    cos2b: int = 0
    cos2l: int = 0

    def __mul__(self, other):
        col = self.column if self.column == other.column else f"{self.column}x{other.column}"
        return Deg(column=col,
                   row_type=f"{self.row_type}*{other.row_type}",
                   coeff=self.coeff * other.coeff,
                   sinb=self.sinb + other.sinb,
                   cosb=self.cosb + other.cosb,
                   sinl=self.sinl + other.sinl,
                   cosl=self.cosl + other.cosl,
                   cos2b=self.cos2b + other.cos2b,
                   cos2l=self.cos2l + other.cos2l)

    def __rmul__(self, other):
        return other * self

    def __pow__(self, power, modulo=None):
        res = self
        for i in range(1, power):
            res *= self
        return res

    def __str__(self):
        # s = f"{self.row_type}"
        s = [str(self.coeff)] if self.coeff != 1 else []
        if self.sinb > 0:
            # s += f" 0"
            return "0"
        # if self.sinb or full_print:
        #     s += f" * sin^{self.sinb} b"
        # if self.cosb or full_print:
        #     s += f" * cos^{self.cosb} b"
        if self.sinl or full_print:
            s.append(f"sin{self.__pow(self.sinl)} l")
        if self.cosl or full_print:
            s.append(f"cos{self.__pow(self.cosl)} l")
        # if self.cos2b or full_print:
        #     s += f" * cos^{self.cos2b} 2b"
        if self.cos2l or full_print:
            s.append(f"cos{self.__pow(self.cos2l)} 2l")
        return " * ".join(s)

    def __pow(self, n):
        return "" if n == 1 else f"^{n}"

    def key(self):
        if self.sinb > 0:
            return (
                self.column,
                0,
                0,
                0
            )
        # if self.cosb
        return (
            self.column,
            # self.sinb,
            # self.cosb,
            self.sinl,
            self.cosl,
            # self.cos2b,
            self.cos2l
        )

    def with_name(self, row_type, column):
        return Deg(column=column,
                   row_type=row_type,
                   coeff=self.coeff,
                   sinb=self.sinb,
                   cosb=self.cosb,
                   sinl=self.sinl,
                   cosl=self.cosl,
                   cos2b=self.cos2b,
                   cos2l=self.cos2l)

    def apply(self, b, l):
        res = self.coeff * sin(l) ** self.sinl * cos(l) ** self.cosl * cos(2 * l) ** self.cos2l

        # hui = sin(b) ** self.sinb * cos(b) ** self.cosb * cos(2 * b) ** self.cos2b
        # res *= hui

        if self.sinb > 0:
            return 0

        return res

    def __eq__(self, other):
        return self.key() == other.key()


R1 = Deg('R1', coeff=1, sinb=0, cosb=1, sinl=0, cosl=1)
R2 = Deg('R2', coeff=1, sinb=0, cosb=1, sinl=1, cosl=0)
R3 = Deg('R3', coeff=1, sinb=1, cosb=0, sinl=0, cosl=0)

la1 = Deg('la1', coeff=1, sinb=1, cosb=0, sinl=0, cosl=1)
la2 = Deg('la2', coeff=1, sinb=1, cosb=0, sinl=1, cosl=0)
la3 = Deg('la3', coeff=1, sinb=0, cosb=1, sinl=0, cosl=0)
ta1 = Deg('ta1', coeff=1, sinb=0, cosb=0, sinl=1, cosl=0)
ta2 = Deg('ta2', coeff=1, sinb=0, cosb=0, sinl=0, cosl=1)

p11 = Deg('p11', coeff=-1, sinb=0, cosb=1, sinl=1, cosl=1)
p12 = Deg('p12', coeff=1, sinb=0, cosb=1, sinl=0, cosl=0, cos2l=1)
# p12_cos2 = Deg('p12_cos2', sinb=0, cosb=1, sinl=0, cosl=2)
# p12_1 = Deg('p12_1', sinb=0, cosb=1, sinl=0, cosl=0)
p13 = Deg('p13', coeff=-1, sinb=1, cosb=0, sinl=1, cosl=0)
p22 = Deg('p22', coeff=1, sinb=0, cosb=1, sinl=1, cosl=1)
p23 = Deg('p23', coeff=1, sinb=1, cosb=0, sinl=0, cosl=1)

q11 = Deg('q11', coeff=-1, sinb=1, cosb=1, sinl=0, cosl=2)
q12 = Deg('q12', coeff=-1, sinb=1, cosb=1, sinl=0, cosl=0, cos2l=1)
# q12_cos2 = Deg('q12_cos2', sinb=1, cosb=1, sinl=0, cosl=2)
# q12_1 = Deg('q12_1', sinb=1, cosb=1, sinl=0, cosl=0)
q13 = Deg('q13', coeff=1, sinb=0, cosb=0, sinl=0, cosl=1, cos2b=1)
# q13_cos2 = Deg('q13_cos2', sinb=0, cosb=2, sinl=0, cosl=1)
# q13_1 = Deg('q13_1', sinb=0, cosb=0, sinl=0, cosl=1)
q22 = Deg('q22', coeff=-1, sinb=1, cosb=1, sinl=2, cosl=0)
q23 = Deg('q23', coeff=1, sinb=0, cosb=0, sinl=1, cosl=0, cos2b=1)
# q23_cos2 = Deg('q23_cos2', sinb=0, cosb=2, sinl=1, cosl=0)
# q23_1 = Deg('q23_1', sinb=0, cosb=0, sinl=1, cosl=0)
q33 = Deg('q33', coeff=1, sinb=1, cosb=1, sinl=0, cosl=0)

mul_dw1dr1 = (R1 * la1).with_name('mul', 'dw1dr1')
mul_dw1dr2 = (R2 * la1).with_name('mul', 'dw1dr2')
mul_dw1dr3 = (R3 * la1).with_name('mul', 'dw1dr3')
mul_dw2dr1 = (R1 * la2).with_name('mul', 'dw2dr1')
mul_dw2dr2 = (R2 * la2).with_name('mul', 'dw2dr2')
mul_dw2dr3 = (R3 * la2).with_name('mul', 'dw2dr3')
mul_dw3dr1 = (R1 * la3).with_name('mul', 'dw3dr1')
mul_dw3dr2 = (R2 * la3).with_name('mul', 'dw3dr2')
mul_dw3dr3 = (R3 * la3).with_name('mul', 'dw3dr3')

mub_dw1dr1 = (R1 * ta1).with_name('mub', 'dw1dr1')
mub_dw1dr2 = (R2 * ta1).with_name('mub', 'dw1dr2')
mub_dw1dr3 = (R3 * ta1).with_name('mub', 'dw1dr3')
mub_dw2dr1 = (R1 * ta2).with_name('mub', 'dw2dr1')
mub_dw2dr2 = (R2 * ta2).with_name('mub', 'dw2dr2')
mub_dw2dr3 = (R3 * ta2).with_name('mub', 'dw2dr3')

mub_dM11dr1 = (R1 * q11).with_name('mub', 'dM11dr1')
mub_dM11dr2 = (R2 * q11).with_name('mub', 'dM11dr2')
mub_dM11dr3 = (R3 * q11).with_name('mub', 'dM11dr3')
mub_dM12dr1 = (R1 * q12).with_name('mub', 'dM12dr1')
mub_dM12dr2 = (R2 * q12).with_name('mub', 'dM12dr2')
mub_dM12dr3 = (R3 * q12).with_name('mub', 'dM12dr3')
# mub_dM12dr1_cos2 = (R1 * q12_cos2).with_name('mub_dM12dr1_cos2')
# mub_dM12dr2_cos2 = (R2 * q12_cos2).with_name('mub_dM12dr2_cos2')
# mub_dM12dr3_cos2 = (R3 * q12_cos2).with_name('mub_dM12dr3_cos2')
# mub_dM12dr1_1 = (R1 * q12_1).with_name('mub_dM12dr1_1')
# mub_dM12dr2_1 = (R2 * q12_1).with_name('mub_dM12dr2_1')
# mub_dM12dr3_1 = (R3 * q12_1).with_name('mub_dM12dr3_1')
mub_dM13dr1 = (R1 * q13).with_name('mub', 'dM13dr1')
mub_dM13dr2 = (R2 * q13).with_name('mub', 'dM13dr2')
mub_dM13dr3 = (R3 * q13).with_name('mub', 'dM13dr3')
mub_dM22dr1 = (R1 * q22).with_name('mub', 'dM22dr1')
mub_dM22dr2 = (R2 * q22).with_name('mub', 'dM22dr2')
mub_dM22dr3 = (R3 * q22).with_name('mub', 'dM22dr3')
mub_dM23dr1 = (R1 * q23).with_name('mub', 'dM23dr1')
mub_dM23dr2 = (R2 * q23).with_name('mub', 'dM23dr2')
mub_dM23dr3 = (R3 * q23).with_name('mub', 'dM23dr3')
mub_dM33dr1 = (R1 * q33).with_name('mub', 'dM33dr1')
mub_dM33dr2 = (R2 * q33).with_name('mub', 'dM33dr2')
mub_dM33dr3 = (R3 * q33).with_name('mub', 'dM33dr3')

mul_dM11dr1 = (R1 * p11).with_name('mul', 'dM11dr1')
mul_dM11dr2 = (R2 * p11).with_name('mul', 'dM11dr2')
mul_dM11dr3 = (R3 * p11).with_name('mul', 'dM11dr3')
mul_dM12dr1 = (R1 * p12).with_name('mul', 'dM12dr1')
mul_dM12dr2 = (R2 * p12).with_name('mul', 'dM12dr2')
mul_dM12dr3 = (R3 * p12).with_name('mul', 'dM12dr3')
# mul_dM12dr1_cos2 = (R1 * p12_cos2).with_name('mul_dM12dr1_cos2')
# mul_dM12dr2_cos2 = (R2 * p12_cos2).with_name('mul_dM12dr2_cos2')
# mul_dM12dr3_cos2 = (R3 * p12_cos2).with_name('mul_dM12dr3_cos2')
# mul_dM12dr1_1 = (R1 * p12_1).with_name('mul_dM12dr1_1')
# mul_dM12dr2_1 = (R2 * p12_1).with_name('mul_dM12dr2_1')
# mul_dM12dr3_1 = (R3 * p12_1).with_name('mul_dM12dr3_1')
mul_dM13dr1 = (R1 * p13).with_name('mul', 'dM13dr1')
mul_dM13dr2 = (R2 * p13).with_name('mul', 'dM13dr2')
mul_dM13dr3 = (R3 * p13).with_name('mul', 'dM13dr3')
mul_dM22dr1 = (R1 * p22).with_name('mul', 'dM22dr1')
mul_dM22dr2 = (R2 * p22).with_name('mul', 'dM22dr2')
mul_dM22dr3 = (R3 * p22).with_name('mul', 'dM22dr3')
mul_dM23dr1 = (R1 * p23).with_name('mul', 'dM23dr1')
mul_dM23dr2 = (R2 * p23).with_name('mul', 'dM23dr2')
mul_dM23dr3 = (R3 * p23).with_name('mul', 'dM23dr3')

vr_dM11dr1 = (R1 * R1 ** 2).with_name('vr', 'dM11dr1')
vr_dM11dr2 = (R2 * R1 ** 2).with_name('vr', 'dM11dr2')
vr_dM11dr3 = (R3 * R1 ** 2).with_name('vr', 'dM11dr3')
vr_dM22dr1 = (R1 * R2 ** 2).with_name('vr', 'dM22dr1')
vr_dM22dr2 = (R2 * R2 ** 2).with_name('vr', 'dM22dr2')
vr_dM22dr3 = (R3 * R2 ** 2).with_name('vr', 'dM22dr3')
vr_dM33dr1 = (R1 * R3 ** 2).with_name('vr', 'dM33dr1')
vr_dM33dr2 = (R2 * R3 ** 2).with_name('vr', 'dM33dr2')
vr_dM33dr3 = (R3 * R3 ** 2).with_name('vr', 'dM33dr3')
vr_dM12dr1 = (R1 * R1 * R2).with_name('vr', 'dM12dr1')
vr_dM12dr2 = (R2 * R1 * R2).with_name('vr', 'dM12dr2')
vr_dM12dr3 = (R3 * R1 * R2).with_name('vr', 'dM12dr3')
vr_dM13dr1 = (R1 * R1 * R3).with_name('vr', 'dM13dr1')
vr_dM13dr2 = (R2 * R1 * R3).with_name('vr', 'dM13dr2')
vr_dM13dr3 = (R3 * R1 * R3).with_name('vr', 'dM13dr3')
vr_dM23dr1 = (R1 * R2 * R3).with_name('vr', 'dM23dr1')
vr_dM23dr2 = (R2 * R2 * R3).with_name('vr', 'dM23dr2')
vr_dM23dr3 = (R3 * R2 * R3).with_name('vr', 'dM23dr3')


def main():
    all = [
        mul_dw1dr1,
        mul_dw1dr2,
        mul_dw1dr3,

        mul_dw2dr1,
        mul_dw2dr2,
        mul_dw2dr3,

        mul_dw3dr1,
        mul_dw3dr2,
        mul_dw3dr3,

        mub_dw1dr1,
        mub_dw1dr2,
        mub_dw1dr3,

        mub_dw2dr1,
        mub_dw2dr2,
        mub_dw2dr3,

        mul_dM11dr1,
        mul_dM11dr2,
        mul_dM11dr3,

        mul_dM12dr1,
        mul_dM12dr2,
        mul_dM12dr3,

        mul_dM13dr1,
        mul_dM13dr2,
        mul_dM13dr3,

        mul_dM22dr1,
        mul_dM22dr2,
        mul_dM22dr3,

        mul_dM23dr1,
        mul_dM23dr2,
        mul_dM23dr3,

        mub_dM11dr1,
        mub_dM11dr2,
        mub_dM11dr3,

        mub_dM12dr1,
        mub_dM12dr2,
        mub_dM12dr3,

        mub_dM13dr1,
        mub_dM13dr2,
        mub_dM13dr3,

        mub_dM22dr1,
        mub_dM22dr2,
        mub_dM22dr3,

        mub_dM23dr1,
        mub_dM23dr2,
        mub_dM23dr3,

        mub_dM33dr1,
        mub_dM33dr2,
        mub_dM33dr3,

        vr_dM11dr1,
        vr_dM11dr2,
        vr_dM11dr3,

        vr_dM22dr1,
        vr_dM22dr2,
        vr_dM22dr3,

        vr_dM33dr1,
        vr_dM33dr2,
        vr_dM33dr3,

        vr_dM12dr1,
        vr_dM12dr2,
        vr_dM12dr3,

        vr_dM13dr1,
        vr_dM13dr2,
        vr_dM13dr3,

        vr_dM23dr1,
        vr_dM23dr2,
        vr_dM23dr3,
    ]

    col2eqs = {}
    all_row_types = OrderedSet()
    for el in all:
        if el.column not in col2eqs:
            col2eqs[el.column] = {}
        row = col2eqs[el.column]
        row[el.row_type] = str(el)
        if el.row_type not in row:
            row[el.row_type] = 0
        all_row_types.add(el.row_type)

    col2eqs_str = {}
    for col, eqs in col2eqs.items():

        s = []
        for row_type in all_row_types:
            eq = eqs.get(row_type, "0")
            s.append(eq)
        col2eqs_str[col] = "\t".join(s)
        # print()

    all_cols = list(col2eqs_str.keys())
    all_cols.sort(key=lambda col: col2eqs_str[col])



    for row_type in all_row_types:
        print(row_type, end='\t')
    print()
    for col in all_cols:
        print(f"{col}\t{col2eqs_str[col]}")


    # all.sort(key=lambda el: el.key())
    #
    # for i in range(len(all) - 1):
    #     print(all[i])
    #     if all[i] != all[i + 1]:
    #         print()


if __name__ == '__main__':
    main()
