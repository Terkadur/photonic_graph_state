from sympy import I, pi, sqrt, sin, cos, exp, sign, simplify, latex, Expr, Function, Add, Mul
from sympy.physics.quantum import qapply
from sympy.physics.quantum.qubit import Qubit, QubitBra, qubit_to_matrix
from sympy.physics.quantum.gate import IdentityGate, X, Y, Z, H, S, T, CNOT, CPHASE


def dephase(state: Expr) -> Expr:
    res = qubit_to_matrix(state)
    glob_phase = 1
    for elem in res:
        if elem != 0:
            glob_phase = sign(elem)
            break

    return state / glob_phase
    

def apply(state: Expr, operation: Expr = None, dp=False) -> Expr:
    if operation is None:
        if dp:
            return simplify(dephase(qapply((state).doit()))).expand()
        else:
            return simplify(qapply((state).doit())).expand()
    else:
        if dp:
            return simplify(dephase(qapply((operation * state).doit()))).expand()
        else:
            return simplify(qapply((operation * state).doit())).expand()


def tensor(state1, state2) -> Expr:
    term1 = apply(state1)
    term2 = apply(state2)

    if isinstance(term1, Qubit):
        if isinstance(term2, Qubit):
            bits = term1.args + term2.args
            string = ''
            for bit in bits:
                string += str(bit)
            return Qubit(string)
        elif isinstance(term2, Add):
            result = tensor(term1, term2.args[0])
            for i in range(1, len(term2.args)):
                result += tensor(term1, term2.args[i])
            return result
        elif isinstance(term2, Mul) and isinstance(term2.args[-1], Qubit):
            result = term2.args[0]
            for i in range(1, len(term2.args) - 1):
                result *= term2.args[i]
            return result * tensor(term1, term2.args[-1])

    elif isinstance(term1, Add):
        result = tensor(term1.args[0], term2)
        for i in range(1, len(term1.args)):
            result += tensor(term1.args[i], term2)
        return result

    elif isinstance(term1, Mul) and isinstance(term1.args[-1], Qubit):
        result = term1.args[0]
        for i in range(1, len(term1.args) - 1):
            result *= term1.args[i]
        return result * tensor(term1.args[-1], term2)

    else:
        raise TypeError("Unexpected input")


class Rx(Function):
    def doit(self, **kwargs):
        return cos(self.args[1]/2)*IdentityGate(self.args[0]) - I*sin(self.args[1]/2)*X(self.args[0])

    def _latex(self, printer):
        return r"R_{x, %s}\left(%s\right)" % (latex(self.args[0]), latex(self.args[1]))


class Ry(Function):
    def doit(self, **kwargs):
        return cos(self.args[1]/2)*IdentityGate(self.args[0]) - I*sin(self.args[1]/2)*Y(self.args[0])

    def _latex(self, printer):
        return r"R_{y, %s}\left(%s\right)" % (latex(self.args[0]), latex(self.args[1]))


class Rz(Function):
    def doit(self, **kwargs):
        return cos(self.args[1]/2)*IdentityGate(self.args[0]) - I*sin(self.args[1]/2)*Z(self.args[0])

    def _latex(self, printer):
        return r"R_{z, %s}\left(%s\right)" % (latex(self.args[0]), latex(self.args[1]))


