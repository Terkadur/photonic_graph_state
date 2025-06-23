from sympy import I, pi, sqrt, sin, cos, exp, sign, simplify, latex, Expr, Function, Add, Mul
from sympy.physics.quantum import qapply
from sympy.physics.quantum.qubit import Qubit, QubitBra, qubit_to_matrix
from sympy.physics.quantum.gate import IdentityGate, X, Y, Z, H, S, T, CNOT, CPHASE


def dephase(state: Expr) -> Expr:
    '''
    Removes the global phase from a SymPy statevector by dividing by the phase of the first component

    :param state: The SymPy statevector to be dephased
    :type state: Expr

    :return: state without the global phase
    :rtype: Expr
    '''
    try:
        res = qubit_to_matrix(state)
        glob_phase = 1
        for elem in res:
            if elem != 0:
                glob_phase = sign(elem)
                break

        return state / glob_phase
    except:
        return state


def apply(state: Expr, operation: Expr = None, dp: bool = False) -> Expr:
    '''
    Simplifies an expression describing a SymPy statevector

    :param state: The SymPy statevector to be simplified
    :type state: Expr

    :param operation: An operation (usually a quantum gate) to apply to the state, defaults to None
    :type operation: Expr

    :param dp: A boolean determining if the statevector should be dephased, defaults to False
    :type dp: bool

    :return: The simplified version of the SymPy statevector
    :rtype: Expr

    '''
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


def tensor(state1: Expr, state2: Expr) -> Expr:
    '''
    Computes the tensor product of two SymPy statevectors (not fully tested)

    :param state1: The first SymPy statevector in the tensor product
    :type state1: Expr

    :param state2: The second SymPy statevector in the tensor product
    :type state2: Expr

    :return: The SymPy statevector representing the tensor product of state1 and state2
    :rtype: Expr
    '''

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
    '''
    A custom quantum gate describing a qubit rotation about the x axis

    :param arg1: The index of the qubit to apply the rotation to
    :type arg1: int

    :param arg2: The angle by which to rotate the qubit
    :type arg2: float
    '''
    def doit(self, **kwargs):
        return cos(self.args[1]/2)*IdentityGate(self.args[0]) - I*sin(self.args[1]/2)*X(self.args[0])

    def _latex(self, printer):
        return r"R_{x, %s}\left(%s\right)" % (latex(self.args[0]), latex(self.args[1]))


class Ry(Function):
    '''
    A custom quantum gate describing a qubit rotation about the y axis

    :param arg1: The index of the qubit to apply the rotation to
    :type arg1: int

    :param arg2: The angle by which to rotate the qubit
    :type arg2: float
    '''
    def doit(self, **kwargs):
        return cos(self.args[1]/2)*IdentityGate(self.args[0]) - I*sin(self.args[1]/2)*Y(self.args[0])

    def _latex(self, printer):
        return r"R_{y, %s}\left(%s\right)" % (latex(self.args[0]), latex(self.args[1]))


class Rz(Function):
    '''
    A custom quantum gate describing a qubit rotation about the z axis

    :param arg1: The index of the qubit to apply the rotation to
    :type arg1: int

    :param arg2: The angle by which to rotate the qubit
    :type arg2: float
    '''
    def doit(self, **kwargs):
        return cos(self.args[1]/2)*IdentityGate(self.args[0]) - I*sin(self.args[1]/2)*Z(self.args[0])

    def _latex(self, printer):
        return r"R_{z, %s}\left(%s\right)" % (latex(self.args[0]), latex(self.args[1]))
