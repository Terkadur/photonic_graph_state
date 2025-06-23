from stabilizer import *
import math
from sympy.physics.quantum.qubit import measure_partial
from itertools import permutations


def rref(state: Stabilizer) -> None:
    N = state.size
    K = N
    KU = 0
    NL = 0
    while NL < N-1 and KU < K-1:
        zeroitem = True
        oneitem = False
        twoitem = False
        r1 = N
        r2 = N
        for k in range(KU, K):
            if state.tab[k, NL] != 0 or state.tab[k, NL+N] != 0:
                r1 = k
                zeroitem = False
                oneitem = True
                break
        for k in range(r1, K):
            if state.tab[k, NL] != 0 or state.tab[k, NL+N] != 0:
                if state.tab[k, NL] != state.tab[r1, NL] or state.tab[k, NL+N] != state.tab[r1, NL+N]:
                    r2 = k
                    oneitem = False
                    twoitem = True
                    break
        if zeroitem:
            NL += 1
        elif oneitem:
            state.swap(KU, r1)
            for i in range(KU+1, K):
                if state.tab[i, NL] != 0 or state.tab[i, NL+N] != 0:
                    state.row_add(KU, i)
            KU += 1
            NL += 1
        elif twoitem:
            state.swap(KU, r1)
            state.swap(KU+1, r2)
            for i in range(KU+2, K):
                if state.tab[i, NL] != 0 or state.tab[i, NL+N] != 0:
                    if state.tab[i, NL] == state.tab[KU, NL] and state.tab[i, NL+N] == state.tab[KU, NL+N]:
                        state.row_add(KU, i)
                    elif state.tab[i, NL] == state.tab[KU+1, NL] and state.tab[i, NL+N] == state.tab[KU+1, NL+N]:
                        state.row_add(KU+1, i)
                    else:
                        state.row_add(KU, i)
                        state.row_add(KU+1, i)
            NL += 1
            KU += 2


def heightfunction(state: Stabilizer) -> list[int]:
    rref(state)
    N = state.size
    leftmost = []
    for i in range(state.size):
        for j in range(state.size):
            if state.tab[i, j] != 0 or state.tab[i, j+N]:
                leftmost.append(j+1)
                break
    height = []
    for i in range(state.size+1):
        count = sum(j > i for j in leftmost)
        height.append(state.size-i-count)
    return height


def plot_height(state: Stabilizer) -> None:
    height = heightfunction(state)
    x_val = [i for i in range(state.size+1)]
    tickers = range(math.floor(min(height)), math.ceil(max(height))+1)
    plt.grid(color='blue', linewidth=0.5)
    plt.plot(x_val, height, color='blue')
    plt.scatter(x_val, height, color='blue')
    plt.yticks(tickers)
    plt.title('Target Height Function')
    plt.xlabel('x')
    plt.ylabel('h(x)')
    plt.show()


def num_emitters(state: Stabilizer) -> int:
    height = heightfunction(state)
    emitters = max(height)
    return emitters


def photonic_circuit_solver(state: Stabilizer) -> list[tuple]:
    '''
    INITIAL STEP: determine the number of emitters and initialize variables
    '''
    n_e = num_emitters(state)
    n_p = state.size
    N = n_e+n_p
    target_state = Stabilizer(N)
    for i in range(n_p):
        target_state.signvector[i] = state.signvector[i]
        for j in range(n_p):
            target_state.tab[i, j] = state.tab[i, j]
            target_state.tab[i, j+N] = state.tab[i, j+n_p]
    protocol = []

    '''loop through the photons in reverse order'''
    for h in range(n_p, 0, -1):
        '''
        STEP (I): transform the stabilizers into the echelon gauge and compute the height function
        '''
        height = heightfunction(target_state)
        photonindex = h-1
        d = height[h]-height[h-1]

        '''
        STEP (II): if h(j) < h(j-1), do a time reversed measurement
        '''
        if d < 0:
            # find stabilizers that act on no photons
            rows = []
            for i in range(N):
                toggler = True
                for j in range(n_p):
                    if target_state.tab[i, j] != 0 or target_state.tab[i, j+N] != 0:
                        toggler = False
                        break
                if toggler:
                    rows.append(i)

            # find number of emitters acted on for each stabilizer
            sums = []
            for row in rows:
                sum = 0
                for j in range(n_p, N):
                    if target_state.tab[row, j] != 0 or target_state.tab[row, j+N] != 0:
                        sum += 1
                sums.append(sum)

            # find stabilizer (row) which acts on minimum emitters
            row = rows[sums.index(min(sums))]

            # find emitter (emit) in row acted on and turn its pauli to Z
            for i in range(n_p, N):
                if target_state.tab[row, i] != 0 or target_state.tab[row, i+N] != 0:
                    emit = i
                    if target_state.tab[row, i] == 1 and target_state.tab[row, i+N] == 0:
                        protocol.append(['H', i])
                        target_state.clifford('H', i)
                    elif target_state.tab[row, i] == 1 and target_state.tab[row, i+N] == 1:
                        protocol.append(['S', i])
                        target_state.clifford('S', i)
                        protocol.append(['H', i])
                        target_state.clifford('H', i)
                    break

            # turn all other emitter paulis in row to Z and CNOT with emit to remove them
            for i in range(emit+1, N):
                if target_state.tab[row, i] != 0 or target_state.tab[row, i+N] != 0:
                    if target_state.tab[row, i] == 1 and target_state.tab[row, i+N] == 0:
                        protocol.append(['H', i])
                        target_state.clifford('H', i)
                    elif target_state.tab[row, i] == 1 and target_state.tab[row, i+N] == 1:
                        protocol.append(['S', i])
                        target_state.clifford('S', i)
                        protocol.append(['H', i])
                        target_state.clifford('H', i)
                    protocol.append(['CNOT', i, emit])
                    target_state.clifford('CNOT', i, emit)

            # correct if row is negative
            if target_state.signvector[row] == 1:
                target_state.clifford('X', emit)
                protocol.append(['X', emit])

            # perform time reversed measurement
            target_state.clifford('H', emit)
            target_state.clifford('CNOT', emit, photonindex)
            protocol.append(['Measure', emit, photonindex])

        # transform the stabilizers into the echelon gauge
        rref(target_state)

        '''
        STEP (III): do a photon absorption
        '''
        # find stabilizer (row) that acts on no photons but photonindex
        for i in range(N):
            toggler = True
            if target_state.tab[i, photonindex] == 0 and target_state.tab[i, photonindex+N] == 0:
                toggler = False
            if toggler:
                for j in range(photonindex):
                    if target_state.tab[i, j] != 0 or target_state.tab[i, j+N] != 0:
                        toggler = False
                        break
            if toggler:
                row = i
                break

        emit = -1
        # turn photonindex's pauli in row to Z
        if target_state.tab[row, photonindex] == 1 and target_state.tab[row, photonindex+N] == 0:
            protocol.append(['H', photonindex])
            target_state.clifford('H', photonindex)
        elif target_state.tab[row, photonindex] == 1 and target_state.tab[row, photonindex+N] == 1:
            protocol.append(['S', photonindex])
            target_state.clifford('S', photonindex)
            protocol.append(['H', photonindex])
            target_state.clifford('H', photonindex)

        # find emitter (emit) acted on by row and turn its pauli to Z
        for i in range(n_p, N):
            if target_state.tab[row, i] != 0 or target_state.tab[row, i+N] != 0:
                emit = i
                if target_state.tab[row, i] == 1 and target_state.tab[row, i+N] == 0:
                    protocol.append(['H', i])
                    target_state.clifford('H', i)
                elif target_state.tab[row, i] == 1 and target_state.tab[row, i+N] == 1:
                    protocol.append(['S', i])
                    target_state.clifford('S', i)
                    protocol.append(['H', i])
                    target_state.clifford('H', i)
                break

        if emit != -1:
            # turn all other emitter paulis in row to Z and CNOT with emit to remove them
            for i in range(emit+1, N):
                if target_state.tab[row, i] != 0 or target_state.tab[row, i+N] != 0:
                    if target_state.tab[row, i] == 1 and target_state.tab[row, i+N] == 0:
                        protocol.append(['H', i])
                        target_state.clifford('H', i)
                    elif target_state.tab[row, i] == 1 and target_state.tab[row, i+N] == 1:
                        protocol.append(['S', i])
                        target_state.clifford('S', i)
                        protocol.append(['H', i])
                        target_state.clifford('H', i)
                    protocol.append(['CNOT', i, emit])
                    target_state.clifford('CNOT', i, emit)

            # correct if row is negative
            if target_state.signvector[row] == 1:
                target_state.clifford('X', emit)
                protocol.append(['X', emit])

            # perform photon absorption
            target_state.clifford('CNOT', emit, photonindex)
            protocol.append(['Emission', emit, photonindex])
        else:
            # correct if row is negative
            if target_state.signvector[row] == 1:
                target_state.clifford('X', photonindex)
                protocol.append(['X', photonindex])

        # remove photonindex from other stabilizers using row
        for i in range(N):
            if target_state.tab[i, photonindex+N] != 0 or target_state.tab[i, photonindex] != 0:
                if i != row:
                    target_state.row_add(row, i)

    # transform the stabilizers into the echelon gauge
    rref(target_state)

    '''
    STEP (IV): decouple and reset the emitters
    '''
    # for each emitter stabilizer find the number of emitters acted on
    sums = []
    for i in range(n_p, N):
        sum = 0
        for j in range(n_p, N):
            if target_state.tab[i, j] != 0 or target_state.tab[i, j+N] != 0:
                sum += 1
        sums.append(sum)

    # while the emitters are coupled
    if max(sums) == 1:
        decoupled = True
    else:
        decoupled = False
    while not decoupled:
        # find stabilizer (row) that acts on the fewest emitters (at least 2)
        for i in range(2, N):
            if i in sums:
                minimum = i
                break
        row = n_p+sums.index(minimum)

        # find emitter (emit) acted on by row and turn its pauli to Z
        for i in range(n_p, N):
            if target_state.tab[row, i] != 0 or target_state.tab[row, i+N] != 0:
                emit = i
                if target_state.tab[row, i] == 1 and target_state.tab[row, i+N] == 0:
                    protocol.append(['H', emit])
                    target_state.clifford('H', emit)
                elif target_state.tab[row, i] == 1 and target_state.tab[row, i+N] == 1:
                    protocol.append(['S', emit])
                    target_state.clifford('S', emit)
                    protocol.append(['H', emit])
                    target_state.clifford('H', emit)
                break

        # turn all other paulis in row to Z and CNOT with emit to remove them
        for i in range(emit+1, N):
            if target_state.tab[row, i] != 0 or target_state.tab[row, i+N] != 0:
                if target_state.tab[row, i] == 1 and target_state.tab[row, i+N] == 0:
                    protocol.append(['H', i])
                    target_state.clifford('H', i)
                elif target_state.tab[row, i] == 1 and target_state.tab[row, i+N] == 1:
                    protocol.append(['S', i])
                    target_state.clifford('S', i)
                    protocol.append(['H', i])
                    target_state.clifford('H', i)
                target_state.clifford('CNOT', i, emit)
                protocol.append(['CNOT', i, emit])

        # remove emit from other stabilizers using row
        for i in range(n_p, N):
            if target_state.tab[i, emit] != 0 or target_state.tab[i, emit+N] != 0:
                if i != row:
                    target_state.row_add(row, i)

        # for each emitter stabilizer find the number of emitters acted on
        sums = []
        for i in range(n_p, N):
            sum = 0
            for j in range(n_p, N):
                if target_state.tab[i, j] != 0 or target_state.tab[i, j+N] != 0:
                    sum += 1
            sums.append(sum)

        # check if emitters are coupled
        if max(sums) == 1:
            decoupled = True
        else:
            decoupled = False

    rref(target_state)

    # for each emitter stabilizer turn its pauli to Z
    for i in range(n_p, N):
        if target_state.tab[i, i] != 0:
            if target_state.tab[i, i] == 1 and target_state.tab[i, i+N] == 0:
                protocol.append(['H', i])
                target_state.clifford('H', i)
            elif target_state.tab[i, i] == 1 and target_state.tab[i, i+N] == 1:
                protocol.append(['S', i])
                target_state.clifford('S', i)
                protocol.append(['H', i])
                target_state.clifford('H', i)

    # correct any negative stabilizers
    for i in range(n_p, N):
        if target_state.signvector[i] == 1:
            target_state.clifford('X', i)
            protocol.append(['X', i])

    protocol.reverse()

    # confirm that target_state is now fully decoupled and zeroed out
    checker_state = Stabilizer(N)
    if np.array_equal(checker_state.tab, target_state.tab) and np.array_equal(checker_state.signvector, target_state.signvector):
        return protocol
    else:
        print('Something went wrong')
        return None


def emitter_cnot(state: Stabilizer) -> list[int]:
    emitter = num_emitters(state)
    procedure = photonic_circuit_solver(state)
    cnots = 0
    for i in range(len(procedure)):
        if procedure[i][0] == 'CNOT':
            cnots += 1
    data = [emitter, cnots]
    return data


def remove_sign(stabs: list[str]) -> list[str]:
    for i in range(len(stabs)):
        stabs[i] = stabs[i].lstrip('-')
    return stabs


def protocol_to_circuit(protocol: list[tuple]) -> QuantumCircuit:
    n_qubits = 0
    for step in protocol:
        for qubit in range(1, len(step)):
            if step[qubit] > n_qubits:
                n_qubits = step[qubit]
    n_qubits += 1

    qc = QuantumCircuit(n_qubits, 1)

    for step in protocol:
        if step[0].lower() == 'h':
            qc.h(step[1])
        elif step[0].lower() == 'cnot':
            qc.cx(step[1], step[2])
        elif step[0].lower() == 'z':
            qc.z(step[1])
        elif step[0].lower() == 'x':
            qc.x(step[1])
        elif step[0].lower() == 'y':
            qc.y(step[1])
        elif step[0].lower() == 's':
            qc.s(step[1])
        elif step[0].lower() == 'cz':
            qc.cz(step[1], step[2])
        elif step[0].lower() == 'emission':
            qc.barrier()
            qc.cx(step[1], step[2])
        elif step[0].lower() == 'measure':
            qc.measure(step[1], 0)
            with qc.if_test((0, 1)):
                qc.x(step[1])
                qc.x(step[2])
        else:
            raise ValueError("Unexpected operation", step[0].lower)

    return qc


def protocol_to_state(protocol: list[tuple], progress=False, check_measurements=True) -> Expr:
    n_qubits = 0
    for step in protocol:
        for qubit in range(1, len(step)):
            if step[qubit] > n_qubits:
                n_qubits = step[qubit]
    n_qubits += 1

    state = Qubit('0'*n_qubits)

    for i in range(len(protocol)):
        if progress:
            clear_output(wait=True)
            print("Generating state: " + str(round(100*i/len(protocol))) + "%")
        step = protocol[i]
        if step[0].lower() == 'h':
            state = apply(state, operation=H(step[1]))
        elif step[0].lower() == 'cnot':
            state = apply(state, operation=CNOT(step[1], step[2]))
        elif step[0].lower() == 'z':
            state = apply(state, operation=Z(step[1]))
        elif step[0].lower() == 'x':
            state = apply(state, operation=X(step[1]))
        elif step[0].lower() == 'y':
            state = apply(state, operation=Y(step[1]))
        elif step[0].lower() == 's':
            state = apply(state, operation=S(step[1]))
        elif step[0].lower() == 'cz':
            state = apply(state, operation=CPHASE(step[1], step[2]))
        elif step[0].lower() == 'emission':
            state = apply(state, operation=CNOT(step[1], step[2]))
        elif step[0].lower() == 'measure':
            result = measure_partial(
                state, step[1])
            state0 = result[0][0]

            if check_measurements:
                state1 = result[1][0]

                state1 = apply(state1, X(step[1]))
                state1 = apply(state1, X(step[2]))

                if apply(state0) != apply(state1):
                    raise ValueError(
                        "Expected output state to be the same regardless of measurement")

            state = state0
        else:
            raise ValueError("Unexpected operation", step[0].lower)

    if progress:
        clear_output(wait=True)

    return apply(state, dp=True)


def commute_emissions(protocol: list) -> None:
    photons = []
    for i in range(len(protocol)):
        step_iter = protocol[i]
        if step_iter[0] == "Emission":
            photons.append(step_iter[2])
        elif step_iter[0] == "X" and step_iter[1] in photons:
            for j in range(len(protocol)):
                step_search = protocol[j]
                if step_search[0] == "Emission" and step_search[2] == step_iter[1]:
                    step_move = protocol.pop(i)
                    protocol.insert(j, step_move)
                    break


def photon_permuter(graph: Graph, progress=False) -> tuple[list, int]:
    state = Stabilizer(graph=graph)
    min_emitters = num_emitters(state)
    best_edges = graph.edges.copy()

    perms = list(permutations(range(graph.n_qubits)))
    for i in range(len(perms)):
        if progress:
            clear_output(wait=True)
            print("Permuting photons: " + str(round(100*i/len(perms))) + "%")
        perm = perms[i]
        new_edges = []
        for edge in graph.edges:
            new_edge0 = perm[edge[0]]
            new_edge1 = perm[edge[1]]
            new_edges.append((new_edge0, new_edge1))
        new_state = Stabilizer(graph=Graph(new_edges))
        new_emitters = num_emitters(new_state)
        if new_emitters < min_emitters:
            min_emitters = new_emitters
            best_edges = new_edges

    if progress:
        clear_output(wait=True)

    return best_edges, min_emitters
