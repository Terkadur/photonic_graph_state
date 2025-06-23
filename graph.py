import networkx as nx
import matplotlib.pyplot as plt
from state import *
from IPython.display import clear_output

class Graph:
    '''
    This is a class that encodes graphs, and contains a few convenient functions

    :param edges: A list of edges describing the graph
    :type edges: list
    '''

    def __init__(self, edges: list[tuple]) -> None:
        n_qubits = 0
        for edge in edges:
            if edge[0] > n_qubits:
                n_qubits = edge[0]
            if edge[1] > n_qubits:
                n_qubits = edge[1]
        n_qubits += 1

        qubits = range(n_qubits)
        for edge in edges:
            if edge[0] not in qubits:
                raise Exception("Unexpected qubit: " + str(edge[0]))
            if edge[1] not in qubits:
                raise Exception("Unexpected qubit: " + str(edge[1]))

        self.n_qubits = n_qubits
        self.edges = edges


    def draw(self, figsize=(4, 4)) -> None:
        G = nx.Graph()
        G.add_edges_from(self.edges)
        plt.figure(figsize=figsize)
        nx.draw(G, with_labels=True, font_weight='bold', node_size=700,
                node_color="cornflowerblue", width=2.0, linewidths=4, edgecolors="royalblue")
        plt.show()


    def state(self, progress: bool = False) -> Expr:
        '''
        Uses SymPy to calculate the statevector associated with the graph

        :param progress: A boolean option to continuously print out the progress of calculating the state (useful for large graphs), defaults to False
        :type progress: bool

        :return: A SymPy expression representing the statevector
        :rtype: Expr
        '''
        state = Qubit('0'*self.n_qubits)

        photons = []
        for edge in self.edges:
            if edge[0] not in photons:
                photons.append(edge[0])
            if edge[1] not in photons:
                photons.append(edge[1])

        for i in range(len(photons)):
            if progress:
                clear_output(wait=True)
                print("Applying Hadamards: " + str(round(100*i/len(photons))) + "%")
            photon = photons[i]
            state = apply(state, H(photon))

        for i in range(len(self.edges)):
            edge = self.edges[i]
            if progress:
                clear_output(wait=True)
                print("Applying CZs: " + str(round(100*i/len(self.edges))) + "%")
            state = apply(state, CPHASE(edge[0], edge[1]))

        if progress:
            clear_output(wait=True)

        return state