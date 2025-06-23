# **Photonic Graph States**
### **Adapted by Tarek Razzaz from [Nishad Manohar's Repository](https://github.com/nrmanohar/python_photonic_graph_state)**

- Created a new class, Graph in graph.py for organizational purposes.
- Added protocol_to_circuit() in protocol.py to convert protocols to Qiskit QuantumCircuits.
- Added Graph.state() in graph.py and protocol_to_state() in protocol.py (and supporting functions in state.py) which use SymPy to compute the expected and generated statevectors.
- Fixed row_add() in stabilizer.py, it didn't account for the sign of row1 and used depracated np.complex().
- Modified photonic_circuit_solver() in protocol.py so that stabilizer signs would be corrected during step (iii), photon absorption. This way emission is always the first operation for each photon.
- Added photon_permutor() in protocol.py to permute the labelling of photons in a graph state to minimize the number of required emitters.
- Demonstrated and confirmed photon generation protocols in testing.ipynb.