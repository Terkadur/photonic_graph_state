# **Photonic Graph States**
### **Adapted by Tarek Razzaz from [Nishad Manohar's Repository](https://github.com/nrmanohar/python_photonic_graph_state)**

- Created a new class, Graph in graph.py for organizational purposes.
- Added protocol_to_circuit() in protocol.py to convert protocols to Qiskit QuantumCircuits.
- Added Graph.state() in graph.py and protocol_to_state() in protocol.py (and supporting functions in state.py) which use SymPy to compute the expected and generated statevectors.
- Fixed row_add() in stabilizer.py, it didn't account for the sign of row1 and used depracated np.complex().
- Modified photonic_circuit_solver() in protocol.py so that stabilizer signs would be corrected during step (iii), photon absorption. This way emission is always the first operation for each photon.
- Added photon_permutor() in protocol.py to permute the labelling of photons in a graph state to minimize the number of required emitters.
- Demonstrated and confirmed photon generation protocols for the examples in [Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41534-022-00522-6/MediaObjects/41534_2022_522_MOESM1_ESM.pdf) in testing.ipynb.
- Tweaked the order of searching through stabilizers during absorption to hopefully shorten circuits.