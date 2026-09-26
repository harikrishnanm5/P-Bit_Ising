# ⚛️ 4-p-Bit Reconfigurable Ising Optimization Core

This project implements a **4-p-Bit Reconfigurable Ising Optimization Core** in eSim and Ngspice for solving combinatorial optimization (MAX-CUT Graph Partitioning) at room temperature using behavioral stochastic circuits.

---

## 📐 Graph & Coupling Topology (4-Node MAX-CUT)
We map a 4-node benchmark graph with 5 active inhibitory (antiferromagnetic) edges:
* **Active Edges**: `(P1 ↔ P2)`, `(P1 ↔ P3)`, `(P2 ↔ P3)`, `(P2 ↔ P4)`, `(P3 ↔ P4)`
* **Excluded Edge**: `(P1 ↔ P4)` (NO direct connection)

```
        (P1) ----------- (P2)
          |   \       /    |
          |     \   /      |
          |       X        |
          |     /   \      |
          |   /       \    |
        (P3) ----------- (P4)
```

### Optimal Ground State Solution:
To maximize the number of cut edges (MAX-CUT):
* **Partition A**: $\{P1, P4\} = \mathbf{5.0\text{ V} \text{ (Logic 1)}}$
* **Partition B**: $\{P2, P3\} = \mathbf{0.0\text{ V} \text{ (Logic 0)}}$
* **Result**: **4 out of 5 edges cut** (Maximum possible cut for this topology). The complementary state $\{P1=0V, P4=0V, P2=5V, P3=5V\}$ is equally optimal.

---

## 🎛️ Reconfigurable Operating Modes
* **4-p-Bit Mode (`V(MODE) = 5.0V`)**:
  All 4 p-bits (P1..P4) and all 5 synaptic couplings are active.
* **2-p-Bit Mode (`V(MODE) = 0.0V`)**:
  P3 and P4 are disabled ($0\text{V}$) and their coupling terms into P1 and P2 are zeroed. The system operates as an isolated 2-node antiferromagnetically coupled pair ($P1 \leftrightarrow P2$).

---

## 🚀 How to Run the Simulation

### Running the 4-p-Bit Ising Core:
```bash
cd d:\E-Sim-A-Ton\p_bit_ising_machine
ngspice ising_4pbit_core.cir
```

### Running the Empirical Sigmoid Testbench:
```bash
ngspice p_bit_sigmoid_validation.cir
```

---

## 📊 Key Circuit Specifications
| Parameter | Value | Description |
| :--- | :--- | :--- |
| **Supply Voltage ($V_{DD}$)** | $5.0\text{ V}$ | Standard logic rail |
| **State Encoding** | $0.0\text{ V} \leftrightarrow 5.0\text{ V}$ | Normalized bipolar spin $s_i \in [-1, +1]$ |
| **RC Filter Delay ($\tau$)** | $10\text{ }\mu\text{s}$ | $R = 1\text{ k}\Omega, C = 10\text{ nF}$ (breaks algebraic loops) |
| **Noise Sources** | SFFM Multi-Tone | 4 out-of-phase asynchronous quasi-random generators |
| **Simulation Window** | $0 \rightarrow 200\text{ }\mu\text{s}$ | Step size $\Delta t = 100\text{ ns}$ |
