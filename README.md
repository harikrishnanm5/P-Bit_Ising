# ⚛️ Reconfigurable 4-p-Bit CMOS Ising Optimization Accelerator

[![eSim](https://img.shields.io/badge/EDA-eSim%20v2.4-blue.svg)](https://esim.fossee.in/)
[![Ngspice](https://img.shields.io/badge/Simulator-Ngspice%2034%2B-green.svg)](https://ngspice.sourceforge.io/)
[![Technology](https://img.shields.io/badge/Process-180nm%20CMOS-orange.svg)]()
[![License](https://img.shields.io/badge/License-MIT-purple.svg)]()

A room-temperature, quantum-inspired **Reconfigurable 4-p-Bit CMOS Ising Accelerator** designed and verified in **eSim** and **Ngspice** for solving NP-hard combinatorial optimization problems (specifically graph MAX-CUT partitioning) using stochastic CMOS probabilistic bits (p-bits).

---

## 📌 Overview

Combinatorial optimization problems—such as MAX-CUT, traveling salesperson, and protein folding—scale exponentially on traditional Von Neumann computing architectures. While quantum annealers solve these using cryogenic qubits at $-273^\circ\text{C}$ (milliKelvin temperatures), **probabilistic computers (p-computers)** emulate Ising spin physics at **room temperature** ($300\text{ K}$) using fluctuating, thermally-perturbed CMOS devices.

This repository provides complete, SPICE-accurate, transistor-level netlists and modern KiCad 6 schematics for:
1. **Discrete CMOS Probabilistic Bits (p-bits)** biased in the sub-threshold metastable regime.
2. **Tunable Transmission-Gate Synaptic Coupling Network** implementing reconfigurable antiferromagnetic interactions ($J_{ij}$).
3. **Dual Operating Modes**:
   - **4-p-Bit Mode (`V(MODE) = 1.8V`)**: Solves a 4-node benchmark graph with 5 active edges.
   - **2-p-Bit Mode (`V(MODE) = 0.0V`)**: Disables nodes P3 & P4, isolating an antiferromagnetic dimer pair ($P_1 \leftrightarrow P_2$).

---

## 📐 Mathematical Formulation

### 1. The Probabilistic Bit (p-bit)
Each p-bit produces a stochastic digital output $m_i \in \{-1, +1\}$ whose probability of being in state $+1$ is governed by the Boltzmann Sigmoid activation function:

$$P(m_i = +1) = \frac{1}{1 + \exp\left(-\frac{2 I_i}{I_0}\right)} = \frac{1}{2}\left(1 + \tanh\left(\frac{I_i}{I_0}\right)\right)$$

where $I_i$ is the synaptic input potential arriving at the $i$-th p-bit and $I_0$ is the noise/thermal scaling constant.

### 2. Ising Hamiltonian & MAX-CUT Mapping
The network energy corresponds to the classical Ising spin Hamiltonian:

$$H(\vec{s}) = -\sum_{i < j} J_{ij} s_i s_j$$

For a MAX-CUT problem on an undirected graph $G = (V, E)$, setting an antiferromagnetic coupling $J_{ij} < 0$ (e.g., $J = -1$) penalizes identical spins ($s_i = s_j$) and rewards opposing spins ($s_i \neq s_j$), driving the system toward states that maximize cut edges:

$$\text{Cut}(\vec{s}) = \sum_{(i,j) \in E} \frac{1 - s_i s_j}{2}$$

---

## 🕸️ Benchmark Graph Topology

We map a 4-node benchmark graph with 5 active antiferromagnetic edges:
* **Active Edges**: $(P_1 \leftrightarrow P_2)$, $(P_1 \leftrightarrow P_3)$, $(P_2 \leftrightarrow P_3)$, $(P_2 \leftrightarrow P_4)$, $(P_3 \leftrightarrow P_4)$
* **Null Edge**: $(P_1 \leftrightarrow P_4)$ (no direct connection)

```text
        (P1) ----------- (P2)
          |   \       /    |
          |     \   /      |
          |       X        |
          |     /   \      |
          |   /       \    |
        (P3) ----------- (P4)
```

### Optimal Ground State Solution:
To maximize the number of cut edges:
* **Partition A**: $\{P_1, P_4\} = \mathbf{+1} \text{ (High / 1.8V)}$
* **Partition B**: $\{P_2, P_3\} = \mathbf{-1} \text{ (Low / 0.0V)}$
* **Result**: **4 out of 5 edges cut** ($80\%$ cut efficiency, the theoretical optimum for this non-planar graph). The symmetric state $\{P_1=-1, P_4=-1, P_2=+1, P_3=+1\}$ has identical minimal energy.

---

## 🔌 Circuit Architecture

```text
               +-------------------------------------------------------+
               |                    V(MODE) = 1.8V                     |
               +-------------------------------------------------------+
                                           |
                                           v
       +-------------+       +----------------------------+       +-------------+
       |   p-Bit 1   |<=====>| Transmission Gate Synapses |<=====>|   p-Bit 2   |
       |  (Inverter  |       |     (R_c = 500 kOhm)       |       |  (Inverter  |
       |  + Rfb + C) |       +----------------------------+       |  + Rfb + C) |
       +-------------+                     ^                      +-------------+
              ^                            |                             ^
              |                            v                             |
              |              +----------------------------+              |
              +=============>| Reconfigurable Interconnect|<=============+
              |              +----------------------------+              |
              v                            ^                             v
       +-------------+                     |                      +-------------+
       |   p-Bit 3   |<====================+=====================>|   p-Bit 4   |
       +-------------+                                            +-------------+
```

1. **p-Bit Unit Cell (`p_bit.cir`)**:
   - Sub-threshold CMOS inverter pair (`W/L = 2.5u/0.18u` for PMOS, `1.0u/0.18u` for NMOS).
   - High-value resistive feedback ($R_{fb} = 500\text{ k}\Omega$) establishing the metastable point ($V_{in} \approx V_{out} \approx 0.90\text{ V}$).
   - Internal integration capacitor ($C = 15\text{ fF}$) and CMOS output sharpening buffer stage.
2. **Coupling Matrix**:
   - Bilateral CMOS transmission gates ($MN + MP$) gated by `MODE` and complementary `MODE_BAR`.
   - $500\text{ k}\Omega$ coupling resistors delivering current sum $I_i = \sum_{j} J_{ij} s_j$ onto each node.
3. **Noise Perturbation**:
   - Independent, incommensurate dither sources generating continuous quasi-random thermal perturbations.

---

## 📁 Repository Structure

| File | Description |
| :--- | :--- |
| `ising_4pbit_core.cir` | **Master 4-p-bit reconfigurable Ising core** (180nm CMOS netlist) |
| `p_bit_ising_machine.kicad_sch` | Full transistor-level schematic in modern KiCad 6+ / eSim format |
| `p_bit_ising_machine.proj` | eSim project workspace configuration |
| `sym-lib-table` | KiCad symbol library table |
| `get_output_values.py` | Automated transient verification & ground state evaluation script |
| `p_bit.cir` | Standalone single CMOS p-bit cell netlist |
| `p_bit_sigmoid_validation.cir`| Empirical validation testbench for the Boltzmann Sigmoid response |
| `step1_characterize_vth.cir` | DC transfer & threshold voltage ($V_{th}$) characterization netlist |
| `step3_analog_perceptron_and.cir` | Analog perceptron summation netlist (AND configuration) |
| `step5_single_pbit_dither.cir` | Single p-bit transient testbench under noise perturbations |
| `step6_coupled_pbit_pair.cir` | Isolated 2-node coupled p-bit dimer netlist |
| `test_4pbit_core.cir` | Testbench for 4-pbit core transient dynamics |
| `test_mode_4pbit.cir` | Explicit testbench verifying 4-p-bit operational mode |
| `test_mode_2pbit.cir` | Explicit testbench verifying 2-p-bit operational mode |
| `test_analog_perceptron.cir` | Synaptic summation testbench |
| `test_verify.cir` | Ngspice `.meas` verification script for steady-state spin values |
| `.gitignore` | Ignores temporary simulation data dumps (`*.txt`, `temp_*.cir`, logs) |

---

## 🚀 Simulation & Quickstart

### Prerequisites
* **Ngspice** (Version 34 or later recommended)
* **Python 3.8+** with `numpy` (for automated evaluation)
* **eSim v2.4+** (optional, for GUI schematic inspection)

### 1. Run Automated Evaluation via Python
Run the bundled evaluation script to simulate the master core and print the sampled spin vectors and satisfied cuts:
```bash
python get_output_values.py
```

**Expected Console Output**:
```text
=================================================================
      4-p-BIT ISING CORE TRANSIENT OUTPUT VALUES
=================================================================
Time (ms)  | P1 (V)   | P2 (V)   | P3 (V)   | P4 (V)   | Spin Vector (s1,s2,s3,s4) 
-----------------------------------------------------------------
0.10       | 1.782    | 0.015    | 0.018    | 1.780    | (+1, -1, -1, +1)          
0.25       | 1.785    | 0.012    | 0.014    | 1.784    | (+1, -1, -1, +1)          
0.50       | 0.014    | 1.783    | 1.781    | 0.016    | (-1, +1, +1, -1)          
0.75       | 0.013    | 1.784    | 1.782    | 0.015    | (-1, +1, +1, -1)          
1.00       | 0.014    | 1.783    | 1.781    | 0.016    | (-1, +1, +1, -1)          
-----------------------------------------------------------------
Final State at t = 1.00 ms: Lowest Energy Ground State s = (-1, +1, +1, -1)
Satisfied Cut Edges: (1,2), (1,3), (2,4), (3,4) [4 out of 5 edges cut]
```

### 2. Run Direct Ngspice Simulation
To observe raw analog metastable voltages and digital spin transitions:
```bash
# Simulate master 4-p-bit core
ngspice ising_4pbit_core.cir

# Validate single p-bit sigmoid distribution
ngspice p_bit_sigmoid_validation.cir

# Verify 2-p-bit vs 4-p-bit operational modes
ngspice test_mode_2pbit.cir
ngspice test_mode_4pbit.cir
```

### 3. Open in eSim GUI
1. Launch **eSim**.
2. Select **Open Project** and navigate to this repository directory.
3. Open `p_bit_ising_machine.proj` to view and edit `p_bit_ising_machine.kicad_sch` in the schematic capture editor.

---

## 📊 Key Circuit Specifications

| Parameter | Value | Details |
| :--- | :--- | :--- |
| **CMOS Technology** | 180nm Bulk CMOS | $L_{eff} = 180\text{ nm}$ |
| **Supply Rail ($V_{DD}$)** | $1.8\text{ V}$ | Low-voltage operation |
| **DC Bias Operating Point** | $0.90\text{ V}$ | Calibrated inverter switching threshold ($V_{th}$) |
| **Coupling Resistance ($R_c$)** | $500\text{ k}\Omega$ | Symmetric bilateral pass gate loading |
| **Metastable Feedback ($R_{fb}$)**| $500\text{ k}\Omega$ | Promotes rapid sub-threshold noise amplification |
| **Coupling Topology** | 4-Node, 5-Edge | Hardware-reconfigurable via CMOS transmission gates |
| **Transient Convergence** | $< 500\text{ }\mu\text{s}$ | Fast relaxation to minimal Ising energy ground state |

---

## 📜 License
This project is licensed under the [MIT License](LICENSE).
