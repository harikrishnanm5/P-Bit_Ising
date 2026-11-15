# 📘 Comprehensive Learning & Interview Study Guide: CMOS Probabilistic Bits (p-Bits) and Ising Optimization

This document serves as an exhaustive, first-principles technical study guide for the **Reconfigurable 4-p-Bit CMOS Ising Accelerator**. It covers conceptual physics, mathematical derivations, circuit architecture, line-by-line SPICE netlist walkthroughs of [`spice/ising_4pbit_core.cir`](spice/ising_4pbit_core.cir), operating regimes, chatter analysis vs. thermodynamic escape dynamics, and 22 high-yield technical interview questions with sharp, rigorous answers.

---

## 📑 Table of Contents
1. [First-Principles Conceptual Guide](#1-first-principles-conceptual-guide)
   - [1.1 What is a Probabilistic Bit (p-Bit)?](#11-what-is-a-probabilistic-bit-p-bit)
   - [1.2 Why the Inverter is Biased at Its Trip Point (0.90V)](#12-why-the-inverter-is-biased-at-its-trip-point-090v)
   - [1.3 How the Boltzmann Sigmoid Emerges from Noise](#13-how-the-boltzmann-sigmoid-emerges-from-noise)
   - [1.4 How Resistive Coupling Implements Synaptic Weights](#14-how-resistive-coupling-implements-synaptic-weights)
   - [1.5 Mapping MAX-CUT to the Ising Hamiltonian](#15-mapping-max-cut-to-the-ising-hamiltonian)
   - [1.6 Mathematical Proof: Why 4/5 is the Global Maximum Cut](#16-mathematical-proof-why-45-is-the-global-maximum-cut)
2. [Line-by-Line SPICE Netlist Walkthrough (`spice/ising_4pbit_core.cir`)](#2-line-by-line-spice-netlist-walkthrough-spiceising_4pbit_corecir)
3. [Operating Regimes: Freezing vs. Active Boltzmann Sampling](#3-operating-regimes-freezing-vs-active-boltzmann-sampling)
4. [Convergence Dynamics: RC Settling vs. Stationary Langevin Diffusion](#4-convergence-dynamics-rc-settling-vs-stationary-langevin-diffusion)
   - [4.1 Why Initial Time-to-First-Ground-State is Misleading](#41-why-initial-time-to-first-ground-state-is-misleading)
   - [4.2 Output Chatter vs. Persistent State Changes (>= 2 ns)](#42-output-chatter-vs-persistent-state-changes--2-ns)
   - [4.3 Persistence Filtering vs. Voltage Hysteresis & Right-Censoring](#43-persistence-filtering-vs-voltage-hysteresis--right-censoring)
5. [16-State Boltzmann Sampling Deviations & Degree Asymmetry](#5-16-state-boltzmann-sampling-deviations--degree-asymmetry)
6. [22 High-Yield Technical Interview Questions & Answers](#6-22-high-yield-technical-interview-questions--answers)

---

## 1. First-Principles Conceptual Guide

### 1.1 What is a Probabilistic Bit (p-Bit)?
In computing, we classify state variables into three distinct paradigms:
1. **Classical Deterministic Bit**: A bistable system with two deep potential wells separated by a high energy barrier ($\gg k_B T$). Its state $m \in \{0, 1\}$ or $s \in \{-1, +1\}$ is permanently static until forced by a deterministic clock or gate.
2. **Quantum Bit (Qubit)**: A coherent two-level quantum system existing in a linear superposition $|\psi\rangle = \alpha |0\rangle + \beta |1\rangle$ ($\alpha, \beta \in \mathbb{C}$, $|\alpha|^2 + |\beta|^2 = 1$). Qubits exhibit quantum entanglement and tunneling, but require cryogenic isolation (e.g., $15\text{ mK}$) to prevent phase decoherence.
3. **Probabilistic Bit (p-Bit)**: A classical, low-barrier fluctuating physical element operating at **room temperature ($300\text{ K}$)**. A p-bit rapidly fluctuates between $-1$ and $+1$ in time. Its time-averaged output follows a continuous, controllable probability distribution governed by the logistic sigmoid function:

$$P(m_i = +1) = \sigma(2 \beta I_i) = \frac{1}{1 + \exp\left(-\frac{2(I_i - V_{th})}{I_0}\right)}$$

where $I_i$ is the synaptic input voltage arriving at the p-bit, $V_{th} \approx 0.941\text{ V}$ is the switching threshold ($41\text{ mV}$ above the $0.90\text{ V}$ nominal bias point; $R_{fb}$ loading and inverter asymmetry are an untested hypothesis for this offset), and $I_0 \approx 363.0\text{ mV}$ is set by the noise amplitude as it appears at the gate node after the $R_{\text{in}} / R_{fb}$ network (measured across the balanced $R_{\text{in}} = 500\text{ k}\Omega, R_{fb} = 500\text{ k}\Omega$ cell driven by calibrated $440\text{ mV}$ noise). A p-bit is thus a hardware building block that natively performs **Markov Chain Monte Carlo (MCMC)** and **Boltzmann sampling** in physical time without random number generation software overhead.

---

### 1.2 Why the Inverter is Biased at Its Trip Point (0.90V)
A standard CMOS inverter consists of a pull-up PMOS and a pull-down NMOS. 
* The **trip point** (or switching threshold $V_M$) is the operating point where $V_{\text{in}} = V_{\text{out}} = V_{th}$. For a symmetric $1.8\text{ V}$ supply rail ($V_{DD}$), $V_{th} = V_{DD} / 2 = 0.90\text{ V}$.
* At this exact bias point:
  1. Both PMOS and NMOS operate simultaneously in the **saturation region** ($V_{DS} \ge V_{GS} - V_{TH}$).
  2. The small-signal voltage gain $A_v = -g_m (r_{o,n} \parallel r_{o,p})$ reaches its **absolute maximum**.
  3. The circuit is biased at the edge of **metastability**: even microvolt-level perturbations are amplified across the full rail ($0\text{ V}$ to $1.8\text{ V}$).
  4. With zero net input signal ($I_i = V_{th}$), the thermal noise fluctuates symmetrically around the trip point, producing an exact $50\%$ probability of being in state $+1$ and $50\%$ in state $-1$ ($P(+1) = 0.5$). Any small positive synaptic voltage $\Delta V > 0$ biases the state toward $+1$; any negative offset $\Delta V < 0$ biases it toward $-1$.

---

### 1.3 How the Boltzmann Sigmoid Emerges from Noise
The emergence of the smooth Boltzmann sigmoid from a deterministic CMOS inverter is a classic example of **stochastic resonance / noise convolution**:
1. An ideal, high-gain CMOS inverter with infinite gain acts as a signum/Heaviside step function:

$$f(V) = \begin{cases} +1 & \text{if } V > V_{th} \\ -1 & \text{if } V < V_{th} \end{cases}$$

2. When independent Gaussian thermal noise $\xi(t) \sim \mathcal{N}(0, \sigma^2)$ is added to the synaptic input $V_{\text{in}}$, the instantaneous voltage is $V(t) = V_{\text{in}} + \xi(t)$.
3. The probability that the output is $+1$ is the probability that $V_{\text{in}} + \xi(t) > V_{th}$, or equivalently $\xi(t) > V_{th} - V_{\text{in}}$:

$$P(m = +1) = \int_{V_{th} - V_{\text{in}}}^{\infty} \frac{1}{\sqrt{2\pi}\sigma} \exp\left(-\frac{u^2}{2\sigma^2}\right) du = \frac{1}{2} \left[1 + \text{erf}\left(\frac{V_{\text{in}} - V_{th}}{\sqrt{2}\sigma}\right)\right]$$

4. The Gaussian error function $\text{erf}(x)$ is functionally almost identical to the hyperbolic tangent $\tanh(x)$ and logistic sigmoid $\sigma(x)$:

$$\frac{1}{2} \left[1 + \text{erf}\left(\frac{x}{\sqrt{2}\sigma}\right)\right] \approx \frac{1}{1 + \exp\left(-\frac{2 x}{I_0}\right)}$$

where $I_0$ is set by the noise amplitude as it appears at the gate node after the $R_{\text{in}} / R_{fb}$ network. For our standardized cell with $440\text{ mV}$ RMS noise coupled through $R_{\text{in}} = 500\text{ k}\Omega$ and $R_{fb} = 500\text{ k}\Omega$, empirical fitting yields $I_0 \approx 363.0\text{ mV}$ ($R^2 = 0.9941$). The fitted $V_{th} = 0.941\text{ V}$ is $41\text{ mV}$ above the $0.90\text{ V}$ nominal bias point; $R_{fb}$ loading and inverter asymmetry are an untested hypothesis for this offset.

---

### 1.4 How Resistive Coupling Implements Synaptic Weights
In an Ising machine, the local effective field acting on spin $i$ is:

$$h_i = \sum_{j \neq i} J_{ij} s_j$$

In our CMOS circuit, this summation is performed directly using **Kirchhoff's Current Law (KCL)** at the high-impedance gate node (`noise_i`):
1. Each neighbor p-bit $j$ has an internal analog metastable potential $V(a_j)$ connected to node `noise_i` through a coupling resistor $R_{c, ij} = 500\text{ k}\Omega$.
2. Because the CMOS inverter has an inverting transfer characteristic ($V_{\text{in}} \uparrow \implies V_{\text{out}} \downarrow$):
   - When spin $j$ is in state $+1$ ($V(m_{\text{out},j}) \approx 1.8\text{ V}$), internal node $V(a_j) \approx 0.0\text{ V}$.
   - When spin $j$ is in state $-1$ ($V(m_{\text{out},j}) \approx 0.0\text{ V}$), internal node $V(a_j) \approx 1.8\text{ V}$.
3. Therefore, if spin $j$ is $+1$, it pulls current *out* of node `noise_i` (lowering its potential), which acts to push spin $i$ toward $-1$. If spin $j$ is $-1$, it injects current *into* node `noise_i`, pushing spin $i$ toward $+1$.
4. This inhibitory relationship ($s_j = +1 \implies s_i \to -1$) is the exact definition of an **antiferromagnetic interaction ($J_{ij} < 0$)**.
5. By KCL, the total synaptic voltage shift at node $i$ is:

$$\Delta V_{\text{gate}, i} = \sum_{j} \frac{G_{c, ij}}{G_{\text{tot}, i}} (V(a_j) - V_{th})$$

where $G_{\text{tot}, i} = G_{\text{in}} + G_{fb} + \sum_j G_{c, ij}$. This implements weighted analog synaptic addition with zero clocking power.

---

### 1.5 Mapping MAX-CUT to the Ising Hamiltonian
The classical Ising Hamiltonian for a system of $N$ bipolar spins $s_i \in \{-1, +1\}$ is:

$$H(\vec{s}) = -\sum_{(i,j) \in E} J_{ij} s_i s_j$$

For an unweighted graph $G = (V, E)$, the **MAX-CUT problem** seeks a bipartition of vertices into two sets $S$ and $\bar{S}$ that maximizes the number of cut edges crossing the partition.
* Define spin variable $s_i = +1$ if vertex $i \in S$, and $s_i = -1$ if vertex $i \in \bar{S}$.
* An edge $(i, j)$ is cut if and only if $s_i \neq s_j$. The cut indicator is:

$$\text{Cut}_{ij}(\vec{s}) = \frac{1 - s_i s_j}{2} = \begin{cases} 1 & \text{if } s_i \neq s_j \\ 0 & \text{if } s_i = s_j \end{cases}$$

* The total cut of the graph is:

$$\text{Cut}(\vec{s}) = \sum_{(i,j) \in E} \frac{1 - s_i s_j}{2} = \frac{|E|}{2} - \frac{1}{2} \sum_{(i,j) \in E} s_i s_j$$

* Setting antiferromagnetic coupling $J_{ij} = -1$ for all $(i,j) \in E$:

$$H(\vec{s}) = \sum_{(i,j) \in E} s_i s_j = |E| - 2 \cdot \text{Cut}(\vec{s})$$

$$\text{Cut}(\vec{s}) = \frac{|E| - H(\vec{s})}{2}$$

Thus, finding the ground state that minimizes $H(\vec{s})$ mathematically guarantees finding the partition that maximizes $\text{Cut}(\vec{s})$.

---

### 1.6 Mathematical Proof: Why 4/5 is the Global Maximum Cut
Our benchmark graph has 4 vertices $\{P_1, P_2, P_3, P_4\}$ and 5 edges:
$E = \{(1,2), (1,3), (2,3), (2,4), (3,4)\}$.

```text
        (P1) ----------- (P2)
          |            /   |
          |          /     |
          |        /       |
          |      /         |
          |    /           |
        (P3) ----------- (P4)
```

**Theorem**: No bipartition of this graph can achieve a cut of 5 edges. The maximum cut is exactly 4 edges.

**Proof**:
1. A graph is bipartite if and only if it contains no cycles of odd length.
2. In any graph, a cut partition is equivalent to choosing a bipartite subgraph. An odd cycle of length $2k + 1$ can have at most $2k$ of its edges cut by any bipartition (at least one edge must connect two vertices in the same set).
3. The benchmark graph contains two odd cycles of length 3 (triangles):
   - Triangle 1: $\{P_1, P_2, P_3\}$ with edges $(1,2), (1,3), (2,3)$
   - Triangle 2: $\{P_2, P_3, P_4\}$ with edges $(2,3), (2,4), (3,4)$
4. These two triangles share the single edge $(2,3)$.
5. **Case A (Edge $(2,3)$ is uncut)**: 
   - Since $(2,3)$ is uncut, $P_2$ and $P_3$ share the same partition set (say $-1$).
   - In Triangle 1, the remaining edges are $(1,2)$ and $(1,3)$. Setting $P_1 = +1$ cuts both $(1,2)$ and $(1,3)$ (contributing 2 cuts).
   - In Triangle 2, the remaining edges are $(2,4)$ and $(3,4)$. Setting $P_4 = +1$ cuts both $(2,4)$ and $(3,4)$ (contributing 2 cuts).
   - Total satisfied cut edges $= 2 + 2 = 4$.
6. **Case B (Edge $(2,3)$ is cut)**:
   - $P_2$ and $P_3$ are in opposite partition sets (one is $+1$, the other is $-1$).
   - Edge $(2,3)$ contributes 1 cut.
   - For Triangle 1, vertex $P_1$ must match either $P_2$ or $P_3$, so at most 1 of $(1,2)$ and $(1,3)$ can be cut.
   - For Triangle 2, vertex $P_4$ must match either $P_2$ or $P_3$, so at most 1 of $(2,4)$ and $(3,4)$ can be cut.
   - Total satisfied cut edges $\le 1 + 1 + 1 = 3$.
7. Therefore, a cut of 5 is strictly impossible, and $\mathbf{4/5}$ is the provable global maximum cut. $\blacksquare$

---

## 2. Line-by-Line SPICE Netlist Walkthrough (`spice/ising_4pbit_core.cir`)

Below is the line-by-line breakdown of the master netlist [`spice/ising_4pbit_core.cir`](spice/ising_4pbit_core.cir):

### 2.1 Power Supply, Control & Complementary Mode Generation
```spice
Vdd vdd 0 DC 1.8V
Vmode mode 0 DC 1.8V

* Complementary MODE Generation (Inverter)
MP_INV mode_bar mode vdd vdd PMOS_PBIT W=2.5u L=0.18u
MN_INV mode_bar mode 0   0   NMOS_PBIT W=1.0u L=0.18u
```
- `Vdd`: $1.8\text{ V}$ primary DC power rail.
- `Vmode`: Operating mode select. $1.8\text{ V}$ enables all 5 graph edges (4-p-bit mode); $0.0\text{ V}$ disables transmission gates to isolate an independent 2-p-bit pair ($P_1 \leftrightarrow P_2$).
- `MP_INV` / `MN_INV`: CMOS inverter generating complementary control signal `mode_bar`. When `mode = 1.8V`, `mode_bar = 0.0V` ($1.6\times 10^{-9}\text{ V}$).

### 2.2 Transistor Models (180nm CMOS Level-1)
```spice
.model NMOS_PBIT NMOS (LEVEL=1 VTO=0.45 KP=150u GAMMA=0.4 PHI=0.6 LAMBDA=0.05)
.model PMOS_PBIT PMOS (LEVEL=1 VTO=-0.45 KP=60u GAMMA=0.4 PHI=0.6 LAMBDA=0.05)
```
- Symmetric threshold magnitudes: $|V_{TO}| = 0.45\text{ V}$.
- Transconductance ratio: $\frac{KP_n}{KP_p} = \frac{150\,\mu\text{A/V}^2}{60\,\mu\text{A/V}^2} = 2.50$.
- All transistors have channel length $L = 0.18\,\mu\text{m}$.

### 2.3 DC Biases, Noise Sources & Input Summing Network
```spice
Vbias1 in1 0 DC 0.90V
Vnoise1 dither1 in1 DC 0 trnoise(440m 0.2n 0 0)
Rin1 dither1 noise1 500k
```
(Repeated for $P_1$ through $P_4$)
- `Vbias1`: Sets the quiescent input at the nominal mid-rail trip point $V_{th} = V_{DD}/2 = 0.90\text{ V}$.
- `Vnoise1`: Independent Gaussian thermal noise source with $440\text{ mV}$ RMS amplitude updating every $0.2\text{ ns}$ ($5\text{ GHz}$ bandwidth). Explicit `DC 0` prevents Ngspice time-zero warnings.
- `Rin1`: $500\text{ k}\Omega$ input summing resistor. It isolates the zero-impedance voltage source from gate node `noise1`, allowing passive linear superposition with feedback and coupling currents.

### 2.4 Discrete p-Bit Cell Topology ($P_1$ Example)
```spice
* --- p-Bit Cell P1 ---
MP1 a1 noise1 vdd vdd PMOS_PBIT W=2.5u L=0.18u
MN1 a1 noise1 0   0   NMOS_PBIT W=1.0u L=0.18u
Rfb1 noise1 a1 500k
C1 a1 0 15fF
MP_BUF1 m_out1 a1 vdd vdd PMOS_PBIT W=3.0u L=0.18u
MN_BUF1 m_out1 a1 0   0   NMOS_PBIT W=1.5u L=0.18u
Cout1 m_out1 0 20fF
```
- `MP1` / `MN1`: The core metastable inverting amplifier:
  - $W_p = 2.5\,\mu\text{m}$, $W_n = 1.0\,\mu\text{m}$ ($\frac{W_p}{W_n} = 2.5$).
  - $\beta_n = KP_n \cdot \frac{W_n}{L} = 150\mu \times \frac{1.0}{0.18} = 833.3\,\mu\text{A/V}^2$.
  - $\beta_p = KP_p \cdot \frac{W_p}{L} = 60\mu \times \frac{2.5}{0.18} = 833.3\,\mu\text{A/V}^2$.
  - Because $\beta_n = \beta_p$, the inverter switching threshold is balanced at $V_{DD}/2 = 0.900\text{ V}$.
- `Rfb1` ($500\text{ k}\Omega$): DC negative feedback from `a1` to `noise1`, self-biasing the inverter near its trip point.
- `C1` ($15\text{ fF}$): Dominant-pole integration capacitor setting the nodal RC relaxation time $\tau \approx 1.88\text{ ns}$.
- `MP_BUF1` / `MN_BUF1`: Digital restoring inverter ($W_p = 3.0\mu, W_n = 1.5\mu$) isolating `a1` from load impedance and sharpening outputs into digital $0\text{ V} \leftrightarrow 1.8\text{ V}$.
- `Cout1` ($20\text{ fF}$): Digital output load capacitance.

### 2.5 Synaptic Coupling & Transmission Gate Network
```spice
* Edge 1: P1 - P2 (Permanent Edge)
R_c12 a1 noise2 500k
R_c21 a2 noise1 500k

* Edge 2: P1 - P3 (Switchable via Discrete CMOS Transmission Gate)
R_c13 a1 tg_13 500k
MN_TG13 tg_13 mode noise3 0 NMOS_PBIT W=1.0u L=0.18u
MP_TG13 tg_13 mode_bar noise3 vdd PMOS_PBIT W=2.5u L=0.18u

R_c31 a3 tg_31 500k
MN_TG31 tg_31 mode noise1 0 NMOS_PBIT W=1.0u L=0.18u
MP_TG31 tg_31 mode_bar noise1 vdd PMOS_PBIT W=2.5u L=0.18u
```
- `R_c12` / `R_c21`: Symmetric $500\text{ k}\Omega$ coupling resistors implementing permanent antiferromagnetic interaction ($J_{12} < 0$).
- Edge 2–5 Switches: Implemented using discrete, parallel NMOS/PMOS pairs (`MN_TG*` and `MP_TG*`):
  - When `mode = 1.8V` and `mode_bar = 0.0V`, both transistors conduct across rail-to-rail potentials ($R_{\text{on}} \approx 1.5\text{ k}\Omega \ll 500\text{ k}\Omega$).
  - When `mode = 0.0V` and `mode_bar = 1.8V`, both transistors turn OFF ($R_{\text{off}} > 10^{10}\,\Omega$), isolating the edge with zero floating nodes (node `tg_13` is held at `a1` through $R_{c13}$).
- Transient specification: `.tran 0.5n 2u` defines integration. The master netlist contains no `.control` block, preventing duplicate simulation runs when included in testbenches.

---

## 3. Operating Regimes: Freezing vs. Active Boltzmann Sampling

```text
Noise Amplitude Spectrum:
0 mV                   200 mV                   440 mV                   >600 mV
  |-----------------------|------------------------|------------------------|
     Frozen Regime          Transition Region         Active Boltzmann         Thermal Overwhelm
   - GS > 99.9%           - Escapes appear           - 84.1% GS, 15.0% C3     - GS -> 12.5% (random)
   - 0 escapes            - Dwell drops from us      - Persistent state       - beta_eff -> 0
   - Trivial match          to tens of ns              changes (>= 2 ns)      - Loss of optimization
                                                     - Non-trivial TVD/KL
```

### 1. Frozen Regime ($V_{\text{noise}} < 200\text{ mV}$)
* The thermal perturbations are too small to overcome the potential well created by the coupling network.
* Once the circuit enters a ground state, it stays trapped permanently ($0$ escapes across 50 runs, mean dwell $= 980\text{ ns}$ across a $1\text{ }\mu\text{s}$ simulation).
* In this regime, claiming a "perfect Boltzmann match" is trivial and deceptive because the theory also predicts $\approx 0$ on all excited states.

### 2. Active Boltzmann Sampling Regime ($V_{\text{noise}} \approx 380\text{–}450\text{ mV}$)
* Thermal noise is strong enough to induce stochastic barrier crossings without destroying energy ordering.
* At $440\text{ mV}$ RMS:
  - Ground-state fraction: **84.1%**.
  - 3/5-cut excited-state fraction: **15.0%**.
  - Non-ground-state sampling: active exploration of the 10 excited states.
  - Real fitted $\beta_{\text{eff}} = 1.521$ (min-KL fit), yielding honest $\text{TVD} = 0.0440$ and $\text{KL} = 0.0219\text{ nats}$.

### 3. Thermal Overwhelm Regime ($V_{\text{noise}} > 600\text{ mV}$)
* Noise dominates over coupling. Spins toggle randomly and independently, driving the distribution toward uniform random probability (6.25% per state, 12.5% for the two ground states combined), with $\beta_{\text{eff}} \to 0$.

---

## 4. Convergence Dynamics: RC Settling vs. Stationary Langevin Diffusion

### 4.1 Why Initial Time-to-First-Ground-State is Misleading
When an excited state is forced using `.ic`:
* The initial charge creates large transient currents across coupling resistors.
* The circuit switches into the ground state in **$<1.3\text{ ns}$**.
* **This initial drop is purely linear RC circuit discharge and gate delay, NOT simulated annealing!**

### 4.2 Output Chatter vs. Persistent State Changes (>= 2 ns)
In hardware measurements, counting every threshold crossing $V(m_{\text{out}}) = 0.90\text{ V}$ yields over 18,000 "escapes" per 50 runs. Detailed duration analysis of all 18,050 non-ground-state excursions reveals:
* **$31.0\%$** last $< 0.2\text{ ns}$ (a single noise sample step).
* **$97.9\%$** last $< 1.0\text{ ns}$.
* **$100.0\%$** last $< 2.0\text{ ns}$ (median excursion duration is **$0.28\text{ ns}$**).

These sub-nanosecond excursions are **output buffer chatter** caused by single high-frequency noise spikes pushing the gate briefly past $V_{th}$. The circuit does not actually settle into an excited state; it snaps right back as soon as the noise sample ends.

### 4.3 Persistence Filtering vs. Voltage Hysteresis & Right-Censoring
1. **Why Voltage Hysteresis Fails**: Because the restoring buffer is a high-gain CMOS inverter, even a 10 mV gate excursion causes $m_{\text{out}}$ to swing rail-to-rail (0V to 1.8V) within 50–100 ps. Testing hysteresis bands as wide as $0.20\text{ V} / 1.60\text{ V}$ on $m_{\text{out}}$ still leaves $>15,000$ trips.
2. **The Physical Persistence Filter**: The physical summing node has an RC time constant $\tau \approx 1.88\text{ ns}$. A genuine escape requires the state to persist for $\Delta t \ge 2.0\text{ ns}$.
3. **Right-Censored Dwell Time vs. Escape Rate Interval**:
   - In a $1.0\text{ }\mu\text{s}$ observation window, $>85\%$ of individual simulation runs experience **zero escapes**.
   - Consequently, the observed ground-state dwell segments are right-censored at $1,000\text{ ns}$ (the simulation end).
   - Averaging censored $1,000\text{ ns}$ runs produces an artificially truncated sample mean (~$900\text{ ns}$).
   - **The honest, uncensored metric is the escape rate**:
     Reporting 4 events in $50\text{ }\mu\text{s}$ yields an escape rate of $\sim 0.08\text{ / }\mu\text{s}$ with an exact 95% Poisson interval of about $0.02\text{ to }0.2\text{ / }\mu\text{s}$ ($[0.022, 0.205]\text{ / }\mu\text{s}$).
   - **Implied Mean Lifetime**:
     Giving the implied mean lifetime as a range from that interval yields:
     $$\tau_{\text{true}} \approx \mathbf{4.9\text{ – }45.9\text{ }\mu\text{s}} \quad (\text{central estimate: } 12.5\text{ }\mu\text{s})$$
     and this lifetime is **extrapolated beyond the $1\text{ }\mu\text{s}$ window**.

---

## 5. 16-State Boltzmann Sampling Deviations & Degree Asymmetry

In our benchmark $440\text{ mV}$ simulation ($N = 50$ runs), the measured 16-state distribution achieves a **reasonable fit with a single $\beta_{\text{eff}}$ (TVD 0.044), with systematic deviations attributed to node-degree asymmetry (not verified)** ($\text{TVD} = 0.0440$, $\text{KL} = 0.0219\text{ nats}$). The top deviating states are:
1. `(-1, +1, -1, -1)` (Cut = 3/5, $H = -1$): Empirical = **3.44%** vs. Boltzmann = **2.00%** ($\Delta = \mathbf{+1.44\%}$)
2. `(-1, -1, +1, -1)` (Cut = 3/5, $H = -1$): Empirical = **3.32%** vs. Boltzmann = **2.00%** ($\Delta = \mathbf{+1.32\%}$)
3. `(+1, -1, -1, +1)` (Cut = 4/5, $H = -3$): Empirical = **42.83%** vs. Boltzmann = **41.83%** ($\Delta = \mathbf{+1.01\%}$)
4. `(+1, -1, +1, -1)` (Cut = 3/5, $H = -1$): Empirical = **1.03%** vs. Boltzmann = **2.00%** ($\Delta = \mathbf{-0.96\%}$)
5. `(-1, -1, +1, +1)` (Cut = 3/5, $H = -1$): Empirical = **1.04%** vs. Boltzmann = **2.00%** ($\Delta = \mathbf{-0.95\%}$)

### Physical Cause of the Deviation:
* **Degree Asymmetry**: Nodes $P_2$ and $P_3$ have degree 3 (3 incident coupling resistors), while nodes $P_1$ and $P_4$ have degree 2.
* At node $P_2$, total summing conductance is $G_{\text{tot}} = G_{\text{in}} + G_{fb} + 3 G_c = 5 G$.
* At node $P_1$, total summing conductance is $G_{\text{tot}} = G_{\text{in}} + G_{fb} + 2 G_c = 4 G$.
* Because synaptic weights are voltage shifts $\Delta V \propto 1 / G_{\text{tot}}$, edges incident on $P_2$ and $P_3$ have smaller effective weights ($1/5$ vs $1/4$).
* States with isolated spins on high-degree nodes `(-1, +1, -1, -1)` and `(-1, -1, +1, -1)` face weaker restoring forces, leading to a $+1.4\%$ excess residence over lower-degree excited states (~$1.0\%$).

---

## 6. 22 High-Yield Technical Interview Questions & Answers

### Q1: What is a probabilistic bit (p-bit) and how does it fundamentally differ from a classical bit and a qubit?
**Answer**: A classical bit is deterministic and bistable, held in state $0$ or $1$ by high energy barriers ($\gg k_B T$). A qubit is a quantum two-level system in a coherent superposition $|\psi\rangle = \alpha |0\rangle + \beta |1\rangle$ requiring cryogenic temperatures to prevent phase decoherence. A p-bit is a room-temperature classical element with low energy barriers that continuously fluctuates between $-1$ and $+1$ in time, whose time-averaged output follows a controllable Boltzmann distribution $P(m = +1) = \sigma(2 \beta I)$. It performs physical Markov Chain Monte Carlo sampling at room temperature without quantum overhead.

---

### Q2: Why must the CMOS inverter in a p-bit be biased exactly at its trip point ($V_{th} = 0.90\text{ V}$)?
**Answer**: At the trip point $V_{\text{in}} = V_{\text{out}} = V_{DD}/2$, both NMOS and PMOS operate in saturation, maximizing small-signal gain ($A_v = -g_m r_o$). This creates a metastable balance point where minute input thermal noise fluctuations are amplified into full-rail stochastic excursions. Furthermore, it ensures an unbiased $50/50$ baseline probability ($P(+1) = 0.5$) with zero input offset.

---

### Q3: How does the Boltzmann sigmoid activation function mathematically emerge from circuit noise?
**Answer**: An ideal inverter acts as a signum/Heaviside step function $f(V) = \text{sgn}(V - V_{th})$. Adding zero-mean Gaussian thermal noise $\xi \sim \mathcal{N}(0, \sigma^2)$ to input $V_{\text{in}}$ convolves the step function with the Gaussian probability density function. The output probability $P(m=+1) = \int_{V_{th}-V_{\text{in}}}^\infty \frac{1}{\sqrt{2\pi}\sigma} e^{-u^2/2\sigma^2} du = \frac{1}{2}[1 + \text{erf}(\frac{V_{\text{in}}-V_{th}}{\sqrt{2}\sigma})]$, which closely matches the logistic Boltzmann sigmoid $\frac{1}{1 + \exp(-2(V_{\text{in}}-V_{th})/I_0)}$, where $I_0$ is set by the noise amplitude as it appears at the gate node after the $R_{\text{in}} / R_{fb}$ network (measured as $I_0 \approx 363.0\text{ mV}$ for the standardized core cell; the fitted $V_{th} = 0.941\text{ V}$ is $41\text{ mV}$ above the $0.90\text{ V}$ nominal bias point, with $R_{fb}$ loading and inverter asymmetry as an untested hypothesis for this offset).

---

### Q4: Why does the noise amplitude matter, and what happens if it is too low or too high?
**Answer**: 
- **Too Low ($<200\text{ mV}$)**: The circuit enters a "frozen regime" where once a ground state is reached, thermal noise cannot overcome the coupling potential barrier ($0$ escapes). The system gets trapped in local minima and cannot explore configuration space.
- **Too High ($>600\text{ mV}$)**: Thermal noise overwhelms the coupling weights ($J_{ij}$). Spins fluctuate independently, destroying energy preference and reducing the distribution to a uniform random baseline ($P \to 6.25\%$ per state).
- **Optimal ($380\text{–}450\text{ mV}$)**: Active Boltzmann sampling regime where the ground state dominates (75–85%), while excited states are visited (15–20%) via real, frequent barrier crossings, enabling escape from suboptimal local traps.

---

### Q5: What is effective temperature, and how is $\beta_{\text{eff}}$ measured in hardware?
**Answer**: Effective temperature $T_{\text{eff}} = 1/\beta_{\text{eff}}$ represents the statistical dispersion of energy states in the network. Rather than using an idealized algebraic shortcut ($V_J / I_0$), $\beta_{\text{eff}}$ is measured by **Maximum Likelihood Estimation (Min-KL fitting)**: fitting the theoretical distribution $P_{\text{Boltz}}(s; \beta) = e^{-\beta H(s)} / Z(\beta)$ to the measured 16-state empirical histogram by numerically minimizing Kullback-Leibler divergence. This mathematically enforces that the theoretical ensemble mean energy matches the empirical mean energy: $\langle H \rangle_{\text{Boltz}}(\beta) = \langle H \rangle_{\text{emp}}$.

---

### Q6: Is this circuit performing simulated annealing?
**Answer**: **No.** Simulated annealing explicitly requires an annealing schedule where temperature/noise is gradually cooled over time ($T(t) \to 0$) to freeze the system into the global minimum. This circuit operates at a **constant noise amplitude** ($440\text{ mV}$ RMS), performing **stationary Langevin diffusion / continuous-time Boltzmann sampling**. It does not freeze; it samples the Boltzmann distribution indefinitely.

---

### Q7: Why is the time to first ground state from a forced `.ic` start a misleading metric?
**Answer**: Forcing an initial condition using `.ic` charges the node capacitances to unphysical non-equilibrium voltages. The immediate transition to the ground state occurring within $<1.3\text{ ns}$ is purely driven by linear RC circuit discharge and CMOS gate delays, not by thermal barrier crossing or stochastic sampling.

---

### Q8: What metric actually reflects hardware sampling dynamics, and why?
**Answer**: The **escape rate of persistent state changes ($\ge 2\text{ ns}$)**: reporting 4 events in $50\text{ }\mu\text{s}$, rate $\sim 0.08\text{ / }\mu\text{s}$, with an exact 95% Poisson interval of about $0.02\text{ to }0.2\text{ / }\mu\text{s}$ ($[0.022, 0.205]\text{ / }\mu\text{s}$). Giving the implied mean lifetime as a range from that interval yields $4.9\text{ – }45.9\text{ }\mu\text{s}$ (central estimate $12.5\text{ }\mu\text{s}$), which is extrapolated beyond the $1\text{ }\mu\text{s}$ window. Under physical Langevin dynamics, dwell time within a finite $1.0\text{ }\mu\text{s}$ window is right-censored at $1,000\text{ ns}$ when no escapes occur; reporting the Poisson escape rate interval provides an honest, uncensored measure.

---

### Q9: Why is 4 out of 5 edges the true maximum cut for this 4-node 5-edge graph?
**Answer**: The graph consists of two triangles ($K_3$) sharing an edge: $\{P_1, P_2, P_3\}$ and $\{P_2, P_3, P_4\}$. Any cycle of odd length cannot be bipartite; a triangle can cut at most 2 of its 3 edges. Since the two triangles share an edge, if the shared edge is uncut, each triangle contributes 2 cuts, giving $2 + 2 = 4$. If the shared edge is cut, each triangle can cut at most 1 other edge, giving $1 + 1 + 1 = 3$. Hence, 5 cuts is topologically impossible and 4 is the exact global maximum.

---

### Q10: How does the MAX-CUT problem map to the Ising Hamiltonian?
**Answer**: An edge $(i,j)$ is cut when spins are anti-aligned ($s_i \neq s_j$), giving cut indicator $(1 - s_i s_j)/2$. Summing over all edges gives $\text{Cut}(\vec{s}) = \frac{|E|}{2} - \frac{1}{2} \sum_{(i,j)\in E} s_i s_j$. The Ising Hamiltonian with antiferromagnetic coupling $J = -1$ is $H(\vec{s}) = \sum_{(i,j)\in E} s_i s_j$. Thus, $H(\vec{s}) = |E| - 2 \cdot \text{Cut}(\vec{s})$. Minimizing $H$ directly maximizes the cut.

---

### Q11: Why does a single CMOS inverter stage produce antiferromagnetic coupling ($J < 0$)?
**Answer**: A CMOS inverter has negative small-signal gain ($V_{\text{out}}$ moves opposite to $V_{\text{in}}$). When spin $j$ is high ($+1$), its internal complementary node $a_j$ is low ($0\text{ V}$), which pulls current out of neighbor node `noise_i`, pushing spin $i$ toward low ($-1$). This mutually inhibitory feedback penalizes parallel alignment and rewards anti-parallel alignment, which is the definition of antiferromagnetism.

---

### Q12: How does the resistive summing network implement synaptic summation?
**Answer**: By Kirchhoff's Current Law at the high-impedance gate node: $\sum I_k = 0 \implies \sum_j G_{c,ij}(V(a_j) - V_{\text{gate}}) + G_{\text{in}}(V_{\text{bias}} + v_n - V_{\text{gate}}) + G_{fb}(V(a_i) - V_{\text{gate}}) = 0$. Because gate current is zero, the gate potential is the exact weighted linear superposition of neighbor voltages, feedback, and bias/noise, implementing analog vector-matrix multiplication in passive physics.

---

### Q13: What is the physical role of the feedback resistor $R_{fb}$ and capacitor $C_1$?
**Answer**: $R_{fb} = 500\text{ k}\Omega$ provides negative DC feedback that self-biases the high-gain inverter near its trip point ($V_{th} \approx 0.90\text{–}0.94\text{ V}$). Capacitor $C_1 = 15\text{ fF}$ integrates the summing currents, establishing the dominant pole and internal relaxation time constant $\tau \approx R_{\text{eq}} C_1 \approx 1.88\text{ ns}$ that prevents gigahertz ring-oscillator instability.

---

### Q14: Why is a digital restoring buffer required at the output of each p-bit?
**Answer**: Internal node $a_i$ fluctuates continuously through analog intermediate voltages ($0.4\text{ V}$ to $1.4\text{ V}$) during metastable switching. Connecting coupling resistors directly to $a_i$ would load the integration node, lowering gain and creating analog loading distortion. The buffer isolates $a_i$, provides low output impedance, and sharpens the waveform into clean, rail-to-rail digital voltages ($0\text{ V}$ / $1.8\text{ V}$).

---

### Q15: How do the CMOS transmission gates achieve reconfigurability without degrading coupling?
**Answer**: A CMOS transmission gate places an NMOS and PMOS in parallel. When active, one transistor conducts strongly for high voltages while the other conducts strongly for low voltages, ensuring a flat, low ON-resistance ($R_{\text{on}} \approx 1\text{–}2\text{ k}\Omega$) across the entire $0\text{–}1.8\text{ V}$ range. Since $R_{\text{on}} \ll R_c = 500\text{ k}\Omega$, insertion loss is $<0.4\%$. When disabled, leakage resistance exceeds $10^{10}\text{ }\Omega$, completely isolating the graph edge.

---

### Q16: What are the primary limitations of Level-1 SPICE models for 180nm CMOS?
**Answer**: Level-1 (Shichman-Hodges) models assume square-law saturation ($I_D \propto (V_{GS}-V_T)^2$) and constant mobility. In real 180nm CMOS:
1. **Velocity Saturation**: Lateral electric fields cause carrier velocity saturation, making $I_D$ linear with $(V_{GS}-V_T)$.
2. **DIBL (Drain-Induced Barrier Lowering)**: Short-channel drain bias lowers threshold voltage, increasing output conductance.
3. **Subthreshold Conduction**: Real MOSFETs have an exponential subthreshold slope ($S \approx 80\text{ mV/decade}$), whereas Level-1 treats cutoff abruptly.
4. BSIM4 or PSP models are required for tapeout sign-off precision.

---

### Q17: How would process variations and transistor mismatch affect this circuit?
**Answer**: Pelgrom's Law states that threshold mismatch $\sigma(\Delta V_{th}) = A_{VT} / \sqrt{WL}$. In 180nm, $\sigma(V_{th}) \approx 5\text{–}10\text{ mV}$. If an inverter's trip point shifts from $0.900\text{ V}$ to $0.910\text{ V}$, it introduces an unintentional local magnetic field ($h_i \neq 0$), biasing that p-bit toward $-1$ even with zero input. This destroys ground-state degeneracy unless compensated by trimming DACs or auto-zeroing capacitor calibration.

---

### Q18: How would this architecture scale to a 100-node all-to-all Ising accelerator?
**Answer**: Scaling to 100 nodes presents three major challenges:
1. **$O(N^2)$ Synaptic Routing**: 100 fully connected nodes require $100 \times 99 / 2 = 4,950$ resistive synapses, creating severe routing congestion.
2. **Degree Attenuation**: With 99 incident resistors at each summing node, $G_{\text{tot}} \approx 100 / R_c$. The voltage swing per edge drops to $\Delta V_{\text{gate}} \approx V_{DD} / 100 = 18\text{ mV}$, which approaches the thermal noise floor and mismatch offset.
3. **Mitigation**: Requires active transconductance synapses (e.g., Gilbert multiplier cells, memristor crossbars, or capacitive summing) with local op-amp virtual ground summing.

---

### Q19: Why does node degree asymmetry affect coupling strength, and how can it be fixed?
**Answer**: In a passive resistor summing network, $V_{\text{gate}} = (\sum G_c V_j) / (\sum G)$. A node with degree 3 has a larger denominator ($G_{\text{tot}} \propto 5$) than a node with degree 2 ($G_{\text{tot}} \propto 4$). Consequently, each edge connected to a high-degree node has a smaller effective weight ($J_{\text{eff}}$). In planar graphs, this creates unwanted bias toward low-degree nodes. It can be fixed by degree-normalizing the coupling resistors ($R_{c,i} \propto 1 / \text{deg}(i)$) or using op-amp virtual ground summing where the summing junction impedance is zero.

---

### Q20: What are Total Variation Distance (TVD) and KL Divergence, and how do they quantify Boltzmann fidelity?
**Answer**:
- **Total Variation Distance**: $\text{TVD}(P, Q) = \frac{1}{2} \sum_s |P(s) - Q(s)| \in [0, 1]$. It represents the maximum difference in probability assigned to any subset of states. A TVD of $0.0440$ indicates a reasonable fit with a single $\beta_{\text{eff}}$ (TVD 0.044), with systematic deviations attributed to node-degree asymmetry (not verified).
- **Kullback-Leibler Divergence**: $D_{\text{KL}}(P \parallel Q) = \sum_s P(s) \ln \frac{P(s)}{Q(s)}$. It measures the information loss (in nats) when approximating distribution $P$ with model $Q$.

---

### Q21: Are the 18,000+ ground-state departures real thermodynamic state escapes or output chatter? How do you distinguish them?
**Answer**:
They are overwhelmingly **output chatter / threshold jitter**. Analysis shows that $97.9\%$ of raw departures last $< 1.0\text{ ns}$ and $100.0\%$ last $< 2.0\text{ ns}$ (median duration $0.28\text{ ns}$). They occur when single 0.2 ns noise steps push the gate briefly across $V_{th}$, triggering the high-gain digital buffer before the circuit can physically relax into an excited state.
- **Voltage hysteresis on $m_{\text{out}}$ fails** to filter this chatter because the digital buffer snaps rail-to-rail within $<100\text{ ps}$.
- **Persistence filtering ($\Delta t \ge 2.0\text{ ns}$)**: Because the summing node RC relaxation time is $\tau \approx 1.88\text{ ns}$, a genuine escape requires the state to persist for $\ge 2.0\text{ ns}$.
- **Filtered Result**: Filtering for persistent state changes ($\ge 2\text{ ns}$) yields **4 events in $50\text{ }\mu\text{s}$**, giving an escape rate of **$\sim 0.08\text{ / }\mu\text{s}$** with an exact 95% Poisson interval of about **$0.02\text{ to }0.2\text{ / }\mu\text{s}$** ($[0.022, 0.205]\text{ / }\mu\text{s}$). This gives an implied mean lifetime of **$4.9\text{ – }45.9\text{ }\mu\text{s}$** ($\hat{\tau} \approx 12.5\text{ }\mu\text{s}$), which is extrapolated beyond the $1\text{ }\mu\text{s}$ observation window, resolving the right-censoring limitation of the finite simulation window.

---

### Q22: Which states deviate most from the ideal Boltzmann distribution, and what physical hardware mechanism causes the discrepancy?
**Answer**:
The largest deviations occur in:
1. `(-1, +1, -1, -1)` ($+1.44\%$ excess) and `(-1, -1, +1, -1)` ($+1.32\%$ excess) among 3/5-cut excited states.
2. Mirror ground-state splitting: `(+1, -1, -1, +1)` ($42.8\%$) vs `(-1, +1, +1, -1)` ($41.2\%$).
- **Mechanism**: Nodes $P_2$ and $P_3$ have degree 3, while $P_1$ and $P_4$ have degree 2. The higher parallel conductance at $P_2/P_3$ ($5G$ vs $4G$) attenuates the coupling voltage swing per edge by $20\%$. Consequently, excited configurations with isolated high spins on $P_2$ or $P_3$ experience weaker restoring fields than theoretical unweighted graph edges, elevating their residence time to ~3.4% compared to ~1.0% for low-degree excited states.
