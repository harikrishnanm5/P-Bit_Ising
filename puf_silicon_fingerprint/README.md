# PROJECT 1: CMOS SRAM / Cross-Coupled Silicon PUF (Physically Unclonable Function)

## 1. Executive Summary & Pitch to Judges
Modern IoT edge devices are vulnerable to firmware cloning, reverse engineering, and physical theft. Traditional cryptographic key storage (eFLASH / EEPROM) consumes static power and is vulnerable to side-channel probing.

This project presents a **Zero-Power CMOS Silicon PUF (Physically Unclonable Function)** designed using standard 180nm / SkyWater 130nm CMOS technology. By exploiting unavoidable microscopic manufacturing process variations (threshold voltage $\Delta V_{th}$ and channel length $\Delta L$), each fabricated silicon die settles into a unique, unclonable, and reproducible digital fingerprint upon power-up.

---

## 2. Circuit Architecture
The core PUF cell consists of a cross-coupled inverter latch with an active power-gating PMOS switch (`M_ENABLE`):

* **Inverter 1**: `M_P1` (PMOS) and `M_N1` (NMOS) driving node `Q`
* **Inverter 2**: `M_P2` (PMOS) and `M_N2` (NMOS) driving node `Q_BAR`
* **Power-Gating Header**: `M_EN` (PMOS) controlled by `V_EN` (Active LOW power-up pulse)

### Working Mechanism:
1. **Reset Phase ($V_{EN} = V_{DD}$)**: The cell is unpowered; nodes `Q` and `Q_BAR` are pulled low.
2. **Evaluation Phase ($V_{EN} = 0\text{ V}$)**: Supply voltage ramps up. Both inverters enter the **metastable operating point** ($V_Q \approx V_{Q\_BAR} \approx V_{DD}/2$).
3. **Mismatch Amplification**: The microscopic asymmetry ($\Delta V_{th} = V_{th1} - V_{th2}$) creates an infinitesimal current difference $\Delta I$. The positive feedback loop rapidly amplifies this difference, pulling one node to $V_{DD}$ (Logic 1) and the other to $GND$ (Logic 0).
4. **Result**: The output bit is 100% deterministic for a single chip, but 100% random across different chips (Inter-die Hamming Distance $\approx 50\%$).

---

## 3. SPICE Netlist (`puf_cell.cir`)
Below is the verified Ngspice netlist including a Monte Carlo process-variation testbench.
