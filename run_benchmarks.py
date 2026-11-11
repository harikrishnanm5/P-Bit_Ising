"""Automated Ngspice benchmark & waveform plotting suite for the 4-p-bit CMOS Ising core.

Generates publication-grade plots in results/:
  1. ising_16state_histogram.png  - 100-run 16-state spin probability distribution (Main Result)
  2. pbit_sigmoid_curve.png       - Single p-bit Boltzmann sigmoid & stochastic waveform
  3. coupled_pair_transient.png   - Coupled 2-p-bit antiferromagnetic dimer transient
  4. ising_4pbit_transient.png    - 4-p-bit Ising core transient waveforms & MAX-CUT trajectory
  5. convergence_analysis.png     - Excited-start ensemble relaxation & <500 us convergence CDF
"""

from collections import Counter
import os
import re
import subprocess
import tempfile
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

GRAPH_EDGES = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
VDD = 1.8
VTH = 0.9


def configure_plot_style():
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "#fafafa",
            "axes.edgecolor": "#cbd5e1",
            "axes.grid": True,
            "grid.color": "#e2e8f0",
            "grid.linestyle": "--",
            "grid.linewidth": 0.7,
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
        }
    )


def run_ngspice_wrdata(netlist_text, signals, workdir, tag):
    """Inject wrdata + quit into a SPICE netlist, execute Ngspice, and return loaded array."""
    out_txt = os.path.join(workdir, f"{tag}.txt").replace(os.sep, "/")
    cir_path = os.path.join(workdir, f"{tag}.cir")

    idx = netlist_text.rfind(".control")
    if idx == -1:
        raise ValueError(f"No .control block found in {tag}")

    ctrl = (
        ".control\n"
        "  run\n"
        f"  wrdata {out_txt} {' '.join(signals)}\n"
        "  quit\n"
        ".endc\n"
        ".end\n"
    )
    patched = netlist_text[:idx] + ctrl
    with open(cir_path, "w") as f:
        f.write(patched)

    subprocess.run(["ngspice", "-b", cir_path], capture_output=True, check=False)
    if not os.path.exists(out_txt):
        raise RuntimeError(f"Ngspice failed to generate {out_txt}")
    return np.loadtxt(out_txt)


def compute_cuts(spins):
    """Compute vector of satisfied MAX-CUT edges (0..5) for shape (N, 4) spin array."""
    return np.sum(
        [(spins[:, i] != spins[:, j]).astype(int) for i, j in GRAPH_EDGES], axis=0
    )


def state_to_cut(state_tuple):
    return sum(1 for i, j in GRAPH_EDGES if state_tuple[i] != state_tuple[j])


def make_core_variant(base_netlist, phases, init_spins=None, tstep="500n", tstop="1m", uic=False):
    """Return a copy of ising_4pbit_core.cir with custom dither phases and optional .ic."""
    text = base_netlist
    freqs = ["5.0k", "7.1k", "11.3k", "13.7k"]
    for i, (freq, ph) in enumerate(zip(freqs, phases), start=1):
        pattern = rf"(Vnoise{i}\s+dither{i}\s+in{i}\s+SIN\([^\)]+)\)"
        text = re.sub(pattern, rf"\g<1> {ph:.2f})", text)

    tran_line = f".tran {tstep} {tstop}" + (" uic" if uic else "")
    if init_spins is not None:
        ic_vals = []
        for i, s in enumerate(init_spins, start=1):
            va = 0.0 if s == 1 else VDD
            vm = VDD if s == 1 else 0.0
            ic_vals.append(f"v(a{i})={va:.1f} v(m_out{i})={vm:.1f}")
        tran_line = f".ic {' '.join(ic_vals)}\n{tran_line}"

    text = re.sub(r"\.tran\s+\S+\s+\S+", tran_line, text)
    return text


def generate_sigmoid_plot(spice_dir, results_dir, workdir):
    print("[1/5] Simulating p-bit sigmoid response & single-cell fluctuation...")
    with open(os.path.join(spice_dir, "p_bit_sigmoid_validation.cir")) as f:
        sig_base = f.read()

    vins = np.round(np.arange(0.84, 0.965, 0.01), 3)
    probs = []
    for v in vins:
        mod = sig_base.replace("Vin in_node 0 DC 0.9V", f"Vin in_node 0 DC {v:.3f}")
        mod = re.sub(r"\.control[\s\S]*", ".tran 0.1n 400n\n.control\n.endc\n.end\n", mod)
        data = run_ngspice_wrdata(mod, ["v(m_out)"], workdir, f"sig_{v:.3f}")
        mask = data[:, 0] >= 10e-9
        probs.append(float(np.mean(data[mask, 1] > VTH)))

    probs = np.array(probs)
    # Fit Boltzmann sigmoid P(m=+1) = 1 / (1 + exp(-2*(Vin - V0)/I0)) via logit regression on interior points
    interior = (probs > 0.005) & (probs < 0.995)
    x_fit = vins[interior]
    y_logit = np.log(probs[interior] / (1.0 - probs[interior]))
    slope, intercept = np.polyfit(x_fit, y_logit, 1)
    i0_fit = 2.0 / slope
    v0_fit = -intercept / slope
    v_dense = np.linspace(0.84, 0.96, 300)
    p_fit = 1.0 / (1.0 + np.exp(-2.0 * (v_dense - v0_fit) / i0_fit))
    ss_res = np.sum((probs - 1.0 / (1.0 + np.exp(-2.0 * (vins - v0_fit) / i0_fit))) ** 2)
    ss_tot = np.sum((probs - np.mean(probs)) ** 2)
    r2 = 1.0 - ss_res / ss_tot

    with open(os.path.join(spice_dir, "p_bit.cir")) as f:
        pbit_net = f.read()
    cell_data = run_ngspice_wrdata(
        pbit_net, ["v(analog_out)", "v(m_out)"], workdir, "pbit_cell"
    )
    t_ns = cell_data[:, 0] * 1e9
    win = t_ns <= 40.0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6), dpi=300)

    ax1.plot(v_dense, p_fit, color="#0284c7", lw=2.2, label=f"Boltzmann Fit ($V_{{th}}={v0_fit:.3f}$ V, $I_0={i0_fit*1e3:.1f}$ mV, $R^2={r2:.4f}$)")
    ax1.scatter(vins, probs, color="#0f172a", s=45, zorder=5, label="Ngspice Empirical $P(m=+1)$")
    ax1.axhline(0.5, color="#94a3b8", ls=":", lw=1)
    ax1.axvline(0.90, color="#94a3b8", ls=":", lw=1)
    ax1.set_title("Single CMOS p-Bit Boltzmann Sigmoid Activation")
    ax1.set_xlabel("Synaptic Input Bias $V_{in}$ (V)")
    ax1.set_ylabel("Output Probability $P(m = +1)$")
    ax1.set_xlim(0.835, 0.965)
    ax1.set_ylim(-0.04, 1.06)
    ax1.legend(loc="upper left", fontsize=8.5, frameon=True)

    ax2.plot(t_ns[win], cell_data[win, 1], color="#f59e0b", lw=1.1, alpha=0.85, label="$V(\\mathrm{analog\\_out})$")
    ax2.plot(t_ns[win], cell_data[win, 3], color="#0f766e", lw=1.5, label="$V(m_{\\mathrm{out}})$ Digital Spin")
    ax2.axhline(VTH, color="#64748b", ls="--", lw=0.9, label="$V_{th} = 0.90$ V")
    ax2.set_title("Metastable Stochastic Switching at $V_{in} = 0.90$ V")
    ax2.set_xlabel("Time (ns)")
    ax2.set_ylabel("Voltage (V)")
    ax2.set_xlim(0, 40)
    ax2.set_ylim(-0.1, 1.95)
    ax2.legend(loc="upper right", fontsize=8.5, frameon=True)

    plt.tight_layout()
    out_path = os.path.join(results_dir, "pbit_sigmoid_curve.png")
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"      Saved {out_path} (Vth={v0_fit:.3f} V, I0={i0_fit*1e3:.1f} mV, R^2={r2:.4f})")
    return v0_fit, i0_fit, r2


def generate_coupled_pair_plot(spice_dir, results_dir, workdir):
    print("[2/5] Simulating coupled 2-p-bit antiferromagnetic pair...")
    with open(os.path.join(spice_dir, "step6_coupled_pbit_pair.cir")) as f:
        net = f.read()

    data = run_ngspice_wrdata(net, ["v(m_out1)", "v(m_out2)"], workdir, "pair")
    t_ms = data[:, 0] * 1e3
    v1, v2 = data[:, 1], data[:, 3]
    s1 = np.where(v1 > VTH, 1, -1)
    s2 = np.where(v2 > VTH, 1, -1)
    anti_pct = 100.0 * float(np.mean(s1 != s2))
    p10 = 100.0 * float(np.mean((s1 == 1) & (s2 == -1)))
    p01 = 100.0 * float(np.mean((s1 == -1) & (s2 == 1)))

    fig, axes = plt.subplots(3, 1, figsize=(11, 5.8), sharex=True, dpi=300, gridspec_kw={"height_ratios": [1, 1, 0.55]})
    axes[0].plot(t_ms, v1, color="#0284c7", lw=1.2)
    axes[0].set_ylabel("$V(m_{\\mathrm{out1}})$ (V)")
    axes[0].set_ylim(-0.15, 1.98)
    axes[0].set_title(
        f"Coupled 2-p-Bit Antiferromagnetic Pair ($R_c = 500\\text{{ k}}\\Omega$): {anti_pct:.1f}% Anti-Correlated"
    )

    axes[1].plot(t_ms, v2, color="#d97706", lw=1.2)
    axes[1].set_ylabel("$V(m_{\\mathrm{out2}})$ (V)")
    axes[1].set_ylim(-0.15, 1.98)

    anti = (s1 != s2).astype(int)
    axes[2].fill_between(t_ms, 0, anti, step="post", color="#0f766e", alpha=0.75, label=f"Anti-correlated $s_1 \\neq s_2$ ({anti_pct:.1f}%: (+1,-1) {p10:.1f}%, (-1,+1) {p01:.1f}%)")
    axes[2].set_ylabel("Cut")
    axes[2].set_yticks([0, 1])
    axes[2].set_yticklabels(["0/1", "1/1"])
    axes[2].set_ylim(-0.1, 1.25)
    axes[2].set_xlabel("Time (ms)")
    axes[2].set_xlim(0, 1.0)
    axes[2].legend(loc="lower right", fontsize=8.5, frameon=True)

    plt.tight_layout()
    out_path = os.path.join(results_dir, "coupled_pair_transient.png")
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"      Saved {out_path} (anti-correlated = {anti_pct:.1f}%)")
    return anti_pct


def generate_4pbit_transient_plot(spice_dir, results_dir, workdir):
    print("[3/5] Simulating nominal 4-p-bit Ising core transient...")
    with open(os.path.join(spice_dir, "ising_4pbit_core.cir")) as f:
        net = f.read()

    data = run_ngspice_wrdata(
        net, ["v(m_out1)", "v(m_out2)", "v(m_out3)", "v(m_out4)"], workdir, "core_nom"
    )
    t_ms = data[:, 0] * 1e3
    outs = data[:, [1, 3, 5, 7]]
    spins = np.where(outs > VTH, 1, -1)
    cuts = compute_cuts(spins)
    gs_pct = 100.0 * float(np.mean(cuts == 4))

    colors = ["#0284c7", "#0f766e", "#d97706", "#7c3aed"]
    labels = ["$P_1$ ($m_{\\mathrm{out1}}$)", "$P_2$ ($m_{\\mathrm{out2}}$)", "$P_3$ ($m_{\\mathrm{out3}}$)", "$P_4$ ($m_{\\mathrm{out4}}$)"]

    fig, axes = plt.subplots(5, 1, figsize=(11, 7.6), sharex=True, dpi=300, gridspec_kw={"height_ratios": [1, 1, 1, 1, 0.85]})
    for idx in range(4):
        axes[idx].plot(t_ms, outs[:, idx], color=colors[idx], lw=1.15)
        axes[idx].set_ylabel(f"{labels[idx]}\n(V)")
        axes[idx].set_ylim(-0.15, 2.0)

    axes[0].set_title(
        f"4-p-Bit CMOS Ising Core Transient ($V_{{\\mathrm{{MODE}}}} = 1.8\\text{{ V}}$, Ground-State Residence = {gs_pct:.1f}%)"
    )

    axes[4].step(t_ms, cuts, where="post", color="#0f766e", lw=1.2)
    axes[4].fill_between(t_ms, 0, cuts, step="post", where=(cuts == 4), color="#0f766e", alpha=0.22, label="MAX-CUT Ground State (4/5 edges)")
    axes[4].axhline(4, color="#0f766e", ls="--", lw=0.9)
    axes[4].set_ylabel("Cut\nEdges")
    axes[4].set_yticks([0, 2, 3, 4])
    axes[4].set_ylim(0, 4.6)
    axes[4].set_xlabel("Time (ms)")
    axes[4].set_xlim(0, 1.0)
    axes[4].legend(loc="lower right", fontsize=8.5, frameon=True)

    plt.tight_layout()
    out_path = os.path.join(results_dir, "ising_4pbit_transient.png")
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"      Saved {out_path} (ground-state residence = {gs_pct:.1f}%)")
    return gs_pct


def generate_16state_histogram(spice_dir, results_dir, workdir, n_runs=100):
    print(f"[4/5] Running {n_runs}-trial Monte Carlo ensemble for 16-state histogram...")
    with open(os.path.join(spice_dir, "ising_4pbit_core.cir")) as f:
        base_net = f.read()

    rng = np.random.default_rng(2026)
    all_16_states = [
        tuple(1 if (i >> b) & 1 else -1 for b in [3, 2, 1, 0]) for i in range(16)
    ]

    time_counts = Counter()
    final_counts = Counter()
    t_grid = np.linspace(0.0, 1000.0, 1001)  # uniform 1 us grid across 1 ms

    for r in range(n_runs):
        phases = rng.uniform(0.0, 360.0, size=4)
        init_state = all_16_states[r % 16]
        net = make_core_variant(
            base_net, phases=phases, init_spins=init_state, tstep="500n", tstop="1m", uic=True
        )
        data = run_ngspice_wrdata(
            net, ["v(m_out1)", "v(m_out2)", "v(m_out3)", "v(m_out4)"], workdir, f"mc_{r}"
        )
        t_us = data[:, 0] * 1e6
        spins = np.where(data[:, [1, 3, 5, 7]] > VTH, 1, -1)

        idx = np.clip(np.searchsorted(t_us, t_grid, side="right") - 1, 0, len(t_us) - 1)
        for row in spins[idx]:
            time_counts[tuple(row)] += 1
        final_counts[tuple(spins[-1])] += 1

    total_samples = sum(time_counts.values())
    probs = [100.0 * time_counts[st] / total_samples for st in all_16_states]
    cuts = [state_to_cut(st) for st in all_16_states]

    gs_total = sum(p for p, c in zip(probs, cuts) if c == 4)
    c3_total = sum(p for p, c in zip(probs, cuts) if c == 3)
    c_low_total = sum(p for p, c in zip(probs, cuts) if c < 3)
    final_gs = sum(v for st, v in final_counts.items() if state_to_cut(st) == 4)

    bar_colors = []
    for c in cuts:
        if c == 4:
            bar_colors.append("#0f766e")
        elif c == 3:
            bar_colors.append("#d97706")
        else:
            bar_colors.append("#94a3b8")

    state_labels = [
        "(" + ",".join("+" if s == 1 else "-" for s in st) + ")"
        for st in all_16_states
    ]

    fig, ax = plt.subplots(figsize=(13, 5.6), dpi=300)
    x = np.arange(16)
    bars = ax.bar(x, probs, color=bar_colors, edgecolor="#1e293b", linewidth=0.6, width=0.72)

    for bar, p, c in zip(bars, probs, cuts):
        if c == 4 or p >= 1.5:
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                bar.get_height() + 0.8,
                f"{p:.1f}%",
                ha="center",
                va="bottom",
                fontsize=8.5,
                fontweight="bold" if c == 4 else "normal",
                color="#0f766e" if c == 4 else "#334155",
            )

    ax.axhline(100.0 / 16.0, color="#64748b", ls="--", lw=1.0, label="Uniform Random Baseline (6.25%)")
    ax.set_xticks(x)
    ax.set_xticklabels(state_labels, rotation=40, ha="right", fontsize=9)
    ax.set_ylabel("State Occupation Probability (%)")
    ax.set_xlabel("4-p-Bit Spin Configuration $(s_1, s_2, s_3, s_4)$")
    ax.set_ylim(0, max(probs) + 8.0)
    ax.set_title(
        f"16-State Boltzmann Sampling Distribution Across {n_runs} Ngspice Runs (1 ms Each)"
    )

    from matplotlib.patches import Patch

    legend_handles = [
        Patch(facecolor="#0f766e", edgecolor="#1e293b", label=f"4/5 Cut Ground States: {gs_total:.1f}%"),
        Patch(facecolor="#d97706", edgecolor="#1e293b", label=f"3/5 Cut Excited States: {c3_total:.1f}%"),
        Patch(facecolor="#94a3b8", edgecolor="#1e293b", label=f"0/5 & 2/5 Cut High-Energy States: {c_low_total:.2f}%"),
        plt.Line2D([0], [0], color="#64748b", ls="--", lw=1.0, label="Uniform Baseline (6.25%)"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=9, frameon=True)

    summary_txt = (
        f"Ensemble Summary (N = {n_runs} runs)\n"
        f"• (-1, +1, +1, -1): {probs[6]:.1f}%\n"
        f"• (+1, -1, -1, +1): {probs[9]:.1f}%\n"
        f"• Combined 4/5 Ground States: {gs_total:.1f}%\n"
        f"• Final-State (t=1 ms) Ground States: {final_gs}/{n_runs} ({100*final_gs/n_runs:.0f}%)"
    )
    ax.text(
        0.985,
        0.96,
        summary_txt,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.8,
        bbox={"boxstyle": "round,pad=0.5", "facecolor": "white", "edgecolor": "#cbd5e1"},
    )

    plt.tight_layout()
    out_path = os.path.join(results_dir, "ising_16state_histogram.png")
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(
        f"      Saved {out_path} | GS=({probs[6]:.2f}% + {probs[9]:.2f}% = {gs_total:.2f}%), final={final_gs}/{n_runs}"
    )
    return probs[6], probs[9], gs_total, final_gs


def generate_convergence_plot(spice_dir, results_dir, workdir, n_runs=100):
    print(f"[5/5] Running {n_runs}-trial excited-start convergence benchmark...")
    with open(os.path.join(spice_dir, "ising_4pbit_core.cir")) as f:
        base_net = f.read()

    rng = np.random.default_rng(777)
    non_gs_states = [
        st
        for i in range(16)
        for st in [tuple(1 if (i >> b) & 1 else -1 for b in [3, 2, 1, 0])]
        if state_to_cut(st) < 4
    ]

    t_grid = np.linspace(0.0, 500.0, 1001)  # 0 to 500 us
    cuts_matrix = []
    first_hit_us = []

    for r in range(n_runs):
        phases = rng.uniform(0.0, 360.0, size=4)
        init_state = non_gs_states[r % len(non_gs_states)]
        net = make_core_variant(
            base_net, phases=phases, init_spins=init_state, tstep="100n", tstop="500u", uic=True
        )
        data = run_ngspice_wrdata(
            net, ["v(m_out1)", "v(m_out2)", "v(m_out3)", "v(m_out4)"], workdir, f"conv_{r}"
        )
        t_us = data[:, 0] * 1e6
        spins = np.where(data[:, [1, 3, 5, 7]] > VTH, 1, -1)
        cuts = compute_cuts(spins)

        hit = np.where(cuts == 4)[0]
        if len(hit) > 0:
            first_hit_us.append(float(t_us[hit[0]]))
        else:
            first_hit_us.append(np.nan)

        idx = np.clip(np.searchsorted(t_us, t_grid, side="right") - 1, 0, len(t_us) - 1)
        cuts_matrix.append(cuts[idx])

    cuts_matrix = np.array(cuts_matrix)
    first_hit_us = np.array(first_hit_us)
    valid_hits = first_hit_us[~np.isnan(first_hit_us)]
    hit_count = len(valid_hits)
    med_us = float(np.median(valid_hits))
    p90_us = float(np.percentile(valid_hits, 90))
    max_us = float(np.max(valid_hits))
    mean_cut = np.mean(cuts_matrix, axis=0)
    gs_frac = 100.0 * np.mean(cuts_matrix == 4, axis=0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.7), dpi=300)

    ax1.plot(t_grid, mean_cut, color="#0f766e", lw=2.0, label="Ensemble Mean Cut Edges (out of 5)")
    ax1.axhline(4.0, color="#0f766e", ls="--", lw=1.0, label="Theoretical Max Cut (4/5)")
    ax1_r = ax1.twinx()
    ax1_r.plot(t_grid, gs_frac, color="#0284c7", lw=1.6, alpha=0.85, label="Instantaneous Ground-State Fraction (%)")
    ax1_r.set_ylabel("Ground-State Fraction (%)", color="#0284c7")
    ax1_r.set_ylim(0, 105)
    ax1.set_title(f"Ensemble Relaxation from Excited Initial States (N = {n_runs})")
    ax1.set_xlabel("Time ($\\mu$s)")
    ax1.set_ylabel("Mean Satisfied Cut Edges", color="#0f766e")
    ax1.set_xlim(0, 500)
    ax1.set_ylim(2.5, 4.25)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1_r.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="lower right", fontsize=8.5, frameon=True)

    sorted_t = np.sort(valid_hits)
    cdf = 100.0 * np.arange(1, len(sorted_t) + 1) / n_runs
    ax2.step(sorted_t, cdf, where="post", color="#0f766e", lw=2.2, label="Empirical Convergence CDF")
    ax2.axvline(500.0, color="#dc2626", ls="--", lw=1.5, label="500 $\\mu$s Target Spec")
    ax2.axvline(max_us, color="#d97706", ls=":", lw=1.4, label=f"Observed Worst Case ({max_us:.2f} $\\mu$s)")
    ax2.set_xscale("log")
    ax2.set_xlim(5e-4, 1000.0)
    ax2.set_ylim(0, 105)
    ax2.set_title("Time to First MAX-CUT Ground State (Log Time Scale)")
    ax2.set_xlabel("Convergence Time $t_{\\mathrm{GS}}$ ($\\mu$s)")
    ax2.set_ylabel("Cumulative Runs Reaching 4/5 Cut (%)")
    ax2.legend(loc="lower right", fontsize=8.5, frameon=True)

    stats_txt = (
        f"Converged < 500 $\\mu$s: {hit_count}/{n_runs} ({100*hit_count/n_runs:.0f}%)\n"
        f"Median $t_{{\\mathrm{{GS}}}}$: {med_us:.3f} $\\mu$s\n"
        f"90th Pct $t_{{\\mathrm{{GS}}}}$: {p90_us:.3f} $\\mu$s\n"
        f"Max $t_{{\\mathrm{{GS}}}}$: {max_us:.2f} $\\mu$s"
    )
    ax2.text(
        0.50,
        0.52,
        stats_txt,
        transform=ax2.transAxes,
        ha="center",
        va="center",
        fontsize=8.8,
        bbox={"boxstyle": "round,pad=0.5", "facecolor": "white", "edgecolor": "#cbd5e1"},
    )

    plt.tight_layout()
    out_path = os.path.join(results_dir, "convergence_analysis.png")
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(
        f"      Saved {out_path} | {hit_count}/{n_runs} < 500us | med={med_us:.3f}us, p90={p90_us:.2f}us, max={max_us:.2f}us"
    )
    return hit_count, med_us, p90_us, max_us


def main():
    configure_plot_style()
    repo_root = os.path.dirname(os.path.abspath(__file__))
    spice_dir = os.path.join(repo_root, "spice")
    results_dir = os.path.join(repo_root, "results")
    os.makedirs(results_dir, exist_ok=True)

    with tempfile.TemporaryDirectory() as workdir:
        generate_sigmoid_plot(spice_dir, results_dir, workdir)
        generate_coupled_pair_plot(spice_dir, results_dir, workdir)
        generate_4pbit_transient_plot(spice_dir, results_dir, workdir)
        generate_16state_histogram(spice_dir, results_dir, workdir, n_runs=100)
        generate_convergence_plot(spice_dir, results_dir, workdir, n_runs=100)

    print("All benchmark plots generated in results/.")


if __name__ == "__main__":
    main()
