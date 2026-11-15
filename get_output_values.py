import os
import subprocess
import numpy as np

GRAPH_EDGES = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]


def compute_cut(spins):
    satisfied = [(i + 1, j + 1) for i, j in GRAPH_EDGES if spins[i] != spins[j]]
    return len(satisfied), satisfied


def get_4pbit_values():
    print("=" * 77)
    print("               4-p-BIT ISING CORE TRANSIENT OUTPUT VALUES")
    print("           (Active Boltzmann Regime: 440mV RMS trnoise Noise)")
    print("=" * 77)

    repo_root = os.path.dirname(os.path.abspath(__file__))
    cir_path = os.path.join(repo_root, "spice", "ising_4pbit_core.cir")
    if not os.path.exists(cir_path):
        print(f"Error: {cir_path} not found.")
        return

    with open(cir_path, "r") as f:
        content = f.read()

    tmp_cir = os.path.join(repo_root, "temp_sim.cir")
    tmp_txt = os.path.join(repo_root, "ising_results.txt")

    idx = content.rfind(".control")
    if idx != -1:
        base = content[:idx]
    else:
        end_idx = content.rfind(".end")
        base = content[:end_idx] if end_idx != -1 else content

    sim_content = base + (
        f".control\n"
        f"  set rndseed = 42\n"
        f"  run\n"
        f"  wrdata {tmp_txt.replace(os.sep, '/')} v(m_out1) v(m_out2) v(m_out3) v(m_out4) v(a1) v(a2) v(a3) v(a4)\n"
        f"  quit\n"
        f".endc\n"
        f".end\n"
    )

    with open(tmp_cir, "w") as f:
        f.write(sim_content)

    try:
        subprocess.run(["ngspice", "-b", tmp_cir], capture_output=True, check=False)

        if not os.path.exists(tmp_txt):
            print("Simulation failed to produce results.")
            return

        data = np.loadtxt(tmp_txt)
        t_us = data[:, 0] * 1e6
        o1, o2, o3, o4 = data[:, 1], data[:, 3], data[:, 5], data[:, 7]

        t_max = t_us[-1]
        sample_points = np.linspace(t_max * 0.1, t_max, 5)

        print(
            f"{'Time (us)':<10} | {'P1 (V)':<7} | {'P2 (V)':<7} | {'P3 (V)':<7} | {'P4 (V)':<7} | {'Spin Vector (s1,s2,s3,s4)':<22} | {'Cut'}"
        )
        print("-" * 77)

        for pt in sample_points:
            idx = np.searchsorted(t_us, pt)
            if idx >= len(t_us):
                idx = len(t_us) - 1

            v1, v2, v3, v4 = o1[idx], o2[idx], o3[idx], o4[idx]
            s1 = 1 if v1 > 0.9 else -1
            s2 = 1 if v2 > 0.9 else -1
            s3 = 1 if v3 > 0.9 else -1
            s4 = 1 if v4 > 0.9 else -1

            s_vec = [s1, s2, s3, s4]
            cut_val, _ = compute_cut(s_vec)
            s_str = f"({'+1' if s1==1 else '-1'}, {'+1' if s2==1 else '-1'}, {'+1' if s3==1 else '-1'}, {'+1' if s4==1 else '-1'})"

            print(
                f"{t_us[idx]:<10.2f} | {v1:<7.3f} | {v2:<7.3f} | {v3:<7.3f} | {v4:<7.3f} | {s_str:<22} | {cut_val}/5"
            )

        print("-" * 77)
        # Final state check
        v1_f, v2_f, v3_f, v4_f = o1[-1], o2[-1], o3[-1], o4[-1]
        s_final = [
            1 if v1_f > 0.9 else -1,
            1 if v2_f > 0.9 else -1,
            1 if v3_f > 0.9 else -1,
            1 if v4_f > 0.9 else -1,
        ]
        cut_f, edges_f = compute_cut(s_final)
        final_str = f"({'+1' if s_final[0]==1 else '-1'}, {'+1' if s_final[1]==1 else '-1'}, {'+1' if s_final[2]==1 else '-1'}, {'+1' if s_final[3]==1 else '-1'})"

        is_gs = cut_f == 4
        gs_tag = " (Ground State)" if is_gs else ""
        print(f"Final State at t = {t_us[-1]:.2f} us: s = {final_str}{gs_tag}")
        edge_str = ", ".join(f"({u},{v})" for u, v in edges_f)
        print(f"Satisfied Cut Edges: {edge_str} [{cut_f} out of 5 edges cut]")

        # Dynamic trajectory stats
        all_spins = np.where(data[:, [1, 3, 5, 7]] > 0.9, 1, -1)
        flips = int(np.sum(np.any(all_spins[1:] != all_spins[:-1], axis=1)))
        cuts_all = [compute_cut(row)[0] for row in all_spins]
        gs_pct = np.mean(np.array(cuts_all) == 4) * 100.0
        c3_pct = np.mean(np.array(cuts_all) == 3) * 100.0
        print(f"Trajectory Summary: {flips} state transitions observed over {t_us[-1]:.2f} us")
        print(f"Occupation: Ground State (Cut = 4/5) = {gs_pct:.1f}% | Excited (Cut = 3/5) = {c3_pct:.1f}%\n")

    finally:
        for fpath in [tmp_cir, tmp_txt]:
            if os.path.exists(fpath):
                try:
                    os.remove(fpath)
                except OSError:
                    pass


if __name__ == "__main__":
    get_4pbit_values()
