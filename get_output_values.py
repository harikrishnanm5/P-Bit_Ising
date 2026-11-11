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

    sim_content = content.replace(
        ".control",
        f".control\n  run\n  wrdata {tmp_txt.replace(os.sep, '/')} v(m_out1) v(m_out2) v(m_out3) v(m_out4) v(a1) v(a2) v(a3) v(a4)\n  quit",
    )

    with open(tmp_cir, "w") as f:
        f.write(sim_content)

    try:
        subprocess.run(["ngspice", "-b", tmp_cir], capture_output=True, check=False)

        if not os.path.exists(tmp_txt):
            print("Simulation failed to produce results.")
            return

        data = np.loadtxt(tmp_txt)
        t_ms = data[:, 0] * 1e3
        o1, o2, o3, o4 = data[:, 1], data[:, 3], data[:, 5], data[:, 7]

        sample_points = [0.10, 0.25, 0.50, 0.75, 1.00]

        print(
            f"{'Time (ms)':<10} | {'P1 (V)':<7} | {'P2 (V)':<7} | {'P3 (V)':<7} | {'P4 (V)':<7} | {'Spin Vector (s1,s2,s3,s4)':<22} | {'Cut'}"
        )
        print("-" * 77)

        final_spins = None
        for sp in sample_points:
            idx = int(np.argmin(np.abs(t_ms - sp)))
            voltages = [o1[idx], o2[idx], o3[idx], o4[idx]]
            spins = tuple(1 if v > 0.9 else -1 for v in voltages)
            cut_cnt, _ = compute_cut(spins)
            spin_str = "(" + ", ".join(f"{s:+d}" for s in spins) + ")"
            print(
                f"{t_ms[idx]:<10.2f} | {voltages[0]:<7.3f} | {voltages[1]:<7.3f} | {voltages[2]:<7.3f} | {voltages[3]:<7.3f} | {spin_str:<22} | {cut_cnt}/5"
            )
            final_spins = spins

        cut_cnt, satisfied = compute_cut(final_spins)
        spin_str = "(" + ", ".join(f"{s:+d}" for s in final_spins) + ")"
        edges_str = ", ".join(f"({i},{j})" for i, j in satisfied)

        all_spins = np.where(np.column_stack([o1, o2, o3, o4]) > 0.9, 1, -1)
        all_cuts = np.sum(
            [(all_spins[:, i] != all_spins[:, j]).astype(int) for i, j in GRAPH_EDGES],
            axis=0,
        )
        gs_pct = 100.0 * float(np.mean(all_cuts == 4))

        print("-" * 77)
        print(f"Final State at t = {t_ms[-1]:.2f} ms: s = {spin_str} ({'Ground State' if cut_cnt == 4 else 'Excited State'})")
        print(f"Satisfied Cut Edges: {edges_str} [{cut_cnt} out of 5 edges cut]")
        print(f"Ground-State Residence (Cut = 4/5 over 1.00 ms): {gs_pct:.1f}%\n")
    finally:
        for tmp_file in [tmp_cir, tmp_txt]:
            if os.path.exists(tmp_file):
                try:
                    os.remove(tmp_file)
                except Exception:
                    pass


if __name__ == "__main__":
    get_4pbit_values()
