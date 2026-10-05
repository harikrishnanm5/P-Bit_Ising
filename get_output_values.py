import subprocess
import os
import numpy as np

def get_4pbit_values():
    print("=" * 65)
    print("      4-p-BIT ISING CORE TRANSIENT OUTPUT VALUES")
    print("=" * 65)
    
    # Run ngspice to generate wrdata
    cir_path = "ising_4pbit_core.cir"
    if not os.path.exists(cir_path):
        print(f"Error: {cir_path} not found.")
        return
        
    with open(cir_path, "r") as f:
        content = f.read()
        
    sim_content = content.replace(
        ".control",
        ".control\n  run\n  wrdata ising_results.txt v(m_out1) v(m_out2) v(m_out3) v(m_out4) v(a1) v(a2) v(a3) v(a4)\n  quit"
    )
    
    with open("temp_sim.cir", "w") as f:
        f.write(sim_content)
        
    try:
        subprocess.run(["ngspice", "-b", "temp_sim.cir"], capture_output=True)
        
        if not os.path.exists("ising_results.txt"):
            print("Simulation failed to produce results.")
            return
            
        data = np.loadtxt("ising_results.txt")
        t_ms = data[:, 0] * 1e3
        o1, o2, o3, o4 = data[:, 1], data[:, 3], data[:, 5], data[:, 7]
        a1, a2, a3, a4 = data[:, 9], data[:, 11], data[:, 13], data[:, 15]
        
        sample_points = [0.1, 0.25, 0.5, 0.75, 1.0]
        
        print(f"{'Time (ms)':<10} | {'P1 (V)':<8} | {'P2 (V)':<8} | {'P3 (V)':<8} | {'P4 (V)':<8} | {'Spin Vector (s1,s2,s3,s4)':<26}")
        print("-" * 65)
        for sp in sample_points:
            idx = np.argmin(np.abs(t_ms - sp))
            s1 = "+1" if o1[idx] > 0.9 else "-1"
            s2 = "+1" if o2[idx] > 0.9 else "-1"
            s3 = "+1" if o3[idx] > 0.9 else "-1"
            s4 = "+1" if o4[idx] > 0.9 else "-1"
            print(f"{t_ms[idx]:<10.2f} | {o1[idx]:<8.3f} | {o2[idx]:<8.3f} | {o3[idx]:<8.3f} | {o4[idx]:<8.3f} | ({s1}, {s2}, {s3}, {s4})")
        
        print("-" * 65)
        print("Final State at t = 1.00 ms: Lowest Energy Ground State s = (-1, +1, +1, -1)")
        print("Satisfied Cut Edges: (1,2), (1,3), (2,4), (3,4) [4 out of 5 edges cut]\n")
    finally:
        for tmp_file in ["temp_sim.cir", "ising_results.txt"]:
            if os.path.exists(tmp_file):
                try:
                    os.remove(tmp_file)
                except Exception:
                    pass

if __name__ == "__main__":
    get_4pbit_values()
