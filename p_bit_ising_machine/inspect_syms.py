import subprocess

res = subprocess.run(['git', 'show', 'b87a2e9:p_bit_ising_machine/p_bit_ising_machine.kicad_sch'], capture_output=True, text=True)
lines = res.stdout.splitlines()

for i, line in enumerate(lines):
    if line.strip().startswith('(symbol (lib_id'):
        print(f"--- Line {i} ---")
        for j in range(i, min(i+15, len(lines))):
            print(lines[j])
            if lines[j].strip() == ')':
                break
