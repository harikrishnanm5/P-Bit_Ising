import re

with open('p_bit_ising_machine.kicad_sch', 'r') as f:
    text = f.read()

# find lib_symbols
lib_idx = text.find('(lib_symbols')
end_lib = text.find('(symbol (lib_id', lib_idx)
lib_text = text[lib_idx:end_lib]

# split by (symbol "
syms = lib_text.split('(symbol "')
for s in syms[1:]:
    name = s.split('"')[0]
    print(f"=== {name} ===")
    lines = s.splitlines()
    for l in lines:
        if '(pin ' in l:
            print(" ", l.strip())
