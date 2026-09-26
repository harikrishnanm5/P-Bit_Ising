import uuid
import re

def gen_uuid():
    return str(uuid.uuid4())

# Read existing lib_symbols
with open('p_bit_ising_machine.kicad_sch', 'r') as f:
    orig = f.read()

lib_start = orig.find('(lib_symbols')
# find matching closing bracket for lib_symbols
depth = 0
lib_end = -1
for i in range(lib_start, len(orig)):
    if orig[i] == '(':
        depth += 1
    elif orig[i] == ')':
        depth -= 1
        if depth == 0:
            lib_end = i + 1
            break

lib_symbols_block = orig[lib_start:lib_end]

sch_lines = [
    '(kicad_sch (version 20211123) (generator eeschema)',
    f'  (uuid {gen_uuid()})',
    '  (paper "A2")',
    '',
    lib_symbols_block,
    ''
]

symbols = []
wires = []
labels = []

def add_wire(p1, p2):
    wires.append(f'  (wire (pts (xy {p1[0]} {p1[1]}) (xy {p2[0]} {p2[1]}))\n    (stroke (width 0) (type default) (color 0 0 0 0))\n    (uuid {gen_uuid()})\n  )')

def add_glabel(text, pos, shape="bidirectional", angle=0):
    labels.append(f'  (global_label "{text}" (shape {shape}) (at {pos[0]} {pos[1]} {angle}) (fields_autoplaced)\n    (effects (font (size 1.27 1.27)) (justify left))\n    (uuid {gen_uuid()})\n  )')

def add_symbol(lib_id, ref, val, pos, angle=0, mirror=False, unit=1):
    u = gen_uuid()
    p1 = gen_uuid()
    p2 = gen_uuid()
    p3 = gen_uuid()
    p4 = gen_uuid()
    
    # mirror/rotation effects if any
    rot_str = f'{angle}'
    sym_str = f'  (symbol (lib_id "{lib_id}") (at {pos[0]} {pos[1]} {rot_str}) (unit {unit})\n'
    sym_str += f'    (in_bom yes) (on_board yes) (fields_autoplaced)\n'
    sym_str += f'    (uuid {u})\n'
    sym_str += f'    (property "Reference" "{ref}" (id 0) (at {pos[0]} {pos[1]-2.54} 0) (effects (font (size 1.27 1.27))))\n'
    sym_str += f'    (property "Value" "{val}" (id 1) (at {pos[0]} {pos[1]+2.54} 0) (effects (font (size 1.27 1.27))))\n'
    sym_str += f'    (property "Footprint" "" (id 2) (at {pos[0]} {pos[1]} 0) (effects (font (size 1.27 1.27)) hide))\n'
    sym_str += f'    (property "Datasheet" "" (id 3) (at {pos[0]} {pos[1]} 0) (effects (font (size 1.27 1.27)) hide))\n'
    sym_str += f'    (pin "1" (uuid {p1}))\n'
    sym_str += f'    (pin "2" (uuid {p2}))\n'
    if 'mosfet' in lib_id:
        sym_str += f'    (pin "3" (uuid {p3}))\n'
        sym_str += f'    (pin "4" (uuid {p4}))\n'
    sym_str += '  )'
    symbols.append(sym_str)

# Build 4 p-bits
for i in range(1, 5):
    y_off = 50 + (i - 1) * 75
    
    # Vin (DC source at 38.1, y_off + 25.4)
    add_symbol("eSim_Sources:dc", f"Vin{i}", "0.9V", (38.1, y_off + 25.4))
    add_symbol("eSim_Power:eSim_GND", f"#PWR_GND_in{i}", "eSim_GND", (38.1, y_off + 35.56))
    add_wire((38.1, y_off + 30.48), (38.1, y_off + 35.56))
    
    # Vnoise (Pulse / AC source at 63.5, y_off + 25.4)
    add_symbol("eSim_Sources:pulse", f"Vnoise{i}", f"Noise_{i}", (63.5, y_off + 25.4))
    
    # Connect Vin to Vnoise
    add_wire((38.1, y_off + 20.32), (38.1, y_off + 15.24))
    add_wire((38.1, y_off + 15.24), (63.5, y_off + 15.24))
    add_wire((63.5, y_off + 15.24), (63.5, y_off + 20.32))
    
    # Noise output node
    add_wire((63.5, y_off + 30.48), (88.9, y_off + 30.48))
    add_glabel(f"noise{i}", (75.0, y_off + 30.48))
    
    # Stage 1: MP1 / MN1 + Rfb + C1
    # Rfb (114.3, y_off - 10)
    add_symbol("Device:R", f"Rfb{i}", "500k", (114.3, y_off - 10), angle=90)
    add_wire((88.9, y_off + 30.48), (88.9, y_off - 10))
    add_wire((88.9, y_off - 10), (106.68, y_off - 10))
    add_wire((121.92, y_off - 10), (139.7, y_off - 10))
    add_wire((139.7, y_off - 10), (139.7, y_off + 30.48))
    
    # MP1 (PMOS at 114.3, y_off + 10)
    add_symbol("eSim_Devices:mosfet_p", f"MP{i}", "W=2.5u L=0.18u", (114.3, y_off + 10))
    add_symbol("eSim_Power:eSim_GND", f"#PWR_VDD1_{i}", "eSim_GND", (114.3, y_off - 2.54)) # VDD label instead
    add_glabel("VDD", (114.3, y_off - 2.54), angle=90)
    add_wire((114.3, y_off + 4.92), (114.3, y_off - 2.54))
    
    # MN1 (NMOS at 114.3, y_off + 45)
    add_symbol("eSim_Devices:mosfet_n", f"MN{i}", "W=1.0u L=0.18u", (114.3, y_off + 45))
    add_symbol("eSim_Power:eSim_GND", f"#PWR_GND_mn{i}", "eSim_GND", (114.3, y_off + 53.34))
    add_wire((114.3, y_off + 50.08), (114.3, y_off + 53.34))
    
    # Gate connections
    add_wire((88.9, y_off + 30.48), (88.9, y_off + 10))
    add_wire((88.9, y_off + 10), (106.68, y_off + 10))
    add_wire((88.9, y_off + 30.48), (88.9, y_off + 45))
    add_wire((88.9, y_off + 45), (106.68, y_off + 45))
    
    # Drain connections to a_i (analog_out)
    add_wire((114.3, y_off + 15.08), (114.3, y_off + 30.48))
    add_wire((114.3, y_off + 39.92), (114.3, y_off + 30.48))
    add_wire((114.3, y_off + 30.48), (139.7, y_off + 30.48))
    
    # C1 (Capacitor to GND at 139.7, y_off + 45)
    add_symbol("Device:C", f"C{i}", "15fF", (139.7, y_off + 45))
    add_symbol("eSim_Power:eSim_GND", f"#PWR_GND_c{i}", "eSim_GND", (139.7, y_off + 53.34))
    add_wire((139.7, y_off + 30.48), (139.7, y_off + 41.18))
    add_wire((139.7, y_off + 48.82), (139.7, y_off + 53.34))
    
    add_glabel(f"a{i}", (145.0, y_off + 30.48))
    add_wire((139.7, y_off + 30.48), (165.1, y_off + 30.48))
    
    # Stage 2: Buffer MP_BUF / MN_BUF
    add_symbol("eSim_Devices:mosfet_p", f"MP_BUF{i}", "W=3.0u L=0.18u", (190.5, y_off + 10))
    add_glabel("VDD", (190.5, y_off - 2.54), angle=90)
    add_wire((190.5, y_off + 4.92), (190.5, y_off - 2.54))
    
    add_symbol("eSim_Devices:mosfet_n", f"MN_BUF{i}", "W=1.5u L=0.18u", (190.5, y_off + 45))
    add_symbol("eSim_Power:eSim_GND", f"#PWR_GND_buf{i}", "eSim_GND", (190.5, y_off + 53.34))
    add_wire((190.5, y_off + 50.08), (190.5, y_off + 53.34))
    
    # Buffer gates
    add_wire((165.1, y_off + 30.48), (165.1, y_off + 10))
    add_wire((165.1, y_off + 10), (182.88, y_off + 10))
    add_wire((165.1, y_off + 30.48), (165.1, y_off + 45))
    add_wire((165.1, y_off + 45), (182.88, y_off + 45))
    
    # Buffer output
    add_wire((190.5, y_off + 15.08), (190.5, y_off + 30.48))
    add_wire((190.5, y_off + 39.92), (190.5, y_off + 30.48))
    add_wire((190.5, y_off + 30.48), (215.9, y_off + 30.48))
    
    # Cout (215.9, y_off + 45)
    add_symbol("Device:C", f"Cout{i}", "20fF", (215.9, y_off + 45))
    add_symbol("eSim_Power:eSim_GND", f"#PWR_GND_cout{i}", "eSim_GND", (215.9, y_off + 53.34))
    add_wire((215.9, y_off + 30.48), (215.9, y_off + 41.18))
    add_wire((215.9, y_off + 48.82), (215.9, y_off + 53.34))
    
    # Digital Out Pin
    add_wire((215.9, y_off + 30.48), (235.0, y_off + 30.48))
    add_glabel(f"m_out{i}", (235.0, y_off + 30.48), shape="output")

# --- Mode Controller Block (Top Right: X=270, Y=30) ---
add_symbol("eSim_Sources:dc", "Vmode", "1.8V", (270.0, 30.0))
add_symbol("eSim_Power:eSim_GND", "#PWR_GND_vmode", "eSim_GND", (270.0, 40.16))
add_wire((270.0, 35.08), (270.0, 40.16))
add_wire((270.0, 24.92), (270.0, 20.0))
add_glabel("mode", (270.0, 20.0))

# MODE Inverter
add_symbol("eSim_Devices:mosfet_p", "MP_INV", "W=2.5u L=0.18u", (295.0, 20.0))
add_glabel("VDD", (295.0, 10.0), angle=90)
add_wire((295.0, 14.92), (295.0, 10.0))

add_symbol("eSim_Devices:mosfet_n", "MN_INV", "W=1.0u L=0.18u", (295.0, 45.0))
add_symbol("eSim_Power:eSim_GND", "#PWR_GND_minv", "eSim_GND", (295.0, 53.34))
add_wire((295.0, 50.08), (295.0, 53.34))

add_wire((285.0, 30.0), (285.0, 20.0))
add_wire((285.0, 20.0), (287.38, 20.0))
add_wire((285.0, 30.0), (285.0, 45.0))
add_wire((285.0, 45.0), (287.38, 45.0))
add_glabel("mode", (280.0, 30.0))

add_wire((295.0, 25.08), (295.0, 32.5))
add_wire((295.0, 39.92), (295.0, 32.5))
add_wire((295.0, 32.5), (305.0, 32.5))
add_glabel("mode_bar", (305.0, 32.5))

# --- Coupling Network Blocks (Right Side: X=270..380) ---
# Edge 1: P1-P2 (Always ON Resistors)
add_symbol("Device:R", "R_c12", "500k", (280.0, 80.0), angle=90)
add_glabel("a1", (268.0, 80.0))
add_wire((268.0, 80.0), (272.38, 80.0))
add_wire((287.62, 80.0), (295.0, 80.0))
add_glabel("noise2", (295.0, 80.0))

add_symbol("Device:R", "R_c21", "500k", (280.0, 100.0), angle=90)
add_glabel("a2", (268.0, 100.0))
add_wire((268.0, 100.0), (272.38, 100.0))
add_wire((287.62, 100.0), (295.0, 100.0))
add_glabel("noise1", (295.0, 100.0))

# Switchable Edges with Transmission Gates
switch_edges = [
    ("13", "a1", "noise3", "31", "a3", "noise1", 130),
    ("23", "a2", "noise3", "32", "a3", "noise2", 180),
    ("24", "a2", "noise4", "42", "a4", "noise2", 230),
    ("34", "a3", "noise4", "43", "a4", "noise3", 280),
]

for tag_f, src_f, dst_f, tag_r, src_r, dst_r, y_pos in switch_edges:
    # Forward direction
    add_symbol("Device:R", f"R_c{tag_f}", "500k", (270.0, y_pos), angle=90)
    add_glabel(src_f, (258.0, y_pos))
    add_wire((258.0, y_pos), (262.38, y_pos))
    add_wire((277.62, y_pos), (285.0, y_pos))
    
    # TG NMOS + PMOS
    add_symbol("eSim_Devices:mosfet_n", f"MN_TG{tag_f}", "W=1.0u", (295.0, y_pos + 8))
    add_symbol("eSim_Devices:mosfet_p", f"MP_TG{tag_f}", "W=2.5u", (295.0, y_pos - 8))
    add_glabel("mode", (285.0, y_pos + 8))
    add_glabel("mode_bar", (285.0, y_pos - 8))
    
    add_wire((285.0, y_pos), (285.0, y_pos - 3))
    add_wire((285.0, y_pos), (285.0, y_pos + 13))
    add_wire((285.0, y_pos - 3), (287.38, y_pos - 3))
    add_wire((285.0, y_pos + 13), (287.38, y_pos + 13))
    
    add_wire((305.0, y_pos), (305.0, y_pos - 3))
    add_wire((305.0, y_pos), (305.0, y_pos + 13))
    add_wire((302.62, y_pos - 3), (305.0, y_pos - 3))
    add_wire((302.62, y_pos + 13), (305.0, y_pos + 13))
    add_glabel(dst_f, (305.0, y_pos))
    
    # Reverse direction
    add_symbol("Device:R", f"R_c{tag_r}", "500k", (335.0, y_pos), angle=90)
    add_glabel(src_r, (323.0, y_pos))
    add_wire((323.0, y_pos), (327.38, y_pos))
    add_wire((342.62, y_pos), (350.0, y_pos))
    
    # TG NMOS + PMOS
    add_symbol("eSim_Devices:mosfet_n", f"MN_TG{tag_r}", "W=1.0u", (360.0, y_pos + 8))
    add_symbol("eSim_Devices:mosfet_p", f"MP_TG{tag_r}", "W=2.5u", (360.0, y_pos - 8))
    add_glabel("mode", (350.0, y_pos + 8))
    add_glabel("mode_bar", (350.0, y_pos - 8))
    
    add_wire((350.0, y_pos), (350.0, y_pos - 3))
    add_wire((350.0, y_pos), (350.0, y_pos + 13))
    add_wire((350.0, y_pos - 3), (352.38, y_pos - 3))
    add_wire((350.0, y_pos + 13), (352.38, y_pos + 13))
    
    add_wire((370.0, y_pos), (370.0, y_pos - 3))
    add_wire((370.0, y_pos), (370.0, y_pos + 13))
    add_wire((367.62, y_pos - 3), (370.0, y_pos - 3))
    add_wire((367.62, y_pos + 13), (370.0, y_pos + 13))
    add_glabel(dst_r, (370.0, y_pos))

sch_lines.extend(symbols)
sch_lines.extend(wires)
sch_lines.extend(labels)
sch_lines.append('  (sheet_instances\n    (path "/" (page "1"))\n  )\n)')

full_content = '\n'.join(sch_lines)

with open('p_bit_ising_machine.kicad_sch', 'w') as f:
    f.write(full_content)

print('Successfully generated complete 4-p-bit schematic! Total size:', len(full_content))
