import uuid

def gen_uuid():
    return str(uuid.uuid4())

with open('p_bit_ising_machine.kicad_sch', 'r') as f:
    orig = f.read()

lib_start = orig.find('(lib_symbols')
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

# Build each P-bit cell i (1..4) with EXACT offset dy
# Baseline Y for P1: y_base = 0
for i in range(1, 5):
    dy = (i - 1) * 120.0  # 120mm spacing between rows
    
    # Vin DC source at (38.1, 76.2 + dy)
    u_vin = gen_uuid()
    symbols.append(f'''  (symbol (lib_id "eSim_Sources:dc") (at 38.1 {76.2 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {u_vin})
    (property "Reference" "Vin{i}" (id 0) (at 45.72 {73.66 + dy} 0) (effects (font (size 1.524 1.524)) (justify left)))
    (property "Value" "0.9V" (id 1) (at 45.72 {76.2 + dy} 0) (effects (font (size 1.524 1.524)) (justify left)))
    (property "Footprint" "" (id 2) (at 38.1 {76.2 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 38.1 {76.2 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
    (pin "2" (uuid {gen_uuid()}))
  )''')

    # GND for Vin at (38.1, 88.9 + dy)
    symbols.append(f'''  (symbol (lib_id "eSim_Power:eSim_GND") (at 38.1 {88.9 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "#PWR_GND_vin{i}" (id 0) (at 38.1 {95.25 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Value" "eSim_GND" (id 1) (at 38.1 {93.98 + dy} 0))
    (property "Footprint" "" (id 2) (at 38.1 {88.9 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 38.1 {88.9 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
  )''')

    # Vnoise pulse source at (63.5, 76.2 + dy)
    symbols.append(f'''  (symbol (lib_id "eSim_Sources:pulse") (at 63.5 {76.2 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "Vnoise{i}" (id 0) (at 71.12 {73.66 + dy} 0) (effects (font (size 1.524 1.524)) (justify left)))
    (property "Value" "Noise_{i}" (id 1) (at 71.12 {76.2 + dy} 0) (effects (font (size 1.524 1.524)) (justify left)))
    (property "Footprint" "" (id 2) (at 63.5 {76.2 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 63.5 {76.2 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
    (pin "2" (uuid {gen_uuid()}))
  )''')

    # Rfb at (101.6, 38.1 + dy, 90 deg)
    symbols.append(f'''  (symbol (lib_id "Device:R") (at 101.6 {38.1 + dy} 90) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "Rfb{i}" (id 0) (at 101.6 {30.48 + dy} 90) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Value" "500k" (id 1) (at 101.6 {33.02 + dy} 90) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Footprint" "" (id 2) (at 101.6 {38.1 + dy} 90) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 101.6 {38.1 + dy} 90) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
    (pin "2" (uuid {gen_uuid()}))
  )''')

    # MP1 at (110.49, 55.88 + dy)
    symbols.append(f'''  (symbol (lib_id "eSim_Devices:mosfet_p") (at 110.49 {55.88 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "MP{i}" (id 0) (at 118.11 {53.34 + dy} 0) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Value" "W=2.5u L=0.18u" (id 1) (at 118.11 {55.88 + dy} 0) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Footprint" "" (id 2) (at 110.49 {55.88 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 110.49 {55.88 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
    (pin "2" (uuid {gen_uuid()}))
    (pin "3" (uuid {gen_uuid()}))
    (pin "4" (uuid {gen_uuid()}))
  )''')

    # MN1 at (109.22, 121.92 + dy)
    symbols.append(f'''  (symbol (lib_id "eSim_Devices:mosfet_n") (at 109.22 {121.92 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "MN{i}" (id 0) (at 118.11 {120.65 + dy} 0) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Value" "W=1.0u L=0.18u" (id 1) (at 118.11 {123.19 + dy} 0) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Footprint" "" (id 2) (at 109.22 {121.92 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 109.22 {121.92 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
    (pin "2" (uuid {gen_uuid()}))
    (pin "3" (uuid {gen_uuid()}))
    (pin "4" (uuid {gen_uuid()}))
  )''')

    # C1 at (139.7, 101.6 + dy)
    symbols.append(f'''  (symbol (lib_id "Device:C") (at 139.7 {101.6 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "C{i}" (id 0) (at 146.05 {99.06 + dy} 0) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Value" "15fF" (id 1) (at 146.05 {101.6 + dy} 0) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Footprint" "" (id 2) (at 139.7 {101.6 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 139.7 {101.6 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
    (pin "2" (uuid {gen_uuid()}))
  )''')

    # GND for C1 at (139.7, 109.22 + dy)
    symbols.append(f'''  (symbol (lib_id "eSim_Power:eSim_GND") (at 139.7 {109.22 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "#PWR_GND_c{i}" (id 0) (at 139.7 {115.57 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Value" "eSim_GND" (id 1) (at 139.7 {114.3 + dy} 0))
    (property "Footprint" "" (id 2) (at 139.7 {109.22 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 139.7 {109.22 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
  )''')

    # GND for MN1 at (114.3, 129.54 + dy)
    symbols.append(f'''  (symbol (lib_id "eSim_Power:eSim_GND") (at 114.3 {129.54 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "#PWR_GND_mn{i}" (id 0) (at 114.3 {135.89 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Value" "eSim_GND" (id 1) (at 114.3 {134.62 + dy} 0))
    (property "Footprint" "" (id 2) (at 114.3 {129.54 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 114.3 {129.54 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
  )''')

    # MP_BUF at (186.69, 55.88 + dy)
    symbols.append(f'''  (symbol (lib_id "eSim_Devices:mosfet_p") (at 186.69 {55.88 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "MP_BUF{i}" (id 0) (at 194.31 {53.34 + dy} 0) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Value" "W=3.0u L=0.18u" (id 1) (at 194.31 {55.88 + dy} 0) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Footprint" "" (id 2) (at 186.69 {55.88 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 186.69 {55.88 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
    (pin "2" (uuid {gen_uuid()}))
    (pin "3" (uuid {gen_uuid()}))
    (pin "4" (uuid {gen_uuid()}))
  )''')

    # MN_BUF at (185.42, 121.92 + dy)
    symbols.append(f'''  (symbol (lib_id "eSim_Devices:mosfet_n") (at 185.42 {121.92 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "MN_BUF{i}" (id 0) (at 194.31 {120.65 + dy} 0) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Value" "W=1.5u L=0.18u" (id 1) (at 194.31 {123.19 + dy} 0) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Footprint" "" (id 2) (at 185.42 {121.92 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 185.42 {121.92 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
    (pin "2" (uuid {gen_uuid()}))
    (pin "3" (uuid {gen_uuid()}))
    (pin "4" (uuid {gen_uuid()}))
  )''')

    # GND for MN_BUF at (190.5, 129.54 + dy)
    symbols.append(f'''  (symbol (lib_id "eSim_Power:eSim_GND") (at 190.5 {129.54 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "#PWR_GND_buf{i}" (id 0) (at 190.5 {135.89 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Value" "eSim_GND" (id 1) (at 190.5 {134.62 + dy} 0))
    (property "Footprint" "" (id 2) (at 190.5 {129.54 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 190.5 {129.54 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
  )''')

    # Cout at (203.2, 101.6 + dy)
    symbols.append(f'''  (symbol (lib_id "Device:C") (at 203.2 {101.6 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "Cout{i}" (id 0) (at 209.55 {99.06 + dy} 0) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Value" "20fF" (id 1) (at 209.55 {101.6 + dy} 0) (effects (font (size 1.27 1.27)) (justify left)))
    (property "Footprint" "" (id 2) (at 203.2 {101.6 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 203.2 {101.6 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
    (pin "2" (uuid {gen_uuid()}))
  )''')

    # GND for Cout at (203.2, 109.22 + dy)
    symbols.append(f'''  (symbol (lib_id "eSim_Power:eSim_GND") (at 203.2 {109.22 + dy} 0) (unit 1)
    (in_bom yes) (on_board yes) (fields_autoplaced)
    (uuid {gen_uuid()})
    (property "Reference" "#PWR_GND_cout{i}" (id 0) (at 203.2 {115.57 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Value" "eSim_GND" (id 1) (at 203.2 {114.3 + dy} 0))
    (property "Footprint" "" (id 2) (at 203.2 {109.22 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (property "Datasheet" "" (id 3) (at 203.2 {109.22 + dy} 0) (effects (font (size 1.27 1.27)) hide))
    (pin "1" (uuid {gen_uuid()}))
  )''')

    # Perfect Wires for Cell i
    add_wire((38.1, 88.9 + dy), (38.1, 76.2 + dy))
    add_wire((38.1, 63.5 + dy), (63.5, 63.5 + dy))
    add_wire((63.5, 63.5 + dy), (63.5, 76.2 + dy))
    add_wire((63.5, 88.9 + dy), (88.9, 88.9 + dy))
    add_wire((88.9, 88.9 + dy), (88.9, 38.1 + dy))
    add_wire((88.9, 38.1 + dy), (97.79, 38.1 + dy))
    add_wire((88.9, 88.9 + dy), (88.9, 55.88 + dy))
    add_wire((88.9, 55.88 + dy), (106.68, 55.88 + dy))
    add_wire((88.9, 88.9 + dy), (88.9, 121.92 + dy))
    add_wire((88.9, 121.92 + dy), (106.68, 121.92 + dy))
    add_wire((105.41, 38.1 + dy), (127.0, 38.1 + dy))
    add_wire((127.0, 38.1 + dy), (127.0, 88.9 + dy))
    add_wire((114.3, 50.8 + dy), (114.3, 40.64 + dy))
    add_wire((114.3, 60.96 + dy), (114.3, 88.9 + dy))
    add_wire((114.3, 88.9 + dy), (114.3, 121.92 + dy))
    add_wire((114.3, 88.9 + dy), (139.7, 88.9 + dy))
    add_wire((139.7, 88.9 + dy), (139.7, 97.79 + dy))
    add_wire((139.7, 105.41 + dy), (139.7, 109.22 + dy))
    add_wire((139.7, 88.9 + dy), (165.1, 88.9 + dy))
    add_wire((165.1, 88.9 + dy), (165.1, 55.88 + dy))
    add_wire((165.1, 55.88 + dy), (182.88, 55.88 + dy))
    add_wire((165.1, 88.9 + dy), (165.1, 121.92 + dy))
    add_wire((165.1, 121.92 + dy), (182.88, 121.92 + dy))
    add_wire((190.5, 50.8 + dy), (190.5, 40.64 + dy))
    add_wire((190.5, 60.96 + dy), (190.5, 88.9 + dy))
    add_wire((190.5, 88.9 + dy), (190.5, 121.92 + dy))
    add_wire((190.5, 88.9 + dy), (203.2, 88.9 + dy))
    add_wire((203.2, 88.9 + dy), (203.2, 97.79 + dy))
    add_wire((203.2, 105.41 + dy), (203.2, 109.22 + dy))
    add_wire((203.2, 88.9 + dy), (228.6, 88.9 + dy))

    # Global Labels for Cell i
    add_glabel(f"in{i}", (50.8, 63.5 + dy), shape="input")
    add_glabel(f"noise{i}", (76.2, 88.9 + dy), shape="input")
    add_glabel(f"a{i}", (133.35, 88.9 + dy), shape="bidirectional")
    add_glabel(f"m_out{i}", (228.6, 88.9 + dy), shape="output")
    add_glabel("VDD", (114.3, 40.64 + dy), shape="input", angle=90)
    add_glabel("VDD", (190.5, 40.64 + dy), shape="input", angle=90)

sch_lines.extend(symbols)
sch_lines.extend(wires)
sch_lines.extend(labels)
sch_lines.append('  (sheet_instances\n    (path "/" (page "1"))\n  )\n)')

full_content = '\n'.join(sch_lines)

with open('p_bit_ising_machine.kicad_sch', 'w') as f:
    f.write(full_content)

print('Generated clean 4-p-bit schematic!')
