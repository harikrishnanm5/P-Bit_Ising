import re

with open(r'D:\E-Sim\mosfet\mosfet.kicad_sch', 'r', encoding='utf-8') as f:
    c = f.read()

syms = re.findall(r'\(symbol "([^"]+)"', c)
print('Symbols in mosfet.kicad_sch:', syms)

match_r = re.search(r'\(symbol "eSim_Devices:r"[\s\S]*?\n    \)', c)
if match_r:
    print('Found eSim_Devices:r, length:', len(match_r.group(0)))

match_c = re.search(r'\(symbol "eSim_Devices:c"[\s\S]*?\n    \)', c)
if match_c:
    print('Found eSim_Devices:c, length:', len(match_c.group(0)))
