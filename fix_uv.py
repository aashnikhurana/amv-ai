import re

with open("app.py", "r") as f:
    content = f.read()

old = '''def interpret_uv(hbd, hba, smiles):
    smiles = smiles or ""
    aromatic = any(c in smiles for c in ['c', 'n', 'o', 's'] if c == c.lower() and c in smiles)
    has_carbonyl = "C=O" in smiles or "C(=O)" in smiles
    if aromatic:
        return "Aromatic ring detected — start at 254 nm"
    elif has_carbonyl:
        return "Carbonyl group detected — try 210 nm (low wavelength, more background noise)"
    else:
        return "⚠️ No UV chromophore detected — UV detection may not work. Consider MS, ELSD, or CAD detector."'''

new = '''def interpret_uv(hbd, hba, smiles):
    smiles = smiles or ""

    # Check for aromatic rings (lowercase letters in SMILES = aromatic)
    aromatic = any(c in smiles for c in ["c", "n", "o", "s"] if c == c.lower())

    # Check for carbonyl
    has_carbonyl = "C=O" in smiles or "C(=O)" in smiles

    # Check for sulfur-sulfur bonds (disulfide, trisulfide)
    # SS appears in SMILES for S-S linkages
    has_sulfur_bond = "SS" in smiles or smiles.count("S") >= 2

    # Check for nitro groups
    has_nitro = "N(=O)" in smiles or "[N+](=O)" in smiles

    # Check for conjugated systems (C=C-C=O or similar)
    has_conjugated = "C=CC=O" in smiles or "C=Cc" in smiles

    if aromatic:
        return "Aromatic ring detected — start at 254 nm. Also try 220 nm for higher sensitivity."
    elif has_nitro:
        return "Nitro group detected — absorbs at 254 nm and 330 nm. Try both."
    elif has_conjugated:
        return "Conjugated system detected — try 280–320 nm range."
    elif has_carbonyl:
        return "Carbonyl group detected — try 210–215 nm. Expect some background noise at low wavelength."
    elif has_sulfur_bond:
        return "Sulfur-sulfur bond detected (disulfide/trisulfide) — absorbs at 210–220 nm. Your compound should be detectable around 215 nm."
    else:
        return "⚠️ No common UV chromophore detected — UV detection may not work reliably. Consider MS, ELSD, or CAD detector."'''

content = content.replace(old, new)

with open("app.py", "w") as f:
    f.write(content)

print("Fixed.")
