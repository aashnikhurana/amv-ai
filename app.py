import streamlit as st
import requests
import numpy as np

st.set_page_config(page_title="AMV-AI", page_icon="🔬", layout="centered")

st.markdown("""
<style>
    .stApp { background-color: #F6FAFD; }
    .stApp, .stMarkdown, p, label { color: #1A3D63 !important; }
    h1, h2, h3 { color: #0A1931 !important; }
    h1 { font-size: 2.2rem !important; }
    .stButton > button, .stDownloadButton > button {
        background-color: #B3CFE5 !important;
        color: #1A3D63 !important;
        border: none !important;
        border-radius: 8px !important;
        padding: 0.5rem 2rem !important;
        font-size: 1rem !important;
        transition: background-color 0.3s ease, color 0.3s ease;
    }
    .stButton > button p, .stDownloadButton > button p {
        color: #1A3D63 !important;
        transition: color 0.3s ease;
    }
    .stButton > button:hover, .stDownloadButton > button:hover { 
        background-color: #4A7FA7 !important; 
    }
    .stButton > button:hover p, .stDownloadButton > button:hover p {
        color: white !important;
    }
    .stRadio > label { color: #1A3D63 !important; font-weight: 500 !important; }
    hr { border-color: #B3CFE5 !important; }
    .stProgress > div > div { background-color: #1A3D63 !important; }
    [data-testid="stMetricValue"] { 
        color: #0A1931 !important; 
        font-weight: 700 !important; 
        justify-content: center;
    }
    [data-testid="stMetricLabel"] { 
        color: #1A3D63 !important; 
        justify-content: center;
    }
    [data-testid="stMetric"] {
        display: flex;
        flex-direction: column;
        align-items: center;
    }
    
    .stage-card, .prop-card, .method-card {
        background-color: transparent;
        border-radius: 8px;
    }
    .stage-card {
        border: 1px solid #B3CFE5;
        border-left: 4px solid #1A3D63;
        padding: 1.2rem 1.5rem;
        margin: 1rem 0;
    }
    .prop-card {
        padding: 1rem 1.5rem;
        margin: 0.5rem 0;
        border: 1px solid #B3CFE5;
    }
    .method-card {
        padding: 0.8rem 1rem;
        margin: -0.5rem 0 0.5rem 0;
    }
    .method-card small {
        color: #1A3D63;
    }
    
    [data-testid="stExpander"] {
        background-color: #F6FAFD !important;
        border: 1px solid #B3CFE5 !important;
        border-radius: 8px !important;
    }
    [data-testid="stExpander"] summary {
        background-color: #E9F1F6 !important;
    }
    [data-testid="stExpander"] summary p, [data-testid="stExpander"] summary span {
        color: #0A1931 !important;
        font-weight: 600 !important;
    }
    [data-testid="stExpander"] div[role="region"] {
        background-color: #F6FAFD !important;
    }
</style>
""", unsafe_allow_html=True)

for key, default in {
    "stage": 0,
    "qualification": None,
    "method_type": None,
    "compound": None,
    "properties": None
}.items():
    if key not in st.session_state:
        st.session_state[key] = default

st.title("🔬 AMV-AI")
st.caption("Analytical Method Lifecycle Tool")
st.progress(st.session_state.stage / 8)
st.caption(f"Stage {st.session_state.stage} of 8")

st.markdown("""
<div style="background-color:#FFF3CD; border:1px solid #FFC107; border-radius:6px; 
padding:0.6rem 1rem; margin-bottom:0.5rem; font-size:0.85rem; color:#856404;">
⚠️ <b>Research & Guidance Tool Only.</b> AMV-AI is not a validated analytical system. 
All outputs must be reviewed by a qualified analytical chemist. Not suitable for 
direct regulatory submission without independent experimental verification.
</div>
""", unsafe_allow_html=True)

st.divider()

def go_next():
    st.session_state.stage += 1

def go_back():
    st.session_state.stage -= 1

def fetch_compound(name):
    try:
        props = "MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,IUPACName,MolecularFormula,IsomericSMILES"
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/{props}/JSON"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            p = r.json()["PropertyTable"]["Properties"][0]
            return {
                "name": name,
                "iupac": p.get("IUPACName", name),
                "formula": p.get("MolecularFormula", "N/A"),
                "mw": p.get("MolecularWeight", "N/A"),
                "logp": p.get("XLogP", "N/A"),
                "tpsa": p.get("TPSA", "N/A"),
                "hbd": p.get("HBondDonorCount", "N/A"),
                "hba": p.get("HBondAcceptorCount", "N/A"),
                "smiles": p.get("IsomericSMILES", ""),
            }
        return None
    except:
        return None

import re

from rules import CHROMOPHORE_RULES, COLUMN_RULES, SAMPLE_PREP_RULES, VALIDATION_SST_RULES

def interpret_uv(smiles, formula="", extra_hints=None):
    """
    Analyzes chemical properties to identify chromophores using the knowledge base.
    Returns:
        html_output (str): Formatted HTML bullet points for the UI.
        features (dict): Normalized boolean flags for downstream logic.
    """
    smiles = smiles or ""
    formula = formula or ""
    hints = (extra_hints or "").lower()

    # 1. Normalize features
    features = {
        "aromatic_ring": False,
        "heteroaromatic_ring": False,
        "carbonyl": False,
        "disulfide": False,
        "nitro": False,
        "azo": False,
        "nucleic_acid": False,
        "peptide_or_protein": False,
        "conjugated_alkene": False,
        "alkene": False,
        "alkyne": False,
    }

    # Heteroaromatic: n, o, s in lowercase inside rings, or keywords
    if any(c in smiles for c in ["n", "o", "s"] if c.islower()) or \
       any(kw in hints for kw in ["pyridine", "indole", "heteroaromatic", "imidazole", "thiazole"]):
        features["heteroaromatic_ring"] = True
        
    # Aromatic: c in lowercase inside rings, or keywords (if heteroaromatic is true, we might still set aromatic to true, but distinguish)
    if ("c" in smiles) or any(kw in hints for kw in ["benzene", "phenyl", "aromatic"]):
        features["aromatic_ring"] = True

    # Carbonyl
    has_carbonyl_base = "C(=O)" in smiles or "C=O" in smiles or "carbonyl" in hints
    has_sulfoxide = "S(=O)" in smiles
    has_phosphate = "P(=O)" in smiles
    if has_carbonyl_base and not has_sulfoxide and not has_phosphate:
        features["carbonyl"] = True

    # Disulfide / Sulfur
    if "SS" in smiles or "SS" in formula or re.search(r'S[2-9]', formula) or \
       any(kw in hints for kw in ["ss bond", "s-s bond", "disulfide", "sulfide"]):
        features["disulfide"] = True

    # Nitro
    if "N(=O)" in smiles or "[N+](=O)" in smiles or "NO2" in formula or "nitro" in hints:
        features["nitro"] = True

    # Azo
    if "N=N" in smiles or "azo" in hints:
        features["azo"] = True
        
    # Conjugated Alkene (heuristic: multiple double bonds close to each other, or hint)
    if "=C-C=" in smiles or "C=C-C=C" in smiles or "conjugat" in hints:
        features["conjugated_alkene"] = True
    elif "=" in smiles and not features["carbonyl"]:
        features["alkene"] = True
        
    # Alkyne
    if "#" in smiles or "alkyne" in hints:
        features["alkyne"] = True

    # Nucleic Acid
    if any(kw in hints for kw in ["dna", "rna", "nucleic", "oligonucleotide"]):
        features["nucleic_acid"] = True

    # Peptide / Protein
    if any(kw in hints for kw in ["peptide", "protein", "mab", "antibody"]):
        features["peptide_or_protein"] = True

    # 2. Map features to recommendations using CHROMOPHORE_RULES
    results = []
    has_strong = features["aromatic_ring"] or features["heteroaromatic_ring"] or features["conjugated_alkene"]

    # Evaluate rules
    if features["heteroaromatic_ring"]:
        rule = CHROMOPHORE_RULES["heteroaromatic_ring"]
        results.append(f"• <b>{rule['description']}</b> → {rule['recommendation']}")
    elif features["aromatic_ring"]: # Else-if to avoid redundant aromatic warnings if heteroaromatic is present
        rule = CHROMOPHORE_RULES["aromatic_ring"]
        results.append(f"• <b>{rule['description']}</b> → {rule['recommendation']}")

    if features["conjugated_alkene"]:
        rule = CHROMOPHORE_RULES["conjugated_alkene"]
        results.append(f"• <b>{rule['description']}</b> → {rule['recommendation']}")

    if features["carbonyl"]:
        rule = CHROMOPHORE_RULES["carbonyl"]
        results.append(f"• <b>{rule['description']}</b> → {rule['recommendation']}")

    if features["disulfide"]:
        rule = CHROMOPHORE_RULES["disulfide"]
        results.append(f"• <b>{rule['description']}</b> → {rule['recommendation']}")

    if features["nitro"]:
        rule = CHROMOPHORE_RULES["nitro"]
        if has_strong:
            results.append(f"• <b>{rule['description']}</b> → Useful bands at 254 nm and 330 nm. Monitor 330 nm if low UV background is poor.")
        else:
            results.append(f"• <b>{rule['description']}</b> → Strong at 254 nm. Recommend monitoring this channel.")

    if features["azo"]:
        rule = CHROMOPHORE_RULES["azo"]
        results.append(f"• <b>{rule['description']}</b> → {rule['recommendation']}")
        
    if features["nucleic_acid"]:
        rule = CHROMOPHORE_RULES["nucleic_acid"]
        results.append(f"• <b>{rule['description']}</b> → {rule['recommendation']}")
        
    if features["peptide_or_protein"]:
        rule = CHROMOPHORE_RULES["peptide_or_protein"]
        results.append(f"• <b>{rule['description']}</b> → {rule['recommendation']}")

    if not results:
        # Check deep UV options
        if features["alkene"]:
            rule = CHROMOPHORE_RULES["alkene"]
            results.append(f"• <b>{rule['description']}</b> → {rule['recommendation']}")
        elif features["alkyne"]:
            rule = CHROMOPHORE_RULES["alkyne"]
            results.append(f"• <b>{rule['description']}</b> → {rule['recommendation']}")
        else:
            results.append("• ⚠️ <b>No strong chromophore</b> → UV response is predicted to be weak; consider PDA scanning to confirm signal and be prepared to use MS, ELSD/CAD, or derivatization if sensitivity is inadequate.")

    html_output = "<br>".join(results)
    
    if any(str(wl) in html_output for wl in [190, 195, 200, 205, 210, 214, 215, 220]):
        html_output += "<br><br>⚠️ <i>Note: Recommended λ values are approximate detection windows. At deep-UV (< 220 nm), organic solvents and buffers may have significant absorbance and cause baseline noise/drift. Acetonitrile has a lower UV cut-off than methanol and is usually preferred.</i>"

    return html_output, features

def assess_additional_properties(smiles, formula, hints, logp, tpsa):
    """
    Evaluates ionization, solubility expectations, and chemical stability based on SMILES and hints.
    """
    smiles = smiles or ""
    formula = formula or ""
    hints = (hints or "").lower()
    
    try: logp_val = float(logp)
    except: logp_val = None
    try: tpsa_val = float(tpsa)
    except: tpsa_val = None

    flags = {
        "acidic": False,
        "basic": False,
        "amphoteric": False,
        "hydrolysis_risk": [],
        "oxidation_risk": []
    }

    # 1. Ionization & pKa awareness
    if any(m in smiles for m in ["C(=O)O", "S(=O)(=O)O", "P(=O)(O)O"]) or \
       any(kw in hints for kw in ["carboxylic acid", "benzoic acid", "acetic acid", "sulfonic acid", "phosphate", "phosphonate", "phenol", "tyrosine"]):
        flags["acidic"] = True
        
    # 'N' not adjacent to C=O (very basic heuristic), or explicit keywords
    if ("N" in smiles and "C(=O)N" not in smiles) or \
       any(kw in hints for kw in ["amine", "guanidinium", "biguanide", "piperidine", "piperazine", "pyrrolidine", "pyridine", "imidazole"]):
        flags["basic"] = True
        
    if any(kw in hints for kw in ["amino acid", "peptide", "betaine"]):
        flags["acidic"] = True
        flags["basic"] = True

    if flags["acidic"] and flags["basic"]:
        flags["amphoteric"] = True
        flags["ionization_warning"] = "Retention and peak shape will be highly pH-dependent due to amphoteric/zwitterionic ionizable groups. Predicted column choice and logP-based retention are less reliable without knowing actual pKa values. Experimental pKa scouting is required before finalizing mobile-phase pH."
    elif flags["acidic"]:
        flags["ionization_warning"] = "Retention and peak shape will be highly pH-dependent due to acidic ionizable groups. Predicted column choice and logP-based retention are less reliable without knowing actual pKa values. Experimental pKa scouting is required before finalizing mobile-phase pH."
    elif flags["basic"]:
        flags["ionization_warning"] = "Retention and peak shape will be highly pH-dependent due to basic ionizable groups. Predicted column choice and logP-based retention are less reliable without knowing actual pKa values. Experimental pKa scouting is required before finalizing mobile-phase pH."
    else:
        flags["ionization_warning"] = None

    # 2. Chemical Stability Flags
    if "C(=O)O" in smiles and "C(=O)OH" not in smiles:
        flags["hydrolysis_risk"].append("ester")
    if any(kw in hints for kw in ["ester", "lactam", "acetal", "ketal", "imine", "schiff"]):
        flags["hydrolysis_risk"].append("hydrolyzable group")
        
    if "C=O" in smiles and "C(=O)H" in smiles:
        flags["oxidation_risk"].append("aldehyde")
    if "c1ccc(O)cc1" in smiles or "phenol" in hints:
        flags["oxidation_risk"].append("phenol")
    if "S" in smiles or any(kw in hints for kw in ["thiol", "disulfide", "sulfide"]):
        flags["oxidation_risk"].append("thiol/sulfur")
        
    stability_warnings = []
    if flags["hydrolysis_risk"]:
        stability_warnings.append(f"Contains {', '.join(set(flags['hydrolysis_risk']))} functionality; verify solution stability over the planned pH range and autosampler time (hydrolysis can bias assay results).")
    if flags["oxidation_risk"]:
        stability_warnings.append(f"Contains {', '.join(set(flags['oxidation_risk']))} functionality; possible oxidation. Recommend minimizing oxygen/light exposure or adding antioxidants where appropriate.")

    flags["stability_warnings"] = stability_warnings

    # 3. Solubility & Diluent Compatibility
    if logp_val is not None:
        if logp_val < 0 and tpsa_val and tpsa_val > 60:
            flags["solubility_statement"] = "Analyte is very polar (logP < 0, high TPSA). Likely good water/buffer solubility; may have limited solubility in high organic solvents."
        elif logp_val > 4:
            flags["solubility_statement"] = "Analyte is very hydrophobic (logP > 4). Likely limited water solubility, better in organic solvents (ACN/MeOH)."
        else:
            flags["solubility_statement"] = "Analyte has moderate lipophilicity. Confirm solubility experimentally in both aqueous buffers and organics."
    else:
        flags["solubility_statement"] = "Solubility expectations cannot be inferred without logP/TPSA."

    flags["solubility_warning"] = (
        "Ensure the analyte is fully soluble in the chosen diluent and in the starting mobile phase; "
        "precipitation or phase separation will invalidate the method. "
        "For RP-HPLC, the injection solvent should not be significantly stronger than the starting mobile phase "
        "(e.g., injecting in 100% ACN onto a 5% organic starting condition can cause fronting/distorted peaks). "
        "Encourage matching or slightly weaker organic content in the diluent."
    )

    return flags

def recommend_column(logp, tpsa, mw, features, method_type, matrix):
    """
    Recommends a column chemistry based on physicochemical properties, method type, and normalized features.
    """
    # Safely cast inputs
    try: logp_val = float(logp)
    except: logp_val = None
    try: mw_val = float(mw)
    except: mw_val = None
    try: tpsa_val = float(tpsa)
    except: tpsa_val = None

    matrix = (matrix or "").lower()
    method_type = (method_type or "").lower()
    is_peptide = features.get("peptide_or_protein", False)
    has_aromatic = features.get("aromatic_ring", False) or features.get("heteroaromatic_ring", False)

    rule = None

    # 1. Peptides and large molecules
    if is_peptide or (mw_val and mw_val > 1500):
        rule = COLUMN_RULES["peptide_large"]

    # 2. Very polar analytes (HILIC) and check for k' < 1
    elif logp_val is not None and logp_val < 0:
        if tpsa_val and tpsa_val > 80:
            rule = COLUMN_RULES["very_polar"]
        else:
            rule = COLUMN_RULES["polar"]

    # 3. Very hydrophobic analytes
    elif logp_val is not None and logp_val > 4.0:
        rule = COLUMN_RULES["very_hydrophobic"]

    # 4. Secondary Check: Aromatic + Hydrophobic / Polar Aromatic
    elif logp_val is not None and has_aromatic:
        if logp_val > 3.5 and mw_val and mw_val > 400:
            rule = COLUMN_RULES["aromatic_bulky"]
        elif tpsa_val and tpsa_val > 60:
            rule = COLUMN_RULES["aromatic_polar"]

    # 5. Lipid-rich / Bioanalytical matrices
    if not rule and any(m in matrix for m in ["plasma", "serum", "oil", "lipid"]):
        rule = COLUMN_RULES["lipid_matrix"]
        
    # 6. Method Type Overrides
    if not rule and ("stability" in method_type or "impurity" in method_type):
        rule = COLUMN_RULES["stability_indicating"]

    # Default fallback
    if not rule:
        rule = COLUMN_RULES["default"]
        
    rec, expl = rule["recommendation"], rule["explanation"]

    return rec, expl + "<br><br><i>Note: This selection is a structure-based starting point; final column choice must be confirmed experimentally (retention, peak shape, resolution, robustness).</i>"

def recommend_sample_prep(matrix, features, method_type):
    """
    Recommends sample preparation steps based on knowledge base rules.
    """
    matrix = (matrix or "").lower()
    method_type = (method_type or "").lower()
    
    if any(m in matrix for m in ["plasma", "serum", "blood", "csf", "tissue"]):
        prep_info = dict(SAMPLE_PREP_RULES["plasma_serum_blood"])
        # Method-type adjustment
        if "bioanalytical" in method_type:
            prep_info["warnings"].append("For bioanalysis, strict matrix effect evaluation (IS normalization) is required.")
        if features.get("peptide_or_protein"):
            prep_info["warnings"].append("For peptides/proteins, avoid harsh PPT. Consider ultrafiltration or milder SPE.")
        return prep_info
        
    if "urine" in matrix:
        return dict(SAMPLE_PREP_RULES["urine"])
        
    if "cell culture" in matrix or "media" in matrix:
        return dict(SAMPLE_PREP_RULES["cell_culture"])

    if any(m in matrix for m in ["tablet", "capsule", "solid"]):
        prep_info = dict(SAMPLE_PREP_RULES["tablet_capsule"])
        if "impurity" in method_type:
            prep_info["warnings"].append("For impurity profiling, avoid aggressive extraction that might degrade the API artificially.")
        return prep_info

    if any(m in matrix for m in ["water", "environmental"]):
        return dict(SAMPLE_PREP_RULES["water_env"])

    if any(m in matrix for m in ["oil", "fat", "food", "lipid"]):
        return dict(SAMPLE_PREP_RULES["oil_lipid"])

    return dict(SAMPLE_PREP_RULES["default"])

def recommend_validation_checks(method_type):
    """
    Returns relevant validation and SST checks based on method type.
    """
    method_type = (method_type or "").lower()
    checks = [VALIDATION_SST_RULES["usp_621"]]
    
    if "stability" in method_type or "impurity" in method_type:
        checks.append(VALIDATION_SST_RULES["stability_indicating"])
    
    if "bioanalytical" in method_type:
        checks.append(VALIDATION_SST_RULES["bioanalytical"])
        
    checks.append(VALIDATION_SST_RULES["robustness"])
    return checks

def collect_warnings(matrix, method_type, features, column_type, logp, uv_html, chem_flags):
    """
    Aggregates warnings across all domains.
    """
    warnings = []
    
    # Sample Prep / Matrix Warnings
    matrix = (matrix or "").lower()
    if any(m in matrix for m in ["plasma", "serum", "tissue"]) and "Do not inject un-centrifuged" not in str(warnings):
        warnings.append("Do not inject un-centrifuged biological samples.")
        
    # Column/Polarity Warnings
    try: logp_val = float(logp)
    except: logp_val = None
    
    if logp_val is not None and logp_val < 0 and "C18" in column_type:
        warnings.append("Polar analyte on standard C18 – risk of k' < 1 and co-elution with void volume.")
        
    if any(m in matrix for m in ["oil", "lipid", "fat"]) and "C18" in column_type:
        warnings.append("Lipid matrix on C18 column – very high risk of irreversible fouling.")
        
    # LC-MS
    if "ms" in method_type.lower() or "bioanalytical" in method_type.lower() or "weak" in uv_html.lower():
        warnings.append("If using LC-MS detection: non-volatile buffers (phosphate, sulfate) are unsuitable; use volatile buffers (e.g. ammonium formate/acetate).")
        
    # Chemistry flags
    if chem_flags and chem_flags.get("ionization_warning"):
        warnings.append("⚛️ " + chem_flags["ionization_warning"])
    if chem_flags and chem_flags.get("stability_warnings"):
        for w in chem_flags["stability_warnings"]:
            warnings.append("🧪 " + w)
        
    return warnings

METHOD_OPTIONS = [
    {
        "label": "Potency / Assay",
        "description": "Measure the amount of API in a product",
        "impact": "Focuses on accuracy and linearity. Calibration range and standard curve guidance provided.",
        "key": "Potency Assay",
        "icon": "⚗️"
    },
    {
        "label": "Stability-Indicating",
        "description": "Separate and quantify API + all degradants",
        "impact": "Specificity is #1. Forced degradation study required. Column selectivity critical.",
        "key": "Stability-Indicating Assay",
        "icon": "🔬"
    },
    {
        "label": "Impurity Profiling",
        "description": "Detect and measure related substances",
        "impact": "High sensitivity required. LOQ guidance and low-UV detection will be flagged.",
        "key": "Impurity/Related Substances",
        "icon": "🔍"
    },
    {
        "label": "Dissolution",
        "description": "How fast the drug releases from the dosage form",
        "impact": "Dissolution medium IS the sample matrix. No protein precipitation needed.",
        "key": "Dissolution",
        "icon": "💊"
    },
    {
        "label": "Identification Test",
        "description": "Confirm presence — not measuring how much",
        "impact": "Retention time and UV spectrum match only. No quantitation needed.",
        "key": "Identification Test",
        "icon": "✅"
    },
    {
        "label": "Bioanalytical",
        "description": "Measure drug in plasma, urine, or tissue",
        "impact": "Matrix extraction is mandatory. LC-MS compatibility flagged. Internal standard required.",
        "key": "Bioanalytical",
        "icon": "🧬"
    },
    {
        "label": "Content Uniformity",
        "description": "Dose-to-dose consistency across tablets",
        "impact": "Individual tablet prep guidance (n=10). Precision is the primary concern.",
        "key": "Content Uniformity",
        "icon": "📊"
    }
]

METHOD_CONTEXT = {
    "Potency Assay": {
        "priority": "Accuracy & Linearity",
        "sensitivity_note": "Typical working range: 50–150% of nominal concentration.",
        "key_warning": "Ensure the reference standard is accurately weighed and the calibration range brackets all expected sample concentrations.",
        "lc_ms_needed": False,
        "extraction_required": False,
        "forced_degradation_required": False,
    },
    "Stability-Indicating Assay": {
        "priority": "Specificity — separation from all degradants",
        "sensitivity_note": "Must resolve API from all forced degradation products (acid, base, oxidative, thermal, photolytic).",
        "key_warning": "Forced degradation study required before method scouting. Column selectivity is the #1 parameter to optimize.",
        "lc_ms_needed": False,
        "extraction_required": False,
        "forced_degradation_required": True,
    },
    "Impurity/Related Substances": {
        "priority": "Sensitivity — LOQ typically 0.05–0.1% of API",
        "sensitivity_note": "Detection limits must cover ICH Q3A/Q3B thresholds. Low UV or MS may be needed for weak chromophores.",
        "key_warning": "If analyte has weak UV, impurity profiling at trace levels may require LC-MS, ELSD/CAD, or derivatization.",
        "lc_ms_needed": False,
        "extraction_required": False,
        "forced_degradation_required": False,
    },
    "Dissolution": {
        "priority": "Sample simplicity — dissolution medium IS the matrix",
        "sensitivity_note": "Medium pH directly affects mobile phase design. No protein precipitation needed.",
        "key_warning": "Confirm analyte is UV-active at expected dissolution concentrations. Filter (0.45 μm) the medium directly; no other prep usually required.",
        "lc_ms_needed": False,
        "extraction_required": False,
        "forced_degradation_required": False,
    },
    "Identification Test": {
        "priority": "Retention time match + UV spectral confirmation",
        "sensitivity_note": "No calibration curve needed. Method confirms presence/absence only.",
        "key_warning": "PDA is strongly recommended so UV spectrum can be compared to a reference standard in the same run.",
        "lc_ms_needed": False,
        "extraction_required": False,
        "forced_degradation_required": False,
    },
    "Bioanalytical": {
        "priority": "Sensitivity + selectivity in biological matrix",
        "sensitivity_note": "Typically pg/mL to ng/mL range; LC-MS/MS is often the required detector.",
        "key_warning": "Matrix extraction (PPT/SPE/LLE) is mandatory. An internal standard is required for quantitation. Validate per FDA/EMA bioanalytical guidance.",
        "lc_ms_needed": True,
        "extraction_required": True,
        "forced_degradation_required": False,
    },
    "Content Uniformity": {
        "priority": "Precision across individual dosage units",
        "sensitivity_note": "USP <905> or Ph.Eur. 2.9.40 applies. Typically n=10 individual tablets.",
        "key_warning": "Prepare each tablet individually (no pooling). Extraction efficiency must be >98% and consistent unit-to-unit.",
        "lc_ms_needed": False,
        "extraction_required": False,
        "forced_degradation_required": False,
    }
}

# STAGE 0
if st.session_state.stage == 0:
    st.subheader("Stage 0 — Instrument Qualification")
    st.write("Has your HPLC instrument been formally qualified?")
    st.caption("📖 [Read the official FDA Guidance on Analytical Procedures and Methods Validation](https://www.fda.gov/media/87801/download) (outlines equipment qualification requirements)")

    q = st.radio(
        "Select your situation:",
        [
            "Yes — all four phases complete (DQ, IQ, OQ, PQ)",
            "Partially — some phases done, some not",
            "No — instrument has not been qualified",
            "I don't know what instrument qualification is"
        ],
        key="q0"
    )

    if st.button("Continue →"):
        if q == "No — instrument has not been qualified":
            st.error("❌ Qualification required. No regulatory body will accept data from an unqualified instrument.")
        elif q == "I don't know what instrument qualification is":
            st.info("📚 Instrument qualification is a four-step process (DQ, IQ, OQ, PQ) that proves your HPLC works correctly. Complete this with your instrument vendor before continuing.")
        else:
            st.session_state.qualification = q
            go_next()
            st.rerun()



# STAGE 1
elif st.session_state.stage == 1:
    st.subheader("Stage 1 — Method Type")
    st.write("What are you trying to do with this method?")

    for opt in METHOD_OPTIONS:
        is_selected = st.session_state.get("method_type") == opt["key"]
        border = "2px solid #0A1931" if is_selected else "1px solid #B3CFE5"
        if st.button(f"{opt['icon']} **{opt['label']}**", key=f"btn_{opt['key']}"):
            st.session_state.method_type = opt["key"]
            go_next()
            st.rerun()
        st.markdown(f"""
        <div class="method-card" style="border:{border};">
            <small>{opt['description']}</small><br>
            <small>→ {opt['impact']}</small>
        </div>
        """, unsafe_allow_html=True)

    st.divider()
    if st.button("← Back"):
        go_back()
        st.rerun()

# STAGE 2
elif st.session_state.stage == 2:
    st.subheader("Stage 2 — Analyte Properties")
    st.write("Enter your compound name and we will fetch its molecular properties from PubChem.")

    compound_name = st.text_input(
        "Compound name:",
        placeholder="e.g. ibuprofen, dimethyl trisulfide, metformin"
    )

    extra_hints = st.text_input(
        "Compound class or hints (optional):",
        placeholder="e.g. peptide, DNA, impurity, high-fat matrix..."
    )

    if st.button("🔍 Fetch Properties"):
        if not compound_name.strip():
            st.warning("Please enter a compound name.")
        else:
            with st.spinner(f"Searching PubChem for {compound_name}..."):
                props = fetch_compound(compound_name.strip())
            if props:
                st.session_state.compound = compound_name.strip()
                st.session_state.properties = props
            else:
                st.error("❌ Not found in PubChem. Check spelling or try the chemical name.")

    if st.session_state.properties:
        p = st.session_state.properties
        smiles = p.get("smiles", "")

        st.success(f"✅ Found: {p['iupac']}")
        
        ctx = METHOD_CONTEXT.get(st.session_state.method_type, {})
        if ctx:
            st.markdown(f"""
            <div class="stage-card">
                <b>🎯 Method Priority:</b> {ctx['priority']}<br><br>
                <b>📏 Sensitivity context:</b> {ctx['sensitivity_note']}<br><br>
                <b>⚠️ Key requirement:</b> {ctx['key_warning']}
            </div>
            """, unsafe_allow_html=True)

            if ctx.get("lc_ms_needed"):
                st.warning("🔬 Bioanalytical methods typically require LC-MS/MS. Ensure mobile phase uses volatile buffers only (ammonium formate/acetate).")
            if ctx.get("extraction_required"):
                st.warning("🧪 Biological matrix detected: protein precipitation, centrifugation, and filtration are mandatory before column injection.")
            if ctx.get("forced_degradation_required"):
                st.warning("⚗️ Stability-indicating method: forced degradation study must precede method scouting to identify potential degradants.")
                
        st.markdown("### Molecular Properties")

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Molecular Weight", f"{p['mw']} g/mol")
        with col2:
            st.metric("LogP", p['logp'])
        with col3:
            st.metric("TPSA", f"{p['tpsa']} Ų")

        col4, col5 = st.columns(2)
        with col4:
            st.metric("H-bond Donors", p['hbd'])
        with col5:
            st.metric("H-bond Acceptors", p['hba'])

        st.warning("⚠️ LogP shown is calculated (XLogP3). For ionizable compounds and salt forms, experimental LogP may differ significantly. Verify against literature before finalizing mobile phase organic content.")
        
        iupac_lower = p.get('iupac', '').lower()
        salt_keywords = ["hydrochloride", "mesylate", "sulfate", "sodium", "potassium", "maleate", "tartrate"]
        if any(kw in iupac_lower for kw in salt_keywords):
            st.warning("⚠️ Salt form detected. LogP for the salt form differs from the free base/acid. Method development should use free base/acid LogP where possible.")

        st.markdown("### pKa (Required for ionizable compounds)")
        st.write("PubChem does not reliably provide pKa. Enter it manually if known.")

        pka_known = st.radio("Do you know the pKa of your compound?", 
            ["No — compound is neutral (no ionizable groups)",
             "Yes — acidic compound (carboxylic acid, phenol, sulfonamide)",
             "Yes — basic compound (amine, imidazole, pyridine)",
             "Yes — amphoteric (has both acidic and basic groups)"])

        pka_data = {"type": "neutral", "acidic": None, "basic": None, "recommendation": ""}

        if pka_known == "Yes — acidic compound (carboxylic acid, phenol, sulfonamide)":
            pka_data["type"] = "acidic"
            pka_data["acidic"] = st.number_input("Acidic pKa", value=4.0, step=0.1)
            pka_data["recommendation"] = f"Run mobile phase at pH {pka_data['acidic']-2:.1f} to {pka_data['acidic']-3:.1f}. This ensures fully protonated (neutral) form for reproducible retention."
            st.info(pka_data["recommendation"])
        elif pka_known == "Yes — basic compound (amine, imidazole, pyridine)":
            pka_data["type"] = "basic"
            pka_data["basic"] = st.number_input("Basic pKa", value=9.0, step=0.1)
            pka_data["recommendation"] = f"Run mobile phase at pH {pka_data['basic']-2:.1f} or lower. This fully protonates the amine and suppresses silanol interactions."
            st.info(pka_data["recommendation"])
        elif pka_known == "Yes — amphoteric (has both acidic and basic groups)":
            pka_data["type"] = "amphoteric"
            c1, c2 = st.columns(2)
            with c1:
                pka_data["acidic"] = st.number_input("Acidic pKa", value=4.0, step=0.1)
            with c2:
                pka_data["basic"] = st.number_input("Basic pKa", value=9.0, step=0.1)
            pka_data["recommendation"] = "Run between the two pKa values or test both pH extremes. This compound class requires experimental pH scouting."
            st.info(pka_data["recommendation"])
        else:
            pka_data["recommendation"] = "pH has minimal effect on retention. Use 0.1% formic acid as default."
            st.info(pka_data["recommendation"])

        st.caption("Don't know your pKa? Calculate it free at chemicalize.com using your SMILES string.")
        st.session_state.pka_data = pka_data

        st.markdown("### What this means for your method")

        uv_result, features = interpret_uv(smiles, p.get('formula', ''), extra_hints)
        chem_flags = assess_additional_properties(smiles, p.get('formula', ''), extra_hints, p.get('logp'), p.get('tpsa'))
        
        st.session_state.features = features  # Save for downstream stages
        
        col_headline, col_explanation = recommend_column(
            logp=p.get('logp'),
            tpsa=p.get('tpsa'),
            mw=p.get('mw'),
            features=features,
            method_type=st.session_state.method_type,
            matrix=extra_hints # Proxy for matrix until stage 3
        )

        warnings = collect_warnings(
            matrix=extra_hints,
            method_type=st.session_state.method_type,
            features=features,
            column_type=col_headline,
            logp=p.get('logp'),
            uv_html=uv_result,
            chem_flags=chem_flags
        )

        validation_checks = recommend_validation_checks(st.session_state.method_type)

        st.markdown(f"""
        <div class="prop-card">
            <b>🧪 Column recommendation:</b><br>
            <b>{col_headline}</b><br>
            {col_explanation}
        </div>
        <div class="prop-card">
            <b>💡 Detection wavelength:</b><br>
            {uv_result}
        </div>
        <div class="prop-card">
            <b>🔬 Formula:</b> {p['formula']}<br>
            <b>📋 IUPAC Name:</b> {p['iupac']}<br>
            <b>🧬 SMILES:</b> {smiles}
        </div>
        """, unsafe_allow_html=True)

        if warnings:
            st.markdown("### ⚠️ Method Warnings")
            for w in warnings:
                st.warning(w)

        st.markdown("### 📋 Validation & System Suitability Checks")
        for check in validation_checks:
            with st.expander(check["description"]):
                for item in check["items"]:
                    st.write(f"- {item}")

        st.divider()
        st.info("ℹ️ **PubChem Identity Caution:** These recommendations assume the PubChem match (IUPAC name, formula, SMILES) corresponds to your exact analyte and form (isomer / salt / hydrate). Please confirm before using the guidance.")
        st.info(f"💧 **Solubility Caution:** {chem_flags['solubility_statement']} {chem_flags['solubility_warning']}")

    st.divider()
    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("← Back"):
            go_back()
            st.rerun()
    with col2:
        if st.button("Continue →"):
            if not st.session_state.properties:
                st.warning("Please fetch compound properties before continuing.")
            else:
                go_next()
                st.rerun()

# STAGE 3
elif st.session_state.stage == 3:
    st.subheader("Stage 3 — Matrix & Sample Preparation")
    st.write("What kind of sample matrix will you be analyzing?")

    matrix_options = [
        "Plasma / Serum / Blood",
        "Urine",
        "Cell Culture Media",
        "Tablet / Capsule / Solid Dosage",
        "Water / Environmental",
        "Oil / Lipid / Fat",
        "Other / Standard Solution"
    ]
    
    # Try to pre-select based on existing hints, otherwise default to 'Other'
    default_idx = len(matrix_options) - 1
    
    selected_matrix = st.selectbox("Select Matrix:", matrix_options, index=default_idx)

    st.markdown("### Recommended Sample Preparation")
    
    features = st.session_state.get("features", {})
    method_type = st.session_state.get("method_type", "")
    
    prep_info = recommend_sample_prep(selected_matrix, features, method_type)
    
    st.info(f"**Overview:** {prep_info['summary']}")
    
    st.markdown("#### Steps:")
    for i, step in enumerate(prep_info['steps'], 1):
        st.markdown(f"{i}. {step}")
        
    if prep_info['warnings']:
        st.markdown("#### ⚠️ Matrix Warnings:")
        for w in prep_info['warnings']:
            st.warning(w)

    st.divider()
    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("← Back"):
            go_back()
            st.rerun()
    with col2:
        if st.button("Continue →"):
            st.session_state.matrix = selected_matrix
            go_next()
            st.rerun()

# STAGE 4 - METHOD SCOUTING
elif st.session_state.stage == 4:
    import anthropic
    import os
    from dotenv import load_dotenv
    load_dotenv()

    st.subheader("Stage 4 — Method Scouting")
    st.write("Based on your compound properties, we will now generate a complete scouting plan.")

    p = st.session_state.properties or {}
    method_type = st.session_state.method_type or "Assay"
    matrix = st.session_state.get("matrix", "Not specified")

    st.markdown(f"""
    <div class="stage-card">
        <b>Compound:</b> {st.session_state.compound}<br>
        <b>Method type:</b> {method_type}<br>
        <b>Matrix:</b> {matrix}<br>
        <b>LogP:</b> {p.get("logp", "N/A")} &nbsp;|&nbsp;
        <b>MW:</b> {p.get("mw", "N/A")} g/mol &nbsp;|&nbsp;
        <b>TPSA:</b> {p.get("tpsa", "N/A")} Ų
    </div>
    """, unsafe_allow_html=True)

    st.divider()

    st.markdown("### Before we scout — a few quick questions")

    num_analytes = st.radio(
        "How many analytes does your method need to separate?",
        [
            "1 — just the API, no impurities needed",
            "2–5 — API plus a few known impurities",
            "6–15 — API plus many impurities or degradants",
            "More than 15 — complex mixture"
        ],
        key="num_analytes"
    )

    instruments = st.multiselect(
        "Which instruments do you have available?",
        ["HPLC-UV", "HPLC-DAD/PDA", "HPLC-Fluorescence", "HPLC-ELSD"],
        default=["HPLC-UV"],
        key="instruments"
    )

    prior_knowledge = st.text_area(
        "Any prior knowledge about this compound? (optional)",
        placeholder="e.g. known to tail on C18, elutes around 5 min, sensitive to pH changes, previously used ACN/water gradient...",
        key="prior_knowledge"
    )

    forced_deg = "Not required"
    if method_type == "Stability-Indicating Assay":
        forced_deg = st.radio(
            "Have forced degradation samples been prepared?",
            [
                "Yes — I have stressed samples ready",
                "No — generate the protocol for me"
            ],
            key="forced_deg"
        )

    st.divider()

    if st.button("🧠 Generate Scouting Plan"):

        api_key = os.getenv("ANTHROPIC_API_KEY")
        if not api_key:
            st.error("API key not found. Check your .env file.")
        else:
            SCOUTING_SYSTEM_PROMPT = """You are a senior analytical chemist with 20+ years of HPLC method development experience in pharmaceutical laboratories. You specialize in ICH Q2(R2) compliant methods.

Your job is to generate a practical, immediately actionable HPLC method scouting plan based on the compound profile and study parameters provided.

DECISION RULES YOU APPLY:

COLUMN SELECTION (based on LogP):
- LogP < 0: HILIC primary. Compound too polar for standard RP.
- LogP 0–1: C8 or mixed-mode. Limited RP retention expected.
- LogP 1–4: C18 primary (e.g. Waters XBridge C18, Phenomenex Luna C18, Agilent Zorbax Eclipse Plus C18). This is the ideal RP range.
- LogP 4–5: C18 with high organic start, or phenyl-hexyl for aromatic compounds.
- LogP > 5: C4/C8 preferred. C18 risks excessive retention.

MOBILE PHASE pH (based on ionization):
- Acidic compound (pKa 3–5): Run at pH 2.0–3.0. Use 0.1% formic acid or 10 mM ammonium formate pH 3.
- Basic compound (pKa 8–10): Run at pH 2.0–3.0 to fully protonate. Low pH suppresses silanol activity.
- Neutral compound: pH less critical. Use 0.1% formic acid as default.
- RULE: Always run at least 2 pH units from the pKa.

ISOCRATIC vs GRADIENT:
- Single analyte, LogP 1–3: Try isocratic first at estimated organic %
- Multiple analytes or wide polarity range: Use gradient
- Complex mixtures: Gradient from 5% organic, ramp to 95% over 20 minutes

FORCED DEGRADATION (stability-indicating methods only):
- Acid hydrolysis: 0.1M HCl, 60°C, 1 hour
- Base hydrolysis: 0.1M NaOH, 60°C, 1 hour
- Oxidation: 3% H2O2, room temperature, 24 hours
- Photolysis: ICH Q1B light exposure, 24 hours
- Thermal: 60°C dry heat, 1 week
- Target: 10–30% degradation of API

YOUR OUTPUT FORMAT:
Provide a complete scouting plan with these exact sections:

1. SCOUTING STRATEGY
Brief explanation of your approach based on the compound properties.

2. COLUMN CANDIDATES
List 2–3 columns to scout in priority order. For each: column name, dimensions, particle size, reason for selection.

3. MOBILE PHASE CONDITIONS
For each column, specify 2 mobile phase conditions to test. Include exact compositions.

4. STARTING GRADIENT
Exact gradient table (time vs %B) for scouting runs.

5. INSTRUMENT SETTINGS
Flow rate, column temperature, injection volume, detection wavelength.

6. SCOUTING EXPERIMENT MATRIX
A clear table showing which experiments to run: Column x Mobile Phase combinations.

7. WHAT TO LOOK FOR
How to interpret the scouting chromatograms. What constitutes a good result vs a poor result.

8. IF STABILITY-INDICATING
Forced degradation protocol and how to use stressed samples in scouting.

9. NEXT STEPS
What to do with the scouting results before proceeding to optimization.

Be specific. Give real numbers. A scientist should be able to walk to their instrument and start immediately."""

            user_message = f"""Generate a complete HPLC method scouting plan for this compound:

COMPOUND PROFILE:
- Name: {st.session_state.compound}
- IUPAC: {p.get("iupac", "N/A")}
- Formula: {p.get("formula", "N/A")}
- Molecular Weight: {p.get("mw", "N/A")} g/mol
- LogP: {p.get("logp", "N/A")}
- TPSA: {p.get("tpsa", "N/A")} Angstrom squared
- H-bond Donors: {p.get("hbd", "N/A")}
- H-bond Acceptors: {p.get("hba", "N/A")}
- SMILES: {p.get("smiles", "N/A")}
- pKa Data: {st.session_state.get('pka_data', {}).get('type', 'N/A')} (Acidic: {st.session_state.get('pka_data', {}).get('acidic', 'N/A')}, Basic: {st.session_state.get('pka_data', {}).get('basic', 'N/A')})

STUDY PARAMETERS:
- Method type: {method_type}
- Sample matrix: {matrix}
- Number of analytes to separate: {num_analytes}
- Available instruments: {", ".join(instruments)}
- Forced degradation status: {forced_deg}

PRIOR KNOWLEDGE:
{prior_knowledge if prior_knowledge else "None provided."}

Generate the complete scouting plan now."""

            with st.spinner("We are generating your scouting plan... (this takes 15–30 seconds)"):
                try:
                    client = anthropic.Anthropic(api_key=api_key)
                    message = client.messages.create(
                        model="claude-sonnet-4-6",
                        max_tokens=9000,
                        system=SCOUTING_SYSTEM_PROMPT,
                        messages=[{"role": "user", "content": user_message}]
                    )
                    scouting_plan = message.content[0].text
                    st.session_state.scouting_plan = scouting_plan
                except Exception as e:
                    st.error(f"Claude API error: {str(e)}")

    if st.session_state.get("scouting_plan"):
        st.divider()
        st.markdown("### Your Scouting Plan")
        st.markdown(st.session_state.scouting_plan)

        st.divider()
        st.download_button(
            label="📄 Download Scouting Plan",
            data=st.session_state.scouting_plan,
            file_name=f"scouting_plan_{st.session_state.compound}.txt",
            mime="text/plain"
        )

        st.divider()
        st.markdown("### Record Your Selected Method Conditions")
        st.write("After running your scouting experiments, enter the conditions you selected. These will carry forward to robustness testing.")

        with st.expander("Enter final selected conditions (required before continuing)"):
            selected_column = st.text_input("Selected column (name, dimensions, particle size)", 
                placeholder="e.g. Waters XBridge C18, 150x4.6mm, 3.5μm", key="sel_column")
            selected_mpa = st.text_input("Mobile Phase A", 
                placeholder="e.g. 0.1% Formic acid in water", key="sel_mpa")
            selected_mpb = st.text_input("Mobile Phase B", 
                placeholder="e.g. Acetonitrile", key="sel_mpb")
            selected_gradient = st.text_area("Gradient or isocratic conditions", 
                placeholder="e.g. 5% B to 60% B over 15 min, hold 3 min, re-equilibrate 5 min", key="sel_gradient")
            selected_flow = st.number_input("Flow rate (mL/min)", 0.1, 5.0, 1.0, 0.1, key="sel_flow")
            selected_temp = st.number_input("Column temperature (°C)", 15.0, 60.0, 30.0, 1.0, key="sel_temp")
            selected_wavelength = st.number_input("Detection wavelength (nm)", 190.0, 400.0, 254.0, 1.0, key="sel_wavelength")
            observed_rt = st.number_input("Observed analyte retention time (min)", 0.0, 60.0, 0.0, 0.1, key="sel_rt")
            observed_tailing = st.number_input("Observed tailing factor", 0.0, 10.0, 0.0, 0.01, key="sel_tailing")
            observed_plates = st.number_input("Observed theoretical plates (N)", 0, 100000, 0, 100, key="sel_plates")
            
            st.session_state.final_method = {
                "column": selected_column,
                "mpa": selected_mpa,
                "mpb": selected_mpb,
                "gradient": selected_gradient,
                "flow": selected_flow,
                "temp": selected_temp,
                "wavelength": selected_wavelength,
                "rt": observed_rt,
                "tailing": observed_tailing,
                "plates": observed_plates
            }

    st.divider()
    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("← Back"):
            go_back()
            st.rerun()
    with col2:
        if st.button("Continue →"):
            if not st.session_state.get("scouting_plan"):
                st.warning("Please generate a scouting plan before continuing.")
            elif not st.session_state.get("final_method", {}).get("column"):
                st.warning("Please record your selected method conditions before continuing.")
            else:
                go_next()
                st.rerun()


# STAGE 5 - ROBUSTNESS TESTING
elif st.session_state.stage == 5:
    import anthropic
    import os
    from dotenv import load_dotenv
    load_dotenv()

    st.subheader("Stage 5 — Robustness Testing")
    st.write("Test how sensitive your method is to small deliberate changes in conditions.")
    st.info("Vary one parameter at a time. Run 3 replicates at each condition. Enter the results below.")

    p = st.session_state.properties or {}
    fm = st.session_state.get("final_method", {})

    if fm:
        st.markdown(f"""
        <div class="method-card">
            <b>Selected Column:</b> {fm.get('column', 'N/A')}<br>
            <b>Mobile Phase A:</b> {fm.get('mpa', 'N/A')}<br>
            <b>Mobile Phase B:</b> {fm.get('mpb', 'N/A')}<br>
            <b>Gradient:</b> {fm.get('gradient', 'N/A')}
        </div>
        """, unsafe_allow_html=True)

    st.markdown("### Your Final Method Conditions")
    col1, col2 = st.columns(2)
    with col1:
        nominal_ph = st.number_input("Nominal mobile phase pH", min_value=1.0, max_value=13.0, value=3.0, step=0.1)
        nominal_organic = st.number_input("Nominal organic % (e.g. 40 for 40% ACN)", min_value=0.0, max_value=100.0, value=40.0, step=1.0)
        nominal_flow = st.number_input("Nominal flow rate (mL/min)", min_value=0.1, max_value=5.0, value=float(fm.get("flow", 1.0)), step=0.1)
    with col2:
        nominal_temp = st.number_input("Nominal column temperature (°C)", min_value=15.0, max_value=60.0, value=float(fm.get("temp", 30.0)), step=1.0)
        nominal_wavelength = st.number_input("Nominal detection wavelength (nm)", min_value=190.0, max_value=400.0, value=float(fm.get("wavelength", 210.0)), step=1.0)
        nominal_assay = st.number_input("Nominal assay result at standard conditions (%)", min_value=0.0, max_value=200.0, value=100.0, step=0.1)

    st.divider()
    st.divider()
    st.markdown("### Robustness Results")
    st.write("Enter the 3 individual replicate assay results (%) at each varied condition.")
    st.write("Leave blank (0.0) if you have not tested that condition yet.")

    parameters = [
        ("pH minus", f"pH {round(nominal_ph - 0.2, 1)} (nominal - 0.2)"),
        ("pH plus", f"pH {round(nominal_ph + 0.2, 1)} (nominal + 0.2)"),
        ("Organic minus", f"{round(nominal_organic - 2.0, 1)}% organic (nominal - 2%)"),
        ("Organic plus", f"{round(nominal_organic + 2.0, 1)}% organic (nominal + 2%)"),
        ("Flow minus", f"{round(nominal_flow - 0.1, 1)} mL/min (nominal - 10%)"),
        ("Flow plus", f"{round(nominal_flow + 0.1, 1)} mL/min (nominal + 10%)"),
        ("Temp minus", f"{round(nominal_temp - 5.0, 0):.0f}°C (nominal - 5°C)"),
        ("Temp plus", f"{round(nominal_temp + 5.0, 0):.0f}°C (nominal + 5°C)"),
        ("Wavelength minus", f"{round(nominal_wavelength - 2.0, 0):.0f} nm (nominal - 2 nm)"),
        ("Wavelength plus", f"{round(nominal_wavelength + 2.0, 0):.0f} nm (nominal + 2 nm)"),
        ("Column lot B", "Different column lot (same spec)"),
    ]

    results = {}
    for key, label in parameters:
        with st.expander(f"📊 {label}"):
            c1, c2, c3 = st.columns(3)
            with c1:
                rep1 = st.number_input("Replicate 1 (%)", min_value=0.0, max_value=200.0, value=0.0, step=0.1, key=f"r1_{key}")
            with c2:
                rep2 = st.number_input("Replicate 2 (%)", min_value=0.0, max_value=200.0, value=0.0, step=0.1, key=f"r2_{key}")
            with c3:
                rep3 = st.number_input("Replicate 3 (%)", min_value=0.0, max_value=200.0, value=0.0, step=0.1, key=f"r3_{key}")
            
            if rep1 > 0 and rep2 > 0 and rep3 > 0:
                vals = [rep1, rep2, rep3]
                mean_val = np.mean(vals)
                rsd_val = (np.std(vals, ddof=1) / mean_val * 100) if mean_val > 0 else 0.0
                diff = abs(mean_val - nominal_assay)
                
                st.caption(f"**Mean:** {mean_val:.2f}% | **%RSD:** {rsd_val:.2f}% | **Diff from Nominal:** {diff:.2f}%")
                results[key] = {"label": label, "mean": round(mean_val, 2), "rsd": round(rsd_val, 2)}
            elif rep1 > 0 or rep2 > 0 or rep3 > 0:
                st.warning("⚠️ Minimum 3 replicates required per ICH Q2(R2)")

    st.divider()

    if st.button("🧠 Assess Robustness"):
        if not results:
            st.warning("Please enter at least one robustness result before assessing.")
        else:
            api_key = os.getenv("ANTHROPIC_API_KEY")

            ROBUSTNESS_PROMPT = """You are a senior pharmaceutical analytical chemist assessing HPLC method robustness per ICH Q2(R2) guidelines.

ACCEPTANCE CRITERIA:
- %RSD at each changed condition must be ≤ 2.0%
- Difference between nominal assay and mean assay at changed condition must be ≤ 2.0%
- A parameter is CRITICAL if changing it causes either criterion to fail
- A parameter is NON-CRITICAL if both criteria pass

YOUR OUTPUT:
1. PARAMETER-BY-PARAMETER ASSESSMENT
For each tested parameter: PASS or FAIL, the values, and why.

2. CRITICAL PARAMETERS
List any parameters that failed. These must be tightly controlled in the method specification.

3. NON-CRITICAL PARAMETERS
List parameters that passed. Normal lab variation is acceptable for these.

4. METHOD OPERABLE DESIGN REGION (MODR)
State the proven acceptable ranges for each non-critical parameter.

5. OVERALL ROBUSTNESS CONCLUSION
Is the method robust? What must be added to the method specification?

Be specific. Cite exact values."""

            results_text = f"Nominal assay result: {nominal_assay}%\n\n"
            for key, data in results.items():
                diff = abs(data["mean"] - nominal_assay)
                results_text += f"- {data['label']}: Mean = {data['mean']}%, %RSD = {data['rsd']}%, Difference from nominal = {diff:.1f}%\n"

            user_msg = f"""Assess the robustness of this HPLC method:

COMPOUND: {st.session_state.compound}
METHOD TYPE: {st.session_state.method_type}

NOMINAL CONDITIONS:
- pH: {nominal_ph}
- Organic %: {nominal_organic}%
- Flow rate: {nominal_flow} mL/min
- Column temperature: {nominal_temp}°C
- Detection wavelength: {nominal_wavelength} nm

ROBUSTNESS RESULTS:
{results_text}

Assess each parameter and provide the complete robustness report."""

            with st.spinner("We are assessing your robustness data..."):
                try:
                    client = anthropic.Anthropic(api_key=api_key)
                    message = client.messages.create(
                        model="claude-sonnet-4-6",
                        max_tokens=4096,
                        system=ROBUSTNESS_PROMPT,
                        messages=[{"role": "user", "content": user_msg}]
                    )
                    robustness_report = message.content[0].text
                    st.session_state.robustness_report = robustness_report
                except Exception as e:
                    st.error(f"Claude API error: {str(e)}")

    if st.session_state.get("robustness_report"):
        st.divider()
        st.markdown("### Robustness Assessment")
        st.markdown(st.session_state.robustness_report)
        st.download_button(
            label="📄 Download Robustness Report",
            data=st.session_state.robustness_report,
            file_name=f"robustness_{st.session_state.compound}.txt",
            mime="text/plain"
        )

    st.divider()
    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("← Back"):
            go_back()
            st.rerun()
    with col2:
        if st.button("Continue →"):
            go_next()
            st.rerun()


# STAGE 6 - PREREQUISITES
elif st.session_state.stage == 6:

    st.subheader("Stage 6 — Prerequisites Confirmation")
    st.write("Before formal validation begins, confirm all prerequisites are in place.")
    st.warning("All items must be confirmed before validation experiments start. Missing prerequisites will be flagged in any regulatory audit.")

    st.markdown("### Reagents & Standards")
    p1 = st.checkbox("All reagents and solvents are available and within expiry date")
    p2 = st.checkbox("Reference standard has a valid certificate of analysis (CoA)")
    p3 = st.checkbox("Reference standard potency/purity is assigned and documented")
    p4 = st.checkbox("Working standard has been prepared and concentration assigned")

    st.markdown("### Instrument & Column")
    p5 = st.checkbox("HPLC column is available (new column recommended for validation)")
    p6 = st.checkbox("Instrument has a current calibration certificate")
    p7 = st.checkbox("Calibration is not overdue")

    st.markdown("### Personnel & Documentation")
    p8 = st.checkbox("Analyst has documented training on HPLC operation")
    p9 = st.checkbox("Analyst has documented training on this method type")
    p10 = st.checkbox("Method SOP draft has been written and is available")

    all_checked = all([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10])

    st.divider()

    if all_checked:
        st.success("✅ All prerequisites confirmed. You are ready to begin formal validation.")
    else:
        missing = sum([not p1, not p2, not p3, not p4, not p5, not p6, not p7, not p8, not p9, not p10])
        st.error(f"❌ {missing} prerequisite(s) not yet confirmed. Resolve these before starting validation.")

    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("← Back"):
            go_back()
            st.rerun()
    with col2:
        if st.button("Continue →"):
            if not all_checked:
                st.error("All prerequisites must be confirmed before proceeding to validation.")
            else:
                go_next()
                st.rerun()

# STAGE 7 - FORMAL VALIDATION
elif st.session_state.stage == 7:
    import anthropic
    import os
    from dotenv import load_dotenv
    load_dotenv()

    st.subheader("Stage 7 — Formal Method Validation")
    st.write("Enter your validation data. We will assess each parameter against ICH Q2(R2) acceptance criteria.")

    p = st.session_state.properties or {}

    st.markdown("### Specificity / Selectivity")
    st.write("ICH Q2(R2) requires this as the first validation parameter.")

    blank_result = st.radio(
        "Blank injection — is there any peak at the analyte retention time?",
        ["No peak detected — blank is clean ✅",
         "Peak detected — interference present ❌",
         "Not yet tested"])

    placebo_result = st.radio(
        "Placebo injection (formulation without API) — any co-eluting peak with analyte?",
        ["No interference detected ✅",
         "Interference detected ❌", 
         "Not applicable — pure API, no excipients",
         "Not yet tested"])

    peak_purity = st.number_input(
        "Peak purity index (PDA/DAD only — leave 0 if UV detector only)",
        min_value=0.0, max_value=1.0, value=0.0, step=0.001, format="%.3f",
        key="peak_purity")

    degradant_resolution = "Not evaluated"
    if "Stability" in st.session_state.method_type:
        degradant_resolution = st.radio(
            "Are all forced degradation products baseline resolved from the API peak?",
            ["Yes — all degradants resolved (Rs ≥ 2.0) ✅",
             "No — one or more degradants co-elute with API ❌",
             "Not yet tested"])

    st.markdown("### System Suitability (SST)")
    st.write("Enter results from 6 replicate injections of working standard.")
    col1, col2, col3 = st.columns(3)
    with col1:
        sst_rsd = st.number_input("Peak area %RSD (6 injections)", 0.0, 100.0, 0.0, 0.01, key="sst_rsd")
    with col2:
        sst_tailing = st.number_input("Tailing factor (mean)", 0.0, 10.0, 0.0, 0.01, key="sst_tailing")
    with col3:
        sst_plates = st.number_input("Theoretical plates N (mean)", 0, 100000, 0, 100, key="sst_plates")

    st.markdown("### Linearity")
    st.write("Minimum 5 concentration levels across 80–120% of target concentration.")
    col1, col2 = st.columns(2)
    with col1:
        linearity_r2 = st.number_input("Correlation coefficient R²", 0.0, 1.0, 0.999, 0.0001, format="%.4f", key="lin_r2")
    with col2:
        linearity_levels = st.number_input("Number of concentration levels tested", 0, 20, 5, 1, key="lin_levels")

    st.markdown("### Accuracy (Recovery)")
    st.write("3 levels × 3 replicates: 80%, 100%, 120% of target concentration.")
    col1, col2, col3 = st.columns(3)
    with col1:
        acc_80 = st.number_input("Mean recovery at 80% level (%)", 0.0, 200.0, 0.0, 0.1, key="acc_80")
    with col2:
        acc_100 = st.number_input("Mean recovery at 100% level (%)", 0.0, 200.0, 0.0, 0.1, key="acc_100")
    with col3:
        acc_120 = st.number_input("Mean recovery at 120% level (%)", 0.0, 200.0, 0.0, 0.1, key="acc_120")
    acc_rsd = st.number_input("Overall %RSD of recoveries", 0.0, 100.0, 0.0, 0.01, key="acc_rsd")

    st.markdown("### Precision")
    col1, col2 = st.columns(2)
    with col1:
        rep_rsd = st.number_input("Repeatability %RSD (6 replicates, same day)", 0.0, 100.0, 0.0, 0.01, key="rep_rsd")
    with col2:
        rep_mean = st.number_input("Repeatability mean assay (%)", 0.0, 200.0, 100.0, 0.1, key="rep_mean")

    st.markdown("#### Intermediate Precision (Day 2 / Different Analyst)")
    st.write("Enter the 6 individual assay results from the second day/analyst:")

    inter_values = []
    cols = st.columns(3)
    for i in range(6):
        with cols[i % 3]:
            val = st.number_input(f"IP Result {i+1} (%)", 0.0, 200.0, 0.0, 0.1, key=f"ip_{i}")
            inter_values.append(val)
            
    inter_valid = [v for v in inter_values if v > 0]
    inter_rsd_calc = 0.0
    inter_mean = 0.0
    if len(inter_valid) >= 6:
        inter_mean = np.mean(inter_valid)
        inter_rsd_calc = (np.std(inter_valid, ddof=1) / inter_mean * 100) if inter_mean > 0 else 0.0
        diff_means = abs(inter_mean - rep_mean)
        st.caption(f"**Calculated IP %RSD:** {inter_rsd_calc:.2f}%")
        st.caption(f"**Difference between Day 1 and Day 2 means:** {diff_means:.2f}% (must be ≤ 2.0%)")
    elif any(v > 0 for v in inter_values):
        st.warning("Minimum 6 determinations required")

    st.markdown("### LOD / LOQ")
    sn_method = st.radio(
        "How was signal-to-noise ratio measured?",
        ["Automatic calculation by chromatography software (recommended)",
         "Manual measurement — peak height divided by peak-to-peak baseline noise",
         "Standard deviation method — from blank injections",
         "Not yet measured"])

    if sn_method == "Manual measurement — peak height divided by peak-to-peak baseline noise":
        st.warning("⚠️ Manual S/N measurement is least reproducible. Document the exact noise region used (start and end retention time). Ensure noise window is at least 20× the peak width and free of peaks.")

    if sn_method == "Not yet measured":
        st.info("LOD and LOQ will be reported as MISSING DATA in the validation assessment.")

    col1, col2 = st.columns(2)
    with col1:
        lod_sn = st.number_input("Signal-to-noise at LOD", 0.0, 1000.0, 0.0, 0.1, key="lod_sn")
    with col2:
        loq_sn = st.number_input("Signal-to-noise at LOQ", 0.0, 1000.0, 0.0, 0.1, key="loq_sn")
    loq_rsd = st.number_input("LOQ %RSD (6 replicates)", 0.0, 100.0, 0.0, 0.01, key="loq_rsd")

    st.markdown("### Solution Stability")
    col1, col2 = st.columns(2)
    with col1:
        stab_24 = st.number_input("% change from T0 at 24 hours", 0.0, 100.0, 0.0, 0.1, key="stab_24")
    with col2:
        stab_48 = st.number_input("% change from T0 at 48 hours", 0.0, 100.0, 0.0, 0.1, key="stab_48")

    st.divider()

    if st.button("✅ Run Full Validation Assessment"):
        api_key = os.getenv("ANTHROPIC_API_KEY")

        VALIDATION_PROMPT = """You are a senior pharmaceutical analytical chemist and regulatory affairs specialist.

Assess the following HPLC method validation data against ICH Q2(R2) acceptance criteria.

ICH Q2(R2) ACCEPTANCE CRITERIA — apply these exactly:

SPECIFICITY:
- Blank: no peak at analyte retention time
- Placebo: no excipient co-eluting with analyte  
- Peak purity index ≥ 0.99 where PDA is available
- Stability-indicating: all degradants resolved from API (Rs ≥ 2.0)

SYSTEM SUITABILITY:
- Peak area %RSD ≤ 2.0% (6 injections)
- Tailing factor ≤ 2.0
- Theoretical plates N ≥ 2000

LINEARITY:
- R² ≥ 0.999 for assay methods
- Minimum 5 concentration levels required

ACCURACY:
- Mean recovery 98.0–102.0% at each level
- Overall %RSD ≤ 2.0%

PRECISION:
- Repeatability %RSD ≤ 2.0%
- Intermediate precision %RSD ≤ 3.0%

LOD:
- Signal-to-noise ≥ 3:1

LOQ:
- Signal-to-noise ≥ 10:1
- %RSD at LOQ ≤ 10%

SOLUTION STABILITY:
- ≤ 2.0% change from T0 at each time point

YOUR OUTPUT FORMAT:

OVERALL VALIDATION STATUS: [VALIDATED / NOT VALIDATED / INCOMPLETE]

Then for each parameter:
PARAMETER NAME
- Value submitted: [value]
- Acceptance criterion: [criterion]
- Status: PASS / FAIL / MISSING DATA
- Notes: [explanation if fail or borderline]

End with:
PRIORITY ACTIONS — numbered list of what must be done before this method can be submitted."""

        validation_data = f"""
COMPOUND: {st.session_state.compound}
METHOD TYPE: {st.session_state.method_type}
REGULATORY FRAMEWORK: ICH Q2(R2)

VALIDATION DATA:

Specificity:
- Blank injection: {blank_result}
- Placebo injection: {placebo_result}
- Peak purity: {peak_purity}
- Degradant resolution: {degradant_resolution}

System Suitability:
- Peak area %RSD: {sst_rsd}%
- Tailing factor: {sst_tailing}
- Theoretical plates: {sst_plates}

Linearity:
- R²: {linearity_r2}
- Number of levels: {linearity_levels}

Accuracy:
- Recovery at 80%: {acc_80}%
- Recovery at 100%: {acc_100}%
- Recovery at 120%: {acc_120}%
- Overall %RSD: {acc_rsd}%

Precision:
- Repeatability mean: {rep_mean}%
- Repeatability %RSD: {rep_rsd}%
- Intermediate precision mean: {inter_mean:.2f}%
- Intermediate precision %RSD: {inter_rsd_calc:.2f}%

LOD/LOQ:
- S/N Method: {sn_method}
- S/N at LOD: {lod_sn}
- S/N at LOQ: {loq_sn}
- LOQ %RSD: {loq_rsd}%

Solution Stability:
- % change at 24h: {stab_24}%
- % change at 48h: {stab_48}%
"""

        with st.spinner("We are assessing your validation data against ICH Q2(R2) criteria..."):
            try:
                client = anthropic.Anthropic(api_key=api_key)
                message = client.messages.create(
                    model="claude-sonnet-4-6",
                    max_tokens=4096,
                    system=VALIDATION_PROMPT,
                    messages=[{"role": "user", "content": validation_data}]
                )
                validation_report = message.content[0].text
                st.session_state.validation_report = validation_report
            except Exception as e:
                st.error(f"Claude API error: {str(e)}")

    if st.session_state.get("validation_report"):
        st.divider()
        st.markdown("### Validation Assessment Report")
        st.markdown(st.session_state.validation_report)

        st.divider()
        st.success("🎉 Validation assessment complete. Download your report below.")
        st.download_button(
            label="📄 Download Validation Report",
            data=st.session_state.validation_report,
            file_name=f"validation_report_{st.session_state.compound}.txt",
            mime="text/plain"
        )

    col1, _ = st.columns([1, 4])
    with col1:
        if st.button("← Back"):
            go_back()
            st.rerun()

# STAGE 8+
elif st.session_state.stage >= 8:
    st.subheader(f"Stage {st.session_state.stage} — Coming soon")
    p = st.session_state.properties or {}
    st.markdown(f"""
    <div class="stage-card">
        <b>Method type:</b> {st.session_state.method_type}<br>
        <b>Matrix:</b> {st.session_state.get('matrix', 'N/A')}<br>
        <b>Compound:</b> {st.session_state.compound}<br>
        <b>LogP:</b> {p.get('logp', 'N/A')}<br>
        <b>Molecular Weight:</b> {p.get('mw', 'N/A')} g/mol
    </div>
    """, unsafe_allow_html=True)
    st.info(f"✅ Stages 0, 1, 2, and 3 complete. Next we build Stage {st.session_state.stage} — Method Scouting.")
    col1, _ = st.columns([1, 4])
    with col1:
        if st.button("← Back"):
            go_back()
            st.rerun()
