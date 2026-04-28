import streamlit as st
import requests

st.set_page_config(page_title="AMV-AI", page_icon="🔬", layout="centered")

st.markdown("""
<style>
    .stApp { background-color: #ffffff; }
    .stApp, .stMarkdown, p, label { color: #1a1a2e !important; }
    h1 { color: #0077b6 !important; font-size: 2.2rem !important; }
    h2, h3 { color: #0096c7 !important; }
    .stButton > button {
        background-color: #0077b6 !important;
        color: white !important;
        border: none !important;
        border-radius: 8px !important;
        padding: 0.5rem 2rem !important;
        font-size: 1rem !important;
    }
    .stButton > button:hover { background-color: #0096c7 !important; }
    .stRadio > label { color: #1a1a2e !important; font-weight: 500 !important; }
    hr { border-color: #90e0ef !important; }
    .stProgress > div > div { background-color: #0077b6 !important; }
    [data-testid="stMetricValue"] { color: #0077b6 !important; font-weight: 700 !important; }
    [data-testid="stMetricLabel"] { color: #1a1a2e !important; }
    .stage-card {
        background-color: #f0f9ff;
        border-left: 4px solid #0077b6;
        border-radius: 8px;
        padding: 1.2rem 1.5rem;
        margin: 1rem 0;
    }
    .prop-card {
        background-color: #f0f9ff;
        border-radius: 8px;
        padding: 1rem 1.5rem;
        margin: 0.5rem 0;
        border: 1px solid #90e0ef;
        color: #1a1a2e !important;
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
st.progress(st.session_state.stage / 7)
st.caption(f"Stage {st.session_state.stage} of 7")
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

def interpret_logp(logp):
    try:
        v = float(logp)
        if v < 0:
            return "Very polar — consider HILIC column"
        elif v < 1:
            return "Polar — C8 or mixed-mode column"
        elif v < 4:
            return "Moderate — C18 reverse phase (ideal range)"
        elif v < 5:
            return "Lipophilic — C18 with high organic start"
        else:
            return "Very lipophilic — C18 or C4, high organic"
    except:
        return "N/A"

def interpret_uv(smiles):
    smiles = smiles or ""
    st.write("DEBUG SMILES:", smiles)

    # Aromatic ring — lowercase letters in SMILES
    has_aromatic = any(c in smiles for c in ["c", "n", "o", "s"] if c.islower())

    # Carbonyl
    has_carbonyl = "C=O" in smiles or "C(=O)" in smiles

    # Sulfur-sulfur bond — CSSSC, CSSC etc
    has_sulfur_bond = "SS" in smiles

    # Nitro
    has_nitro = "N(=O)" in smiles or "[N+](=O)" in smiles

    if has_aromatic:
        return "Aromatic ring detected — start at 254 nm. Also try 220 nm for higher sensitivity."
    elif has_nitro:
        return "Nitro group detected — absorbs at 254 nm and 330 nm. Try both."
    elif has_carbonyl:
        return "Carbonyl group detected — try 210–215 nm. Expect some background noise."
    elif has_sulfur_bond:
        return "Sulfur-sulfur bond detected (disulfide/trisulfide) — absorbs at 210–220 nm. Your compound should be detectable around 215 nm."
    else:
        return "⚠️ No common UV chromophore detected — UV detection may not work reliably. Consider MS, ELSD, or CAD."

# STAGE 0
if st.session_state.stage == 0:
    st.subheader("Stage 0 — Instrument Qualification")
    st.write("Has your HPLC instrument been formally qualified?")

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

    method_type = st.radio(
        "Select the purpose of your method:",
        [
            "Measure the amount of active drug (API) in a product",
            "Separate and measure the drug AND all its impurities/degradants",
            "Measure impurities only (not the API content)",
            "Test how fast the drug releases from the dosage form",
            "Confirm the drug is present — not measuring how much",
            "Measure the drug in a biological sample (plasma, urine, tissue)",
            "Check dose-to-dose consistency across tablets or capsules"
        ],
        key="q1"
    )

    method_map = {
        "Measure the amount of active drug (API) in a product": "Potency Assay",
        "Separate and measure the drug AND all its impurities/degradants": "Stability-Indicating Assay",
        "Measure impurities only (not the API content)": "Impurity/Related Substances",
        "Test how fast the drug releases from the dosage form": "Dissolution",
        "Confirm the drug is present — not measuring how much": "Identification Test",
        "Measure the drug in a biological sample (plasma, urine, tissue)": "Bioanalytical",
        "Check dose-to-dose consistency across tablets or capsules": "Content Uniformity"
    }

    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("← Back"):
            go_back()
            st.rerun()
    with col2:
        if st.button("Continue →"):
            st.session_state.method_type = method_map[method_type]
            go_next()
            st.rerun()

# STAGE 2
elif st.session_state.stage == 2:
    st.subheader("Stage 2 — Analyte Properties")
    st.write("Enter your compound name and we will fetch its molecular properties from PubChem.")

    compound_name = st.text_input(
        "Compound name:",
        placeholder="e.g. ibuprofen, dimethyl trisulfide, metformin"
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

        st.markdown("### What this means for your method")

        uv_result = interpret_uv(smiles)

        st.markdown(f"""
        <div class="prop-card">
            <b>🧪 Column recommendation:</b><br>
            LogP = {p['logp']} → {interpret_logp(p['logp'])}
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

# STAGE 3+
elif st.session_state.stage >= 3:
    st.subheader(f"Stage {st.session_state.stage} — Coming soon")
    p = st.session_state.properties or {}
    st.markdown(f"""
    <div class="stage-card">
        <b>Method type:</b> {st.session_state.method_type}<br>
        <b>Compound:</b> {st.session_state.compound}<br>
        <b>LogP:</b> {p.get('logp', 'N/A')}<br>
        <b>Molecular Weight:</b> {p.get('mw', 'N/A')} g/mol
    </div>
    """, unsafe_allow_html=True)
    st.info("✅ Stages 0, 1, and 2 complete. Next we build Stage 3 — Method Scouting.")
    col1, _ = st.columns([1, 4])
    with col1:
        if st.button("← Back"):
            go_back()
            st.rerun()
