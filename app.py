import streamlit as st

st.set_page_config(page_title="AMV-AI", page_icon="🔬", layout="centered")

# ── STYLING ──────────────────────────────────────
st.markdown("""
<style>
    /* White background */
    .stApp {
        background-color: #ffffff;
    }

    /* Main text dark */
    .stApp, .stMarkdown, p, label {
        color: #1a1a2e !important;
    }

    /* Title */
    h1 {
        color: #0077b6 !important;
        font-size: 2.2rem !important;
    }

    /* Subheaders */
    h2, h3 {
        color: #0096c7 !important;
    }

    /* Primary button */
    .stButton > button {
        background-color: #0077b6 !important;
        color: white !important;
        border: none !important;
        border-radius: 8px !important;
        padding: 0.5rem 2rem !important;
        font-size: 1rem !important;
    }

    .stButton > button:hover {
        background-color: #0096c7 !important;
    }

    /* Radio buttons */
    .stRadio > label {
        color: #1a1a2e !important;
        font-weight: 500 !important;
    }

    /* Divider */
    hr {
        border-color: #90e0ef !important;
    }

    /* Success box */
    .stSuccess {
        background-color: #e0f7fa !important;
        color: #006064 !important;
    }

    /* Info box */
    .stInfo {
        background-color: #e1f5fe !important;
        color: #01579b !important;
    }

    /* Warning box */
    .stWarning {
        background-color: #fff8e1 !important;
    }

    /* Error box */
    .stError {
        background-color: #fce4ec !important;
    }

    /* Progress bar */
    .stProgress > div > div {
        background-color: #0077b6 !important;
    }

    /* Caption */
    .stCaption {
        color: #0096c7 !important;
    }

    /* Card style for stage content */
    .stage-card {
        background-color: #f0f9ff;
        border-left: 4px solid #0077b6;
        border-radius: 8px;
        padding: 1.2rem 1.5rem;
        margin: 1rem 0;
    }
</style>
""", unsafe_allow_html=True)

# ── SESSION STATE SETUP ───────────────────────────
if "stage" not in st.session_state:
    st.session_state.stage = 0
if "qualification" not in st.session_state:
    st.session_state.qualification = None
if "method_type" not in st.session_state:
    st.session_state.method_type = None

# ── HEADER ────────────────────────────────────────
st.title("🔬 AMV-AI")
st.caption("Analytical Method Lifecycle Tool")

# ── PROGRESS BAR ──────────────────────────────────
total_stages = 7
progress = st.session_state.stage / total_stages
st.progress(progress)
st.caption(f"Stage {st.session_state.stage} of {total_stages}")

st.divider()

# ── NAVIGATION HELPERS ────────────────────────────
def go_next():
    st.session_state.stage += 1

def go_back():
    st.session_state.stage -= 1

# ── STAGE 0 ───────────────────────────────────────
if st.session_state.stage == 0:

    st.subheader("Stage 0 — Instrument Qualification")
    st.write("Before we begin — has your HPLC instrument been formally qualified?")
    st.write("This is required before any method development or validation can begin.")

    qualification = st.radio(
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
        if qualification == "No — instrument has not been qualified":
            st.error("❌ Instrument qualification is required before proceeding. Without DQ/IQ/OQ/PQ documentation, no regulatory body will accept your data.")
        elif qualification == "I don't know what instrument qualification is":
            st.info("📚 Instrument qualification is a four-step process (DQ, IQ, OQ, PQ) that proves your HPLC is working correctly before any regulated analysis. Please complete this with your instrument vendor before continuing.")
        else:
            st.session_state.qualification = qualification
            go_next()
            st.rerun()

# ── STAGE 1 ───────────────────────────────────────
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

# ── PLACEHOLDER FOR NEXT STAGES ───────────────────
elif st.session_state.stage >= 2:

    st.subheader(f"Stage {st.session_state.stage} — Coming next")

    st.markdown(f"""
    <div class="stage-card">
        <b>Method type selected:</b> {st.session_state.method_type}<br>
        <b>Instrument status:</b> {st.session_state.qualification}
    </div>
    """, unsafe_allow_html=True)

    st.info("✅ Stages 0 and 1 complete. Next we build Stage 2 — Analyte Properties.")

    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("← Back"):
            go_back()
            st.rerun()

