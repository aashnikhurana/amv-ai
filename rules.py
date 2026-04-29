# rules.py

CHROMOPHORE_RULES = {
    "aromatic_ring": {
        "description": "Aromatic ring(s) detected",
        "recommendation": (
            "Strong UV around 200–210 nm; practical choices 220–254 nm depending on "
            "mobile-phase cut-off. Start at 254 nm and also monitor 220 nm for sensitivity."
        ),
    },
    "heteroaromatic_ring": {
        "description": "Heteroaromatic ring(s) detected",
        "recommendation": (
            "Strong UV around 200–230 nm; typical choices 230–260 nm. "
            "Start at 254 nm and adjust based on matrix interference."
        ),
    },
    "conjugated_alkene": {
        "description": "Conjugated C=C system",
        "recommendation": (
            "π→π* transitions give useful UV response around 200–230 nm; "
            "consider 210–220 nm as primary wavelength."
        ),
    },
    "carbonyl": {
        "description": "Carbonyl chromophore",
        "recommendation": (
            "Strong π→π* band near 190–200 nm and weaker n→π* near 270–300 nm; "
            "use 210–220 nm if background permits, or 260–280 nm for selectivity."
        ),
    },
    "disulfide": {
        "description": "Sulfur/disulfide chromophore",
        "recommendation": (
            "Strong absorption at 210–215 nm. Use ACN and low-UV buffers; "
            "avoid MeOH at very low wavelengths."
        ),
    },
    "nitro": {
        "description": "Nitro group",
        "recommendation": (
            "Strong bands in the 250–280 nm range and often a higher-wavelength band; "
            "start with 254 nm and optionally 330 nm if baseline allows."
        ),
    },
    "azo": {
        "description": "Azo (–N=N–) chromophore",
        "recommendation": (
            "Absorbs strongly in 300–500 nm region depending on substitution; "
            "use reported λmax when known, otherwise screen 254–400 nm."
        ),
    },
    "nucleic_acid": {
        "description": "DNA/RNA / nucleic acid",
        "recommendation": (
            "λmax ≈ 260 nm. Pure DNA typically has A260/A280 ≈ 1.8; pure RNA ≈ 2.0; "
            "A260/A230 ≈ 2.1–2.5 for clean nucleic acid."
        ),
    },
    "peptide_or_protein": {
        "description": "Peptide/protein",
        "recommendation": (
            "Aromatic residues absorb at ~280 nm (Trp/Tyr); peptide backbone absorbs "
            "strongly at 190–230 nm (commonly monitored at 214–230 nm)."
        ),
    },
    "alkene": {
        "description": "Isolated C=C (alkene)",
        "recommendation": (
            "Weak chromophore with λmax around 180–190 nm; may be monitored at ~200 nm "
            "if mobile-phase background is acceptable."
        ),
    },
    "alkyne": {
        "description": "C≡C (alkyne)",
        "recommendation": (
            "Very weak UV chromophore; any signal will be in deep UV (<200 nm). "
            "Consider MS or other detectors for robust quantitation."
        ),
    },
}

COLUMN_RULES = {
    "peptide_large": {
        "recommendation": "Wide-pore RP (C8 or C18, e.g., 300 Å) — optimized for peptides/proteins.",
        "explanation": (
            "Use a wide-pore C8 or C18 column with gradient elution and volatile buffers. "
            "Bias toward C8 when very hydrophobic to avoid irreversible binding."
        ),
    },
    "very_polar": {
        "recommendation": "HILIC or Ion-Exchange — analyte is highly polar.",
        "explanation": (
            "LogP < 0 and high TPSA suggest k' < 1 on RP. Use HILIC or ion-exchange. "
            "If RP is mandatory, use Aqueous C18 at very low organic."
        ),
    },
    "polar": {
        "recommendation": "HILIC or Aqueous C18 (AQ) — polar analyte.",
        "explanation": (
            "RP retention may be limited. Start with low organic on an AQ-stable C18 "
            "or switch to HILIC for better retention."
        ),
    },
    "very_hydrophobic": {
        "recommendation": "C8 Reversed-Phase — very hydrophobic analyte.",
        "explanation": (
            "LogP > 4 indicates strong hydrophobicity. C18 may retain too strongly; "
            "C8 (or C4) shortens run time and improves elution, with strong organic washes."
        ),
    },
    "aromatic_bulky": {
        "recommendation": "C8 or Phenyl-Hexyl — bulky aromatic analyte.",
        "explanation": (
            "Bulky hydrophobic aromatics often elute slowly on C18. C8 reduces retention; "
            "phenyl phases can improve π–π selectivity for closely related impurities."
        ),
    },
    "aromatic_polar": {
        "recommendation": "Phenyl or Embedded-Polar C18 — polar aromatic analyte.",
        "explanation": (
            "Polar aromatics may show improved peak shape and selectivity on phenyl or "
            "embedded-polar RP phases compared with standard C18."
        ),
    },
    "lipid_matrix": {
        "recommendation": "C8 with Guard Column — lipid-rich matrix.",
        "explanation": (
            "Lipid-rich samples can foul C18 columns. Use C8 with a guard column and strong "
            "organic wash steps to minimize irreversible binding."
        ),
    },
    "stability_indicating": {
        "recommendation": "C18 Reversed-Phase — robust, general-purpose selectivity.",
        "explanation": (
            "For stability-indicating and impurity methods, start with a well-characterized C18. "
            "If critical pairs are not resolved, screen C8, phenyl, and embedded-polar phases."
        ),
    },
    "default": {
        "recommendation": "C18 Reversed-Phase — default small-molecule choice.",
        "explanation": (
            "For typical small molecules with moderate logP, C18 offers reliable retention, "
            "wide availability, and well-understood behavior."
        ),
    },
}

SAMPLE_PREP_RULES = {
    "plasma_serum_blood": {
        "summary": "Biological samples require protein removal and clarification before HPLC injection.",
        "steps": [
            "Perform protein precipitation (PPT) with 3–4× volume of cold ACN or MeOH (or acidified ACN/MeOH).",
            "Vortex thoroughly and incubate on ice if needed.",
            "Centrifuge at ≥10,000 × g for at least 5–10 minutes.",
            "Optionally apply SPE or LLE to further clean up the extract and reduce matrix effects.",
            "Filter the supernatant through a 0.2–0.45 μm filter before injection.",
        ],
        "warnings": [
            "Do not inject un-centrifuged or unfiltered biological samples; they will clog the column and increase backpressure.",
            "Strongly protein-bound analytes may co-precipitate during PPT; assess recovery and consider SPE/LLE if low.",
            "Hemolyzed or lipemic specimens may require additional cleanup or dilution to maintain robustness.",
        ],
    },
    "urine": {
        "summary": "Urine has lower protein content but may still require cleanup.",
        "steps": [
            "Centrifuge to remove particulates.",
            "Filter (0.2–0.45 μm) prior to injection.",
            "For trace-level analytes, use SPE or LLE for preconcentration and cleanup.",
        ],
        "warnings": [
            "Salts and endogenous metabolites can still affect column performance; use a guard column for long sequences.",
        ],
    },
    "cell_culture": {
        "summary": "Cell culture media contain proteins, salts, and surfactants that need removal.",
        "steps": [
            "Centrifuge to remove cells and debris.",
            "Perform PPT or ultrafiltration to remove proteins.",
            "Consider SPE to remove surfactants and concentrate analyte.",
            "Filter before injection.",
        ],
        "warnings": [
            "Residual proteins/surfactants can cause severe tailing and fouling; assess recovery for the chosen cleanup strategy.",
        ],
    },
    "tablet_capsule": {
        "summary": "Ensure complete extraction of API from solid dosage forms.",
        "steps": [
            "Crush or disperse the tablet/capsule contents.",
            "Dissolve or extract in a suitable solvent system (e.g., buffer/organic mixture) based on analyte solubility.",
            "Sonicate or mechanically shake to ensure complete extraction.",
            "Centrifuge or filter (0.45 μm or 0.2 μm) to remove undissolved excipients.",
        ],
        "warnings": [
            "Incomplete dissolution or extraction leads to low assay values and poor precision.",
            "Do not inject suspensions containing visible particles.",
        ],
    },
    "water_env": {
        "summary": "Environmental water samples often require filtration and enrichment.",
        "steps": [
            "Filter the sample (0.2–0.45 μm) to remove particulates.",
            "For low-level analytes, perform SPE or LLE to concentrate and clean up before injection.",
        ],
        "warnings": [
            "Natural organic matter can slowly foul columns; use a guard column and periodic strong washes.",
        ],
    },
    "oil_lipid": {
        "summary": "Oil and high-fat matrices contain large amounts of lipids that can foul RP columns.",
        "steps": [
            "If possible, perform a defatting step with a non-polar solvent (e.g., hexane).",
            "Use SPE or LLE to isolate the analyte from the lipid matrix.",
            "Include a guard column and apply strong organic washes between runs.",
            "Filter the final extract before injection.",
        ],
        "warnings": [
            "Injecting crude lipid extracts onto analytical columns can cause irreversible fouling, especially on C18.",
            "Consider sacrificing inexpensive guard columns regularly when working with very fatty matrices.",
        ],
    },
    "default": {
        "summary": (
            "Basic cleanup (centrifugation and filtration) protects the column and improves robustness."
        ),
        "steps": [
            "Centrifuge the sample if particulates are present.",
            "Filter through a 0.45 μm or 0.2 μm filter before injection.",
        ],
        "warnings": [
            "Particulates will block frits and increase backpressure; always clarify samples where feasible.",
        ],
    },
}

VALIDATION_SST_RULES = {
    "usp_621": {
        "description": "System Suitability (USP <621>)",
        "items": [
            "Establish system suitability criteria: retention time window, theoretical plates (N), tailing factor/asymmetry, and resolution for critical pairs.",
            "Evaluate %RSD of replicate standard injections (typically n ≥ 5).",
            "Ensure SST passes before reporting any sample results; failed SST invalidates the analytical run.",
        ],
    },
    "stability_indicating": {
        "description": "Stability-Indicating Method Requirements",
        "items": [
            "Perform forced degradation (acid/base, oxidation, heat, light) to generate relevant degradants.",
            "Demonstrate baseline separation between the main analyte peak and its degradants.",
            "Show that the method can quantify the intact analyte in the presence of degradation products.",
        ],
    },
    "bioanalytical": {
        "description": "Bioanalytical Validation Considerations",
        "items": [
            "Assess matrix effects using post-column or post-extraction addition with internal standard normalization.",
            "Validate extraction recovery, carryover, and dilution integrity across the calibration range.",
            "Include QC samples at low, mid, and high concentrations in every analytical batch.",
        ],
    },
    "robustness": {
        "description": "Robustness and Ruggedness",
        "items": [
            "Challenge the method by small deliberate changes (pH ±0.2, organic ±2–5 %, temperature and flow variations).",
            "Evaluate impact of column lot changes, minor mobile-phase composition variations, and different analysts/instruments.",
            "Confirm that system suitability and key performance metrics remain within acceptable limits under these variations.",
        ],
    },
    "dissolution": {
        "description": "Dissolution-specific checks",
        "items": [
            "Confirm filtration does not adsorb analyte.",
            "Verify solution stability in dissolution medium over the full test window.",
            "Demonstrate acceptable precision at the target sampling time points.",
            "Confirm the method is suitable for the dissolution medium composition and pH."
        ],
    },
    "identification": {
        "description": "Identification-focused checks",
        "items": [
            "Confirm analyte identity by retention time match to a reference standard.",
            "Use PDA or spectral confirmation when available.",
            "Demonstrate specificity against excipients or matrix peaks."
        ],
    },
    "content_uniformity": {
        "description": "Content uniformity checks",
        "items": [
            "Prepare individual dosage units separately; do not pool samples.",
            "Demonstrate consistent extraction efficiency across units.",
            "Verify method precision is appropriate for unit-to-unit comparison.",
            "Ensure assay calculations use the correct analyte form (salt / hydrate / free base)."
        ],
    },
}
