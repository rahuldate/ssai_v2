import streamlit as st
from modules.structure import render_structure_module
from modules.woe import render_woe_module
from modules.read_across import render_read_across_module
from modules.potency import render_potency_module
from modules.qra import render_qra_module
from modules.dossier import render_dossier_module

st.set_page_config(page_title="Skin Sensitization Decision Support", page_icon="🧬", layout="wide")

st.markdown(
    '<div style="font-size:32px;font-weight:800;color:#0d6efd;margin-bottom:4px;">'
    '🧬 Skin Sensitization Decision-Support Tool</div>',
    unsafe_allow_html=True,
)
st.caption(
    "A calculation aid for structural screening, assay-based hazard classification, "
    "and dermal-exposure risk assessment. Every number you see is computed from what "
    "you enter — there is no pre-filled or simulated 'demo' output. Not a substitute "
    "for expert toxicological judgment or a validated regulatory tool."
)

page = st.sidebar.radio(
    "Workflow step",
    [
        "1. Structure & Alerts",
        "2. Assay-Based WoE",
        "3. Read-Across Prediction",
        "4. Potency (1A/1B)",
        "5. QRA (Exposure & Risk)",
        "6. Export Summary",
    ],
)
st.sidebar.markdown("---")
st.sidebar.caption(
    "Steps build on each other in this session: structure → hazard → read-across → "
    "potency → risk → export. You can also use any step on its own."
)

if page.startswith("1"):
    render_structure_module()
elif page.startswith("2"):
    render_woe_module()
elif page.startswith("3"):
    render_read_across_module()
elif page.startswith("4"):
    render_potency_module()
elif page.startswith("5"):
    render_qra_module()
else:
    render_dossier_module()

st.markdown(
    '<div style="text-align:center;font-size:12px;color:#6c757d;margin-top:40px;'
    'border-top:1px solid #dee2e6;padding-top:12px;">'
    'Decision-support calculator — outputs require review by a qualified toxicologist.'
    '</div>',
    unsafe_allow_html=True,
)
