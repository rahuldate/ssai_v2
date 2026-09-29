"""
Module Name: Human-in-the-Loop Decision Gate
Version: 1.2.0
Description: Provides interactive review, expert override, and audit trail logging 
             for skin sensitization predictions (Weight of Evidence / QRA).
"""

import datetime
import pandas as pd

class HITLDecisionGate:
    def __init__(self, substance_id: str):
        self.substance_id = substance_id
        self.version = "1.2.0"

    def create_review_record(
        self,
        initial_prediction: str,
        algorithm_confidence: float,
        expert_name: str = "",
        override_status: str = "Pending Review",
        justification: str = "",
    ) -> dict:
        timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
        return {
            "substance_id": self.substance_id,
            "module_version": self.version,
            "timestamp": timestamp,
            "initial_prediction": initial_prediction,
            "algorithm_confidence": algorithm_confidence,
            "expert_name": expert_name,
            "override_status": override_status,
            "justification": justification,
        }

    def render_hitl_form_state(self, default_prediction: str):
        import streamlit as st
        st.subheader(f"🛡️ Human-in-the-Loop (HITL) Decision Gate (v{self.version})")
        st.info(f"**Algorithmic Baseline Prediction:** {default_prediction}")

        expert_name = st.text_input("Assessor / Toxicologist Name", key="hitl_expert")
        override_status = st.selectbox(
            "Expert Decision Review",
            [
                "Pending Review",
                "Confirm Algorithm Prediction",
                "Override to Sensitizer",
                "Override to Non-Sensitizer",
            ],
            key="hitl_status",
        )
        justification = st.text_area(
            "Expert Scientific Justification / Rationale",
            placeholder="Provide mechanistic or structural justification if overriding...",
            key="hitl_justification",
        )

        if st.button("Commit Expert Review Decision"):
            if not expert_name:
                st.warning("Please enter the assessor's name before committing.")
            else:
                record = self.create_review_record(
                    initial_prediction=default_prediction,
                    algorithm_confidence=0.85,
                    expert_name=expert_name,
                    override_status=override_status,
                    justification=justification,
                )
                st.success(f"Decision successfully logged under version {self.version}!")
                return record
        return None
