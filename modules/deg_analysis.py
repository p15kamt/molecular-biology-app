"""
Ολοκληρωμένο Module για Differential Gene Expression Analysis

Αυτό το module περιλαμβάνει όλες τις λειτουργίες που απαιτούνται για την
ανάλυση διαφορικής γονιδιακής έκφρασης, συμπεριλαμβανομένου στατιστικού
ελέγχου, correction για multiple testing, και οπτικοποίηση αποτελεσμάτων.


Ημερομηνία: 2025
"""

import streamlit as st
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats, sparse
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import io
import tempfile
import os
import sys
from pathlib import Path
import warnings

# Προσθήκη utils directory στο path
sys.path.append(str(Path(__file__).parent.parent / "utils"))

try:
    from memory_utils import memory_monitor, display_memory_info
    from advanced_memory import memory_manager
except ImportError:
    def memory_monitor(func):
        return func
    def display_memory_info():
        pass

class DEGAnalysisPageComplete:
    """Ολοκληρωμένη κλάση για ανάλυση διαφορικής γονιδιακής έκφρασης"""
    
    def __init__(self):
        self.adata = None
        self.deg_results = None
        self.comparison_groups = {}
        self.selected_method = "wilcoxon"
        self.max_genes = 50000  # PRODUCTION DEFAULT: όλα τα γονίδια
        
    def render(self):
        """Κεντρική μέθοδος για την εμφάνιση της σελίδας"""
        
        st.title("🧬 Ανάλυση Διαφορικής Έκφρασης")
        st.markdown("### Στατιστική ανάλυση γονιδιακής έκφρασης μεταξύ ομάδων")
        
        # Memory information
        display_memory_info()
        
        # Tabs για διαφορετικές λειτουργίες
        tab1, tab2, tab3, tab4 = st.tabs([
            "📊 Δεδομένα & Ομάδες", 
            "⚙️ Παράμετροι Ανάλυσης",
            "🔄 Εκτέλεση DEG",
            "📈 Αποτελέσματα & Plots"
        ])
        
        with tab1:
            self.render_data_selection()
            
        with tab2:
            self.render_analysis_parameters()
            
        with tab3:
            self.render_deg_execution()
            
        with tab4:
            self.render_results_visualization()
    
    def render_data_selection(self):
        """Επιλογή δεδομένων και ορισμός ομάδων"""
        
        st.header("📊 Επιλογή Δεδομένων")
        
        # Έλεγχος για διαθέσιμα δεδομένα
        data_sources = []
        
        if 'adata' in st.session_state and st.session_state.adata is not None:
            data_sources.append("Προεπεξεργασμένα Δεδομένα")
        
        if 'integrated_data' in st.session_state and st.session_state.integrated_data is not None:
            data_sources.append("Integrated Δεδομένα")
        
        if not data_sources:
            st.warning("⚠️ Δεν βρέθηκαν διαθέσιμα δεδομένα. Πηγαίνετε πρώτα στο Preprocessing ή Integration module.")
            return
        
        # Επιλογή πηγής δεδομένων
        selected_source = st.selectbox(
            "Επιλέξτε πηγή δεδομένων:",
            data_sources
        )
        
        if selected_source == "Προεπεξεργασμένα Δεδομένα":
            self.adata = st.session_state.adata
        else:
            self.adata = st.session_state.integrated_data
        
        if self.adata is not None:
            # Εμφάνιση πληροφοριών dataset
            st.subheader("📋 Πληροφορίες Dataset")
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Κύτταρα", f"{self.adata.n_obs:,}")
            with col2:
                st.metric("Γονίδια", f"{self.adata.n_vars:,}")
            with col3:
                memory_size = memory_manager.estimate_memory_usage(self.adata)
                st.metric("Μέγεθος", f"{memory_size:.1f} MB")
            
            # Εμφάνιση διαθέσιμων metadata columns
            st.subheader("🏷️ Διαθέσιμα Metadata")
            
            categorical_columns = []
            for col in self.adata.obs.columns:
                if self.adata.obs[col].dtype == 'object' or self.adata.obs[col].dtype.name == 'category':
                    categorical_columns.append(col)
            
            if categorical_columns:
                st.write("**Διαθέσιμες κατηγορικές μεταβλητές:**")
                for col in categorical_columns:
                    unique_values = self.adata.obs[col].nunique()
                    st.write(f"- **{col}**: {unique_values} unique values")
                
                # Επιλογή στήλης για σύγκριση
                self.comparison_column = st.selectbox(
                    "Επιλέξτε στήλη για σύγκριση:",
                    categorical_columns,
                    help="Η στήλη που περιέχει τις ομάδες για σύγκριση"
                )
                
                if self.comparison_column:
                    self.render_group_selection()
            else:
                st.warning("⚠️ Δεν βρέθηκαν κατηγορικές μεταβλητές για σύγκριση")
    
    def render_group_selection(self):
        """Επιλογή ομάδων για σύγκριση"""
        
        st.subheader("👥 Επιλογή Ομάδων για Σύγκριση")
        
        # Εμφάνιση διαθέσιμων ομάδων
        available_groups = self.adata.obs[self.comparison_column].unique()
        
        st.write(f"**Διαθέσιμες ομάδες στη στήλη '{self.comparison_column}':**")
        group_counts = self.adata.obs[self.comparison_column].value_counts()
        
        for group, count in group_counts.items():
            st.write(f"- **{group}**: {count:,} κύτταρα")
        
        # Επιλογή τύπου σύγκρισης
        comparison_type = st.radio(
            "Τύπος σύγκρισης:",
            ["Ένα προς Όλα (One vs All)", "Ένα προς Ένα (Pairwise)"],
            help="Επιλέξτε τον τύπο σύγκρισης που θέλετε να εκτελέσετε"
        )
        
        if comparison_type == "Ένα προς Όλα (One vs All)":
            self.comparison_groups['type'] = 'one_vs_all'
            self.comparison_groups['reference_group'] = st.selectbox(
                "Επιλέξτε ομάδα αναφοράς:",
                available_groups,
                help="Η ομάδα που θα συγκριθεί με όλες τις άλλες"
            )
        
        else:  # Pairwise
            self.comparison_groups['type'] = 'pairwise'
            
            col1, col2 = st.columns(2)
            with col1:
                self.comparison_groups['group1'] = st.selectbox(
                    "Ομάδα 1:",
                    available_groups,
                    key="group1_select"
                )
            
            with col2:
                available_group2 = [g for g in available_groups if g != self.comparison_groups.get('group1')]
                if available_group2:
                    self.comparison_groups['group2'] = st.selectbox(
                        "Ομάδα 2:",
                        available_group2,
                        key="group2_select"
                    )
    
    def render_analysis_parameters(self):
        """Ρύθμιση παραμετρών ανάλυσης"""
        
        st.header("⚙️ Παράμετροι Ανάλυσης")
        
        if self.adata is None:
            st.warning("⚠️ Επιλέξτε πρώτα δεδομένα")
            return
        
        # Επιλογή στατιστικής μεθόδου
        col1, col2 = st.columns(2)
        
        with col1:
            self.selected_method = st.selectbox(
                "Στατιστική μέθοδος:",
                ["wilcoxon", "t-test", "logreg"],
                help="Επιλέξτε τη στατιστική μέθοδο για τον έλεγχο"
            )
        
        with col2:
            self.correction_method = st.selectbox(
                "Multiple testing correction:",
                ["benjamini-hochberg", "bonferroni", "none"],
                help="Μέθοδος διόρθωσης για πολλαπλούς ελέγχους"
            )
        
        # Παράμετροι filtering
        st.subheader("🔍 Παράμετροι Filtering")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            self.min_pct = st.slider(
                "Minimum % κυττάρων:",
                min_value=0.0, max_value=1.0, value=0.1, step=0.05,
                help="Ελάχιστο ποσοστό κυττάρων που εκφράζουν το γονίδιο"
            )
        
        with col2:
            self.logfc_threshold = st.slider(
                "Log fold change threshold:",
                min_value=0.0, max_value=2.0, value=0.25, step=0.05,
                help="Ελάχιστο log fold change για σημαντικότητα"
            )
        
        with col3:
            # PRODUCTION MODE: Ανάλυση ΌΛΩΝ των γονιδίων
            self.max_genes = st.number_input(
                "Μέγιστος αριθμός γονιδίων:",
                min_value=100, max_value=50000, value=50000, step=1000,
                help="Μέγιστος αριθμός γονιδίων προς ανάλυση (50K = όλα τα γονίδια)"
            )
            
            # Προειδοποίηση για production use
            if self.max_genes < 10000:
                st.warning("⚠️ Για production analysis, συνιστάται >10K γονίδια")
        
        # Προηγμένες επιλογές
        with st.expander("🔧 Προηγμένες Επιλογές"):
            
            self.use_raw = st.checkbox(
                "Χρήση raw counts",
                value=False,
                help="Χρήση raw counts αντί για normalized"
            )
            
            self.filter_highly_variable = st.checkbox(
                "Φιλτράρισμα για highly variable genes",
                value=False,  # PRODUCTION: Default σε όλα τα γονίδια
                help="Ανάλυση μόνο των highly variable genes (συνιστάται για πολύ μεγάλα datasets)"
            )
            
            if self.selected_method == "wilcoxon":
                self.tie_correct = st.checkbox(
                    "Tie correction",
                    value=True,
                    help="Διόρθωση για tied values στο Wilcoxon test"
                )
    
    @memory_monitor
    def render_deg_execution(self):
        """Εκτέλεση DEG ανάλυσης"""
        
        st.header("🔄 Εκτέλεση DEG Analysis")
        
        if self.adata is None or not self.comparison_groups:
            st.warning("⚠️ Ολοκληρώστε την επιλογή δεδομένων και ομάδων")
            return
        
        # Έλεγχος αν έχουν επιλεγεί οι απαραίτητες ομάδες
        if self.comparison_groups.get('type') == 'pairwise':
            if not self.comparison_groups.get('group1') or not self.comparison_groups.get('group2'):
                st.warning("⚠️ Επιλέξτε και τις δύο ομάδες για pairwise σύγκριση")
                return
        elif self.comparison_groups.get('type') == 'one_vs_all':
            if not self.comparison_groups.get('reference_group'):
                st.warning("⚠️ Επιλέξτε reference group για one-vs-all σύγκριση")
                return
        
        # Προεπισκόπηση ανάλυσης
        st.subheader("📋 Προεπισκόπηση Ανάλυσης")
        
        # DEBUG: Ας δούμε τι υπάρχει στη στήλη
        if self.comparison_column in self.adata.obs.columns:
            unique_groups = self.adata.obs[self.comparison_column].unique()
            st.info(f"🔍 Διαθέσιμες ομάδες στη στήλη '{self.comparison_column}': {list(unique_groups)}")
            st.info(f"🔍 Αριθμός ομάδων: {len(unique_groups)}")
            
            # Εκτίμηση αριθμού comparisons
            if self.comparison_groups['type'] == 'one_vs_all':
                n_comparisons = max(0, len(unique_groups) - 1)
                comparison_desc = f"'{self.comparison_groups.get('reference_group', 'N/A')}' vs άλλες ομάδες"
            else:
                n_comparisons = 1 if len(unique_groups) >= 2 else 0
                group1 = self.comparison_groups.get('group1', 'N/A')
                group2 = self.comparison_groups.get('group2', 'N/A')
                comparison_desc = f"'{group1}' vs '{group2}'"
        else:
            st.error(f"❌ Η στήλη '{self.comparison_column}' δεν υπάρχει στα δεδομένα!")
            n_comparisons = 0
            comparison_desc = "Μη έγκυρη στήλη"
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Συγκρίσεις", n_comparisons)
        with col2:
            n_genes = min(self.max_genes, self.adata.n_vars)
            st.metric("Γονίδια προς ανάλυση", f"{n_genes:,}")
        with col3:
            # Απλή εκτίμηση χρόνου
            estimated_time = self.estimate_analysis_time()
            st.metric("Εκτιμώμενος χρόνος", f"{estimated_time:.1f}min")
        
        st.info(f"**Σύγκριση:** {comparison_desc}")
        st.info(f"**Μέθοδος:** {self.selected_method} με {self.correction_method} correction")
        
        # Απλή εκτίμηση χρόνου
        estimated_time = self.estimate_analysis_time()
        
        with st.expander("⏱️ Analysis Estimate"):
            st.metric("Estimated Time", f"{estimated_time:.1f} minutes")
            st.info(f"Based on {self.adata.n_obs:,} cells and {self.adata.n_vars:,} genes")
        
        # Απλό κουμπί εκτέλεσης
        if st.button("🚀 Εκτέλεση DEG Analysis", type="primary"):
            try:
                # Debug info
                st.info(f"🔍 Debug: Comparison groups = {self.comparison_groups}")
                st.info(f"🔍 Debug: Dataset shape = {self.adata.shape}")
                
                with st.spinner("Εκτέλεση DEG analysis..."):
                    self.deg_results = self.perform_deg_analysis()
                
                st.info(f"🔍 Debug: Results type = {type(self.deg_results)}")
                if self.deg_results is not None:
                    st.info(f"🔍 Debug: Results shape = {self.deg_results.shape if hasattr(self.deg_results, 'shape') else 'No shape'}")
                
                if self.deg_results is not None and not self.deg_results.empty:
                    st.success("✅ DEG Analysis ολοκληρώθηκε επιτυχώς!")
                    
                    # Αποθήκευση στο session state
                    st.session_state['deg_results'] = self.deg_results
                    
                    # Εμφάνιση περίληψης
                    self.display_analysis_summary()
                else:
                    st.error("❌ Η ανάλυση δεν παρήγαγε αποτελέσματα")
                    if self.deg_results is None:
                        st.error("DEG results is None")
                    elif self.deg_results.empty:
                        st.error("DEG results is empty DataFrame")
                
            except Exception as e:
                st.error(f"❌ Σφάλμα κατά την ανάλυση: {str(e)}")
                import traceback
                st.error(f"Full traceback: {traceback.format_exc()}")
    
    def _run_production_analysis(self):
        """Production analysis - πλήρη ακρίβεια, όλα τα δεδομένα"""
        
        st.info("🏭 PRODUCTION MODE: Επεξεργασία όλων των δεδομένων για μέγιστη ακρίβεια")
        
        try:
            with st.spinner("Production analysis σε εξέλιξη..."):
                
                # Production recommendations
                recommendations = production_engine.get_production_recommendations(self.adata, self.comparison_column)
                
                with st.expander("📋 Production Analysis Plan"):
                    st.json(recommendations)
                
                # Εκτέλεση production analysis
                if self.comparison_groups.get('type') == 'pairwise':
                    group1 = self.comparison_groups.get('group1')
                    group2 = self.comparison_groups.get('group2')
                    
                    if group1 and group2:
                        self.deg_results = production_engine.production_deg_analysis(
                            self.adata,
                            self.comparison_column,
                            group1,
                            group2,
                            self.selected_method,
                            self.correction_method
                        )
                else:
                    st.warning("⚠️ Production one-vs-all δεν είναι ακόμα διαθέσιμο - χρήση standard method")
                    self.deg_results = self.perform_deg_analysis()
                
                if self.deg_results is not None and not self.deg_results.empty:
                    st.success("✅ Production Analysis ολοκληρώθηκε με πλήρη ακρίβεια!")
                    
                    # Production quality metrics
                    if 'analysis_metadata' in self.deg_results:
                        attrs = self.deg_results['analysis_metadata']
                        with st.expander("🏆 Production Quality Metrics"):
                            col1, col2, col3 = st.columns(3)
                            with col1:
                                st.metric("Cells Analyzed", f"{attrs.get('total_cells_analyzed', 0):,}")
                            with col2:
                                st.metric("Genes Analyzed", f"{attrs.get('total_genes_analyzed', 0):,}")
                            with col3:
                                st.metric("Data Integrity", attrs.get('data_integrity', 'unknown'))
                            
                            st.info(f"Analysis Time: {attrs['elapsed_time']/60:.1f} minutes")
                            st.info(f"Method: {attrs['method']} with {attrs['correction_method']} correction")
                    
                    # Αποθήκευση
                    st.session_state['deg_results'] = self.deg_results
                    self.display_analysis_summary()
                    
        except Exception as e:
            st.error(f"❌ Production analysis error: {str(e)}")
    
    def _run_fast_analysis(self):
        """Fast analysis - βελτιστοποιημένη ταχύτητα με διατήρηση λειτουργικότητας"""
        
        st.info("⚡ FAST MODE: Βελτιστοποιημένη ανάλυση με smart algorithms")
        
        try:
            with st.spinner("Εκτέλεση βελτιστοποιημένης ανάλυσης..."):
                
                # Χρήση optimized engine για fast mode
                if self.comparison_groups.get('type') == 'pairwise':
                    group1 = self.comparison_groups.get('group1')
                    group2 = self.comparison_groups.get('group2')
                    
                    if group1 and group2:
                        # Use standard analysis
                        self.deg_results = self.perform_deg_analysis()
                else:
                    # One vs all με optimized engine
                    reference_group = self.comparison_groups.get('reference_group')
                    if reference_group:
                        # Use standard analysis
                        self.deg_results = self.perform_deg_analysis()
                
                if self.deg_results is not None and not self.deg_results.empty:
                    # Multiple testing correction
                    if self.correction_method != 'none' and 'analysis_metadata' not in self.deg_results:
                        corrected_pvals = self.apply_multiple_testing_correction(
                            self.deg_results['pvals'].values
                        )
                        self.deg_results['pvals_adj'] = corrected_pvals
                        
                        # Significance calculation
                        self.deg_results['significant'] = (
                            (self.deg_results['pvals_adj'] < 0.05) & 
                            (np.abs(self.deg_results['logfoldchanges']) > self.logfc_threshold)
                        )
                    
                    st.success("✅ Fast Analysis ολοκληρώθηκε!")
                    
                    # Performance metrics
                    if isinstance(self.deg_results, dict) and 'analysis_metadata' in self.deg_results:
                        metadata = self.deg_results['analysis_metadata']
                        with st.expander("📊 Performance Metrics"):
                            col1, col2, col3 = st.columns(3)
                            with col1:
                                st.metric("Χρόνος Εκτέλεσης", f"{metadata['elapsed_time']:.1f}s")
                            with col2:
                                st.metric("Κύτταρα Analyzed", f"{metadata['analysis_cells']:,}")
                            with col3:
                                st.metric("Γονίδια Analyzed", f"{metadata['analysis_genes']:,}")
                    
                    # Αποθήκευση
                    st.session_state['deg_results'] = self.deg_results
                    self.display_analysis_summary()
                
        except Exception as e:
            st.error(f"❌ Fast analysis error: {str(e)}")
            st.info("🔄 Fallback σε standard method...")
            try:
                self.deg_results = self.perform_deg_analysis()
                if self.deg_results is not None:
                    st.session_state['deg_results'] = self.deg_results
                    self.display_analysis_summary()
            except Exception as e2:
                st.error(f"❌ Fallback method απέτυχε: {str(e2)}")
    
    def estimate_analysis_time(self):
        """Ρεαλιστική εκτίμηση χρόνου ανάλυσης"""
        
        # Εκτίμηση βάσει του τι θα αναλυθεί πραγματικά
        if self.adata.n_obs > 50000:
            # Θα γίνει subsampling
            estimated_cells = min(20000, self.adata.n_obs)
        else:
            estimated_cells = self.adata.n_obs
        
        n_genes = min(self.max_genes, self.adata.n_vars)
        
        # ΔΙΟΡΘΩΜΕΝΗ εκτίμηση - το πρόβλημα ήταν στον υπολογισμό, όχι στην ανάλυση
        base_time = (estimated_cells * n_genes) / 100_000_000  # 100M operations per second (πιο realistic)
        
        # Πολλαπλασιαστής βάσει μεθόδου
        method_multiplier = {
            'wilcoxon': 1.0,
            't-test': 0.8,
            'logreg': 2.0
        }
        
        return base_time * method_multiplier.get(self.selected_method, 1.0)
    
    def perform_deg_analysis(self):
        """Εκτέλεση της πραγματικής DEG ανάλυσης"""
        
        # Προετοιμασία δεδομένων
        adata_analysis = self.prepare_data_for_analysis()
        
        if adata_analysis is None:
            return None
        
        # Εκτέλεση ανάλυσης βάσει τύπου σύγκρισης
        if self.comparison_groups['type'] == 'one_vs_all':
            return self.perform_one_vs_all_analysis(adata_analysis)
        else:
            return self.perform_pairwise_analysis(adata_analysis)
    
    def prepare_data_for_analysis(self):
        """Προετοιμασία δεδομένων για ανάλυση με intelligent subsampling"""
        
        # ΚΡΙΣΙΜΗ ΒΕΛΤΙΩΣΗ: Intelligent subsampling για μεγάλα datasets
        if self.adata.n_obs > 50000:  # >50K κύτταρα
            st.warning(f"⚠️ Πολύ μεγάλο dataset ({self.adata.n_obs:,} κύτταρα)")
            
            # Στρατηγικό subsampling που διατηρεί την ποικιλία
            max_cells_for_deg = 20000  # Reasonable για DEG analysis
            
            if self.comparison_column in self.adata.obs.columns:
                # Stratified sampling - διατηρεί αναλογίες ομάδων
                st.info(f"🎯 Stratified sampling: {max_cells_for_deg:,} κύτταρα διατηρώντας αναλογίες ομάδων")
                
                groups = self.adata.obs[self.comparison_column].unique()
                sampled_indices = []
                
                for group in groups:
                    group_mask = self.adata.obs[self.comparison_column] == group
                    group_indices = np.where(group_mask)[0]
                    
                    # Proportional sampling
                    n_group_cells = len(group_indices)
                    n_sample_group = int((n_group_cells / self.adata.n_obs) * max_cells_for_deg)
                    n_sample_group = max(100, min(n_sample_group, n_group_cells))  # Min 100, max available
                    
                    if n_sample_group < n_group_cells:
                        selected_indices = np.random.choice(group_indices, n_sample_group, replace=False)
                    else:
                        selected_indices = group_indices
                    
                    sampled_indices.extend(selected_indices)
                    st.info(f"  - {group}: {len(selected_indices):,} κύτταρα")
                
                # Δημιουργία subsampled dataset
                adata_analysis = self.adata[sampled_indices, :].copy()
                st.success(f"✅ Stratified sampling: {adata_analysis.n_obs:,} κύτταρα, {adata_analysis.n_vars:,} γονίδια")
            else:
                # Random sampling αν δεν υπάρχουν ομάδες
                st.info(f"🎲 Random sampling: {max_cells_for_deg:,} κύτταρα")
                indices = np.random.choice(self.adata.n_obs, max_cells_for_deg, replace=False)
                adata_analysis = self.adata[indices, :].copy()
        else:
            # Μικρό dataset - χρήση όλων των κυττάρων
            adata_analysis = self.adata.copy()
            st.info(f"📊 Χρήση όλων των {adata_analysis.n_obs:,} κυττάρων")
        
        # Φιλτράρισμα για highly variable genes αν επιλέχθηκε
        if self.filter_highly_variable:
            if 'highly_variable' not in adata_analysis.var.columns:
                st.info("🔍 Υπολογισμός highly variable genes...")
                sc.pp.highly_variable_genes(adata_analysis, min_mean=0.0125, max_mean=3, min_disp=0.5)
            
            adata_analysis = adata_analysis[:, adata_analysis.var.highly_variable]
        
        # PRODUCTION MODE: Χρήση όλων των γονιδίων εκτός αν ο χρήστης επιλέξει διαφορετικά
        if adata_analysis.n_vars > self.max_genes and self.max_genes < 50000:
            st.warning(f"⚠️ Περιορισμός από {adata_analysis.n_vars:,} σε {self.max_genes:,} γονίδια")
            
            # Επιλογή των πιο μεταβλητών γονιδίων
            if 'highly_variable_rank' in adata_analysis.var.columns:
                top_genes = adata_analysis.var.nsmallest(self.max_genes, 'highly_variable_rank').index
                st.info("🧬 Χρήση top highly variable genes")
            elif 'highly_variable' in adata_analysis.var.columns:
                hvg_genes = adata_analysis.var_names[adata_analysis.var['highly_variable']]
                if len(hvg_genes) > self.max_genes:
                    top_genes = hvg_genes[:self.max_genes]
                else:
                    top_genes = hvg_genes
                st.info(f"🧬 Χρήση {len(top_genes)} highly variable genes")
            else:
                # Fallback: τυχαία επιλογή
                top_genes = np.random.choice(adata_analysis.var_names, self.max_genes, replace=False)
                st.warning("🎲 Τυχαία επιλογή γονιδίων (δεν υπάρχουν HVG rankings)")
            
            adata_analysis = adata_analysis[:, top_genes]
        else:
            st.success(f"✅ Ανάλυση όλων των {adata_analysis.n_vars:,} γονιδίων")
        
        st.info(f"📊 Ανάλυση με {adata_analysis.n_obs:,} κύτταρα και {adata_analysis.n_vars:,} γονίδια")
        
        return adata_analysis
    
    def perform_one_vs_all_analysis(self, adata_analysis):
        """One vs All ανάλυση"""
        
        reference_group = self.comparison_groups['reference_group']
        
        st.info(f"🔬 Εκτέλεση One vs All: '{reference_group}' vs άλλες ομάδες")
        
        # ΚΡΙΣΙΜΗ ΔΙΟΡΘΩΣΗ: Proper data preparation
        st.info("🔄 Προετοιμασία δεδομένων για DEG analysis...")
        
        # Έλεγχος αν τα δεδομένα είναι κατάλληλα
        if sparse.issparse(adata_analysis.X):
            st.info("📊 Sparse matrix detected - converting for analysis...")
            # Ensure proper sparse format
            adata_analysis.X = adata_analysis.X.tocsr()
        
        # ΠΑΝΤΑ κάνουμε normalization για σωστά αποτελέσματα
        st.info("🔧 Normalization για σωστή στατιστική ανάλυση...")
        sc.pp.normalize_total(adata_analysis, target_sum=1e4)
        sc.pp.log1p(adata_analysis)
        
        # Validation των δεδομένων
        if np.any(np.isnan(adata_analysis.X.data if sparse.issparse(adata_analysis.X) else adata_analysis.X)):
            st.error("❌ NaN values detected in data!")
            return pd.DataFrame()
        
        st.success("✅ Δεδομένα έτοιμα για ανάλυση")
        
        # Scanpy DEG analysis
        sc.tl.rank_genes_groups(
            adata_analysis,
            groupby=self.comparison_column,
            groups=[reference_group],
            reference='rest',
            method=self.selected_method,
            use_raw=self.use_raw,
            tie_correct=self.tie_correct if self.selected_method == 'wilcoxon' else True
        )
        
        # Μετατροπή αποτελεσμάτων σε DataFrame
        result_df = sc.get.rank_genes_groups_df(adata_analysis, group=reference_group)
        
        # Validation και safe handling των p-values
        result_df['pvals'] = result_df['pvals'].fillna(1.0)  # NaN → 1.0
        result_df['pvals'] = np.clip(result_df['pvals'], 1e-300, 1.0)  # Clip extreme values
        
        # Multiple testing correction
        if self.correction_method != 'none':
            corrected_pvals = self.apply_multiple_testing_correction(result_df['pvals'].values)
            result_df['pvals_adj'] = corrected_pvals
        else:
            result_df['pvals_adj'] = result_df['pvals']
        
        # Validation των corrected p-values
        result_df['pvals_adj'] = result_df['pvals_adj'].fillna(1.0)
        result_df['pvals_adj'] = np.clip(result_df['pvals_adj'], 1e-300, 1.0)
        
        # Προσθήκη significance flag
        result_df['significant'] = (
            (result_df['pvals_adj'] < 0.05) & 
            (np.abs(result_df['logfoldchanges']) > self.logfc_threshold)
        )
        
        # Ταξινόμηση βάσει significance
        result_df = result_df.sort_values(['significant', 'pvals_adj'], ascending=[False, True])
        
        return result_df
    
    def perform_pairwise_analysis(self, adata_analysis):
        """Pairwise ανάλυση"""
        
        group1 = self.comparison_groups.get('group1')
        group2 = self.comparison_groups.get('group2')
        
        if not group1 or not group2:
            st.error("❌ Δεν έχουν επιλεγεί και οι δύο ομάδες για σύγκριση")
            return pd.DataFrame()
        
        st.info(f"🔬 Εκτέλεση Pairwise: '{group1}' vs '{group2}'")
        
        # Φιλτράρισμα για τις δύο ομάδες
        mask = adata_analysis.obs[self.comparison_column].isin([group1, group2])
        adata_subset = adata_analysis[mask, :].copy()
        
        # Προετοιμασία δεδομένων για DEG analysis
        if not self.use_raw:
            # Log transformation αν δεν υπάρχει ήδη
            if 'log1p' not in adata_subset.uns:
                st.info("🔄 Log normalization για pairwise DEG analysis...")
                sc.pp.normalize_total(adata_subset, target_sum=1e4)
                sc.pp.log1p(adata_subset)
        
        # Scanpy DEG analysis
        sc.tl.rank_genes_groups(
            adata_subset,
            groupby=self.comparison_column,
            groups=[group1],
            reference=group2,
            method=self.selected_method,
            use_raw=self.use_raw,
            tie_correct=self.tie_correct if self.selected_method == 'wilcoxon' else True
        )
        
        # Μετατροπή αποτελεσμάτων σε DataFrame
        result_df = sc.get.rank_genes_groups_df(adata_subset, group=group1)
        
        # Validation και safe handling των p-values
        result_df['pvals'] = result_df['pvals'].fillna(1.0)  # NaN → 1.0
        result_df['pvals'] = np.clip(result_df['pvals'], 1e-300, 1.0)  # Clip extreme values
        
        # Multiple testing correction
        if self.correction_method != 'none':
            corrected_pvals = self.apply_multiple_testing_correction(result_df['pvals'].values)
            result_df['pvals_adj'] = corrected_pvals
        else:
            result_df['pvals_adj'] = result_df['pvals']
        
        # Validation των corrected p-values
        result_df['pvals_adj'] = result_df['pvals_adj'].fillna(1.0)
        result_df['pvals_adj'] = np.clip(result_df['pvals_adj'], 1e-300, 1.0)
        
        # Προσθήκη significance flag
        result_df['significant'] = (
            (result_df['pvals_adj'] < 0.05) & 
            (np.abs(result_df['logfoldchanges']) > self.logfc_threshold)
        )
        
        # Ταξινόμηση βάσει significance
        result_df = result_df.sort_values(['significant', 'pvals_adj'], ascending=[False, True])
        
        return result_df
    
    def apply_multiple_testing_correction(self, pvalues):
        """Εφαρμογή multiple testing correction"""
        
        if self.correction_method == 'benjamini-hochberg':
            rejected, pvals_corrected, _, _ = multipletests(pvalues, method='fdr_bh')
        elif self.correction_method == 'bonferroni':
            rejected, pvals_corrected, _, _ = multipletests(pvalues, method='bonferroni')
        else:
            pvals_corrected = pvalues
        
        return pvals_corrected
    
    def display_analysis_summary(self):
        """Εμφάνιση περίληψης ανάλυσης"""
        
        st.subheader("📊 Περίληψη Αποτελεσμάτων")
        
        if self.deg_results is not None:
            total_genes = len(self.deg_results)
            significant_genes = self.deg_results['significant'].sum()
            upregulated = ((self.deg_results['significant']) & 
                          (self.deg_results['logfoldchanges'] > 0)).sum()
            downregulated = ((self.deg_results['significant']) & 
                            (self.deg_results['logfoldchanges'] < 0)).sum()
            
            # DEBUG: Ας δούμε γιατί δεν βρίσκει σημαντικά γονίδια
            st.subheader("🔍 Debug Analysis")
            
            # P-value distribution
            pval_stats = {
                'Min p-value': self.deg_results['pvals_adj'].min(),
                'Max p-value': self.deg_results['pvals_adj'].max(),
                'Mean p-value': self.deg_results['pvals_adj'].mean(),
                'P-values < 0.05': (self.deg_results['pvals_adj'] < 0.05).sum(),
                'P-values < 0.01': (self.deg_results['pvals_adj'] < 0.01).sum()
            }
            
            # LogFC distribution  
            logfc_stats = {
                'Min LogFC': self.deg_results['logfoldchanges'].min(),
                'Max LogFC': self.deg_results['logfoldchanges'].max(),
                'Mean LogFC': self.deg_results['logfoldchanges'].mean(),
                'LogFC > 0.25': (np.abs(self.deg_results['logfoldchanges']) > 0.25).sum(),
                'LogFC > 0.1': (np.abs(self.deg_results['logfoldchanges']) > 0.1).sum()
            }
            
            col1, col2 = st.columns(2)
            with col1:
                st.write("**P-value Stats:**")
                for key, value in pval_stats.items():
                    st.write(f"- {key}: {value}")
            
            with col2:
                st.write("**LogFC Stats:**")
                for key, value in logfc_stats.items():
                    st.write(f"- {key}: {value:.3f}" if isinstance(value, float) else f"- {key}: {value}")
            
            # Significance criteria debug
            pval_criteria = (self.deg_results['pvals_adj'] < 0.05).sum()
            logfc_criteria = (np.abs(self.deg_results['logfoldchanges']) > self.logfc_threshold).sum()
            both_criteria = ((self.deg_results['pvals_adj'] < 0.05) & 
                           (np.abs(self.deg_results['logfoldchanges']) > self.logfc_threshold)).sum()
            
            st.write(f"**Significance Criteria Debug:**")
            st.write(f"- Genes with p-adj < 0.05: {pval_criteria:,}")
            st.write(f"- Genes with |LogFC| > {self.logfc_threshold}: {logfc_criteria:,}")
            st.write(f"- Genes meeting BOTH criteria: {both_criteria:,}")
            
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("Συνολικά Γονίδια", f"{total_genes:,}")
            with col2:
                st.metric("Σημαντικά Γονίδια", f"{significant_genes:,}")
            with col3:
                st.metric("Upregulated", f"{upregulated:,}")
            with col4:
                st.metric("Downregulated", f"{downregulated:,}")
    
    def render_results_visualization(self):
        """Οπτικοποίηση αποτελεσμάτων"""
        
        st.header("📈 Αποτελέσματα & Οπτικοποίηση")
        
        if self.deg_results is None:
            if 'deg_results' in st.session_state:
                self.deg_results = st.session_state['deg_results']
            else:
                st.warning("⚠️ Εκτελέστε πρώτα DEG analysis")
                return
        
        # Tabs για διαφορετικές οπτικοποιήσεις
        viz_tab1, viz_tab2, viz_tab3, viz_tab4 = st.tabs([
            "📊 Results Table", 
            "🌋 Volcano Plot", 
            "🔥 Heatmap",
            "💾 Export Results"
        ])
        
        with viz_tab1:
            self.render_results_table()
        
        with viz_tab2:
            self.render_volcano_plot()
        
        with viz_tab3:
            self.render_heatmap()
        
        with viz_tab4:
            self.render_export_options()
    
    def render_results_table(self):
        """Εμφάνιση πίνακα αποτελεσμάτων"""
        
        st.subheader("📊 Πίνακας Αποτελεσμάτων")
        
        # Φιλτράρισμα αποτελεσμάτων
        col1, col2, col3 = st.columns(3)
        
        with col1:
            show_only_significant = st.checkbox("Μόνο σημαντικά γονίδια", value=True)
        
        with col2:
            max_pvalue = st.slider("Max p-value:", 0.001, 0.1, 0.05, step=0.001)
        
        with col3:
            min_logfc = st.slider("Min |log FC|:", 0.0, 2.0, 0.25, step=0.05)
        
        # Φιλτράρισμα
        filtered_results = self.deg_results.copy()
        
        if show_only_significant:
            filtered_results = filtered_results[filtered_results['significant']]
        
        filtered_results = filtered_results[
            (filtered_results['pvals_adj'] <= max_pvalue) &
            (np.abs(filtered_results['logfoldchanges']) >= min_logfc)
        ]
        
        st.write(f"**Εμφάνιση {len(filtered_results):,} από {len(self.deg_results):,} γονίδια**")
        
        # Εμφάνιση πίνακα
        if not filtered_results.empty:
            # Μορφοποίηση για καλύτερη εμφάνιση
            display_df = filtered_results.copy()
            display_df['pvals'] = display_df['pvals'].apply(lambda x: f"{x:.2e}")
            display_df['pvals_adj'] = display_df['pvals_adj'].apply(lambda x: f"{x:.2e}")
            display_df['logfoldchanges'] = display_df['logfoldchanges'].round(3)
            
            st.dataframe(
                display_df[['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'significant']],
                use_container_width=True,
                height=400
            )
        else:
            st.info("Δεν βρέθηκαν γονίδια που να πληρούν τα κριτήρια φιλτραρίσματος")
    
    def render_volcano_plot(self):
        """Δημιουργία volcano plot"""
        
        st.subheader("🌋 Volcano Plot")
        
        if self.deg_results is None or self.deg_results.empty:
            st.warning("Δεν υπάρχουν αποτελέσματα για οπτικοποίηση")
            return
        
        # Παράμετροι plot
        col1, col2 = st.columns(2)
        
        with col1:
            pvalue_threshold = st.slider(
                "P-value threshold:", 
                0.001, 0.1, 0.05, step=0.001,
                key="volcano_pval"
            )
        
        with col2:
            logfc_threshold = st.slider(
                "Log FC threshold:", 
                0.0, 2.0, 0.25, step=0.05,
                key="volcano_logfc"
            )
        
        if st.button("🎨 Δημιουργία Volcano Plot"):
            
            # Προετοιμασία δεδομένων για plot
            plot_data = self.deg_results.copy()
            
            # Ασφαλής υπολογισμός -log10 με αποφυγή division by zero
            pvals_adj_safe = plot_data['pvals_adj'].copy()
            
            # Αντικατάσταση 0 values με πολύ μικρή τιμή
            pvals_adj_safe = np.where(pvals_adj_safe <= 0, 1e-300, pvals_adj_safe)
            
            # Clip extreme values για αποφυγή inf
            pvals_adj_safe = np.clip(pvals_adj_safe, 1e-300, 1.0)
            
            plot_data['-log10(p_adj)'] = -np.log10(pvals_adj_safe)
            
            # Κατηγοριοποίηση γονιδίων
            plot_data['category'] = 'Non-significant'
            plot_data.loc[
                (plot_data['pvals_adj'] < pvalue_threshold) & 
                (plot_data['logfoldchanges'] > logfc_threshold), 
                'category'
            ] = 'Upregulated'
            plot_data.loc[
                (plot_data['pvals_adj'] < pvalue_threshold) & 
                (plot_data['logfoldchanges'] < -logfc_threshold), 
                'category'
            ] = 'Downregulated'
            
            # Δημιουργία plotly figure
            fig = px.scatter(
                plot_data,
                x='logfoldchanges',
                y='-log10(p_adj)',
                color='category',
                hover_data=['names'],
                color_discrete_map={
                    'Non-significant': 'lightgray',
                    'Upregulated': 'red',
                    'Downregulated': 'blue'
                },
                title='Volcano Plot - Differential Gene Expression',
                labels={
                    'logfoldchanges': 'Log2 Fold Change',
                    '-log10(p_adj)': '-log10(Adjusted P-value)'
                }
            )
            
            # Προσθήκη threshold lines (με ασφαλή log calculation)
            safe_pvalue_threshold = max(pvalue_threshold, 1e-300)
            fig.add_hline(y=-np.log10(safe_pvalue_threshold), line_dash="dash", line_color="black")
            fig.add_vline(x=logfc_threshold, line_dash="dash", line_color="black")
            fig.add_vline(x=-logfc_threshold, line_dash="dash", line_color="black")
            
            fig.update_layout(height=600)
            st.plotly_chart(fig, use_container_width=True)
            
            # Στατιστικές
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Non-significant", (plot_data['category'] == 'Non-significant').sum())
            with col2:
                st.metric("Upregulated", (plot_data['category'] == 'Upregulated').sum())
            with col3:
                st.metric("Downregulated", (plot_data['category'] == 'Downregulated').sum())
    
    def render_heatmap(self):
        """Δημιουργία heatmap με top genes"""
        
        st.subheader("🔥 Heatmap - Top Differentially Expressed Genes")
        
        if self.adata is None or self.deg_results is None:
            st.warning("Δεν υπάρχουν δεδομένα για heatmap")
            return
        
        # Παράμετροι heatmap
        col1, col2 = st.columns(2)
        
        with col1:
            n_top_genes = st.slider(
                "Αριθμός top genes:", 
                5, 100, 20, step=5
            )
        
        with col2:
            split_up_down = st.checkbox(
                "Χωριστά up/down regulated", 
                value=True
            )
        
        if st.button("🎨 Δημιουργία Heatmap"):
            try:
                # Επιλογή top genes
                significant_genes = self.deg_results[self.deg_results['significant']].copy()
                
                if significant_genes.empty:
                    st.warning("Δεν βρέθηκαν σημαντικά γονίδια για heatmap")
                    return
                
                if split_up_down:
                    # Χωριστά up και down regulated
                    up_genes = significant_genes[
                        significant_genes['logfoldchanges'] > 0
                    ].nsmallest(n_top_genes//2, 'pvals_adj')
                    
                    down_genes = significant_genes[
                        significant_genes['logfoldchanges'] < 0
                    ].nsmallest(n_top_genes//2, 'pvals_adj')
                    
                    selected_genes = pd.concat([up_genes, down_genes])
                else:
                    # Top genes συνολικά
                    selected_genes = significant_genes.nsmallest(n_top_genes, 'pvals_adj')
                
                gene_names = selected_genes['names'].tolist()
                
                # Εξαγωγή expression data
                gene_mask = self.adata.var_names.isin(gene_names)
                expression_data = self.adata[:, gene_mask].X
                
                if hasattr(expression_data, 'toarray'):
                    expression_data = expression_data.toarray()
                
                # Δημιουργία DataFrame για heatmap
                heatmap_df = pd.DataFrame(
                    expression_data.T,
                    index=self.adata.var_names[gene_mask],
                    columns=self.adata.obs_names
                )
                
                # Subsampling αν πολλά κύτταρα
                if heatmap_df.shape[1] > 500:
                    sample_cells = np.random.choice(
                        heatmap_df.columns, 
                        500, 
                        replace=False
                    )
                    heatmap_df = heatmap_df[sample_cells]
                
                # Δημιουργία heatmap με matplotlib
                fig, ax = plt.subplots(figsize=(12, max(6, len(gene_names) * 0.3)))
                
                sns.heatmap(
                    heatmap_df,
                    cmap='RdBu_r',
                    center=0,
                    ax=ax,
                    cbar_kws={'label': 'Expression Level'}
                )
                
                ax.set_title(f'Top {len(gene_names)} Differentially Expressed Genes')
                ax.set_xlabel('Cells (sample)')
                ax.set_ylabel('Genes')
                
                plt.tight_layout()
                st.pyplot(fig)
                
            except Exception as e:
                st.error(f"❌ Σφάλμα στη δημιουργία heatmap: {str(e)}")
    
    def render_export_options(self):
        """Επιλογές export αποτελεσμάτων"""
        
        st.subheader("💾 Export Αποτελεσμάτων")
        
        if self.deg_results is None:
            st.warning("Δεν υπάρχουν αποτελέσματα για export")
            return
        
        # Επιλογές export
        export_format = st.selectbox(
            "Μορφή export:",
            ["CSV", "Excel", "TSV"]
        )
        
        export_options = st.multiselect(
            "Τι να συμπεριληφθεί:",
            ["Όλα τα γονίδια", "Μόνο σημαντικά", "Top 100", "Top 500"],
            default=["Όλα τα γονίδια"]
        )
        
        if st.button("📥 Download Results"):
            try:
                # Προετοιμασία δεδομένων για export
                export_data = self.deg_results.copy()
                
                if "Μόνο σημαντικά" in export_options:
                    export_data = export_data[export_data['significant']]
                elif "Top 100" in export_options:
                    export_data = export_data.nsmallest(100, 'pvals_adj')
                elif "Top 500" in export_options:
                    export_data = export_data.nsmallest(500, 'pvals_adj')
                
                # Export βάσει μορφής
                if export_format == "CSV":
                    csv_buffer = io.StringIO()
                    export_data.to_csv(csv_buffer, index=False)
                    
                    st.download_button(
                        label="💾 Download CSV",
                        data=csv_buffer.getvalue(),
                        file_name="deg_results.csv",
                        mime="text/csv"
                    )
                
                elif export_format == "Excel":
                    excel_buffer = io.BytesIO()
                    with pd.ExcelWriter(excel_buffer, engine='openpyxl') as writer:
                        export_data.to_excel(writer, sheet_name='DEG_Results', index=False)
                        
                        # Προσθήκη summary sheet
                        summary_data = {
                            'Metric': ['Total Genes', 'Significant Genes', 'Upregulated', 'Downregulated'],
                            'Count': [
                                len(self.deg_results),
                                self.deg_results['significant'].sum(),
                                ((self.deg_results['significant']) & (self.deg_results['logfoldchanges'] > 0)).sum(),
                                ((self.deg_results['significant']) & (self.deg_results['logfoldchanges'] < 0)).sum()
                            ]
                        }
                        pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)
                    
                    st.download_button(
                        label="💾 Download Excel",
                        data=excel_buffer.getvalue(),
                        file_name="deg_results.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
                
                else:  # TSV
                    tsv_buffer = io.StringIO()
                    export_data.to_csv(tsv_buffer, sep='\t', index=False)
                    
                    st.download_button(
                        label="💾 Download TSV",
                        data=tsv_buffer.getvalue(),
                        file_name="deg_results.tsv",
                        mime="text/tab-separated-values"
                    )
                
            except Exception as e:
                st.error(f"❌ Σφάλμα κατά το export: {str(e)}")
