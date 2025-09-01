"""
Ολοκληρωμένο Module για την Ενοποίηση Δεδομένων με Scanorama

Αυτό το module περιλαμβάνει όλες τις λειτουργίες που απαιτούνται για την
ενοποίηση πολλαπλών scRNA-seq datasets, συμπεριλαμβανομένου του batch correction
με Scanorama και την αξιολόγηση της ποιότητας της integration.


Ημερομηνία: 2025
"""

import streamlit as st
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import io
import tempfile
import os
import sys
from pathlib import Path

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

# Scanorama import με error handling
try:
    import scanorama
    SCANORAMA_AVAILABLE = True
except ImportError:
    SCANORAMA_AVAILABLE = False

class DataIntegrationPageComplete:
    """Ολοκληρωμένη κλάση για τη σελίδα ενοποίησης δεδομένων"""
    
    def __init__(self):
        self.datasets = {}
        self.integrated_data = None
        self.integration_method = "simple_concatenation"
        
    def render(self):
        """Κεντρική μέθοδος για την εμφάνιση της σελίδας"""
        
        st.title("🔗 Ενοποίηση Δεδομένων")
        st.markdown("### Batch correction και ενοποίηση πολλαπλών scRNA-seq datasets")
        
        if not SCANORAMA_AVAILABLE:
            st.warning("⚠️ Scanorama δεν είναι εγκατεστημένο. Χρήση εναλλακτικής μεθόδου.")
        
        # Memory information
        display_memory_info()
        
        # Tabs για διαφορετικές λειτουργίες
        tab1, tab2, tab3 = st.tabs([
            "📁 Φόρτωση & Παράμετροι", 
            "🔄 Εκτέλεση Integration",
            "📊 Αποτελέσματα"
        ])
        
        with tab1:
            self.render_setup()
            
        with tab2:
            self.render_execution()
            
        with tab3:
            self.render_results()
    
    def render_setup(self):
        """Φόρτωση datasets και ρύθμιση παραμετρών"""
        
        st.header("📁 Φόρτωση Datasets")
        
        # Επιλογή μεθόδου φόρτωσης
        upload_method = st.radio(
            "Επιλέξτε μέθοδο φόρτωσης:",
            ["Χρήση Προεπεξεργασμένων Δεδομένων", "Ανέβασμα Αρχείων"],
            horizontal=True
        )
        
        if upload_method == "Χρήση Προεπεξεργασμένων Δεδομένων":
            # Έλεγχος για δεδομένα στο session state
            if 'adata' in st.session_state and st.session_state.adata is not None:
                st.info("📊 Βρέθηκαν προεπεξεργασμένα δεδομένα από το Preprocessing module")
                
                if st.button("Χρήση Προεπεξεργασμένων Δεδομένων"):
                    self.datasets['preprocessed_data'] = st.session_state.adata.copy()
                    st.success("✅ Προεπεξεργασμένα δεδομένα προστέθηκαν!")
            else:
                st.warning("⚠️ Δεν βρέθηκαν προεπεξεργασμένα δεδομένα.")
        
        else:
            # Multiple file uploader
            uploaded_files = st.file_uploader(
                "Επιλέξτε H5AD αρχεία:",
                type=['h5ad'],
                accept_multiple_files=True
            )
            
            if uploaded_files:
                for uploaded_file in uploaded_files:
                    dataset_name = uploaded_file.name.replace('.h5ad', '')
                    
                    if dataset_name not in self.datasets:
                        try:
                            with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp_file:
                                tmp_file.write(uploaded_file.getvalue())
                                tmp_path = tmp_file.name
                            
                            with st.spinner(f"Φόρτωση {dataset_name}..."):
                                adata = sc.read_h5ad(tmp_path)
                                adata = memory_manager.optimize_adata_for_memory(adata)
                                self.datasets[dataset_name] = adata
                            
                            os.unlink(tmp_path)
                            st.success(f"✅ Dataset '{dataset_name}' φορτώθηκε!")
                            
                        except Exception as e:
                            st.error(f"❌ Σφάλμα φόρτωσης {dataset_name}: {str(e)}")
        
        # Εμφάνιση φορτωμένων datasets
        if self.datasets:
            st.subheader("📊 Φορτωμένα Datasets")
            
            for dataset_name, adata in self.datasets.items():
                with st.expander(f"Dataset: {dataset_name}"):
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.metric("Κύτταρα", f"{adata.n_obs:,}")
                    with col2:
                        st.metric("Γονίδια", f"{adata.n_vars:,}")
                    with col3:
                        memory_size = memory_manager.estimate_memory_usage(adata)
                        st.metric("Μέγεθος", f"{memory_size:.1f} MB")
            
            # Παράμετροι Integration
            st.subheader("⚙️ Παράμετροι Integration")
            
            col1, col2 = st.columns(2)
            
            with col1:
                available_methods = ["simple_concatenation"]
                if SCANORAMA_AVAILABLE and len(self.datasets) > 1:
                    available_methods.insert(0, "scanorama")
                
                self.integration_method = st.selectbox(
                    "Μέθοδος Integration:",
                    available_methods
                )
            
            with col2:
                self.intersection_only = st.checkbox(
                    "Χρήση μόνο κοινών γονιδίων",
                    value=True,
                    help="Διατήρηση μόνο των γονιδίων που υπάρχουν σε όλα τα datasets"
                )
    
    @memory_monitor
    def render_execution(self):
        """Εκτέλεση integration"""
        
        st.header("🔄 Εκτέλεση Integration")
        
        if not self.datasets:
            st.warning("⚠️ Φορτώστε πρώτα datasets")
            return
        
        if len(self.datasets) < 1:
            st.warning("⚠️ Χρειάζεται τουλάχιστον 1 dataset")
            return
        
        # Προεπισκόπηση
        st.subheader("📋 Προεπισκόπηση")
        
        total_cells = sum(adata.n_obs for adata in self.datasets.values())
        common_genes = self.get_common_genes()
        total_genes = len(common_genes) if self.intersection_only and common_genes else max(adata.n_vars for adata in self.datasets.values())
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Συνολικά Κύτταρα", f"{total_cells:,}")
        with col2:
            st.metric("Γονίδια", f"{total_genes:,}")
        with col3:
            estimated_memory = sum(memory_manager.estimate_memory_usage(adata) for adata in self.datasets.values()) * 1.5
            st.metric("Εκτιμώμενη Μνήμη", f"{estimated_memory:.1f} MB")
        
        # Κουμπί εκτέλεσης
        if st.button("🚀 Εκτέλεση Integration", type="primary"):
            try:
                with st.spinner("Εκτέλεση integration..."):
                    self.integrated_data = self.perform_integration()
                
                if self.integrated_data is not None:
                    st.success("✅ Integration ολοκληρώθηκε επιτυχώς!")
                    
                    # Αποθήκευση στο session state
                    st.session_state['integrated_data'] = self.integrated_data
                    
                    # Εμφάνιση στατιστικών
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("Κύτταρα", f"{self.integrated_data.n_obs:,}")
                    with col2:
                        st.metric("Γονίδια", f"{self.integrated_data.n_vars:,}")
                
            except Exception as e:
                st.error(f"❌ Σφάλμα κατά την integration: {str(e)}")
    
    def get_common_genes(self):
        """Εύρεση κοινών γονιδίων"""
        if not self.datasets:
            return []
        
        common_genes = None
        for adata in self.datasets.values():
            if common_genes is None:
                common_genes = set(adata.var_names)
            else:
                common_genes &= set(adata.var_names)
        
        return list(common_genes) if common_genes else []
    
    def perform_integration(self):
        """Εκτέλεση της integration"""
        
        if len(self.datasets) == 1:
            # Αν έχουμε μόνο ένα dataset, απλά το επιστρέφουμε
            dataset_name, adata = next(iter(self.datasets.items()))
            adata_copy = adata.copy()
            adata_copy.obs['dataset'] = dataset_name
            return adata_copy
        
        if self.integration_method == "scanorama" and SCANORAMA_AVAILABLE:
            return self.perform_scanorama_integration()
        else:
            return self.perform_simple_concatenation()
    
    def perform_scanorama_integration(self):
        """Integration με Scanorama"""
        
        st.info("🔬 Εκτέλεση Scanorama integration...")
        
        # Προετοιμασία δεδομένων
        datasets_list = []
        genes_list = []
        
        for dataset_name, adata in self.datasets.items():
            if hasattr(adata.X, 'toarray'):
                X_dense = adata.X.toarray()
            else:
                X_dense = adata.X
            
            datasets_list.append(X_dense)
            genes_list.append(adata.var_names.tolist())
        
        # Εκτέλεση Scanorama
        integrated_data, genes = scanorama.integrate(
            datasets_list,
            genes_list,
            knn=20,
            sigma=15.0,
            alpha=0.1,
            batch_size=5000
        )
        
        # Δημιουργία AnnData object
        integrated_X = np.vstack(integrated_data)
        
        # Συλλογή obs data
        obs_list = []
        for dataset_name, adata in self.datasets.items():
            obs_subset = adata.obs.copy()
            obs_subset['dataset'] = dataset_name
            obs_list.append(obs_subset)
        
        integrated_obs = pd.concat(obs_list, ignore_index=True)
        
        # Δημιουργία var data
        integrated_var = pd.DataFrame(index=genes)
        
        # Δημιουργία AnnData
        integrated_adata = sc.AnnData(
            X=integrated_X,
            obs=integrated_obs,
            var=integrated_var
        )
        
        return memory_manager.optimize_adata_for_memory(integrated_adata)
    
    def perform_simple_concatenation(self):
        """Απλή concatenation"""
        
        st.info("🔗 Εκτέλεση simple concatenation...")
        
        # Εύρεση κοινών γονιδίων αν χρειάζεται
        if self.intersection_only:
            common_genes = self.get_common_genes()
            if not common_genes and len(self.datasets) > 1:
                st.error("❌ Δεν βρέθηκαν κοινά γονίδια")
                return None
        
        # Προετοιμασία datasets
        processed_datasets = []
        
        for dataset_name, adata in self.datasets.items():
            adata_copy = adata.copy()
            
            # Φιλτράρισμα για κοινά γονίδια
            if self.intersection_only and common_genes and len(self.datasets) > 1:
                adata_copy = adata_copy[:, common_genes]
            
            # Προσθήκη dataset information
            adata_copy.obs['dataset'] = dataset_name
            
            processed_datasets.append(adata_copy)
        
        # Concatenation
        if len(processed_datasets) == 1:
            integrated_adata = processed_datasets[0]
        else:
            integrated_adata = sc.concat(processed_datasets, axis=0, join='outer')
        
        return memory_manager.optimize_adata_for_memory(integrated_adata)
    
    def render_results(self):
        """Εμφάνιση αποτελεσμάτων"""
        
        st.header("📊 Αποτελέσματα Integration")
        
        if self.integrated_data is None:
            if 'integrated_data' in st.session_state:
                self.integrated_data = st.session_state['integrated_data']
            else:
                st.warning("⚠️ Εκτελέστε πρώτα integration")
                return
        
        # Βασικές στατιστικές
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Συνολικά Κύτταρα", f"{self.integrated_data.n_obs:,}")
        with col2:
            st.metric("Συνολικά Γονίδια", f"{self.integrated_data.n_vars:,}")
        with col3:
            if 'dataset' in self.integrated_data.obs.columns:
                n_datasets = self.integrated_data.obs['dataset'].nunique()
                st.metric("Αριθμός Datasets", n_datasets)
        
        # Dataset distribution
        if 'dataset' in self.integrated_data.obs.columns:
            st.subheader("📊 Κατανομή Κυττάρων ανά Dataset")
            
            dataset_counts = self.integrated_data.obs['dataset'].value_counts()
            
            fig = px.bar(
                x=dataset_counts.index,
                y=dataset_counts.values,
                labels={'x': 'Dataset', 'y': 'Αριθμός Κυττάρων'},
                title="Κατανομή Κυττάρων ανά Dataset"
            )
            st.plotly_chart(fig, use_container_width=True)
        
        # Visualization
        st.subheader("📈 Οπτικοποίηση")
        
        viz_method = st.selectbox(
            "Μέθοδος Visualization:",
            ["PCA", "UMAP", "t-SNE"]
        )
        
        if st.button(f"Δημιουργία {viz_method} Plot"):
            try:
                with st.spinner(f"Υπολογισμός {viz_method}..."):
                    
                    if viz_method == "PCA":
                        sc.tl.pca(self.integrated_data)
                        basis = 'pca'
                    elif viz_method == "UMAP":
                        sc.tl.pca(self.integrated_data)
                        sc.pp.neighbors(self.integrated_data)
                        sc.tl.umap(self.integrated_data)
                        basis = 'umap'
                    else:  # t-SNE
                        sc.tl.pca(self.integrated_data)
                        sc.tl.tsne(self.integrated_data)
                        basis = 'tsne'
                    
                    # Plot
                    fig, ax = plt.subplots(figsize=(10, 8))
                    
                    if 'dataset' in self.integrated_data.obs.columns:
                        sc.pl.embedding(
                            self.integrated_data, 
                            basis=basis, 
                            color='dataset',
                            ax=ax,
                            show=False,
                            frameon=False
                        )
                        ax.set_title(f'{viz_method} - Colored by Dataset')
                    
                    st.pyplot(fig)
                    
            except Exception as e:
                st.error(f"❌ Σφάλμα στην οπτικοποίηση: {str(e)}")
        
        # Export
        st.subheader("💾 Export Data")
        
        export_format = st.selectbox("Μορφή Export:", ["H5AD", "CSV"])
        
        if st.button("📥 Download Data"):
            try:
                if export_format == "H5AD":
                    with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp_file:
                        self.integrated_data.write_h5ad(tmp_file.name)
                        
                        with open(tmp_file.name, 'rb') as f:
                            st.download_button(
                                label="💾 Download H5AD",
                                data=f.read(),
                                file_name="integrated_data.h5ad",
                                mime="application/octet-stream"
                            )
                        
                        os.unlink(tmp_file.name)
                
                else:  # CSV
                    if hasattr(self.integrated_data.X, 'toarray'):
                        X_dense = self.integrated_data.X.toarray()
                    else:
                        X_dense = self.integrated_data.X
                    
                    df = pd.DataFrame(
                        X_dense,
                        index=self.integrated_data.obs_names,
                        columns=self.integrated_data.var_names
                    )
                    
                    csv_buffer = io.StringIO()
                    df.to_csv(csv_buffer)
                    
                    st.download_button(
                        label="💾 Download CSV",
                        data=csv_buffer.getvalue(),
                        file_name="integrated_data.csv",
                        mime="text/csv"
                    )
                
            except Exception as e:
                st.error(f"❌ Σφάλμα κατά το export: {str(e)}")
