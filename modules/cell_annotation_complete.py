"""
Ολοκληρωμένο Module για Cell Type Annotation

Αυτό το module περιλαμβάνει προηγμένες λειτουργίες για την αυτοματοποιημένη
και χειροκίνητη αναγνώριση τύπων κυττάρων σε scRNA-seq δεδομένα, χρησιμοποιώντας
marker genes, clustering και machine learning τεχνικές.


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
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
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

# Decoupler import με error handling
try:
    import decoupler as dc
    DECOUPLER_AVAILABLE = True
except ImportError:
    DECOUPLER_AVAILABLE = False

class CellAnnotationPageComplete:
    """Ολοκληρωμένη κλάση για cell type annotation"""
    
    def __init__(self):
        self.adata = None
        self.marker_genes = {}
        self.annotation_results = None
        self.cell_types = []
        
    def render(self):
        """Κεντρική μέθοδος για την εμφάνιση της σελίδας"""
        
        st.title("🏷️ Σχολιασμός Κυττάρων")
        st.markdown("### Αυτοματοποιημένη και χειροκίνητη αναγνώριση τύπων κυττάρων")
        
        if not DECOUPLER_AVAILABLE:
            st.warning("⚠️ Decoupler δεν είναι εγκατεστημένο. Θα χρησιμοποιηθούν εναλλακτικές μέθοδοι.")
        
        # Memory information
        display_memory_info()
        
        # Tabs για διαφορετικές λειτουργίες
        tab1, tab2, tab3, tab4, tab5 = st.tabs([
            "📊 Δεδομένα & Clustering", 
            "🧬 Marker Genes",
            "🔍 Automated Annotation", 
            "✏️ Manual Curation",
            "📈 Results & Visualization"
        ])
        
        with tab1:
            self.render_data_clustering()
            
        with tab2:
            self.render_marker_genes()
            
        with tab3:
            self.render_automated_annotation()
            
        with tab4:
            self.render_manual_curation()
            
        with tab5:
            self.render_results_visualization()
    
    def render_data_clustering(self):
        """Δεδομένα και clustering για annotation"""
        
        st.header("📊 Δεδομένα & Clustering")
        
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
            
            # Έλεγχος για υπάρχοντα clusters
            cluster_columns = [col for col in self.adata.obs.columns 
                             if 'cluster' in col.lower() or 'leiden' in col.lower() or 'louvain' in col.lower()]
            
            if cluster_columns:
                st.subheader("🎯 Υπάρχοντα Clusters")
                
                selected_cluster_col = st.selectbox(
                    "Επιλέξτε cluster column:",
                    cluster_columns
                )
                
                if selected_cluster_col:
                    cluster_counts = self.adata.obs[selected_cluster_col].value_counts()
                    st.write(f"**Clusters στη στήλη '{selected_cluster_col}':**")
                    
                    for cluster, count in cluster_counts.items():
                        st.write(f"- Cluster {cluster}: {count:,} κύτταρα")
                    
                    # Visualization των clusters
                    if st.button("📊 Visualize Clusters"):
                        self.visualize_clusters(selected_cluster_col)
            else:
                st.subheader("🎯 Clustering")
                st.info("Δεν βρέθηκαν υπάρχοντα clusters. Θα εκτελεστεί νέο clustering.")
                
                if st.button("🔄 Εκτέλεση Clustering"):
                    self.perform_clustering()
    
    @memory_monitor
    def perform_clustering(self):
        """Εκτέλεση clustering"""
        
        try:
            with st.spinner("Εκτέλεση clustering..."):
                
                # PCA αν δεν υπάρχει
                if 'X_pca' not in self.adata.obsm:
                    sc.tl.pca(self.adata)
                
                # Neighbors
                sc.pp.neighbors(self.adata)
                
                # Leiden clustering
                sc.tl.leiden(self.adata, resolution=0.5)
                
                # UMAP για visualization
                sc.tl.umap(self.adata)
                
                st.success("✅ Clustering ολοκληρώθηκε!")
                
                # Εμφάνιση αποτελεσμάτων
                n_clusters = self.adata.obs['leiden'].nunique()
                st.info(f"📊 Βρέθηκαν {n_clusters} clusters")
                
                # Visualization
                self.visualize_clusters('leiden')
                
        except Exception as e:
            st.error(f"❌ Σφάλμα στο clustering: {str(e)}")
    
    def visualize_clusters(self, cluster_col):
        """Visualization των clusters"""
        
        # Έλεγχος για UMAP
        if 'X_umap' not in self.adata.obsm:
            st.warning("⚠️ UMAP δεν είναι διαθέσιμο. Εκτελώντας UMAP...")
            try:
                if 'X_pca' not in self.adata.obsm:
                    sc.tl.pca(self.adata)
                sc.pp.neighbors(self.adata)
                sc.tl.umap(self.adata)
            except Exception as e:
                st.error(f"❌ Σφάλμα στον υπολογισμό UMAP: {str(e)}")
                return
        
        # UMAP plot με clusters
        umap_coords = self.adata.obsm['X_umap']
        
        plot_df = pd.DataFrame({
            'UMAP_1': umap_coords[:, 0],
            'UMAP_2': umap_coords[:, 1],
            'Cluster': self.adata.obs[cluster_col].astype(str)
        })
        
        fig = px.scatter(
            plot_df,
            x='UMAP_1',
            y='UMAP_2',
            color='Cluster',
            title=f'UMAP - Clusters ({cluster_col})',
            hover_data=['Cluster']
        )
        
        fig.update_layout(height=600)
        st.plotly_chart(fig, use_container_width=True)
    
    def render_marker_genes(self):
        """Διαχείριση marker genes"""
        
        st.header("🧬 Marker Genes")
        
        if self.adata is None:
            st.warning("⚠️ Επιλέξτε πρώτα δεδομένα")
            return
        
        # Tabs για διαφορετικές πηγές marker genes
        marker_tab1, marker_tab2, marker_tab3 = st.tabs([
            "📚 Built-in Markers",
            "📁 Upload Custom Markers", 
            "🔍 Find Cluster Markers"
        ])
        
        with marker_tab1:
            self.render_builtin_markers()
        
        with marker_tab2:
            self.render_custom_markers()
        
        with marker_tab3:
            self.render_cluster_markers()
    
    def render_builtin_markers(self):
        """Built-in marker genes"""
        
        st.subheader("📚 Built-in Marker Databases")
        
        # Προκαθορισμένα marker genes για κοινούς τύπους κυττάρων
        builtin_markers = {
            "Immune Cells": {
                "T Cells": ["CD3D", "CD3E", "CD3G", "CD8A", "CD4"],
                "B Cells": ["CD19", "MS4A1", "CD79A", "CD79B"],
                "NK Cells": ["KLRD1", "KLRF1", "NCR1", "GNLY"],
                "Macrophages": ["CD68", "CD163", "CSF1R", "AIF1"],
                "Dendritic Cells": ["CD1C", "FCER1A", "CLEC9A"],
                "Neutrophils": ["S100A8", "S100A9", "FCGR3B"]
            },
            "Brain Cells": {
                "Neurons": ["SNAP25", "SYT1", "STMN2", "TUBB3"],
                "Astrocytes": ["GFAP", "AQP4", "S100B", "ALDH1L1"],
                "Oligodendrocytes": ["MBP", "MOG", "PLP1", "CNP"],
                "Microglia": ["CX3CR1", "P2RY12", "TMEM119", "AIF1"],
                "Endothelial": ["PECAM1", "VWF", "CDH5", "FLT1"]
            },
            "Tissue Stem Cells": {
                "Mesenchymal": ["THY1", "ENG", "NT5E", "PDGFRB"],
                "Hematopoietic": ["CD34", "KIT", "FLT3", "PROM1"],
                "Epithelial": ["EPCAM", "KRT8", "KRT18", "CDH1"]
            }
        }
        
        # Επιλογή κατηγορίας
        selected_category = st.selectbox(
            "Επιλέξτε κατηγορία:",
            list(builtin_markers.keys())
        )
        
        if selected_category:
            st.write(f"**Διαθέσιμοι τύποι κυττάρων για {selected_category}:**")
            
            # Εμφάνιση marker genes
            for cell_type, markers in builtin_markers[selected_category].items():
                with st.expander(f"{cell_type}"):
                    st.write(f"**Marker genes:** {', '.join(markers)}")
                    
                    # Έλεγχος διαθεσιμότητας στο dataset
                    available_markers = [m for m in markers if m in self.adata.var_names]
                    if available_markers:
                        st.success(f"✅ Διαθέσιμα: {', '.join(available_markers)}")
                        
                        if st.button(f"📊 Visualize {cell_type} Markers", key=f"viz_{cell_type}"):
                            self.visualize_marker_expression(available_markers, cell_type)
                    else:
                        st.warning("⚠️ Κανένα marker gene δεν βρέθηκε στο dataset")
            
            # Επιλογή για χρήση στο annotation
            if st.button(f"✅ Χρήση {selected_category} Markers για Annotation"):
                self.marker_genes = builtin_markers[selected_category]
                st.success(f"✅ Marker genes από {selected_category} φορτώθηκαν!")
    
    def render_custom_markers(self):
        """Upload custom marker genes"""
        
        st.subheader("📁 Upload Custom Marker Genes")
        
        # File uploader
        uploaded_file = st.file_uploader(
            "Επιλέξτε αρχείο με marker genes:",
            type=['csv', 'tsv', 'txt'],
            help="Format: Cell_Type, Gene1, Gene2, Gene3..."
        )
        
        if uploaded_file is not None:
            try:
                # Διαβάστε το αρχείο
                if uploaded_file.name.endswith('.csv'):
                    df = pd.read_csv(uploaded_file)
                else:
                    df = pd.read_csv(uploaded_file, sep='\t')
                
                st.write("**Προεπισκόπηση αρχείου:**")
                st.dataframe(df.head(), use_container_width=True)
                
                # Parsing των marker genes
                if st.button("📊 Parse Marker Genes"):
                    self.parse_custom_markers(df)
                
            except Exception as e:
                st.error(f"❌ Σφάλμα στη φόρτωση αρχείου: {str(e)}")
        
        # Manual entry
        st.subheader("✏️ Manual Entry")
        
        cell_type_name = st.text_input("Όνομα τύπου κυττάρου:")
        marker_genes_text = st.text_area(
            "Marker genes (ένα ανά γραμμή ή διαχωρισμένα με κόμμα):",
            placeholder="CD3D\nCD3E\nCD4"
        )
        
        if st.button("➕ Προσθήκη Custom Markers") and cell_type_name and marker_genes_text:
            # Parse τα genes
            if ',' in marker_genes_text:
                genes = [g.strip().upper() for g in marker_genes_text.split(',')]
            else:
                genes = [g.strip().upper() for g in marker_genes_text.split('\n') if g.strip()]
            
            # Έλεγχος διαθεσιμότητας
            available_genes = [g for g in genes if g in self.adata.var_names]
            
            if available_genes:
                if 'Custom' not in self.marker_genes:
                    self.marker_genes['Custom'] = {}
                
                self.marker_genes['Custom'][cell_type_name] = available_genes
                st.success(f"✅ Προστέθηκαν {len(available_genes)} marker genes για {cell_type_name}")
                st.info(f"Διαθέσιμα genes: {', '.join(available_genes)}")
            else:
                st.warning("⚠️ Κανένα από τα genes δεν βρέθηκε στο dataset")
    
    def parse_custom_markers(self, df):
        """Parse custom marker genes από DataFrame"""
        
        # Υποθέτουμε format: Cell_Type, Gene1, Gene2, ...
        custom_markers = {}
        
        for _, row in df.iterrows():
            cell_type = str(row.iloc[0])
            genes = [str(g).upper().strip() for g in row.iloc[1:] if pd.notna(g) and str(g).strip()]
            
            # Έλεγχος διαθεσιμότητας
            available_genes = [g for g in genes if g in self.adata.var_names]
            
            if available_genes:
                custom_markers[cell_type] = available_genes
        
        if custom_markers:
            self.marker_genes['Custom'] = custom_markers
            st.success(f"✅ Φορτώθηκαν marker genes για {len(custom_markers)} τύπους κυττάρων")
            
            # Εμφάνιση περίληψης
            for cell_type, genes in custom_markers.items():
                st.write(f"**{cell_type}:** {len(genes)} genes ({', '.join(genes[:5])}{'...' if len(genes) > 5 else ''})")
        else:
            st.warning("⚠️ Δεν βρέθηκαν διαθέσιμα marker genes")
    
    def render_cluster_markers(self):
        """Εύρεση marker genes για clusters"""
        
        st.subheader("🔍 Find Cluster Markers")
        
        # Έλεγχος για clusters
        cluster_columns = [col for col in self.adata.obs.columns 
                         if 'cluster' in col.lower() or 'leiden' in col.lower() or 'louvain' in col.lower()]
        
        if not cluster_columns:
            st.warning("⚠️ Δεν βρέθηκαν clusters. Εκτελέστε πρώτα clustering.")
            return
        
        selected_cluster_col = st.selectbox(
            "Επιλέξτε cluster column:",
            cluster_columns,
            key="marker_cluster_col"
        )
        
        # Παράμετροι για marker finding
        col1, col2 = st.columns(2)
        
        with col1:
            min_pct = st.slider(
                "Min percentage:", 
                0.1, 1.0, 0.25, step=0.05,
                help="Ελάχιστο ποσοστό κυττάρων που εκφράζουν το gene"
            )
        
        with col2:
            logfc_threshold = st.slider(
                "Log FC threshold:", 
                0.1, 2.0, 0.25, step=0.05,
                help="Ελάχιστο log fold change"
            )
        
        if st.button("🔍 Find Cluster Markers"):
            try:
                with st.spinner("Αναζήτηση cluster markers..."):
                    
                    # Scanpy rank genes groups
                    sc.tl.rank_genes_groups(
                        self.adata,
                        groupby=selected_cluster_col,
                        method='wilcoxon',
                        min_pct=min_pct,
                        logfc_threshold=logfc_threshold
                    )
                    
                    # Εξαγωγή top markers για κάθε cluster
                    cluster_markers = {}
                    
                    for cluster in self.adata.obs[selected_cluster_col].unique():
                        cluster_df = sc.get.rank_genes_groups_df(self.adata, group=str(cluster))
                        top_genes = cluster_df.head(10)['names'].tolist()
                        cluster_markers[f"Cluster_{cluster}"] = top_genes
                    
                    # Αποθήκευση
                    self.marker_genes['Cluster_Markers'] = cluster_markers
                    
                    st.success("✅ Cluster markers βρέθηκαν!")
                    
                    # Εμφάνιση αποτελεσμάτων
                    for cluster, markers in cluster_markers.items():
                        with st.expander(f"{cluster} - Top Markers"):
                            st.write(f"**Top 10 markers:** {', '.join(markers)}")
                    
            except Exception as e:
                st.error(f"❌ Σφάλμα στην εύρεση markers: {str(e)}")
    
    def visualize_marker_expression(self, genes, cell_type):
        """Visualization της έκφρασης marker genes"""
        
        if 'X_umap' not in self.adata.obsm:
            st.warning("⚠️ UMAP δεν είναι διαθέσιμο για visualization")
            return
        
        # Subplots για πολλαπλά genes
        n_genes = min(len(genes), 6)  # Μέγιστο 6 genes
        cols = min(3, n_genes)
        rows = (n_genes + cols - 1) // cols
        
        fig = make_subplots(
            rows=rows, cols=cols,
            subplot_titles=genes[:n_genes]
        )
        
        umap_coords = self.adata.obsm['X_umap']
        
        for i, gene in enumerate(genes[:n_genes]):
            row = i // cols + 1
            col = i % cols + 1
            
            gene_idx = self.adata.var_names.get_loc(gene)
            
            if hasattr(self.adata.X, 'toarray'):
                expression = self.adata.X[:, gene_idx].toarray().flatten()
            else:
                expression = self.adata.X[:, gene_idx]
            
            fig.add_trace(
                go.Scatter(
                    x=umap_coords[:, 0],
                    y=umap_coords[:, 1],
                    mode='markers',
                    marker=dict(
                        color=expression,
                        colorscale='viridis',
                        showscale=(i == 0),
                        size=3
                    ),
                    showlegend=False,
                    name=gene
                ),
                row=row, col=col
            )
        
        fig.update_layout(
            height=300 * rows,
            title_text=f"{cell_type} Marker Expression"
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    @memory_monitor
    def render_automated_annotation(self):
        """Automated cell type annotation"""
        
        st.header("🔍 Automated Annotation")
        
        if self.adata is None:
            st.warning("⚠️ Επιλέξτε πρώτα δεδομένα")
            return
        
        if not self.marker_genes:
            st.warning("⚠️ Φορτώστε πρώτα marker genes")
            return
        
        # Επιλογή μεθόδου annotation
        annotation_method = st.selectbox(
            "Μέθοδος Annotation:",
            ["Score-based", "Clustering-based", "Machine Learning"] + 
            (["Decoupler"] if DECOUPLER_AVAILABLE else [])
        )
        
        # Παράμετροι annotation
        with st.expander("🔧 Παράμετροι Annotation"):
            
            score_threshold = st.slider(
                "Score threshold:",
                0.0, 1.0, 0.5, step=0.05,
                help="Ελάχιστο score για annotation"
            )
            
            min_genes = st.slider(
                "Minimum genes per cell type:",
                1, 10, 3,
                help="Ελάχιστος αριθμός marker genes που πρέπει να εκφράζονται"
            )
        
        # Κουμπί annotation
        if st.button("🚀 Εκτέλεση Automated Annotation"):
            
            try:
                with st.spinner("Εκτέλεση automated annotation..."):
                    
                    if annotation_method == "Score-based":
                        self.annotation_results = self.perform_score_based_annotation(score_threshold, min_genes)
                    elif annotation_method == "Clustering-based":
                        self.annotation_results = self.perform_clustering_based_annotation()
                    elif annotation_method == "Machine Learning":
                        self.annotation_results = self.perform_ml_annotation()
                    elif annotation_method == "Decoupler" and DECOUPLER_AVAILABLE:
                        self.annotation_results = self.perform_decoupler_annotation()
                
                if self.annotation_results is not None:
                    st.success("✅ Automated annotation ολοκληρώθηκε!")
                    
                    # Αποθήκευση στο adata
                    self.adata.obs['cell_type_auto'] = self.annotation_results['cell_type']
                    self.adata.obs['annotation_score'] = self.annotation_results['score']
                    
                    # Αποθήκευση στο session state
                    st.session_state['annotated_data'] = self.adata
                    
                    # Εμφάνιση περίληψης
                    self.display_annotation_summary()
                
            except Exception as e:
                st.error(f"❌ Σφάλμα στο automated annotation: {str(e)}")
    
    def perform_score_based_annotation(self, score_threshold, min_genes):
        """Score-based annotation"""
        
        st.info("📊 Εκτέλεση score-based annotation...")
        
        # Υπολογισμός scores για κάθε κύτταρο
        cell_scores = []
        
        for cell_idx in range(self.adata.n_obs):
            
            if hasattr(self.adata.X, 'toarray'):
                cell_expression = self.adata.X[cell_idx, :].toarray().flatten()
            else:
                cell_expression = self.adata.X[cell_idx, :]
            
            cell_type_scores = {}
            
            # Υπολογισμός score για κάθε cell type
            for category, cell_types in self.marker_genes.items():
                for cell_type, markers in cell_types.items():
                    
                    # Βρες indices των marker genes
                    marker_indices = [self.adata.var_names.get_loc(gene) 
                                    for gene in markers if gene in self.adata.var_names]
                    
                    if len(marker_indices) >= min_genes:
                        # Mean expression των marker genes
                        marker_expression = cell_expression[marker_indices]
                        score = np.mean(marker_expression)
                        cell_type_scores[cell_type] = score
            
            # Επιλογή καλύτερου cell type
            if cell_type_scores:
                best_cell_type = max(cell_type_scores, key=cell_type_scores.get)
                best_score = cell_type_scores[best_cell_type]
                
                if best_score >= score_threshold:
                    cell_scores.append({
                        'cell_type': best_cell_type,
                        'score': best_score
                    })
                else:
                    cell_scores.append({
                        'cell_type': 'Unknown',
                        'score': best_score
                    })
            else:
                cell_scores.append({
                    'cell_type': 'Unknown',
                    'score': 0.0
                })
        
        return pd.DataFrame(cell_scores)
    
    def perform_clustering_based_annotation(self):
        """Clustering-based annotation"""
        
        st.info("🎯 Εκτέλεση clustering-based annotation...")
        
        # Έλεγχος για clusters
        cluster_columns = [col for col in self.adata.obs.columns 
                         if 'cluster' in col.lower() or 'leiden' in col.lower()]
        
        if not cluster_columns:
            st.error("❌ Δεν βρέθηκαν clusters για annotation")
            return None
        
        cluster_col = cluster_columns[0]  # Χρήση πρώτου cluster column
        
        # Υπολογισμός mean expression ανά cluster
        cluster_annotations = {}
        
        for cluster in self.adata.obs[cluster_col].unique():
            cluster_mask = self.adata.obs[cluster_col] == cluster
            cluster_cells = self.adata[cluster_mask, :]
            
            # Mean expression για το cluster
            if hasattr(cluster_cells.X, 'toarray'):
                cluster_mean = np.mean(cluster_cells.X.toarray(), axis=0)
            else:
                cluster_mean = np.mean(cluster_cells.X, axis=0)
            
            # Βρες καλύτερο cell type για το cluster
            best_cell_type = "Unknown"
            best_score = 0.0
            
            for category, cell_types in self.marker_genes.items():
                for cell_type, markers in cell_types.items():
                    
                    marker_indices = [self.adata.var_names.get_loc(gene) 
                                    for gene in markers if gene in self.adata.var_names]
                    
                    if marker_indices:
                        score = np.mean(cluster_mean[marker_indices])
                        if score > best_score:
                            best_score = score
                            best_cell_type = cell_type
            
            cluster_annotations[cluster] = {
                'cell_type': best_cell_type,
                'score': best_score
            }
        
        # Assign annotations σε κύτταρα
        cell_annotations = []
        for cluster in self.adata.obs[cluster_col]:
            annotation = cluster_annotations[cluster]
            cell_annotations.append(annotation)
        
        return pd.DataFrame(cell_annotations)
    
    def perform_ml_annotation(self):
        """Machine learning-based annotation"""
        
        st.info("🤖 Εκτέλεση ML-based annotation...")
        st.warning("⚠️ ML annotation είναι experimental feature")
        
        # Placeholder για ML implementation
        # Θα χρησιμοποιούσε trained models ή reference datasets
        
        return self.perform_score_based_annotation(0.5, 3)  # Fallback to score-based
    
    def perform_decoupler_annotation(self):
        """Decoupler-based annotation"""
        
        st.info("🔬 Εκτέλεση Decoupler annotation...")
        
        # Placeholder για Decoupler implementation
        # Θα χρησιμοποιούσε decoupler library
        
        return self.perform_score_based_annotation(0.5, 3)  # Fallback to score-based
    
    def display_annotation_summary(self):
        """Εμφάνιση περίληψης annotation"""
        
        st.subheader("📊 Περίληψη Annotation")
        
        if self.annotation_results is not None:
            
            # Κατανομή cell types
            type_counts = self.annotation_results['cell_type'].value_counts()
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.write("**Κατανομή Cell Types:**")
                for cell_type, count in type_counts.items():
                    percentage = (count / len(self.annotation_results)) * 100
                    st.write(f"- {cell_type}: {count:,} ({percentage:.1f}%)")
            
            with col2:
                # Pie chart
                fig = px.pie(
                    values=type_counts.values,
                    names=type_counts.index,
                    title="Cell Type Distribution"
                )
                st.plotly_chart(fig, use_container_width=True)
            
            # Score distribution
            st.subheader("📈 Score Distribution")
            
            fig = px.histogram(
                self.annotation_results,
                x='score',
                title='Annotation Score Distribution',
                nbins=50
            )
            st.plotly_chart(fig, use_container_width=True)
    
    def render_manual_curation(self):
        """Manual curation των annotations"""
        
        st.header("✏️ Manual Curation")
        
        if self.annotation_results is None:
            if 'annotated_data' in st.session_state:
                self.adata = st.session_state['annotated_data']
                if 'cell_type_auto' in self.adata.obs.columns:
                    self.annotation_results = pd.DataFrame({
                        'cell_type': self.adata.obs['cell_type_auto'],
                        'score': self.adata.obs.get('annotation_score', 0.5)
                    })
            else:
                st.warning("⚠️ Εκτελέστε πρώτα automated annotation")
                return
        
        st.subheader("🔧 Manual Corrections")
        
        # Επιλογή clusters για manual curation
        cluster_columns = [col for col in self.adata.obs.columns 
                         if 'cluster' in col.lower() or 'leiden' in col.lower()]
        
        if cluster_columns:
            cluster_col = st.selectbox(
                "Cluster column για curation:",
                cluster_columns
            )
            
            # Εμφάνιση current annotations ανά cluster
            cluster_annotations = {}
            for cluster in self.adata.obs[cluster_col].unique():
                cluster_mask = self.adata.obs[cluster_col] == cluster
                cluster_types = self.adata.obs.loc[cluster_mask, 'cell_type_auto'].value_counts()
                most_common = cluster_types.index[0] if len(cluster_types) > 0 else "Unknown"
                cluster_annotations[cluster] = most_common
            
            st.write("**Current Cluster Annotations:**")
            
            # Manual correction interface
            corrections = {}
            for cluster, current_type in cluster_annotations.items():
                col1, col2 = st.columns([1, 2])
                
                with col1:
                    st.write(f"Cluster {cluster}:")
                
                with col2:
                    new_type = st.text_input(
                        f"Cell type:",
                        value=current_type,
                        key=f"cluster_{cluster}_type"
                    )
                    if new_type != current_type:
                        corrections[cluster] = new_type
            
            # Apply corrections
            if corrections and st.button("✅ Apply Manual Corrections"):
                
                # Δημιουργία corrected annotations
                corrected_annotations = self.adata.obs['cell_type_auto'].copy()
                
                for cluster, new_type in corrections.items():
                    cluster_mask = self.adata.obs[cluster_col] == cluster
                    corrected_annotations.loc[cluster_mask] = new_type
                
                # Αποθήκευση
                self.adata.obs['cell_type_curated'] = corrected_annotations
                
                st.success(f"✅ Εφαρμόστηκαν {len(corrections)} διορθώσεις!")
                
                # Update annotation results
                self.annotation_results['cell_type'] = corrected_annotations
        
        # Bulk corrections
        st.subheader("🔄 Bulk Corrections")
        
        col1, col2 = st.columns(2)
        
        with col1:
            old_type = st.text_input("Παλιός τύπος κυττάρου:")
        
        with col2:
            new_type = st.text_input("Νέος τύπος κυττάρου:")
        
        if st.button("🔄 Replace All") and old_type and new_type:
            
            if 'cell_type_curated' in self.adata.obs.columns:
                target_col = 'cell_type_curated'
            else:
                target_col = 'cell_type_auto'
            
            mask = self.adata.obs[target_col] == old_type
            n_replaced = mask.sum()
            
            if n_replaced > 0:
                self.adata.obs.loc[mask, target_col] = new_type
                st.success(f"✅ Αντικαταστάθηκαν {n_replaced} κύτταρα από '{old_type}' σε '{new_type}'")
            else:
                st.warning(f"⚠️ Δεν βρέθηκαν κύτταρα τύπου '{old_type}'")
    
    def render_results_visualization(self):
        """Visualization των annotation results"""
        
        st.header("📈 Results & Visualization")
        
        if self.adata is None:
            st.warning("⚠️ Δεν υπάρχουν annotated data")
            return
        
        # Έλεγχος για annotation columns
        annotation_cols = [col for col in self.adata.obs.columns 
                         if 'cell_type' in col.lower()]
        
        if not annotation_cols:
            st.warning("⚠️ Δεν βρέθηκαν annotation results")
            return
        
        # Επιλογή annotation column
        selected_annotation = st.selectbox(
            "Επιλέξτε annotation column:",
            annotation_cols
        )
        
        # Visualization options
        viz_type = st.selectbox(
            "Τύπος Visualization:",
            ["UMAP Plot", "Statistics", "Marker Expression", "Export Results"]
        )
        
        if viz_type == "UMAP Plot":
            self.render_annotation_umap(selected_annotation)
        elif viz_type == "Statistics":
            self.render_annotation_statistics(selected_annotation)
        elif viz_type == "Marker Expression":
            self.render_annotation_markers(selected_annotation)
        else:  # Export Results
            self.render_annotation_export(selected_annotation)
    
    def render_annotation_umap(self, annotation_col):
        """UMAP plot με annotations"""
        
        if 'X_umap' not in self.adata.obsm:
            st.warning("⚠️ UMAP δεν είναι διαθέσιμο")
            return
        
        umap_coords = self.adata.obsm['X_umap']
        
        plot_df = pd.DataFrame({
            'UMAP_1': umap_coords[:, 0],
            'UMAP_2': umap_coords[:, 1],
            'Cell_Type': self.adata.obs[annotation_col].astype(str)
        })
        
        fig = px.scatter(
            plot_df,
            x='UMAP_1',
            y='UMAP_2',
            color='Cell_Type',
            title=f'UMAP - Cell Type Annotations ({annotation_col})',
            hover_data=['Cell_Type']
        )
        
        fig.update_layout(height=700)
        st.plotly_chart(fig, use_container_width=True)
    
    def render_annotation_statistics(self, annotation_col):
        """Στατιστικά annotation"""
        
        st.subheader("📊 Annotation Statistics")
        
        # Κατανομή cell types
        type_counts = self.adata.obs[annotation_col].value_counts()
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**Cell Type Counts:**")
            for cell_type, count in type_counts.items():
                percentage = (count / len(self.adata.obs)) * 100
                st.write(f"- {cell_type}: {count:,} ({percentage:.1f}%)")
        
        with col2:
            # Bar chart
            fig = px.bar(
                x=type_counts.values,
                y=type_counts.index,
                orientation='h',
                title="Cell Type Distribution",
                labels={'x': 'Number of Cells', 'y': 'Cell Type'}
            )
            st.plotly_chart(fig, use_container_width=True)
        
        # Score distribution αν υπάρχει
        if 'annotation_score' in self.adata.obs.columns:
            st.subheader("📈 Score Distribution by Cell Type")
            
            fig = px.box(
                x=self.adata.obs[annotation_col],
                y=self.adata.obs['annotation_score'],
                title='Annotation Scores by Cell Type'
            )
            fig.update_xaxes(tickangle=45)
            st.plotly_chart(fig, use_container_width=True)
    
    def render_annotation_markers(self, annotation_col):
        """Marker expression για annotated cell types"""
        
        st.subheader("🧬 Marker Expression by Cell Type")
        
        if not self.marker_genes:
            st.warning("⚠️ Δεν υπάρχουν φορτωμένα marker genes")
            return
        
        # Επιλογή cell type
        cell_types = self.adata.obs[annotation_col].unique()
        selected_cell_type = st.selectbox(
            "Επιλέξτε cell type:",
            cell_types
        )
        
        # Βρες markers για αυτό το cell type
        relevant_markers = []
        for category, types in self.marker_genes.items():
            if selected_cell_type in types:
                relevant_markers.extend(types[selected_cell_type])
        
        if not relevant_markers:
            st.warning(f"⚠️ Δεν βρέθηκαν marker genes για {selected_cell_type}")
            return
        
        # Υπολογισμός mean expression
        cell_type_mask = self.adata.obs[annotation_col] == selected_cell_type
        other_mask = ~cell_type_mask
        
        marker_data = []
        
        for marker in relevant_markers[:10]:  # Top 10 markers
            if marker in self.adata.var_names:
                gene_idx = self.adata.var_names.get_loc(marker)
                
                if hasattr(self.adata.X, 'toarray'):
                    target_exp = np.mean(self.adata.X[cell_type_mask, gene_idx].toarray())
                    other_exp = np.mean(self.adata.X[other_mask, gene_idx].toarray())
                else:
                    target_exp = np.mean(self.adata.X[cell_type_mask, gene_idx])
                    other_exp = np.mean(self.adata.X[other_mask, gene_idx])
                
                marker_data.append({
                    'Gene': marker,
                    f'{selected_cell_type}': target_exp,
                    'Other Cells': other_exp,
                    'Fold_Change': target_exp / (other_exp + 1e-6)
                })
        
        if marker_data:
            marker_df = pd.DataFrame(marker_data)
            
            # Bar plot
            fig = go.Figure()
            
            fig.add_trace(go.Bar(
                name=selected_cell_type,
                x=marker_df['Gene'],
                y=marker_df[selected_cell_type]
            ))
            
            fig.add_trace(go.Bar(
                name='Other Cells',
                x=marker_df['Gene'],
                y=marker_df['Other Cells']
            ))
            
            fig.update_layout(
                title=f'Marker Expression - {selected_cell_type}',
                barmode='group',
                xaxis_tickangle=-45
            )
            
            st.plotly_chart(fig, use_container_width=True)
    
    def render_annotation_export(self, annotation_col):
        """Export annotation results"""
        
        st.subheader("💾 Export Annotation Results")
        
        # Export options
        export_format = st.selectbox(
            "Μορφή export:",
            ["CSV", "Excel", "H5AD"]
        )
        
        include_options = st.multiselect(
            "Τι να συμπεριληφθεί:",
            ["Cell metadata", "Annotation scores", "Cluster information"],
            default=["Cell metadata", "Annotation scores"]
        )
        
        if st.button("📥 Download Annotation Results"):
            
            try:
                # Προετοιμασία δεδομένων
                export_data = pd.DataFrame({
                    'Cell_ID': self.adata.obs_names,
                    'Cell_Type': self.adata.obs[annotation_col]
                })
                
                if "Annotation scores" in include_options and 'annotation_score' in self.adata.obs.columns:
                    export_data['Annotation_Score'] = self.adata.obs['annotation_score']
                
                if "Cell metadata" in include_options:
                    for col in self.adata.obs.columns:
                        if col not in export_data.columns:
                            export_data[col] = self.adata.obs[col]
                
                # Export
                if export_format == "CSV":
                    csv_buffer = io.StringIO()
                    export_data.to_csv(csv_buffer, index=False)
                    
                    st.download_button(
                        label="💾 Download CSV",
                        data=csv_buffer.getvalue(),
                        file_name="cell_annotations.csv",
                        mime="text/csv"
                    )
                
                elif export_format == "Excel":
                    excel_buffer = io.BytesIO()
                    with pd.ExcelWriter(excel_buffer, engine='openpyxl') as writer:
                        export_data.to_excel(writer, sheet_name='Annotations', index=False)
                        
                        # Summary sheet
                        summary = self.adata.obs[annotation_col].value_counts().reset_index()
                        summary.columns = ['Cell_Type', 'Count']
                        summary.to_excel(writer, sheet_name='Summary', index=False)
                    
                    st.download_button(
                        label="💾 Download Excel",
                        data=excel_buffer.getvalue(),
                        file_name="cell_annotations.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
                
                else:  # H5AD
                    with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp_file:
                        self.adata.write_h5ad(tmp_file.name)
                        
                        with open(tmp_file.name, 'rb') as f:
                            st.download_button(
                                label="💾 Download H5AD",
                                data=f.read(),
                                file_name="annotated_data.h5ad",
                                mime="application/octet-stream"
                            )
                        
                        os.unlink(tmp_file.name)
                
            except Exception as e:
                st.error(f"❌ Σφάλμα κατά το export: {str(e)}")
