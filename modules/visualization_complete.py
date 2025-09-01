"""
Ολοκληρωμένο Module για Οπτικοποιήσεις

Αυτό το module περιλαμβάνει εκτενείς δυνατότητες οπτικοποίησης για
scRNA-seq δεδομένα, συμπεριλαμβανομένου dimensionality reduction,
gene expression plots, quality control visualizations και πολλά άλλα.


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
import io
import base64
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

# Matplotlib backend setup
plt.switch_backend('Agg')

class VisualizationPageComplete:
    """Ολοκληρωμένη κλάση για οπτικοποιήσεις scRNA-seq δεδομένων"""
    
    def __init__(self):
        self.adata = None
        self.plot_data = None
        
    def render(self):
        """Κεντρική μέθοδος για την εμφάνιση της σελίδας"""
        
        st.title("📈 Οπτικοποιήσεις")
        st.markdown("### Διαδραστικά plots και γραφήματα για scRNA-seq δεδομένα")
        
        # Memory information
        display_memory_info()
        
        # Tabs για διαφορετικές κατηγορίες οπτικοποιήσεων
        tab1, tab2, tab3, tab4, tab5 = st.tabs([
            "📊 Δεδομένα & Setup", 
            "🗺️ Dimensionality Reduction",
            "🧬 Gene Expression", 
            "📊 Quality Control",
            "🎨 Advanced Plots"
        ])
        
        with tab1:
            self.render_data_setup()
            
        with tab2:
            self.render_dimensionality_reduction()
            
        with tab3:
            self.render_gene_expression()
            
        with tab4:
            self.render_quality_control()
            
        with tab5:
            self.render_advanced_plots()
    
    def render_data_setup(self):
        """Επιλογή δεδομένων και setup"""
        
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
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Κύτταρα", f"{self.adata.n_obs:,}")
            with col2:
                st.metric("Γονίδια", f"{self.adata.n_vars:,}")
            with col3:
                memory_size = memory_manager.estimate_memory_usage(self.adata)
                st.metric("Μέγεθος", f"{memory_size:.1f} MB")
            with col4:
                should_subsample = self.adata.n_obs > 10000
                st.metric("Subsample", "Προτείνεται" if should_subsample else "Όχι")
            
            # Subsample options για μεγάλα datasets
            if should_subsample:
                st.subheader("🎯 Subsample Options")
                
                use_subsample = st.checkbox(
                    "Χρήση subsample για γρηγορότερα plots",
                    value=True,
                    help="Συνιστάται για datasets >10K κύτταρα"
                )
                
                if use_subsample:
                    subsample_size = st.slider(
                        "Μέγεθος subsample:",
                        min_value=1000, 
                        max_value=min(20000, self.adata.n_obs), 
                        value=min(5000, self.adata.n_obs),
                        step=500
                    )
                    
                    if st.button("🎲 Δημιουργία Subsample"):
                        self.plot_data = self.create_subsample(subsample_size)
                        st.success(f"✅ Subsample δημιουργήθηκε με {self.plot_data.n_obs:,} κύτταρα")
                else:
                    self.plot_data = self.adata
            else:
                self.plot_data = self.adata
            
            # Εμφάνιση διαθέσιμων metadata
            st.subheader("🏷️ Διαθέσιμα Metadata")
            
            if self.plot_data is not None:
                categorical_cols = []
                numerical_cols = []
                
                for col in self.plot_data.obs.columns:
                    if self.plot_data.obs[col].dtype == 'object' or self.plot_data.obs[col].dtype.name == 'category':
                        categorical_cols.append(col)
                    else:
                        numerical_cols.append(col)
                
                col1, col2 = st.columns(2)
                
                with col1:
                    if categorical_cols:
                        st.write("**Κατηγορικές μεταβλητές:**")
                        for col in categorical_cols[:10]:  # Εμφάνιση πρώτων 10
                            unique_count = self.plot_data.obs[col].nunique()
                            st.write(f"- {col} ({unique_count} categories)")
                
                with col2:
                    if numerical_cols:
                        st.write("**Αριθμητικές μεταβλητές:**")
                        for col in numerical_cols[:10]:  # Εμφάνιση πρώτων 10
                            st.write(f"- {col}")
    
    def create_subsample(self, n_cells):
        """Δημιουργία subsample του dataset"""
        
        if n_cells >= self.adata.n_obs:
            return self.adata
        
        # Random subsample
        indices = np.random.choice(self.adata.n_obs, n_cells, replace=False)
        return self.adata[indices, :].copy()
    
    @memory_monitor
    def render_dimensionality_reduction(self):
        """Dimensionality reduction visualizations"""
        
        st.header("🗺️ Dimensionality Reduction")
        
        if self.plot_data is None:
            st.warning("⚠️ Επιλέξτε πρώτα δεδομένα")
            return
        
        # Επιλογή μεθόδου
        col1, col2 = st.columns(2)
        
        with col1:
            method = st.selectbox(
                "Μέθοδος Dimensionality Reduction:",
                ["UMAP", "t-SNE", "PCA", "Diffusion Map"]
            )
        
        with col2:
            # Επιλογή χρωματισμού
            color_options = ["None"] + list(self.plot_data.obs.columns)
            color_by = st.selectbox(
                "Χρωματισμός βάσει:",
                color_options
            )
        
        # Παράμετροι για κάθε μέθοδο
        if method in ["UMAP", "t-SNE"]:
            with st.expander("🔧 Παράμετροι"):
                if method == "UMAP":
                    n_neighbors = st.slider("N neighbors:", 5, 100, 15)
                    min_dist = st.slider("Min distance:", 0.1, 1.0, 0.5, step=0.1)
                else:  # t-SNE
                    perplexity = st.slider("Perplexity:", 5, 100, 30)
                    learning_rate = st.slider("Learning rate:", 10, 1000, 200)
        
        # Κουμπί υπολογισμού
        if st.button(f"🚀 Υπολογισμός {method}"):
            try:
                with st.spinner(f"Υπολογισμός {method}..."):
                    
                    if method == "UMAP":
                        if 'X_pca' not in self.plot_data.obsm:
                            sc.tl.pca(self.plot_data)
                        sc.pp.neighbors(self.plot_data, n_neighbors=n_neighbors)
                        sc.tl.umap(self.plot_data, min_dist=min_dist)
                        basis = 'umap'
                        
                    elif method == "t-SNE":
                        if 'X_pca' not in self.plot_data.obsm:
                            sc.tl.pca(self.plot_data)
                        sc.tl.tsne(self.plot_data, perplexity=perplexity, learning_rate=learning_rate)
                        basis = 'tsne'
                        
                    elif method == "PCA":
                        sc.tl.pca(self.plot_data)
                        basis = 'pca'
                        
                    else:  # Diffusion Map
                        if 'X_pca' not in self.plot_data.obsm:
                            sc.tl.pca(self.plot_data)
                        sc.pp.neighbors(self.plot_data)
                        sc.tl.diffmap(self.plot_data)
                        basis = 'diffmap'
                    
                    # Δημιουργία plot
                    self.create_embedding_plot(basis, color_by, method)
                    
            except Exception as e:
                st.error(f"❌ Σφάλμα στον υπολογισμό {method}: {str(e)}")
    
    def create_embedding_plot(self, basis, color_by, method_name):
        """Δημιουργία embedding plot"""
        
        # Εξαγωγή coordinates
        coords = self.plot_data.obsm[f'X_{basis}']
        
        # Δημιουργία DataFrame για plotting
        plot_df = pd.DataFrame({
            f'{method_name}_1': coords[:, 0],
            f'{method_name}_2': coords[:, 1]
        })
        
        # Προσθήκη χρωματισμού
        if color_by != "None":
            plot_df['color'] = self.plot_data.obs[color_by].values
        
        # Plotly figure
        if color_by != "None":
            if self.plot_data.obs[color_by].dtype == 'object' or self.plot_data.obs[color_by].dtype.name == 'category':
                # Categorical coloring
                fig = px.scatter(
                    plot_df,
                    x=f'{method_name}_1',
                    y=f'{method_name}_2',
                    color='color',
                    title=f'{method_name} Plot - Colored by {color_by}',
                    hover_data=['color']
                )
            else:
                # Continuous coloring
                fig = px.scatter(
                    plot_df,
                    x=f'{method_name}_1',
                    y=f'{method_name}_2',
                    color='color',
                    color_continuous_scale='viridis',
                    title=f'{method_name} Plot - Colored by {color_by}',
                    hover_data=['color']
                )
        else:
            # No coloring
            fig = px.scatter(
                plot_df,
                x=f'{method_name}_1',
                y=f'{method_name}_2',
                title=f'{method_name} Plot'
            )
        
        fig.update_layout(height=600)
        st.plotly_chart(fig, use_container_width=True)
        
        # Export option
        self.add_plot_export_option(fig, f"{method_name}_plot")
    
    def render_gene_expression(self):
        """Gene expression visualizations"""
        
        st.header("🧬 Gene Expression Plots")
        
        if self.plot_data is None:
            st.warning("⚠️ Επιλέξτε πρώτα δεδομένα")
            return
        
        # Επιλογή γονιδίων
        st.subheader("🔍 Επιλογή Γονιδίων")
        
        # Search box για γονίδια
        gene_search = st.text_input(
            "Αναζήτηση γονιδίων:",
            placeholder="Πληκτρολογήστε όνομα γονιδίου...",
            help="Μπορείτε να πληκτρολογήσετε μέρος του ονόματος"
        )
        
        if gene_search:
            # Φιλτράρισμα γονιδίων
            matching_genes = [g for g in self.plot_data.var_names if gene_search.upper() in g.upper()]
            
            if matching_genes:
                selected_genes = st.multiselect(
                    f"Βρέθηκαν {len(matching_genes)} γονίδια:",
                    matching_genes[:50],  # Εμφάνιση πρώτων 50
                    help="Επιλέξτε έως 10 γονίδια για visualization"
                )
                
                if len(selected_genes) > 10:
                    st.warning("⚠️ Επιλέξτε έως 10 γονίδια για βέλτιστη απόδοση")
                    selected_genes = selected_genes[:10]
                
                if selected_genes:
                    # Τύπος plot
                    plot_type = st.selectbox(
                        "Τύπος Plot:",
                        ["Dot Plot", "Violin Plot", "Expression on Embedding", "Heatmap"]
                    )
                    
                    # Grouping variable
                    if plot_type in ["Dot Plot", "Violin Plot", "Heatmap"]:
                        categorical_cols = [col for col in self.plot_data.obs.columns 
                                         if self.plot_data.obs[col].dtype == 'object' or 
                                         self.plot_data.obs[col].dtype.name == 'category']
                        
                        if categorical_cols:
                            group_by = st.selectbox(
                                "Ομαδοποίηση βάσει:",
                                categorical_cols
                            )
                        else:
                            st.warning("Δεν βρέθηκαν κατηγορικές μεταβλητές για ομαδοποίηση")
                            return
                    
                    # Δημιουργία plot
                    if st.button("🎨 Δημιουργία Plot"):
                        self.create_gene_expression_plot(selected_genes, plot_type, 
                                                       group_by if plot_type in ["Dot Plot", "Violin Plot", "Heatmap"] else None)
            else:
                st.info(f"Δεν βρέθηκαν γονίδια που να περιέχουν '{gene_search}'")
    
    def create_gene_expression_plot(self, genes, plot_type, group_by=None):
        """Δημιουργία gene expression plot"""
        
        try:
            if plot_type == "Dot Plot":
                self.create_dot_plot(genes, group_by)
            elif plot_type == "Violin Plot":
                self.create_violin_plot(genes, group_by)
            elif plot_type == "Expression on Embedding":
                self.create_expression_embedding(genes)
            elif plot_type == "Heatmap":
                self.create_expression_heatmap(genes, group_by)
                
        except Exception as e:
            st.error(f"❌ Σφάλμα στη δημιουργία {plot_type}: {str(e)}")
    
    def create_dot_plot(self, genes, group_by):
        """Δημιουργία dot plot"""
        
        # Υπολογισμός mean expression και percentage
        groups = self.plot_data.obs[group_by].unique()
        
        plot_data = []
        
        for group in groups:
            group_mask = self.plot_data.obs[group_by] == group
            group_cells = self.plot_data[group_mask, :]
            
            for gene in genes:
                if gene in self.plot_data.var_names:
                    gene_idx = self.plot_data.var_names.get_loc(gene)
                    
                    if hasattr(group_cells.X, 'toarray'):
                        expression = group_cells.X[:, gene_idx].toarray().flatten()
                    else:
                        expression = group_cells.X[:, gene_idx]
                    
                    mean_exp = np.mean(expression)
                    pct_exp = (expression > 0).sum() / len(expression) * 100
                    
                    plot_data.append({
                        'Gene': gene,
                        'Group': group,
                        'Mean_Expression': mean_exp,
                        'Pct_Expressed': pct_exp
                    })
        
        plot_df = pd.DataFrame(plot_data)
        
        # Plotly dot plot
        fig = px.scatter(
            plot_df,
            x='Gene',
            y='Group',
            size='Pct_Expressed',
            color='Mean_Expression',
            color_continuous_scale='Reds',
            title=f'Dot Plot - Gene Expression by {group_by}',
            hover_data=['Mean_Expression', 'Pct_Expressed']
        )
        
        fig.update_layout(height=400 + len(groups) * 30)
        st.plotly_chart(fig, use_container_width=True)
        
        self.add_plot_export_option(fig, "dot_plot")
    
    def create_violin_plot(self, genes, group_by):
        """Δημιουργία violin plot"""
        
        # Subplots για πολλαπλά γονίδια
        n_genes = len(genes)
        cols = min(3, n_genes)
        rows = (n_genes + cols - 1) // cols
        
        fig = make_subplots(
            rows=rows, cols=cols,
            subplot_titles=genes,
            vertical_spacing=0.1
        )
        
        for i, gene in enumerate(genes):
            row = i // cols + 1
            col = i % cols + 1
            
            if gene in self.plot_data.var_names:
                gene_idx = self.plot_data.var_names.get_loc(gene)
                
                if hasattr(self.plot_data.X, 'toarray'):
                    expression = self.plot_data.X[:, gene_idx].toarray().flatten()
                else:
                    expression = self.plot_data.X[:, gene_idx]
                
                # Δεδομένα για violin plot
                violin_data = pd.DataFrame({
                    'Expression': expression,
                    'Group': self.plot_data.obs[group_by].values
                })
                
                for group in violin_data['Group'].unique():
                    group_exp = violin_data[violin_data['Group'] == group]['Expression']
                    
                    fig.add_trace(
                        go.Violin(
                            y=group_exp,
                            name=str(group),
                            showlegend=(i == 0),  # Μόνο στο πρώτο subplot
                            legendgroup=str(group)
                        ),
                        row=row, col=col
                    )
        
        fig.update_layout(
            height=300 * rows,
            title_text=f"Violin Plots - Gene Expression by {group_by}"
        )
        
        st.plotly_chart(fig, use_container_width=True)
        self.add_plot_export_option(fig, "violin_plot")
    
    def create_expression_embedding(self, genes):
        """Expression overlay σε embedding"""
        
        # Έλεγχος για διαθέσιμα embeddings
        available_embeddings = []
        for emb in ['umap', 'tsne', 'pca']:
            if f'X_{emb}' in self.plot_data.obsm:
                available_embeddings.append(emb.upper())
        
        if not available_embeddings:
            st.warning("⚠️ Δεν βρέθηκαν διαθέσιμα embeddings. Εκτελέστε πρώτα dimensionality reduction.")
            return
        
        embedding = st.selectbox("Επιλέξτε embedding:", available_embeddings)
        basis = embedding.lower()
        
        # Subplots για πολλαπλά γονίδια
        n_genes = len(genes)
        cols = min(2, n_genes)
        rows = (n_genes + cols - 1) // cols
        
        fig = make_subplots(
            rows=rows, cols=cols,
            subplot_titles=genes
        )
        
        coords = self.plot_data.obsm[f'X_{basis}']
        
        for i, gene in enumerate(genes):
            row = i // cols + 1
            col = i % cols + 1
            
            if gene in self.plot_data.var_names:
                gene_idx = self.plot_data.var_names.get_loc(gene)
                
                if hasattr(self.plot_data.X, 'toarray'):
                    expression = self.plot_data.X[:, gene_idx].toarray().flatten()
                else:
                    expression = self.plot_data.X[:, gene_idx]
                
                fig.add_trace(
                    go.Scatter(
                        x=coords[:, 0],
                        y=coords[:, 1],
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
            height=400 * rows,
            title_text=f"Gene Expression on {embedding}"
        )
        
        st.plotly_chart(fig, use_container_width=True)
        self.add_plot_export_option(fig, f"expression_{embedding}")
    
    def create_expression_heatmap(self, genes, group_by):
        """Expression heatmap"""
        
        # Υπολογισμός mean expression ανά ομάδα
        groups = self.plot_data.obs[group_by].unique()
        heatmap_data = []
        
        for group in groups:
            group_mask = self.plot_data.obs[group_by] == group
            group_cells = self.plot_data[group_mask, :]
            
            group_means = []
            for gene in genes:
                if gene in self.plot_data.var_names:
                    gene_idx = self.plot_data.var_names.get_loc(gene)
                    
                    if hasattr(group_cells.X, 'toarray'):
                        expression = group_cells.X[:, gene_idx].toarray().flatten()
                    else:
                        expression = group_cells.X[:, gene_idx]
                    
                    group_means.append(np.mean(expression))
                else:
                    group_means.append(0)
            
            heatmap_data.append(group_means)
        
        # Plotly heatmap
        fig = go.Figure(data=go.Heatmap(
            z=heatmap_data,
            x=genes,
            y=groups,
            colorscale='viridis',
            colorbar=dict(title="Mean Expression")
        ))
        
        fig.update_layout(
            title=f'Expression Heatmap by {group_by}',
            xaxis_title='Genes',
            yaxis_title=group_by,
            height=200 + len(groups) * 30
        )
        
        st.plotly_chart(fig, use_container_width=True)
        self.add_plot_export_option(fig, "expression_heatmap")
    
    def render_quality_control(self):
        """Quality control visualizations"""
        
        st.header("📊 Quality Control Plots")
        
        if self.plot_data is None:
            st.warning("⚠️ Επιλέξτε πρώτα δεδομένα")
            return
        
        # QC metrics που θα εμφανιστούν
        qc_metrics = []
        
        # Standard QC metrics
        standard_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
        for metric in standard_metrics:
            if metric in self.plot_data.obs.columns:
                qc_metrics.append(metric)
        
        # Additional metrics
        additional_metrics = [col for col in self.plot_data.obs.columns 
                            if col.startswith(('pct_', 'n_', 'total_')) and col not in qc_metrics]
        qc_metrics.extend(additional_metrics[:5])  # Πρώτα 5 επιπλέον
        
        if not qc_metrics:
            st.warning("⚠️ Δεν βρέθηκαν QC metrics. Εκτελέστε πρώτα preprocessing.")
            return
        
        # Επιλογή metrics για visualization
        selected_metrics = st.multiselect(
            "Επιλέξτε QC metrics:",
            qc_metrics,
            default=qc_metrics[:3]
        )
        
        if selected_metrics:
            # Τύπος plot
            qc_plot_type = st.selectbox(
                "Τύπος QC Plot:",
                ["Histograms", "Violin Plots", "Scatter Plots", "Box Plots"]
            )
            
            if st.button("📊 Δημιουργία QC Plots"):
                self.create_qc_plots(selected_metrics, qc_plot_type)
    
    def create_qc_plots(self, metrics, plot_type):
        """Δημιουργία QC plots"""
        
        if plot_type == "Histograms":
            self.create_qc_histograms(metrics)
        elif plot_type == "Violin Plots":
            self.create_qc_violins(metrics)
        elif plot_type == "Scatter Plots":
            self.create_qc_scatter(metrics)
        else:  # Box Plots
            self.create_qc_boxplots(metrics)
    
    def create_qc_histograms(self, metrics):
        """QC histograms"""
        
        n_metrics = len(metrics)
        cols = min(2, n_metrics)
        rows = (n_metrics + cols - 1) // cols
        
        fig = make_subplots(
            rows=rows, cols=cols,
            subplot_titles=metrics
        )
        
        for i, metric in enumerate(metrics):
            row = i // cols + 1
            col = i % cols + 1
            
            values = self.plot_data.obs[metric].values
            
            fig.add_trace(
                go.Histogram(
                    x=values,
                    name=metric,
                    showlegend=False,
                    nbinsx=50
                ),
                row=row, col=col
            )
        
        fig.update_layout(
            height=300 * rows,
            title_text="QC Metrics - Histograms"
        )
        
        st.plotly_chart(fig, use_container_width=True)
        self.add_plot_export_option(fig, "qc_histograms")
    
    def create_qc_violins(self, metrics):
        """QC violin plots"""
        
        # Δημιουργία long format data
        plot_data = []
        for metric in metrics:
            values = self.plot_data.obs[metric].values
            plot_data.extend([{'Metric': metric, 'Value': val} for val in values])
        
        plot_df = pd.DataFrame(plot_data)
        
        fig = px.violin(
            plot_df,
            x='Metric',
            y='Value',
            title='QC Metrics - Violin Plots'
        )
        
        fig.update_layout(height=500)
        st.plotly_chart(fig, use_container_width=True)
        self.add_plot_export_option(fig, "qc_violins")
    
    def create_qc_scatter(self, metrics):
        """QC scatter plots"""
        
        if len(metrics) < 2:
            st.warning("⚠️ Χρειάζονται τουλάχιστον 2 metrics για scatter plot")
            return
        
        # Επιλογή x και y
        col1, col2 = st.columns(2)
        with col1:
            x_metric = st.selectbox("X-axis:", metrics, index=0)
        with col2:
            y_metric = st.selectbox("Y-axis:", metrics, index=1 if len(metrics) > 1 else 0)
        
        if x_metric != y_metric:
            fig = px.scatter(
                x=self.plot_data.obs[x_metric],
                y=self.plot_data.obs[y_metric],
                labels={'x': x_metric, 'y': y_metric},
                title=f'QC Scatter: {x_metric} vs {y_metric}'
            )
            
            fig.update_layout(height=500)
            st.plotly_chart(fig, use_container_width=True)
            self.add_plot_export_option(fig, "qc_scatter")
    
    def create_qc_boxplots(self, metrics):
        """QC box plots"""
        
        # Δημιουργία long format data
        plot_data = []
        for metric in metrics:
            values = self.plot_data.obs[metric].values
            plot_data.extend([{'Metric': metric, 'Value': val} for val in values])
        
        plot_df = pd.DataFrame(plot_data)
        
        fig = px.box(
            plot_df,
            x='Metric',
            y='Value',
            title='QC Metrics - Box Plots'
        )
        
        fig.update_layout(height=500)
        st.plotly_chart(fig, use_container_width=True)
        self.add_plot_export_option(fig, "qc_boxplots")
    
    def render_advanced_plots(self):
        """Advanced plotting options"""
        
        st.header("🎨 Advanced Plots")
        
        if self.plot_data is None:
            st.warning("⚠️ Επιλέξτε πρώτα δεδομένα")
            return
        
        # Τύποι advanced plots
        advanced_plot_type = st.selectbox(
            "Τύπος Advanced Plot:",
            ["Correlation Heatmap", "PCA Loadings", "Feature Plot", "Density Plot"]
        )
        
        if advanced_plot_type == "Correlation Heatmap":
            self.render_correlation_heatmap()
        elif advanced_plot_type == "PCA Loadings":
            self.render_pca_loadings()
        elif advanced_plot_type == "Feature Plot":
            self.render_feature_plot()
        else:  # Density Plot
            self.render_density_plot()
    
    def render_correlation_heatmap(self):
        """Correlation heatmap των QC metrics"""
        
        # Επιλογή numerical columns
        numerical_cols = [col for col in self.plot_data.obs.columns 
                         if self.plot_data.obs[col].dtype in ['int64', 'float64']]
        
        if len(numerical_cols) < 2:
            st.warning("⚠️ Χρειάζονται τουλάχιστον 2 numerical variables για correlation")
            return
        
        selected_cols = st.multiselect(
            "Επιλέξτε variables για correlation:",
            numerical_cols,
            default=numerical_cols[:min(10, len(numerical_cols))]
        )
        
        if len(selected_cols) >= 2 and st.button("📊 Δημιουργία Correlation Heatmap"):
            
            # Υπολογισμός correlation matrix
            corr_matrix = self.plot_data.obs[selected_cols].corr()
            
            fig = px.imshow(
                corr_matrix,
                title="Correlation Heatmap",
                color_continuous_scale='RdBu',
                aspect='auto'
            )
            
            fig.update_layout(height=600)
            st.plotly_chart(fig, use_container_width=True)
            self.add_plot_export_option(fig, "correlation_heatmap")
    
    def render_pca_loadings(self):
        """PCA loadings plot"""
        
        if 'X_pca' not in self.plot_data.obsm:
            st.warning("⚠️ Εκτελέστε πρώτα PCA")
            return
        
        if 'PCs' not in self.plot_data.varm:
            st.warning("⚠️ PCA loadings δεν είναι διαθέσιμα")
            return
        
        # Επιλογή PC components
        n_pcs = self.plot_data.varm['PCs'].shape[1]
        
        col1, col2 = st.columns(2)
        with col1:
            pc1 = st.selectbox("PC X-axis:", range(1, n_pcs + 1), index=0)
        with col2:
            pc2 = st.selectbox("PC Y-axis:", range(1, n_pcs + 1), index=1)
        
        if st.button("📊 Δημιουργία PCA Loadings Plot"):
            
            loadings = self.plot_data.varm['PCs']
            
            fig = px.scatter(
                x=loadings[:, pc1-1],
                y=loadings[:, pc2-1],
                hover_name=self.plot_data.var_names,
                labels={'x': f'PC{pc1}', 'y': f'PC{pc2}'},
                title=f'PCA Loadings: PC{pc1} vs PC{pc2}'
            )
            
            fig.update_layout(height=600)
            st.plotly_chart(fig, use_container_width=True)
            self.add_plot_export_option(fig, "pca_loadings")
    
    def render_feature_plot(self):
        """Feature plot για συγκεκριμένο feature"""
        
        # Επιλογή feature
        all_features = list(self.plot_data.obs.columns) + list(self.plot_data.var_names)
        
        feature = st.selectbox(
            "Επιλέξτε feature:",
            all_features
        )
        
        if st.button("📊 Δημιουργία Feature Plot"):
            
            if feature in self.plot_data.obs.columns:
                # Metadata feature
                values = self.plot_data.obs[feature].values
            else:
                # Gene feature
                gene_idx = self.plot_data.var_names.get_loc(feature)
                if hasattr(self.plot_data.X, 'toarray'):
                    values = self.plot_data.X[:, gene_idx].toarray().flatten()
                else:
                    values = self.plot_data.X[:, gene_idx]
            
            # Histogram
            fig = px.histogram(
                x=values,
                title=f'Distribution of {feature}',
                labels={'x': feature, 'y': 'Count'}
            )
            
            fig.update_layout(height=400)
            st.plotly_chart(fig, use_container_width=True)
            self.add_plot_export_option(fig, f"feature_{feature}")
    
    def render_density_plot(self):
        """Density plot"""
        
        # Έλεγχος για embeddings
        available_embeddings = []
        for emb in ['umap', 'tsne', 'pca']:
            if f'X_{emb}' in self.plot_data.obsm:
                available_embeddings.append(emb.upper())
        
        if not available_embeddings:
            st.warning("⚠️ Δεν βρέθηκαν διαθέσιμα embeddings για density plot")
            return
        
        embedding = st.selectbox("Επιλέξτε embedding:", available_embeddings)
        
        if st.button("📊 Δημιουργία Density Plot"):
            
            coords = self.plot_data.obsm[f'X_{embedding.lower()}']
            
            fig = px.density_contour(
                x=coords[:, 0],
                y=coords[:, 1],
                title=f'{embedding} Density Plot'
            )
            
            fig.update_layout(height=600)
            st.plotly_chart(fig, use_container_width=True)
            self.add_plot_export_option(fig, f"density_{embedding}")
    
    def add_plot_export_option(self, fig, filename_prefix):
        """Προσθήκη export option για plot"""
        
        if st.button(f"💾 Export {filename_prefix}", key=f"export_{filename_prefix}"):
            
            # Export ως HTML
            html_str = fig.to_html()
            
            st.download_button(
                label="📥 Download HTML",
                data=html_str,
                file_name=f"{filename_prefix}.html",
                mime="text/html"
            )
            
            # Export ως PNG (αν είναι δυνατό)
            try:
                img_bytes = fig.to_image(format="png", width=1200, height=800)
                st.download_button(
                    label="📥 Download PNG",
                    data=img_bytes,
                    file_name=f"{filename_prefix}.png",
                    mime="image/png"
                )
            except Exception:
                st.info("PNG export απαιτεί kaleido package")
