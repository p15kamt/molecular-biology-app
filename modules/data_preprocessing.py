"""
Module για την Προεπεξεργασία Δεδομένων scRNA-seq

Αυτό το module περιλαμβάνει όλες τις λειτουργίες που απαιτούνται για την
προεπεξεργασία δεδομένων single-cell RNA sequencing, συμπεριλαμβανομένου
του quality control, filtering, και normalization.


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
import tempfile
import os
import sys
from pathlib import Path
import warnings

# Προσθήκη utils directory στο path
sys.path.append(str(Path(__file__).parent.parent / "utils"))

try:
    from memory_utils import (
        memory_monitor, display_memory_info, 
        optimize_adata_memory, suggest_optimization,
        cleanup_memory
    )
    from advanced_memory import memory_manager
    from robust_error_handler import safe_execute, memory_safe_operation, emergency_cleanup
except ImportError:
    # Fallback functions αν δεν υπάρχει το memory_utils
    def memory_monitor(func):
        return func
    def display_memory_info():
        pass
    def optimize_adata_memory(adata):
        return adata
    def suggest_optimization(adata):
        return []
    def cleanup_memory():
        pass
    
    # Fallback για robust error handler
    def safe_execute(func_name="Λειτουργία"):
        def decorator(func):
            return func
        return decorator
    
    def memory_safe_operation(operation_name="Λειτουργία"):
        class DummyContext:
            def __enter__(self): return self
            def __exit__(self, *args): return False
        return DummyContext()
    
    def emergency_cleanup():
        pass
    
    # Fallback για advanced memory manager
    class FallbackMemoryManager:
        def get_system_memory_info(self):
            return {'total_gb': 8, 'available_gb': 4, 'process_mb': 500}
        def should_use_streaming(self, adata):
            return False
        def optimize_adata_for_memory(self, adata):
            return adata
        def safe_display_data(self, adata, max_cells=100, max_genes=100):
            return pd.DataFrame({"Info": ["Fallback mode"]})
        def progressive_qc_calculation(self, adata):
            return adata
        def cleanup_session_memory(self):
            pass
    
    memory_manager = FallbackMemoryManager()

# Απόκρυψη warnings για καθαρότερη εμφάνιση
warnings.filterwarnings('ignore')

# Ρύθμιση scanpy
sc.settings.verbosity = 3  # Επίπεδο verbosity
sc.settings.set_figure_params(dpi=80, facecolor='white')

class DataPreprocessingPage:
    """Κλάση για τη σελίδα προεπεξεργασίας δεδομένων"""
    
    def __init__(self):
        self.adata = None
        self.original_adata = None
        
    def render(self):
        """Κεντρική μέθοδος για την εμφάνιση της σελίδας"""
        
        st.title("📊 Προεπεξεργασία Δεδομένων scRNA-seq")
        st.markdown("### Ανεβάστε και προεπεξεργαστείτε τα δεδομένα σας")
        
        # Tabs για οργάνωση της λειτουργικότητας
        tab1, tab2, tab3, tab4 = st.tabs([
            "📁 Ανέβασμα Δεδομένων", 
            "🔍 Quality Control", 
            "🧹 Filtering & Preprocessing",
            "📊 Αποτελέσματα"
        ])
        
        with tab1:
            self._render_upload_section()
            
        with tab2:
            if self.adata is not None:
                self._render_quality_control()
            else:
                st.info("Παρακαλώ ανεβάστε πρώτα δεδομένα στο tab 'Ανέβασμα Δεδομένων'")
                
        with tab3:
            if self.adata is not None:
                self._render_preprocessing_section()
            else:
                st.info("Παρακαλώ ανεβάστε πρώτα δεδομένα στο tab 'Ανέβασμα Δεδομένων'")
                
        with tab4:
            if self.adata is not None:
                self._render_results_section()
            else:
                st.info("Παρακαλώ ανεβάστε πρώτα δεδομένα στο tab 'Ανέβασμα Δεδομένων'")
    
    def _render_upload_section(self):
        """Τμήμα ανεβάσματος αρχείων"""
        
        st.markdown("#### 📤 Ανέβασμα Αρχείου Δεδομένων")
        
        # Οδηγίες για το ανέβασμα
        with st.expander("📖 Υποστηριζόμενα Formats & Οδηγίες"):
            st.markdown("""
            **Υποστηριζόμενα formats:**
            - **H5AD**: AnnData objects (προτιμώμενο για scRNA-seq)
            - **CSV**: Gene expression matrix με cells στις γραμμές, genes στις στήλες
            - **TSV**: Tab-separated values
            - **Excel**: .xlsx files
            
            **Προτεινόμενη δομή δεδομένων:**
            - Γραμμές: Κύτταρα (cells)
            - Στήλες: Γονίδια (genes)
            - Τιμές: Αριθμός reads/counts
            
            **Σημείωση:** Για καλύτερα αποτελέσματα, χρησιμοποιήστε H5AD format που περιέχει metadata.
            """)
        
        # File uploader με υποστήριξη μεγάλων αρχείων
        uploaded_file = st.file_uploader(
            "Επιλέξτε αρχείο δεδομένων",
            type=['h5ad', 'csv', 'tsv', 'xlsx', 'xls'],
            help="Υποστηριζόμενα formats: H5AD, CSV, TSV, Excel. Μέγιστο μέγεθος: 1GB"
        )
        
        # Πληροφορίες για μεγάλα αρχεία και memory management
        col1, col2 = st.columns([3, 1])
        
        with col1:
            with st.expander("💡 Tips για Μεγάλα Αρχεία"):
                st.markdown("""
                **Για αρχεία >100MB συνιστάται:**
                - **H5AD format** για καλύτερη συμπίεση και ταχύτητα
                - **Chunked processing** θα ενεργοποιηθεί αυτόματα
                - **Memory-efficient operations** για μείωση RAM usage
                - **Progressive loading** με πρόοδο στην οθόνη
                
                **Συμβουλές:**
                - Κλείστε άλλες εφαρμογές για περισσότερη μνήμη
                - Χρησιμοποιήστε Subsample option για γρήγορο testing
                - Περιμένετε 2-5 λεπτά για πολύ μεγάλα αρχεία
                """)
        
        with col2:
            st.markdown("#### 🧹 Memory Control")
            if st.button("🗑️ Cleanup Memory", help="Καθαρισμός μνήμης για καλύτερη απόδοση"):
                memory_manager.cleanup_session_memory()
                st.experimental_rerun()
            
            # Εμφάνιση memory status
            memory_info = memory_manager.get_system_memory_info()
            st.metric(
                "RAM Usage", 
                f"{memory_info['process_mb']:.0f}MB",
                help="Τρέχουσα χρήση μνήμης από την εφαρμογή"
            )
        
        if uploaded_file is not None:
            with st.spinner("Φόρτωση δεδομένων..."):
                success = self._load_data(uploaded_file)
                
                if success:
                    st.success("✅ Τα δεδομένα φορτώθηκαν επιτυχώς!")
                    self._display_data_summary()
                else:
                    st.error("❌ Σφάλμα κατά τη φόρτωση των δεδομένων.")
    
    def _load_data(self, uploaded_file):
        """Φόρτωση δεδομένων από διάφορα formats με memory optimization"""
        
        try:
            file_extension = uploaded_file.name.split('.')[-1].lower()
            file_size_mb = uploaded_file.size / (1024 * 1024)
            
            # Εμφάνιση μεγέθους αρχείου
            st.info(f"📁 Μέγεθος αρχείου: {file_size_mb:.1f} MB")
            
            # Memory-efficient file writing για μεγάλα αρχεία
            with tempfile.NamedTemporaryFile(delete=False, suffix=f'.{file_extension}') as tmp_file:
                if file_size_mb > 100:  # Για αρχεία >100MB
                    st.info("🔄 Μεγάλο αρχείο - χρήση chunked loading...")
                    progress_bar = st.progress(0)
                    
                    # Chunked reading
                    chunk_size = 8192  # 8KB chunks
                    bytes_written = 0
                    
                    uploaded_file.seek(0)  # Reset file position
                    while True:
                        chunk = uploaded_file.read(chunk_size)
                        if not chunk:
                            break
                        tmp_file.write(chunk)
                        bytes_written += len(chunk)
                        
                        # Update progress
                        progress = min(bytes_written / uploaded_file.size, 1.0)
                        progress_bar.progress(progress)
                    
                    progress_bar.empty()
                else:
                    # Κανονική φόρτωση για μικρά αρχεία
                    tmp_file.write(uploaded_file.read())
                
                tmp_path = tmp_file.name
            
            if file_extension == 'h5ad':
                self.adata = sc.read_h5ad(tmp_path)
                
            elif file_extension in ['csv', 'tsv']:
                separator = ',' if file_extension == 'csv' else '\t'
                df = pd.read_csv(tmp_path, sep=separator, index_col=0)
                
                # Διαχωρισμός numeric από categorical columns
                numeric_cols = []
                categorical_cols = []
                
                for col in df.columns:
                    try:
                        # Δοκιμή μετατροπής σε numeric
                        pd.to_numeric(df[col], errors='raise')
                        numeric_cols.append(col)
                    except (ValueError, TypeError):
                        categorical_cols.append(col)
                
                # Δημιουργία AnnData object
                if len(numeric_cols) > 0:
                    # Χρήση μόνο των numeric columns για το X matrix
                    X_data = df[numeric_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
                    self.adata = sc.AnnData(X_data)
                    
                    # Προσθήκη categorical data στα obs
                    if len(categorical_cols) > 0:
                        for col in categorical_cols:
                            self.adata.obs[col] = df[col].values
                else:
                    # Αν δεν υπάρχουν numeric columns, μετατροπή όλων
                    X_data = df.apply(pd.to_numeric, errors='coerce').fillna(0)
                    self.adata = sc.AnnData(X_data)
                
            elif file_extension in ['xlsx', 'xls']:
                df = pd.read_excel(tmp_path, index_col=0)
                
                # Εφαρμογή της ίδιας λογικής όπως στα CSV
                numeric_cols = []
                categorical_cols = []
                
                for col in df.columns:
                    try:
                        pd.to_numeric(df[col], errors='raise')
                        numeric_cols.append(col)
                    except (ValueError, TypeError):
                        categorical_cols.append(col)
                
                if len(numeric_cols) > 0:
                    X_data = df[numeric_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
                    self.adata = sc.AnnData(X_data)
                    
                    if len(categorical_cols) > 0:
                        for col in categorical_cols:
                            self.adata.obs[col] = df[col].values
                else:
                    X_data = df.apply(pd.to_numeric, errors='coerce').fillna(0)
                    self.adata = sc.AnnData(X_data)
            
            # Advanced memory optimization
            should_stream = memory_manager.should_use_streaming(self.adata)
            if should_stream:
                st.warning("⚡ Μεγάλο dataset - ενεργοποίηση streaming mode")
            
            # Memory optimization
            self.adata = memory_manager.optimize_adata_for_memory(self.adata)
            
            # Αντιγραφή των αρχικών δεδομένων (μόνο αν δεν είναι πολύ μεγάλο)
            estimated_size = memory_manager.get_system_memory_info()['process_mb']
            if estimated_size < 1000:  # <1GB
                self.original_adata = self.adata.copy()
            else:
                st.info("💾 Παράλειψη backup copy λόγω μεγέθους - χρησιμοποιήστε reload αν χρειάζεται reset")
                self.original_adata = None
            
            # Καθαρισμός προσωρινού αρχείου
            os.unlink(tmp_path)
            
            # Βασικός έλεγχος δεδομένων
            if self.adata.n_obs == 0 or self.adata.n_vars == 0:
                st.error("Τα δεδομένα φαίνονται κενά ή έχουν λάθος format.")
                return False
            
            # Επιλογή subsampling για μεγάλα datasets
            if self.adata.n_obs > 10000:  # Για >10K κύτταρα
                st.warning(f"⚠️ Μεγάλο dataset: {self.adata.n_obs:,} κύτταρα")
                
                use_subsample = st.checkbox(
                    "🎯 Χρήση Subsample για γρηγορότερη ανάλυση",
                    value=True,
                    help="Συνιστάται για γρήγορο testing με μεγάλα datasets"
                )
                
                if use_subsample:
                    sample_size = st.slider(
                        "Αριθμός κυττάρων για subsample:",
                        min_value=1000,
                        max_value=min(50000, self.adata.n_obs),
                        value=min(5000, self.adata.n_obs),
                        step=1000,
                        help="Μικρότερο sample = γρηγορότερη ανάλυση"
                    )
                    
                    if st.button("✂️ Εφαρμογή Subsample"):
                        # Random subsample
                        sc.pp.subsample(self.adata, n_obs=sample_size)
                        st.success(f"✅ Subsample ολοκληρώθηκε: {self.adata.n_obs:,} κύτταρα")
                        st.experimental_rerun()
                
            return True
            
        except Exception as e:
            st.error(f"Σφάλμα κατά τη φόρτωση: {str(e)}")
            return False
    
    def _display_data_summary(self):
        """Εμφάνιση σύνοψης των δεδομένων"""
        
        st.markdown("#### 📋 Σύνοψη Δεδομένων")
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Αριθμός Κυττάρων", f"{self.adata.n_obs:,}")
        with col2:
            st.metric("Αριθμός Γονιδίων", f"{self.adata.n_vars:,}")
        with col3:
            try:
                # Υπολογισμός συνολικών counts με ασφαλή τρόπο
                total_sum = self.adata.X.sum()
                if hasattr(total_sum, 'item'):  # Για numpy arrays/matrices
                    total_counts = int(total_sum.item())
                else:
                    total_counts = int(total_sum)
                display_counts = f"{total_counts:,}"
            except (TypeError, ValueError):
                display_counts = "N/A"
            st.metric("Συνολικά Counts", display_counts)
        with col4:
            try:
                # Υπολογισμός sparsity με ασφαλή τρόπο
                non_zero_count = (self.adata.X > 0).sum()
                if hasattr(non_zero_count, 'item'):
                    non_zero_count = non_zero_count.item()
                total_elements = self.adata.n_obs * self.adata.n_vars
                sparsity = 1 - (non_zero_count / total_elements)
                display_sparsity = f"{sparsity:.1%}"
            except (TypeError, ValueError, ZeroDivisionError):
                display_sparsity = "N/A"
            st.metric("Sparsity", display_sparsity)
        
        # Εμφάνιση sample των δεδομένων
        st.markdown("#### 👀 Preview Δεδομένων")
        
        # Memory-safe εμφάνιση δεδομένων με advanced memory manager
        sample_data = memory_manager.safe_display_data(
            self.adata, 
            max_cells=10, 
            max_genes=10
        )
        
        st.dataframe(sample_data, use_container_width=True)
        
        # Πληροφορίες metadata
        if len(self.adata.obs.columns) > 0:
            st.markdown("#### 🏷️ Διαθέσιμα Metadata (Observations)")
            st.write(list(self.adata.obs.columns))
            
        if len(self.adata.var.columns) > 0:
            st.markdown("#### 🧬 Διαθέσιμα Metadata (Variables/Genes)")
            st.write(list(self.adata.var.columns))
    
    def _render_quality_control(self):
        """Τμήμα quality control"""
        
        st.markdown("#### 🔍 Quality Control Metrics")
        
        # Υπολογισμός QC metrics
        with st.spinner("Υπολογισμός QC metrics..."):
            self._calculate_qc_metrics()
        
        # Εμφάνιση histograms
        self._plot_qc_histograms()
        
        # Εμφάνιση scatter plots
        self._plot_qc_scatter()
        
        # Στατιστικά
        self._display_qc_statistics()
    
    def _calculate_qc_metrics(self):
        """Υπολογισμός quality control metrics με memory protection"""
        
        import numpy as np
        
        # Έλεγχος τύπου dataset
        if self.adata.n_obs < 1000 and self.adata.n_vars < 100:
            st.warning("🔍 Μικρό dataset - χρήση απλοποιημένων QC metrics")
            self._calculate_simple_qc_metrics()
            return
        
        # Έλεγχος μνήμης πριν τον υπολογισμό
        memory_info = memory_manager.get_system_memory_info()
        dataset_size_mb = memory_manager.estimate_memory_usage(self.adata)
        
        st.info(f"💾 Υπολογισμός QC metrics για dataset {dataset_size_mb:.1f}MB...")
        
        # Αν το dataset είναι πολύ μεγάλο, χρήση βελτιωμένης calculation
        if dataset_size_mb > 300 or memory_info['available_gb'] < 3:
            st.info("⚡ Μεγάλο dataset - χρήση optimized QC calculation...")
            self.adata = memory_manager.progressive_qc_calculation(self.adata)
            return
        
        try:
            # Βεβαίωση ότι το X matrix είναι numeric
            if hasattr(self.adata.X, 'astype'):
                try:
                    # Μετατροπή σε numeric αν δεν είναι ήδη
                    if not np.issubdtype(self.adata.X.dtype, np.number):
                        self.adata.X = self.adata.X.astype(np.float32)  # float32 αντί για float64
                except (ValueError, TypeError):
                    # Αν αποτύχει, προσπάθεια με pandas
                    import pandas as pd
                    st.info("🔄 Μετατροπή δεδομένων σε numeric format...")
                    df_temp = pd.DataFrame(self.adata.X, 
                                         index=self.adata.obs_names, 
                                         columns=self.adata.var_names)
                    df_numeric = df_temp.apply(pd.to_numeric, errors='coerce').fillna(0)
                    self.adata.X = df_numeric.values.astype(np.float32)
            
            # Memory-safe υπολογισμοί
            with st.spinner("Υπολογισμός QC metrics..."):
                
                # Αριθμός γονιδίων ανά κύτταρο (memory-safe)
                if hasattr(self.adata.X, 'toarray'):
                    # Sparse matrix - chunked calculation
                    n_genes = np.zeros(self.adata.n_obs)
                    chunk_size = min(1000, self.adata.n_obs)
                    
                    for i in range(0, self.adata.n_obs, chunk_size):
                        end_idx = min(i + chunk_size, self.adata.n_obs)
                        chunk = self.adata.X[i:end_idx, :]
                        n_genes[i:end_idx] = (chunk > 0).sum(axis=1).A1
                    
                    self.adata.obs['n_genes'] = n_genes
                else:
                    # Dense matrix
                    n_genes = (self.adata.X > 0).sum(axis=1)
                    self.adata.obs['n_genes'] = np.array(n_genes).flatten()
                
                # Συνολικά counts ανά κύτταρο (memory-safe)
                if hasattr(self.adata.X, 'toarray'):
                    # Sparse matrix - chunked calculation
                    total_counts = np.zeros(self.adata.n_obs)
                    
                    for i in range(0, self.adata.n_obs, chunk_size):
                        end_idx = min(i + chunk_size, self.adata.n_obs)
                        chunk = self.adata.X[i:end_idx, :]
                        total_counts[i:end_idx] = chunk.sum(axis=1).A1
                    
                    self.adata.obs['total_counts'] = total_counts
                else:
                    # Dense matrix
                    total_counts = self.adata.X.sum(axis=1)
                    self.adata.obs['total_counts'] = np.array(total_counts).flatten()
                
        except MemoryError:
            st.error("❌ Memory Error - χρήση progressive calculation...")
            self.adata = memory_manager.progressive_qc_calculation(self.adata)
            return
        except Exception as e:
            st.error(f"❌ Σφάλμα στον υπολογισμό QC: {str(e)}")
            st.info("🔄 Χρήση εναλλακτικής μεθόδου...")
            self.adata = memory_manager.progressive_qc_calculation(self.adata)
            return
        
        # Υπολογισμός μιτοχονδριακών και ριβοσωμικών γονιδίων με memory protection
        try:
            # Μιτοχονδριακά γονίδια
            mito_genes = self.adata.var_names.str.startswith('MT-') | self.adata.var_names.str.startswith('mt-')
            if mito_genes.any():
                if hasattr(self.adata.X, 'toarray'):
                    # Sparse matrix - memory-safe calculation
                    mito_subset = self.adata[:, mito_genes]
                    mito_counts = np.array(mito_subset.X.sum(axis=1)).flatten()
                else:
                    mito_counts = self.adata[:, mito_genes].X.sum(axis=1)
                    mito_counts = np.array(mito_counts).flatten()
                
                total_counts_safe = np.array(self.adata.obs['total_counts'])
                # Αποφυγή διαίρεσης με 0
                pct_mt = np.where(total_counts_safe > 0, 
                                 (mito_counts / total_counts_safe) * 100, 0)
                self.adata.obs['pct_counts_mt'] = pct_mt
            else:
                self.adata.obs['pct_counts_mt'] = np.zeros(self.adata.n_obs)
            
            # Ριβοσωμικά γονίδια
            ribo_genes = self.adata.var_names.str.startswith('RPS') | self.adata.var_names.str.startswith('RPL')
            if ribo_genes.any():
                if hasattr(self.adata.X, 'toarray'):
                    # Sparse matrix - memory-safe calculation
                    ribo_subset = self.adata[:, ribo_genes]
                    ribo_counts = np.array(ribo_subset.X.sum(axis=1)).flatten()
                else:
                    ribo_counts = self.adata[:, ribo_genes].X.sum(axis=1)
                    ribo_counts = np.array(ribo_counts).flatten()
                
                total_counts_safe = np.array(self.adata.obs['total_counts'])
                # Αποφυγή διαίρεσης με 0
                pct_rb = np.where(total_counts_safe > 0, 
                                 (ribo_counts / total_counts_safe) * 100, 0)
                self.adata.obs['pct_counts_ribo'] = pct_rb
            else:
                self.adata.obs['pct_counts_ribo'] = np.zeros(self.adata.n_obs)
            
            # Αριθμός κυττάρων ανά γονίδιο (memory-safe)
            if hasattr(self.adata.X, 'toarray'):
                # Sparse matrix - avoid full conversion
                n_cells = np.array((self.adata.X > 0).sum(axis=0)).flatten()
            else:
                n_cells = (self.adata.X > 0).sum(axis=0)
                n_cells = np.array(n_cells).flatten()
            
            self.adata.var['n_cells'] = n_cells
            
        except MemoryError:
            st.error("❌ Memory Error στον υπολογισμό secondary metrics - χρήση fallback...")
            # Fallback - βασικά metrics μόνο
            if 'n_genes' not in self.adata.obs.columns:
                self.adata.obs['n_genes'] = 1000  # Default value
            if 'total_counts' not in self.adata.obs.columns:
                self.adata.obs['total_counts'] = 5000  # Default value
            self.adata.obs['pct_counts_mt'] = np.zeros(self.adata.n_obs)
            self.adata.obs['pct_counts_ribo'] = np.zeros(self.adata.n_obs)
    
    def _calculate_simple_qc_metrics(self):
        """Απλοποιημένος υπολογισμός QC για μη-scRNA-seq datasets"""
        
        st.info("🔧 Υπολογισμός απλοποιημένων QC metrics...")
        
        # Για μικρά datasets όπως Iris, χρησιμοποιούμε τον αριθμό features ως "genes"
        self.adata.obs['n_genes'] = self.adata.n_vars  # Όλα τα features
        
        # Συνολικά "counts" = άθροισμα όλων των τιμών ανά sample
        if hasattr(self.adata.X, 'toarray'):
            total_counts = np.array(self.adata.X.sum(axis=1)).flatten()
        else:
            total_counts = self.adata.X.sum(axis=1)
        
        self.adata.obs['total_counts'] = total_counts
        
        # Για μη-scRNA-seq data, δεν υπάρχουν μιτοχονδριακά γονίδια
        self.adata.obs['pct_counts_mt'] = np.zeros(self.adata.n_obs)
        self.adata.obs['pct_counts_ribo'] = np.zeros(self.adata.n_obs)
        
        st.success("✅ Απλοποιημένα QC metrics υπολογίστηκαν!")
        st.info(f"📊 n_genes: {self.adata.obs['n_genes'].iloc[0]} (features)")
        st.info(f"📊 total_counts range: {total_counts.min():.1f} - {total_counts.max():.1f}")
    
    def _plot_qc_histograms(self):
        """Δημιουργία histograms για QC metrics με memory protection"""
        
        st.markdown("##### 📊 Κατανομές QC Metrics")
        
        try:
            # Έλεγχος αν τα QC metrics υπάρχουν
            required_metrics = ['n_genes', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo']
            missing_metrics = [m for m in required_metrics if m not in self.adata.obs.columns]
            
            if missing_metrics:
                st.warning(f"⚠️ Λείπουν metrics: {', '.join(missing_metrics)} - υπολογισμός...")
                self._calculate_qc_metrics()
            
            # Memory-safe plotting με subsample αν χρειάζεται
            plot_data = self.adata.obs
            
            # Αν πολλά κύτταρα, χρήση subsample για plots
            if len(plot_data) > 10000:
                st.info("📊 Subsample για visualization (10K κύτταρα)")
                sample_indices = np.random.choice(len(plot_data), 10000, replace=False)
                plot_data = plot_data.iloc[sample_indices]
            
            # Δημιουργία subplots
            fig = make_subplots(
                rows=2, cols=2,
                subplot_titles=[
                    'Αριθμός Γονιδίων ανά Κύτταρο',
                    'Συνολικά Counts ανά Κύτταρο', 
                    'Ποσοστό Μιτοχονδριακών Γονιδίων',
                    'Ποσοστό Ριβοσωμικών Γονιδίων'
                ]
            )
            
            # Histogram 1: Αριθμός γονιδίων
            fig.add_trace(
                go.Histogram(x=plot_data['n_genes'], nbinsx=50, name='n_genes'),
                row=1, col=1
            )
            
            # Histogram 2: Συνολικά counts
            fig.add_trace(
                go.Histogram(x=plot_data['total_counts'], nbinsx=50, name='total_counts'),
                row=1, col=2
            )
            
            # Histogram 3: Μιτοχονδριακά
            fig.add_trace(
                go.Histogram(x=plot_data['pct_counts_mt'], nbinsx=50, name='pct_mt'),
                row=2, col=1
            )
            
            # Histogram 4: Ριβοσωμικά
            fig.add_trace(
                go.Histogram(x=plot_data['pct_counts_ribo'], nbinsx=50, name='pct_ribo'),
                row=2, col=2
            )
            
            fig.update_layout(height=600, showlegend=False, title_text="Quality Control Metrics")
            st.plotly_chart(fig, use_container_width=True)
            
        except Exception as e:
            st.error(f"❌ Σφάλμα στην οπτικοποίηση QC: {str(e)}")
            st.info("💡 Δοκιμάστε με μικρότερο dataset ή κάντε restart της εφαρμογής")
    
    def _plot_qc_scatter(self):
        """Δημιουργία scatter plots για QC metrics με memory protection"""
        
        st.markdown("##### 🔍 Συσχετίσεις QC Metrics")
        
        try:
            # Memory-safe plotting με subsample αν χρειάζεται
            plot_data = self.adata.obs
            
            # Αν πολλά κύτταρα, χρήση subsample για plots
            if len(plot_data) > 10000:
                st.info("📊 Subsample για scatter plots (10K κύτταρα)")
                sample_indices = np.random.choice(len(plot_data), 10000, replace=False)
                plot_data = plot_data.iloc[sample_indices]
            
            col1, col2 = st.columns(2)
            
            with col1:
                # Scatter plot: total_counts vs n_genes
                fig1 = px.scatter(
                    x=plot_data['total_counts'],
                    y=plot_data['n_genes'],
                    title='Συνολικά Counts vs Αριθμός Γονιδίων',
                    labels={'x': 'Συνολικά Counts', 'y': 'Αριθμός Γονιδίων'}
                )
                st.plotly_chart(fig1, use_container_width=True)
            
            with col2:
                # Scatter plot: total_counts vs pct_mt
                fig2 = px.scatter(
                    x=plot_data['total_counts'],
                    y=plot_data['pct_counts_mt'],
                    title='Συνολικά Counts vs % Μιτοχονδριακών',
                    labels={'x': 'Συνολικά Counts', 'y': '% Μιτοχονδριακών Γονιδίων'}
                )
                st.plotly_chart(fig2, use_container_width=True)
                
        except Exception as e:
            st.error(f"❌ Σφάλμα στα scatter plots: {str(e)}")
            st.info("💡 Χρήση εναλλακτικής απλής οπτικοποίησης...")
            
            # Fallback - απλά box plots
            try:
                fig = go.Figure()
                fig.add_trace(go.Box(y=plot_data['n_genes'], name='N Genes'))
                fig.update_layout(title="Απλοποιημένο QC Plot")
                st.plotly_chart(fig, use_container_width=True)
            except:
                st.warning("⚠️ Δεν είναι δυνατή η οπτικοποίηση - πολύ μεγάλο dataset")
    
    def _display_qc_statistics(self):
        """Εμφάνιση στατιστικών QC"""
        
        st.markdown("##### 📈 Στατιστικά QC Metrics")
        
        qc_stats = self.adata.obs[['n_genes', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo']].describe()
        st.dataframe(qc_stats, use_container_width=True)
    
    def _render_preprocessing_section(self):
        """Τμήμα preprocessing και filtering"""
        
        st.markdown("#### 🧹 Filtering & Preprocessing")
        
        # Παράμετροι filtering στο sidebar
        st.sidebar.markdown("### ⚙️ Παράμετροι Filtering")
        
        # Filtering κυττάρων
        st.sidebar.markdown("#### Filtering Κυττάρων")
        # Smart defaults για διαφορετικούς τύπους datasets
        if self.adata.n_obs < 1000 and self.adata.n_vars < 100:
            st.sidebar.info("🔍 Μικρό dataset - χαλαρές παράμετροι")
            default_min_genes = 1
        else:
            default_min_genes = 200
            
        min_genes = st.sidebar.number_input(
            "Ελάχιστος αριθμός γονιδίων ανά κύτταρο", 
            min_value=0, max_value=10000, value=default_min_genes, step=50
        )
        max_genes = st.sidebar.number_input(
            "Μέγιστος αριθμός γονιδίων ανά κύτταρο", 
            min_value=1000, max_value=20000, value=5000, step=100
        )
        max_counts = st.sidebar.number_input(
            "Μέγιστα συνολικά counts ανά κύτταρο", 
            min_value=1000, max_value=100000, value=30000, step=1000
        )
        max_mt = st.sidebar.slider(
            "Μέγιστο % μιτοχονδριακών γονιδίων", 
            min_value=0.0, max_value=50.0, value=20.0, step=0.5
        )
        
        # Filtering γονιδίων
        st.sidebar.markdown("#### Filtering Γονιδίων")
        min_cells = st.sidebar.number_input(
            "Ελάχιστος αριθμός κυττάρων που εκφράζουν το γονίδιο", 
            min_value=1, max_value=100, value=3, step=1
        )
        
        # Κουμπί εκτέλεσης
        if st.button("🚀 Εκτέλεση Preprocessing", type="primary"):
            with st.spinner("Εκτέλεση preprocessing..."):
                self._apply_preprocessing(min_genes, max_genes, max_counts, max_mt, min_cells)
    
    def _apply_preprocessing(self, min_genes, max_genes, max_counts, max_mt, min_cells):
        """Εφαρμογή preprocessing steps"""
        
        initial_cells = self.adata.n_obs
        initial_genes = self.adata.n_vars
        
        st.markdown("#### 📋 Preprocessing Steps")
        
        # Step 1: Filtering κυττάρων
        st.write("**Βήμα 1: Filtering κυττάρων**")
        
        # Filter by number of genes
        cell_filter = (
            (self.adata.obs['n_genes'] >= min_genes) &
            (self.adata.obs['n_genes'] <= max_genes) &
            (self.adata.obs['total_counts'] <= max_counts) &
            (self.adata.obs['pct_counts_mt'] <= max_mt)
        )
        
        # CRITICAL CHECK before applying filter
        cells_passing = cell_filter.sum()
        st.write(f"- Cells passing filter: {cells_passing:,}/{self.adata.n_obs:,}")
        
        if cells_passing == 0:
            st.error("❌ ΚΡΙΣΙΜΟ: Όλα τα κύτταρα θα φιλτραριστούν!")
            
            # Smart suggestions για διαφορετικούς τύπους datasets
            if self.adata.n_obs < 1000 and self.adata.n_vars < 100:
                st.warning("🔍 Φαίνεται ότι αυτό ΔΕΝ είναι scRNA-seq dataset!")
                st.info("💡 Για μικρά datasets (όπως Iris), χρησιμοποιήστε πολύ χαλαρές παραμέτρους:")
                st.info("   - Min genes: 1-2")
                st.info("   - Max genes: 1000+") 
                st.info("   - Max counts: 100,000+")
                st.info("   - Max MT%: 100%")
            else:
                st.error("Χαλαρώστε τις παραμέτρους φιλτραρίσματος!")
            
            # Εμφάνιση current ranges για guidance
            st.write("**Τρέχουσες τιμές στα δεδομένα:**")
            st.write(f"- n_genes range: {self.adata.obs['n_genes'].min()} - {self.adata.obs['n_genes'].max()}")
            st.write(f"- total_counts range: {self.adata.obs['total_counts'].min():.0f} - {self.adata.obs['total_counts'].max():.0f}")
            st.write(f"- pct_counts_mt range: {self.adata.obs['pct_counts_mt'].min():.1f}% - {self.adata.obs['pct_counts_mt'].max():.1f}%")
            return
        
        self.adata = self.adata[cell_filter, :]
        st.write(f"- Κύτταρα μετά το filtering: {self.adata.n_obs} (από {initial_cells})")
        
        # Step 2: Filtering γονιδίων
        st.write("**Βήμα 2: Filtering γονιδίων**")
        
        # Ξαναϋπολογισμός n_cells μετά το cell filtering
        try:
            n_cells_per_gene = (self.adata.X > 0).sum(axis=0)
            
            # Μετατροπή σε flat array αν χρειάζεται
            if hasattr(n_cells_per_gene, 'A1'):  # sparse matrix result
                n_cells_per_gene = n_cells_per_gene.A1
            elif hasattr(n_cells_per_gene, 'flatten'):  # numpy array
                n_cells_per_gene = n_cells_per_gene.flatten()
            
            # Έλεγχος μεγέθους πριν την εκχώρηση
            if len(n_cells_per_gene) == self.adata.n_vars:
                self.adata.var['n_cells'] = n_cells_per_gene
            else:
                st.warning(f"Mismatch: {len(n_cells_per_gene)} vs {self.adata.n_vars} γονίδια")
                # Fallback: υπολογισμός χωρίς matrix operations
                self.adata.var['n_cells'] = 1  # Default safe value
        except Exception as e:
            st.warning(f"Σφάλμα στον υπολογισμό n_cells: {str(e)}")
            self.adata.var['n_cells'] = 1  # Default safe value
            
        gene_filter = self.adata.var['n_cells'] >= min_cells
        self.adata = self.adata[:, gene_filter]
        st.write(f"- Γονίδια μετά το filtering: {self.adata.n_vars} (από {initial_genes})")
        
        # Step 3: Normalization
        st.write("**Βήμα 3: Normalization**")
        
        # Αποθήκευση raw counts
        self.adata.raw = self.adata
        
        # Normalization σε 10,000 reads per cell
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        
        # Log transformation
        sc.pp.log1p(self.adata)
        
        st.write("- Normalization και log transformation ολοκληρώθηκαν")
        
        # Step 4: Highly variable genes
        st.write("**Βήμα 4: Εντοπισμός highly variable genes**")
        
        try:
            # Έλεγχος αν το dataset είναι αρκετά μεγάλο για HVG analysis
            if self.adata.n_vars >= 10 and self.adata.n_obs >= 10:
                sc.pp.highly_variable_genes(self.adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
                n_hvg = self.adata.var['highly_variable'].sum()
                st.write(f"- Βρέθηκαν {n_hvg} highly variable genes")
            else:
                # Για μικρά datasets, θεωρούμε όλα τα γονίδια ως highly variable
                self.adata.var['highly_variable'] = True
                self.adata.var['dispersions_norm'] = 1.0
                self.adata.var['means'] = self.adata.X.mean(axis=0)
                if hasattr(self.adata.var['means'], 'A1'):
                    self.adata.var['means'] = self.adata.var['means'].A1
                st.write(f"- Dataset μικρό: όλα τα {self.adata.n_vars} γονίδια θεωρούνται highly variable")
        except Exception as e:
            st.warning(f"Πρόβλημα με HVG analysis: {str(e)}")
            # Fallback: όλα τα γονίδια ως highly variable
            self.adata.var['highly_variable'] = True
            self.adata.var['dispersions_norm'] = 1.0
            st.write(f"- Fallback: όλα τα {self.adata.n_vars} γονίδια θεωρούνται highly variable")
        
        st.success("✅ Preprocessing ολοκληρώθηκε επιτυχώς!")
        
        # Αποθήκευση στο session state για χρήση από άλλες σελίδες
        st.session_state['preprocessed_adata'] = self.adata.copy()
        st.session_state['preprocessing_completed'] = True
    
    def _render_results_section(self):
        """Τμήμα αποτελεσμάτων"""
        
        st.markdown("#### 📊 Αποτελέσματα Preprocessing")
        
        if hasattr(st.session_state, 'preprocessed_adata'):
            processed_adata = st.session_state['preprocessed_adata']
            
            # Σύγκριση πριν και μετά
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("##### Πριν το Preprocessing")
                if self.original_adata is not None:
                    st.metric("Κύτταρα", f"{self.original_adata.n_obs:,}")
                    st.metric("Γονίδια", f"{self.original_adata.n_vars:,}")
                else:
                    st.info("💾 Original data δεν διατηρήθηκε λόγω μεγέθους")
                    st.metric("Κύτταρα", "N/A")
                    st.metric("Γονίδια", "N/A")
                
            with col2:
                st.markdown("##### Μετά το Preprocessing")
                st.metric("Κύτταρα", f"{processed_adata.n_obs:,}")
                st.metric("Γονίδια", f"{processed_adata.n_vars:,}")
                if 'highly_variable' in processed_adata.var.columns:
                    hvg_count = processed_adata.var['highly_variable'].sum()
                    st.metric("Highly Variable Genes", f"{hvg_count:,}")
            
            # Δυνατότητα download
            if st.button("💾 Αποθήκευση Προεπεξεργασμένων Δεδομένων"):
                self._save_processed_data(processed_adata)
        else:
            st.info("Δεν υπάρχουν προεπεξεργασμένα δεδομένα. Εκτελέστε πρώτα το preprocessing.")
    
    def _save_processed_data(self, adata):
        """Αποθήκευση προεπεξεργασμένων δεδομένων"""
        
        try:
            # Αποθήκευση στο session state για χρήση από άλλες σελίδες
            st.session_state['preprocessed_adata'] = adata.copy()
            st.session_state['preprocessing_completed'] = True
            
            # Δημιουργία προσωρινού αρχείου με ασφαλή χειρισμό
            import time
            timestamp = int(time.time())
            temp_filename = f"preprocessed_data_{timestamp}.h5ad"
            
            # Δημιουργία temp directory αν δεν υπάρχει
            temp_dir = Path("temp_files")
            temp_dir.mkdir(exist_ok=True)
            
            temp_path = temp_dir / temp_filename
            
            # Αποθήκευση με compression για μικρότερο μέγεθος
            adata.write_h5ad(temp_path, compression='gzip')
            
            # Διάβασμα αρχείου για download
            with open(temp_path, 'rb') as f:
                data = f.read()
            
            st.download_button(
                label="📥 Κατέβασμα H5AD αρχείου",
                data=data,
                file_name="preprocessed_data.h5ad",
                mime="application/octet-stream",
                help="Κατεβάστε το προεπεξεργασμένο dataset"
            )
            
            # Ασφαλής καθαρισμός αρχείου
            try:
                if temp_path.exists():
                    temp_path.unlink()
            except:
                pass  # Αν δεν μπορεί να διαγραφεί, δεν πειράζει
                
            st.success("✅ Δεδομένα αποθηκεύτηκαν επιτυχώς!")
            st.info("💡 Τα προεπεξεργασμένα δεδομένα είναι διαθέσιμα για χρήση στις επόμενες σελίδες")
            
        except Exception as e:
            st.error(f"Σφάλμα κατά την αποθήκευση: {str(e)}")
            # Τουλάχιστον αποθήκευση στο session state
            try:
                st.session_state['preprocessed_adata'] = adata.copy()
                st.session_state['preprocessing_completed'] = True
                st.info("💾 Δεδομένα αποθηκεύτηκαν στη μνήμη της συνεδρίας")
            except Exception:
                st.error("❌ Αποτυχία αποθήκευσης στη μνήμη συνεδρίας")
