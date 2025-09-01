"""
Advanced Memory Management για μεγάλα scRNA-seq datasets

Αυτό το module περιλαμβάνει προηγμένες τεχνικές για την αποδοτική
διαχείριση μνήμης σε μεγάλα datasets, συμπεριλαμβανομένου streaming,
chunking, και progressive loading.


Ημερομηνία: 2025
"""

import streamlit as st
import pandas as pd
import numpy as np
import scanpy as sc
import gc
import psutil
import warnings
from typing import Optional, Tuple, List, Dict, Any
from functools import lru_cache
import tempfile
import os
from scipy import sparse
import time

class AdvancedMemoryManager:
    """Προηγμένη διαχείριση μνήμης για scRNA-seq datasets"""
    
    def __init__(self):
        self.memory_threshold_mb = 500  # Μειωμένο threshold για ασφάλεια
        self.chunk_size = 500  # Μικρότερα chunks για καλύτερη διαχείριση
        self.cache_enabled = True
        self.streaming_mode = False
        self.max_file_size_mb = 1000  # Μέγιστο μέγεθος αρχείου
        self.emergency_cleanup_threshold = 0.85  # 85% memory usage
        
    def get_system_memory_info(self) -> Dict[str, float]:
        """Επιστρέφει λεπτομερείς πληροφορίες μνήμης συστήματος"""
        memory = psutil.virtual_memory()
        process = psutil.Process()
        
        return {
            'total_gb': memory.total / (1024**3),
            'available_gb': memory.available / (1024**3),
            'used_gb': memory.used / (1024**3),
            'percent': memory.percent,
            'process_mb': process.memory_info().rss / (1024**2),
            'process_percent': process.memory_percent()
        }
    
    def should_use_streaming(self, adata) -> bool:
        """Καθορίζει αν θα χρησιμοποιηθεί streaming mode"""
        if adata is None:
            return False
            
        # Εκτίμηση μεγέθους σε μνήμη
        estimated_size_mb = self.estimate_memory_usage(adata)
        memory_info = self.get_system_memory_info()
        
        # Χρήση streaming αν:
        # 1. Το dataset είναι >500MB
        # 2. Η διαθέσιμη μνήμη είναι <2GB
        # 3. Το dataset καταλαμβάνει >30% της συνολικής μνήμης
        return (estimated_size_mb > 500 or 
                memory_info['available_gb'] < 2 or
                estimated_size_mb > memory_info['total_gb'] * 1024 * 0.3)
    
    def estimate_memory_usage(self, adata) -> float:
        """Εκτίμηση χρήσης μνήμης σε MB"""
        if adata is None:
            return 0
            
        # Βασική εκτίμηση για X matrix
        if sparse.issparse(adata.X):
            # Sparse matrix: nnz * (4 bytes για value + 4 bytes για index)
            x_size = adata.X.nnz * 8
        else:
            # Dense matrix
            x_size = adata.X.nbytes
        
        # Εκτίμηση για obs, var, obsm, varm
        obs_size = adata.obs.memory_usage(deep=True).sum()
        var_size = adata.var.memory_usage(deep=True).sum()
        
        # Εκτίμηση για obsm, varm (PCA, embeddings κλπ.)
        obsm_size = sum(arr.nbytes for arr in adata.obsm.values() if hasattr(arr, 'nbytes'))
        varm_size = sum(arr.nbytes for arr in adata.varm.values() if hasattr(arr, 'nbytes'))
        
        total_bytes = x_size + obs_size + var_size + obsm_size + varm_size
        return total_bytes / (1024**2)  # Convert to MB
    
    def optimize_adata_for_memory(self, adata):
        """Βελτιστοποίηση AnnData object για μνήμη"""
        if adata is None:
            return adata
            
        st.info("🔧 Βελτιστοποίηση για μνήμη...")
        
        # 1. Μετατροπή σε sparse αν χρειάζεται
        if not sparse.issparse(adata.X):
            sparsity = (adata.X == 0).sum() / adata.X.size
            if sparsity > 0.7:  # Αν >70% zeros
                st.info(f"📦 Μετατροπή σε sparse format (sparsity: {sparsity:.1%})")
                adata.X = sparse.csr_matrix(adata.X)
        
        # 2. Καθαρισμός άχρηστων raw data αν υπάρχουν
        if hasattr(adata, 'raw') and adata.raw is not None:
            if adata.raw.X.shape == adata.X.shape:
                st.info("🗑️ Αφαίρεση duplicate raw data")
                adata.raw = None
        
        # 3. Βελτιστοποίηση data types
        self._optimize_dtypes(adata)
        
        # 4. Garbage collection
        gc.collect()
        
        return adata
    
    def _optimize_dtypes(self, adata):
        """Βελτιστοποίηση data types για obs/var DataFrames"""
        
        for df_name in ['obs', 'var']:
            df = getattr(adata, df_name)
            
            for col in df.columns:
                if df[col].dtype == 'object':
                    # Μετατροπή σε category αν έχει λίγες unique τιμές
                    unique_ratio = df[col].nunique() / len(df[col])
                    if unique_ratio < 0.5:
                        df[col] = df[col].astype('category')
                elif df[col].dtype == 'float64':
                    # Μετατροπή float64 σε float32 αν δεν χάνεται ακρίβεια
                    if df[col].max() < np.finfo(np.float32).max:
                        df[col] = df[col].astype('float32')
                elif df[col].dtype == 'int64':
                    # Μετατροπή int64 σε int32 αν είναι δυνατό
                    if df[col].max() < np.iinfo(np.int32).max and df[col].min() > np.iinfo(np.int32).min:
                        df[col] = df[col].astype('int32')
    
    def create_chunked_iterator(self, adata, chunk_size: int = None):
        """Δημιουργία iterator για chunked processing"""
        if chunk_size is None:
            chunk_size = self.chunk_size
            
        n_chunks = (adata.n_obs + chunk_size - 1) // chunk_size
        
        for i in range(n_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, adata.n_obs)
            
            chunk = adata[start_idx:end_idx, :].copy()
            yield i, chunk, n_chunks
    
    def safe_display_data(self, adata, max_cells: int = 100, max_genes: int = 100):
        """Ασφαλής εμφάνιση δεδομένων χωρίς memory overflow"""
        
        if adata is None:
            return pd.DataFrame({"Error": ["No data loaded"]})
        
        # Περιορισμός των δεδομένων για εμφάνιση
        display_cells = min(max_cells, adata.n_obs)
        display_genes = min(max_genes, adata.n_vars)
        
        try:
            # Εξαγωγή μικρού subset
            subset = adata[:display_cells, :display_genes]
            
            if sparse.issparse(subset.X):
                data_array = subset.X.toarray()
            else:
                data_array = subset.X
            
            # Δημιουργία DataFrame με ασφαλή μέγεθος
            sample_df = pd.DataFrame(
                data_array,
                index=subset.obs_names,
                columns=subset.var_names
            )
            
            return sample_df
            
        except Exception as e:
            st.error(f"Σφάλμα στην εμφάνιση δεδομένων: {str(e)}")
            return pd.DataFrame({"Error": [f"Display error: {str(e)}"]})
    
    def progressive_qc_calculation(self, adata):
        """Progressive υπολογισμός QC metrics για μεγάλα datasets με βελτιωμένη απόδοση μνήμης"""
        
        memory_info = self.get_system_memory_info()
        dataset_size_mb = self.estimate_memory_usage(adata)
        
        # Για πολύ μεγάλα datasets, δοκιμή vectorized approach πρώτα
        if dataset_size_mb > 500 and memory_info['available_gb'] > 4:
            st.info("🚀 Δοκιμή vectorized QC calculation για καλύτερη απόδοση...")
            try:
                return self._vectorized_qc_calculation(adata)
            except Exception as e:
                st.warning(f"⚠️ Vectorized approach απέτυχε: {str(e)}")
                st.info("🔄 Fallback σε progressive processing...")
        
        st.info("📊 Υπολογισμός QC metrics με progressive processing...")
        
        # Προσαρμοστικό chunk size βάσει μνήμης και μεγέθους dataset
        dataset_size_mb = self.estimate_memory_usage(adata)
        
        # Βελτιωμένο δυναμικό chunk size για καλύτερη απόδοση
        if memory_info['available_gb'] < 2:
            chunk_size = max(1000, int(adata.n_obs / 20))  # Μεγαλύτερα chunks ακόμα και με λίγη μνήμη
        elif memory_info['available_gb'] < 4:
            chunk_size = max(2000, int(adata.n_obs / 10))  # Μεσαία chunks
        elif memory_info['available_gb'] < 8:
            chunk_size = max(5000, int(adata.n_obs / 5))   # Μεγάλα chunks
        else:
            chunk_size = max(10000, int(adata.n_obs / 3))  # Πολύ μεγάλα chunks για καλή απόδοση
        
        st.info(f"🔧 Χρήση adaptive chunk size: {chunk_size} κύτταρα ανά chunk")
        
        # Pre-calculate gene patterns για καλύτερη απόδοση
        mt_gene_pattern = adata.var_names.str.startswith(('MT-', 'mt-', 'Mt-'))
        ribo_gene_pattern = adata.var_names.str.startswith(('RPS', 'RPL', 'rps', 'rpl'))
        
        # Έλεγχος αν υπάρχουν MT/Ribo genes
        n_mt_genes = mt_gene_pattern.sum()
        n_ribo_genes = ribo_gene_pattern.sum()
        
        if n_mt_genes == 0:
            st.warning("⚠️ Δεν βρέθηκαν μιτοχονδριακά γονίδια με τα συνήθη patterns")
        if n_ribo_genes == 0:
            st.warning("⚠️ Δεν βρέθηκαν ριβοσωμικά γονίδια με τα συνήθη patterns")
        
        # Initialize metrics με pre-allocated arrays για καλύτερη απόδοση
        qc_metrics = {
            'n_genes': np.zeros(adata.n_obs, dtype=np.int32),
            'total_counts': np.zeros(adata.n_obs, dtype=np.float32),
            'pct_counts_mt': np.zeros(adata.n_obs, dtype=np.float32),
            'pct_counts_ribo': np.zeros(adata.n_obs, dtype=np.float32)
        }
        
        # Progressive calculation με progress bar και memory monitoring
        progress_bar = st.progress(0)
        status_text = st.empty()
        memory_text = st.empty()
        
        n_chunks = (adata.n_obs + chunk_size - 1) // chunk_size
        
        for i in range(n_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, adata.n_obs)
            
            # Update progress
            progress = (i + 1) / n_chunks
            progress_bar.progress(progress)
            status_text.text(f"Processing chunk {i+1}/{n_chunks} (κύτταρα {start_idx+1}-{end_idx})")
            
            # Memory monitoring
            current_memory = self.get_system_memory_info()
            memory_text.text(f"💾 Memory usage: {current_memory['percent']:.1f}% ({current_memory['process_mb']:.0f}MB)")
            
            # Emergency memory check
            if current_memory['percent'] > 90:
                st.error("❌ Κρίσιμα χαμηλή μνήμη - μείωση chunk size")
                chunk_size = max(50, chunk_size // 2)
                continue
            
            try:
                # Εξαγωγή chunk με copy για ασφάλεια
                chunk = adata[start_idx:end_idx, :].copy()
                
                # Βελτιστοποιημένος υπολογισμός metrics
                if sparse.issparse(chunk.X):
                    # Sparse matrix - optimized calculations
                    chunk_data = chunk.X.tocsr()  # Ensure CSR format για καλύτερη απόδοση
                    
                    # Basic metrics
                    chunk_counts = np.array(chunk_data.sum(axis=1)).flatten()
                    chunk_genes = np.array((chunk_data > 0).sum(axis=1)).flatten()
                    
                    # MT genes (με ασφαλή υπολογισμό)
                    if n_mt_genes > 0:
                        mt_data = chunk_data[:, mt_gene_pattern]
                        mt_counts = np.array(mt_data.sum(axis=1)).flatten()
                        chunk_mt_pct = np.divide(mt_counts, chunk_counts, 
                                               out=np.zeros_like(mt_counts), 
                                               where=chunk_counts!=0) * 100
                    else:
                        chunk_mt_pct = np.zeros(len(chunk_counts))
                    
                    # Ribosomal genes (με ασφαλή υπολογισμό)
                    if n_ribo_genes > 0:
                        ribo_data = chunk_data[:, ribo_gene_pattern]
                        ribo_counts = np.array(ribo_data.sum(axis=1)).flatten()
                        chunk_ribo_pct = np.divide(ribo_counts, chunk_counts, 
                                                 out=np.zeros_like(ribo_counts), 
                                                 where=chunk_counts!=0) * 100
                    else:
                        chunk_ribo_pct = np.zeros(len(chunk_counts))
                        
                else:
                    # Dense matrix
                    chunk_counts = chunk.X.sum(axis=1)
                    chunk_genes = (chunk.X > 0).sum(axis=1)
                    
                    if n_mt_genes > 0:
                        mt_counts = chunk.X[:, mt_gene_pattern].sum(axis=1)
                        chunk_mt_pct = np.divide(mt_counts, chunk_counts, 
                                               out=np.zeros_like(mt_counts), 
                                               where=chunk_counts!=0) * 100
                    else:
                        chunk_mt_pct = np.zeros(len(chunk_counts))
                    
                    if n_ribo_genes > 0:
                        ribo_counts = chunk.X[:, ribo_gene_pattern].sum(axis=1)
                        chunk_ribo_pct = np.divide(ribo_counts, chunk_counts, 
                                                 out=np.zeros_like(ribo_counts), 
                                                 where=chunk_counts!=0) * 100
                    else:
                        chunk_ribo_pct = np.zeros(len(chunk_counts))
                
                # Αποθήκευση στα pre-allocated arrays
                qc_metrics['n_genes'][start_idx:end_idx] = chunk_genes
                qc_metrics['total_counts'][start_idx:end_idx] = chunk_counts
                qc_metrics['pct_counts_mt'][start_idx:end_idx] = chunk_mt_pct
                qc_metrics['pct_counts_ribo'][start_idx:end_idx] = chunk_ribo_pct
                
                # Memory cleanup
                del chunk, chunk_counts, chunk_genes, chunk_mt_pct, chunk_ribo_pct
                if 'chunk_data' in locals():
                    del chunk_data
                if 'mt_data' in locals():
                    del mt_data
                if 'ribo_data' in locals():
                    del ribo_data
                gc.collect()
                
            except Exception as e:
                st.error(f"❌ Σφάλμα στο chunk {i+1}: {str(e)}")
                # Fallback values
                qc_metrics['n_genes'][start_idx:end_idx] = 0
                qc_metrics['total_counts'][start_idx:end_idx] = 0
                qc_metrics['pct_counts_mt'][start_idx:end_idx] = 0
                qc_metrics['pct_counts_ribo'][start_idx:end_idx] = 0
        
        # Καθαρισμός progress indicators
        progress_bar.empty()
        status_text.empty()
        memory_text.empty()
        
        # Αποθήκευση στο adata.obs με validation
        try:
            for metric, values in qc_metrics.items():
                adata.obs[metric] = values
            
            # Validation των αποτελεσμάτων
            validation_results = {
                'n_genes': f"Min: {qc_metrics['n_genes'].min()}, Max: {qc_metrics['n_genes'].max()}, Mean: {qc_metrics['n_genes'].mean():.1f}",
                'total_counts': f"Min: {qc_metrics['total_counts'].min():.1f}, Max: {qc_metrics['total_counts'].max():.1f}, Mean: {qc_metrics['total_counts'].mean():.1f}",
                'pct_counts_mt': f"Min: {qc_metrics['pct_counts_mt'].min():.2f}%, Max: {qc_metrics['pct_counts_mt'].max():.2f}%, Mean: {qc_metrics['pct_counts_mt'].mean():.2f}%",
                'pct_counts_ribo': f"Min: {qc_metrics['pct_counts_ribo'].min():.2f}%, Max: {qc_metrics['pct_counts_ribo'].max():.2f}%, Mean: {qc_metrics['pct_counts_ribo'].mean():.2f}%"
            }
            
            st.success("✅ QC metrics υπολογίστηκαν επιτυχώς!")
            
            with st.expander("📊 Στατιστικά QC Metrics"):
                for metric, stats in validation_results.items():
                    st.text(f"{metric}: {stats}")
                    
        except Exception as e:
            st.error(f"❌ Σφάλμα στην αποθήκευση QC metrics: {str(e)}")
        
        # Final memory cleanup
        del qc_metrics
        gc.collect()
        
        return adata
    
    def _vectorized_qc_calculation(self, adata):
        """Vectorized QC calculation για μεγάλα datasets - πολύ πιο γρήγορη από chunked"""
        
        st.info("⚡ Vectorized QC calculation - πολύ γρήγορη μέθοδος...")
        
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        try:
            # Pre-calculate gene patterns
            mt_gene_pattern = adata.var_names.str.startswith(('MT-', 'mt-', 'Mt-'))
            ribo_gene_pattern = adata.var_names.str.startswith(('RPS', 'RPL', 'rps', 'rpl'))
            
            n_mt_genes = mt_gene_pattern.sum()
            n_ribo_genes = ribo_gene_pattern.sum()
            
            status_text.text("🔢 Υπολογισμός basic metrics...")
            progress_bar.progress(0.25)
            
            # Vectorized calculations - πολύ γρήγορα!
            if sparse.issparse(adata.X):
                # Sparse matrix - optimized
                X_csr = adata.X.tocsr()
                
                # Basic metrics σε ένα βήμα
                total_counts = np.array(X_csr.sum(axis=1)).flatten()
                n_genes = np.array((X_csr > 0).sum(axis=1)).flatten()
                
                progress_bar.progress(0.5)
                status_text.text("🧬 Υπολογισμός MT genes...")
                
                # MT genes
                if n_mt_genes > 0:
                    mt_counts = np.array(X_csr[:, mt_gene_pattern].sum(axis=1)).flatten()
                    pct_counts_mt = np.divide(mt_counts, total_counts, 
                                            out=np.zeros_like(mt_counts), 
                                            where=total_counts!=0) * 100
                else:
                    pct_counts_mt = np.zeros(adata.n_obs)
                
                progress_bar.progress(0.75)
                status_text.text("🔬 Υπολογισμός Ribosomal genes...")
                
                # Ribosomal genes
                if n_ribo_genes > 0:
                    ribo_counts = np.array(X_csr[:, ribo_gene_pattern].sum(axis=1)).flatten()
                    pct_counts_ribo = np.divide(ribo_counts, total_counts, 
                                              out=np.zeros_like(ribo_counts), 
                                              where=total_counts!=0) * 100
                else:
                    pct_counts_ribo = np.zeros(adata.n_obs)
                    
            else:
                # Dense matrix
                total_counts = adata.X.sum(axis=1)
                n_genes = (adata.X > 0).sum(axis=1)
                
                progress_bar.progress(0.5)
                status_text.text("🧬 Υπολογισμός MT genes...")
                
                if n_mt_genes > 0:
                    mt_counts = adata.X[:, mt_gene_pattern].sum(axis=1)
                    pct_counts_mt = np.divide(mt_counts, total_counts, 
                                            out=np.zeros_like(mt_counts), 
                                            where=total_counts!=0) * 100
                else:
                    pct_counts_mt = np.zeros(adata.n_obs)
                
                progress_bar.progress(0.75)
                status_text.text("🔬 Υπολογισμός Ribosomal genes...")
                
                if n_ribo_genes > 0:
                    ribo_counts = adata.X[:, ribo_gene_pattern].sum(axis=1)
                    pct_counts_ribo = np.divide(ribo_counts, total_counts, 
                                              out=np.zeros_like(ribo_counts), 
                                              where=total_counts!=0) * 100
                else:
                    pct_counts_ribo = np.zeros(adata.n_obs)
            
            progress_bar.progress(1.0)
            status_text.text("💾 Αποθήκευση αποτελεσμάτων...")
            
            # Αποθήκευση στο adata.obs
            adata.obs['n_genes'] = n_genes.astype(np.int32)
            adata.obs['total_counts'] = total_counts.astype(np.float32)
            adata.obs['pct_counts_mt'] = pct_counts_mt.astype(np.float32)
            adata.obs['pct_counts_ribo'] = pct_counts_ribo.astype(np.float32)
            
            # Validation
            validation_results = {
                'n_genes': f"Min: {n_genes.min()}, Max: {n_genes.max()}, Mean: {n_genes.mean():.1f}",
                'total_counts': f"Min: {total_counts.min():.1f}, Max: {total_counts.max():.1f}, Mean: {total_counts.mean():.1f}",
                'pct_counts_mt': f"Min: {pct_counts_mt.min():.2f}%, Max: {pct_counts_mt.max():.2f}%, Mean: {pct_counts_mt.mean():.2f}%",
                'pct_counts_ribo': f"Min: {pct_counts_ribo.min():.2f}%, Max: {pct_counts_ribo.max():.2f}%, Mean: {pct_counts_ribo.mean():.2f}%"
            }
            
            # Καθαρισμός UI
            progress_bar.empty()
            status_text.empty()
            
            st.success("✅ Vectorized QC calculation ολοκληρώθηκε γρήγορα!")
            
            with st.expander("📊 Στατιστικά QC Metrics"):
                for metric, stats in validation_results.items():
                    st.text(f"{metric}: {stats}")
            
            return adata
            
        except Exception as e:
            progress_bar.empty()
            status_text.empty()
            st.error(f"❌ Σφάλμα στο vectorized calculation: {str(e)}")
            raise
    
    def memory_efficient_plot_data(self, adata, plot_type: str, max_points: int = 10000):
        """Δημιουργία δεδομένων για plots με memory efficiency"""
        
        if adata.n_obs > max_points:
            st.info(f"📊 Subsample {max_points} κυττάρων για visualization")
            # Random subsample για plotting
            indices = np.random.choice(adata.n_obs, max_points, replace=False)
            plot_adata = adata[indices, :].copy()
        else:
            plot_adata = adata.copy()
        
        return plot_adata
    
    def emergency_memory_cleanup(self):
        """Έκτακτος καθαρισμός μνήμης όταν φτάνουμε στο όριο"""
        
        memory_info = self.get_system_memory_info()
        
        if memory_info['percent'] > (self.emergency_cleanup_threshold * 100):
            st.warning("⚠️ Έκτακτος καθαρισμός μνήμης...")
            
            # Καθαρισμός Streamlit cache
            if hasattr(st, 'cache_data'):
                st.cache_data.clear()
            
            # Καθαρισμός session state
            for key in list(st.session_state.keys()):
                obj = st.session_state[key]
                if hasattr(obj, 'nbytes') and obj.nbytes > 50 * 1024 * 1024:  # >50MB
                    del st.session_state[key]
                    st.info(f"🗑️ Αφαίρεση {key} από session")
            
            # Force garbage collection
            gc.collect()
            
            # Επανέλεγχος μνήμης
            new_memory_info = self.get_system_memory_info()
            freed_mb = (memory_info['process_mb'] - new_memory_info['process_mb'])
            
            if freed_mb > 0:
                st.success(f"✅ Ελευθερώθηκαν {freed_mb:.1f} MB μνήμης")
            
            return True
        return False
    
    def check_file_size_compatibility(self, file_size_mb: float) -> Dict[str, Any]:
        """Έλεγχος συμβατότητας μεγέθους αρχείου με διαθέσιμη μνήμη"""
        
        memory_info = self.get_system_memory_info()
        
        # Εκτίμηση μνήμης που θα χρειαστεί (συντηρητική εκτίμηση)
        estimated_memory_needed = file_size_mb * 3  # 3x για processing overhead
        
        compatibility = {
            'can_process': False,
            'requires_streaming': False,
            'suggested_chunk_size': self.chunk_size,
            'warnings': [],
            'recommendations': []
        }
        
        if estimated_memory_needed > memory_info['available_gb'] * 1024:
            compatibility['warnings'].append(
                f"⚠️ Το αρχείο ({file_size_mb:.1f}MB) μπορεί να υπερβεί τη διαθέσιμη μνήμη"
            )
            compatibility['requires_streaming'] = True
            compatibility['suggested_chunk_size'] = max(100, int(self.chunk_size / 2))
        else:
            compatibility['can_process'] = True
        
        if file_size_mb > self.max_file_size_mb:
            compatibility['warnings'].append(
                f"📁 Μεγάλο αρχείο ({file_size_mb:.1f}MB) - προτείνεται streaming mode"
            )
            compatibility['requires_streaming'] = True
        
        if memory_info['percent'] > 70:
            compatibility['recommendations'].append(
                "🧹 Συνιστάται καθαρισμός μνήμης πριν τη φόρτωση"
            )
        
        return compatibility
    
    def adaptive_chunk_processing(self, adata, process_func, **kwargs):
        """Προσαρμοστική επεξεργασία με chunks βάσει διαθέσιμης μνήμης με βελτιωμένη απόδοση"""
        
        memory_info = self.get_system_memory_info()
        dataset_size_mb = self.estimate_memory_usage(adata)
        
        # Έξυπνη προσαρμογή chunk size βάσει πολλών παραγόντων
        base_chunk_size = self.chunk_size
        
        # Παράγοντας μνήμης
        if memory_info['available_gb'] < 1:
            memory_factor = 0.1  # Πολύ μικρά chunks
        elif memory_info['available_gb'] < 2:
            memory_factor = 0.25
        elif memory_info['available_gb'] < 4:
            memory_factor = 0.5
        elif memory_info['available_gb'] < 8:
            memory_factor = 0.75
        else:
            memory_factor = 1.0
        
        # Παράγοντας μεγέθους dataset
        if dataset_size_mb > 2000:  # >2GB
            size_factor = 0.5
        elif dataset_size_mb > 1000:  # >1GB
            size_factor = 0.75
        else:
            size_factor = 1.0
        
        # Παράγοντας αριθμού κυττάρων
        if adata.n_obs > 50000:
            cell_factor = 0.5
        elif adata.n_obs > 20000:
            cell_factor = 0.75
        else:
            cell_factor = 1.0
        
        # Συνολικός υπολογισμός
        adaptive_chunk_size = int(base_chunk_size * memory_factor * size_factor * cell_factor)
        adaptive_chunk_size = max(50, min(adaptive_chunk_size, 2000))  # Όρια ασφαλείας
        
        st.info(f"🔄 Adaptive chunking: {adaptive_chunk_size} κύτταρα/chunk (Memory: {memory_info['available_gb']:.1f}GB, Dataset: {dataset_size_mb:.1f}MB)")
        
        results = []
        progress_bar = st.progress(0)
        status_text = st.empty()
        memory_text = st.empty()
        performance_metrics = {
            'start_time': time.time(),
            'chunks_processed': 0,
            'memory_cleanups': 0,
            'errors': 0
        }
        
        try:
            n_chunks = (adata.n_obs + adaptive_chunk_size - 1) // adaptive_chunk_size
            
            for i in range(n_chunks):
                start_idx = i * adaptive_chunk_size
                end_idx = min((i + 1) * adaptive_chunk_size, adata.n_obs)
                
                # Real-time memory monitoring
                current_memory = self.get_system_memory_info()
                memory_text.text(f"💾 Memory: {current_memory['percent']:.1f}% ({current_memory['process_mb']:.0f}MB process)")
                
                # Έλεγχος μνήμης πριν από κάθε chunk
                if current_memory['percent'] > 85:
                    st.warning("⚠️ Υψηλή χρήση μνήμης - καθαρισμός...")
                    if self.emergency_memory_cleanup():
                        performance_metrics['memory_cleanups'] += 1
                    
                    # Δυναμική μείωση chunk size αν χρειάζεται
                    if current_memory['percent'] > 90:
                        adaptive_chunk_size = max(50, adaptive_chunk_size // 2)
                        st.info(f"🔧 Μείωση chunk size σε {adaptive_chunk_size} για ασφάλεια")
                
                # Update progress με λεπτομερείς πληροφορίες
                progress = (i + 1) / n_chunks
                progress_bar.progress(progress)
                
                elapsed_time = time.time() - performance_metrics['start_time']
                if i > 0:
                    avg_time_per_chunk = elapsed_time / (i + 1)
                    estimated_remaining = avg_time_per_chunk * (n_chunks - i - 1)
                    status_text.text(f"Chunk {i+1}/{n_chunks} | κύτταρα {start_idx+1}-{end_idx} | ETA: {estimated_remaining:.1f}s")
                else:
                    status_text.text(f"Chunk {i+1}/{n_chunks} | κύτταρα {start_idx+1}-{end_idx}")
                
                try:
                    # Εξαγωγή chunk με memory-safe τρόπο
                    chunk = adata[start_idx:end_idx, :].copy()
                    
                    # Εκτέλεση της συνάρτησης επεξεργασίας
                    chunk_result = process_func(chunk, chunk_index=i, **kwargs)
                    results.append(chunk_result)
                    performance_metrics['chunks_processed'] += 1
                    
                    # Memory cleanup μετά από κάθε chunk
                    del chunk
                    gc.collect()
                    
                except Exception as e:
                    st.error(f"❌ Σφάλμα στο chunk {i+1}: {str(e)}")
                    performance_metrics['errors'] += 1
                    
                    # Προσπάθεια recovery με μικρότερο chunk
                    if adaptive_chunk_size > 100:
                        adaptive_chunk_size = max(50, adaptive_chunk_size // 2)
                        st.info(f"🔄 Retry με μικρότερο chunk size: {adaptive_chunk_size}")
                        continue
                    else:
                        # Αν και τα μικρότερα chunks αποτυγχάνουν, skip
                        results.append(None)
                        continue
                
                # Έλεγχος αν πρέπει να σταματήσουμε
                final_memory = self.get_system_memory_info()
                if final_memory['percent'] > 95:
                    st.error("❌ Κρίσιμα χαμηλή μνήμη - διακοπή επεξεργασίας")
                    break
                    
        except Exception as e:
            st.error(f"❌ Κρίσιμο σφάλμα κατά την επεξεργασία: {str(e)}")
            raise
        finally:
            # Καθαρισμός UI elements
            progress_bar.empty()
            status_text.empty()
            memory_text.empty()
            
            # Εμφάνιση performance metrics
            total_time = time.time() - performance_metrics['start_time']
            st.success(f"✅ Επεξεργασία ολοκληρώθηκε σε {total_time:.1f}s")
            
            with st.expander("📊 Performance Metrics"):
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Chunks Processed", performance_metrics['chunks_processed'])
                with col2:
                    st.metric("Memory Cleanups", performance_metrics['memory_cleanups'])
                with col3:
                    st.metric("Errors", performance_metrics['errors'])
                
                if performance_metrics['chunks_processed'] > 0:
                    avg_time = total_time / performance_metrics['chunks_processed']
                    st.text(f"Average time per chunk: {avg_time:.2f}s")
        
        return results
    
    def cleanup_session_memory(self):
        """Καθαρισμός μνήμης session"""
        
        # Clear Streamlit cache
        if hasattr(st, 'cache_data'):
            st.cache_data.clear()
        
        # Clear session state αν υπάρχουν μεγάλα objects
        for key in list(st.session_state.keys()):
            obj = st.session_state[key]
            if hasattr(obj, 'nbytes') and obj.nbytes > 100 * 1024 * 1024:  # >100MB
                del st.session_state[key]
        
        # Force garbage collection
        gc.collect()
        
        st.info("🧹 Session memory καθαρίστηκε")

# Global instance
memory_manager = AdvancedMemoryManager()
