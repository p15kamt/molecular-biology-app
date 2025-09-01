"""
Utility functions για Memory Management

Αυτό το module περιέχει βοηθητικές συναρτήσεις για τη διαχείριση μνήμης
και την optimization της επίδοσης για μεγάλα datasets.


Ημερομηνία: 2025
"""

import psutil
import gc
import streamlit as st
import numpy as np
import pandas as pd
from functools import wraps

def get_memory_usage():
    """Επιστρέφει την τρέχουσα χρήση μνήμης σε MB"""
    process = psutil.Process()
    memory_info = process.memory_info()
    return memory_info.rss / 1024 / 1024  # Convert to MB

def get_available_memory():
    """Επιστρέφει τη διαθέσιμη μνήμη συστήματος σε MB"""
    memory = psutil.virtual_memory()
    return memory.available / 1024 / 1024  # Convert to MB

def memory_monitor(func):
    """Decorator για monitoring της χρήσης μνήμης"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        # Μνήμη πριν την εκτέλεση
        memory_before = get_memory_usage()
        
        try:
            result = func(*args, **kwargs)
            
            # Μνήμη μετά την εκτέλεση
            memory_after = get_memory_usage()
            memory_diff = memory_after - memory_before
            
            if memory_diff > 100:  # Αν η αύξηση > 100MB
                st.info(f"💾 Memory usage: +{memory_diff:.1f} MB")
            
            return result
            
        except MemoryError:
            st.error("❌ Ανεπαρκής μνήμη! Δοκιμάστε με μικρότερο dataset ή subsample.")
            gc.collect()  # Force garbage collection
            raise
            
        finally:
            # Καθαρισμός μνήμης αν χρειάζεται
            if get_memory_usage() > 2000:  # Αν > 2GB
                gc.collect()
    
    return wrapper

@st.cache_data(show_spinner=False)
def cached_data_load(file_path, file_type):
    """Cache-enabled data loading για αποφυγή επαναφόρτωσης"""
    import scanpy as sc
    
    if file_type == 'h5ad':
        return sc.read_h5ad(file_path)
    elif file_type in ['csv', 'tsv']:
        separator = ',' if file_type == 'csv' else '\t'
        return pd.read_csv(file_path, sep=separator, index_col=0)
    elif file_type in ['xlsx', 'xls']:
        return pd.read_excel(file_path, index_col=0)

def optimize_adata_memory(adata):
    """Βελτιστοποίηση μνήμης για AnnData object"""
    
    # Μετατροπή σε sparse format αν δεν είναι ήδη
    if not hasattr(adata.X, 'toarray'):
        from scipy.sparse import csr_matrix
        adata.X = csr_matrix(adata.X)
    
    # Καθαρισμός άχρηστων metadata
    if 'raw' in adata.__dict__ and adata.raw is not None:
        # Διατήρηση μόνο αν χρειάζεται
        if len(adata.raw.var_names) == len(adata.var_names):
            adata.raw = None
    
    # Memory optimization για obs/var DataFrames
    for df_name in ['obs', 'var']:
        df = getattr(adata, df_name)
        
        # Μετατροπή object columns σε category όπου είναι δυνατό
        for col in df.columns:
            if df[col].dtype == 'object':
                unique_ratio = df[col].nunique() / len(df[col])
                if unique_ratio < 0.5:  # Αν <50% unique values
                    df[col] = df[col].astype('category')
    
    return adata

def estimate_memory_requirements(n_obs, n_vars, data_type='float32'):
    """Εκτίμηση απαιτήσεων μνήμης για dataset"""
    
    if data_type == 'float32':
        bytes_per_element = 4
    elif data_type == 'float64':
        bytes_per_element = 8
    else:
        bytes_per_element = 4  # Default
    
    # Βασική εκτίμηση για dense matrix
    dense_memory_mb = (n_obs * n_vars * bytes_per_element) / (1024 * 1024)
    
    # Εκτίμηση για sparse matrix (υποθέτοντας 10% sparsity)
    sparse_memory_mb = dense_memory_mb * 0.1
    
    return {
        'dense_mb': dense_memory_mb,
        'sparse_mb': sparse_memory_mb,
        'recommended_format': 'sparse' if dense_memory_mb > 1000 else 'dense'
    }

def display_memory_info():
    """Εμφάνιση πληροφοριών μνήμης στο Streamlit"""
    
    current_usage = get_memory_usage()
    available_memory = get_available_memory()
    total_memory = psutil.virtual_memory().total / 1024 / 1024
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric(
            "Χρήση Μνήμης", 
            f"{current_usage:.0f} MB",
            help="Τρέχουσα χρήση μνήμης από την εφαρμογή"
        )
    
    with col2:
        st.metric(
            "Διαθέσιμη Μνήμη", 
            f"{available_memory:.0f} MB",
            help="Διαθέσιμη μνήμη συστήματος"
        )
    
    with col3:
        usage_percent = (current_usage / total_memory) * 100
        st.metric(
            "% Χρήση", 
            f"{usage_percent:.1f}%",
            help="Ποσοστό χρήσης της συνολικής μνήμης"
        )
    
    # Warning αν η μνήμη είναι χαμηλή
    if available_memory < 500:  # <500MB available
        st.warning("⚠️ Χαμηλή διαθέσιμη μνήμη! Συνιστάται subsample ή restart της εφαρμογής.")
    
    return {
        'current_mb': current_usage,
        'available_mb': available_memory,
        'usage_percent': usage_percent
    }

def suggest_optimization(adata):
    """Προτάσεις βελτίωσης για το dataset"""
    
    suggestions = []
    
    # Έλεγχος μεγέθους
    if adata.n_obs > 50000:
        suggestions.append("🎯 Χρησιμοποιήστε subsample για γρηγορότερη ανάλυση")
    
    if adata.n_vars > 20000:
        suggestions.append("🧬 Εφαρμόστε feature selection για μείωση γονιδίων")
    
    # Έλεγχος sparsity
    if hasattr(adata.X, 'toarray'):
        sparsity = 1 - (adata.X.nnz / (adata.n_obs * adata.n_vars))
        if sparsity < 0.5:
            suggestions.append("💾 Dataset δεν είναι αρκετά sparse - ελέγξτε τα δεδομένα")
    
    # Έλεγχος μνήμης
    memory_est = estimate_memory_requirements(adata.n_obs, adata.n_vars)
    if memory_est['dense_mb'] > 1000:
        suggestions.append("🗜️ Χρησιμοποιήστε sparse format για εξοικονόμηση μνήμης")
    
    return suggestions

def cleanup_memory():
    """Καθαρισμός μνήμης και garbage collection"""
    gc.collect()
    
    # Καθαρισμός Streamlit cache αν χρειάζεται
    if get_memory_usage() > 1500:  # Αν > 1.5GB
        st.cache_data.clear()
        st.info("🧹 Εκκαθάριση cache για εξοικονόμηση μνήμης")
