"""
Robust Error Handler για την εφαρμογή Μοριακής Βιολογίας

Αυτό το module περιέχει global error handling και recovery mechanisms
για να εξασφαλίσει ότι η εφαρμογή δεν κρασάρει ποτέ.


Ημερομηνία: 2025
"""

import streamlit as st
import traceback
import gc
import psutil
import numpy as np
from functools import wraps
import sys

class RobustErrorHandler:
    """Κλάση για robust error handling σε όλη την εφαρμογή"""
    
    def __init__(self):
        self.error_count = 0
        self.max_errors_before_reset = 3
        
    def safe_execute(self, func_name="Άγνωστη λειτουργία"):
        """Decorator για ασφαλή εκτέλεση συναρτήσεων"""
        
        def decorator(func):
            @wraps(func)
            def wrapper(*args, **kwargs):
                try:
                    # Έλεγχος μνήμης πριν την εκτέλεση
                    self._check_memory_before_execution()
                    
                    # Εκτέλεση της συνάρτησης
                    result = func(*args, **kwargs)
                    
                    # Reset error count σε επιτυχή εκτέλεση
                    self.error_count = 0
                    
                    return result
                    
                except MemoryError:
                    return self._handle_memory_error(func_name)
                    
                except Exception as e:
                    return self._handle_general_error(e, func_name)
            
            return wrapper
        return decorator
    
    def _check_memory_before_execution(self):
        """Έλεγχος μνήμης πριν την εκτέλεση"""
        
        memory = psutil.virtual_memory()
        
        # Αν η μνήμη είναι <1GB διαθέσιμη, καθαρισμός
        if memory.available < 1024 * 1024 * 1024:  # 1GB
            st.warning("⚠️ Χαμηλή διαθέσιμη μνήμη - εκτέλεση καθαρισμού...")
            self._emergency_cleanup()
    
    def _handle_memory_error(self, func_name):
        """Διαχείριση Memory Error"""
        
        self.error_count += 1
        
        st.error(f"❌ **Memory Error στη λειτουργία: {func_name}**")
        st.error("💾 Ανεπαρκής μνήμη για την εκτέλεση αυτής της λειτουργίας")
        
        # Άμεσος καθαρισμός μνήμης
        self._emergency_cleanup()
        
        # Προτάσεις λύσης
        st.markdown("### 💡 **Προτεινόμενες Λύσεις:**")
        st.markdown("""
        1. **🔄 Κάντε restart της εφαρμογής** (Ctrl+C στο terminal και ξανά `python run_app.py`)
        2. **📉 Χρησιμοποιήστε μικρότερο dataset** ή κάντε subsample
        3. **🧹 Κλείστε άλλες εφαρμογές** για εξοικονόμηση μνήμης
        4. **⚡ Χρησιμοποιήστε το "Progressive Loading"** mode
        """)
        
        # Κουμπί για emergency cleanup
        if st.button("🆘 Emergency Memory Cleanup", type="primary"):
            self._emergency_cleanup()
            st.experimental_rerun()
        
        # Αν πολλά errors, προτείνω restart
        if self.error_count >= self.max_errors_before_reset:
            st.error("🚨 **Πολλαπλά memory errors - συνιστάται restart της εφαρμογής**")
            
        return None
    
    def _handle_general_error(self, error, func_name):
        """Διαχείριση γενικών errors"""
        
        self.error_count += 1
        
        st.error(f"❌ **Σφάλμα στη λειτουργία: {func_name}**")
        st.error(f"**Τύπος σφάλματος:** {type(error).__name__}")
        st.error(f"**Μήνυμα:** {str(error)}")
        
        # Εμφάνιση detailed error σε expander
        with st.expander("🔍 Λεπτομέρειες Σφάλματος (για developers)"):
            st.code(traceback.format_exc())
        
        # Προτάσεις ανάλογα με τον τύπο error
        if "memory" in str(error).lower():
            st.markdown("### 💡 **Πιθανή αιτία:** Πρόβλημα μνήμης")
            st.markdown("**Λύσεις:** Μικρότερο dataset, restart εφαρμογής")
            
        elif "index" in str(error).lower():
            st.markdown("### 💡 **Πιθανή αιτία:** Πρόβλημα με δεδομένα")
            st.markdown("**Λύσεις:** Ελέγξτε τη μορφή των δεδομένων, προεπεξεργασία")
            
        elif "import" in str(error).lower():
            st.markdown("### 💡 **Πιθανή αιτία:** Λείπει dependency")
            st.markdown("**Λύσεις:** Εγκατάσταση dependencies, έλεγχος virtual environment")
            
        else:
            st.markdown("### 💡 **Γενικές Λύσεις:**")
            st.markdown("""
            1. **🔄 Δοκιμάστε ξανά** - μπορεί να είναι προσωρινό πρόβλημα
            2. **🔄 Restart της εφαρμογής**
            3. **📧 Επικοινωνήστε με την ομάδα** αν το πρόβλημα επιμένει
            """)
        
        # Κουμπί retry
        if st.button("🔄 Δοκιμή Ξανά", key=f"retry_{func_name}"):
            st.experimental_rerun()
        
        return None
    
    def _emergency_cleanup(self):
        """Έκτακτος καθαρισμός μνήμης"""
        
        try:
            # Clear Streamlit cache
            if hasattr(st, 'cache_data'):
                st.cache_data.clear()
            
            # Clear session state μεγάλων objects
            for key in list(st.session_state.keys()):
                try:
                    obj = st.session_state[key]
                    if hasattr(obj, 'nbytes') and obj.nbytes > 100 * 1024 * 1024:  # >100MB
                        del st.session_state[key]
                        st.info(f"🗑️ Αφαίρεση {key} από session")
                except:
                    pass
            
            # Force garbage collection
            gc.collect()
            
            # Εμφάνιση νέας κατάστασης μνήμης
            memory = psutil.virtual_memory()
            st.success(f"🧹 Emergency cleanup - Διαθέσιμη μνήμη: {memory.available / (1024**3):.1f}GB")
            
        except Exception as e:
            st.warning(f"⚠️ Πρόβλημα στον καθαρισμό μνήμης: {str(e)}")
    
    def memory_safe_operation(self, operation_name="Λειτουργία"):
        """Context manager για memory-safe operations"""
        
        class MemorySafeContext:
            def __init__(self, handler, op_name):
                self.handler = handler
                self.op_name = op_name
                self.initial_memory = None
            
            def __enter__(self):
                # Καταγραφή αρχικής μνήμης
                self.initial_memory = psutil.virtual_memory().available
                
                # Έλεγχος αν έχουμε αρκετή μνήμη
                if self.initial_memory < 512 * 1024 * 1024:  # <512MB
                    st.warning(f"⚠️ Χαμηλή μνήμη για {self.op_name} - εκτέλεση cleanup...")
                    self.handler._emergency_cleanup()
                
                return self
            
            def __exit__(self, exc_type, exc_val, exc_tb):
                if exc_type is not None:
                    if exc_type == MemoryError:
                        self.handler._handle_memory_error(self.op_name)
                        return True  # Suppress the exception
                    else:
                        self.handler._handle_general_error(exc_val, self.op_name)
                        return True  # Suppress the exception
                
                # Καθαρισμός μετά την εκτέλεση
                gc.collect()
                
                return False
        
        return MemorySafeContext(self, operation_name)
    
    def create_fallback_data(self, original_data, operation="visualization"):
        """Δημιουργία fallback data όταν αποτυγχάνει η κανονική επεξεργασία"""
        
        try:
            if hasattr(original_data, 'n_obs'):
                # AnnData object
                n_cells = min(1000, original_data.n_obs)
                n_genes = min(2000, original_data.n_vars)
                
                # Random subsample
                cell_indices = np.random.choice(original_data.n_obs, n_cells, replace=False)
                gene_indices = np.random.choice(original_data.n_vars, n_genes, replace=False)
                
                fallback_adata = original_data[cell_indices, :][:, gene_indices].copy()
                
                st.info(f"🎯 Fallback mode: {n_cells} κύτταρα × {n_genes} γονίδια")
                
                return fallback_adata
                
            elif isinstance(original_data, np.ndarray):
                # NumPy array
                if original_data.size > 10000:
                    # Subsample
                    flat_data = original_data.flatten()
                    sample_size = min(10000, len(flat_data))
                    sample_indices = np.random.choice(len(flat_data), sample_size, replace=False)
                    return flat_data[sample_indices]
                
                return original_data
                
            else:
                # Άλλοι τύποι δεδομένων
                return original_data
                
        except Exception as e:
            st.error(f"❌ Αδυναμία δημιουργίας fallback data: {str(e)}")
            return None

# Global instance
error_handler = RobustErrorHandler()

# Convenience functions
def safe_execute(func_name="Λειτουργία"):
    """Shorthand για safe execution"""
    return error_handler.safe_execute(func_name)

def memory_safe_operation(operation_name="Λειτουργία"):
    """Shorthand για memory safe operation"""
    return error_handler.memory_safe_operation(operation_name)

def emergency_cleanup():
    """Shorthand για emergency cleanup"""
    return error_handler._emergency_cleanup()

def create_fallback_data(original_data, operation="visualization"):
    """Shorthand για fallback data creation"""
    return error_handler.create_fallback_data(original_data, operation)
