"""
Utils package για τη διαδραστική εφαρμογή μοριακής βιολογίας

Αυτό το package περιέχει βοηθητικές συναρτήσεις και utilities
για τη βελτιστοποίηση της απόδοσης της εφαρμογής.


Ημερομηνία: 2025
"""

# Import των κυρίων utilities
try:
    from .memory_utils import (
        memory_monitor, 
        display_memory_info, 
        optimize_adata_memory, 
        suggest_optimization,
        cleanup_memory,
        get_memory_usage,
        get_available_memory
    )
    
    __all__ = [
        'memory_monitor',
        'display_memory_info', 
        'optimize_adata_memory',
        'suggest_optimization',
        'cleanup_memory',
        'get_memory_usage',
        'get_available_memory'
    ]
    
except ImportError as e:
    # Fallback functions αν υπάρχει πρόβλημα με imports
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
    def get_memory_usage():
        return 0
    def get_available_memory():
        return 1000
    
    __all__ = [
        'memory_monitor',
        'display_memory_info', 
        'optimize_adata_memory',
        'suggest_optimization',
        'cleanup_memory',
        'get_memory_usage',
        'get_available_memory'
    ]
