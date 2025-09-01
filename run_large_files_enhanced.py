#!/usr/bin/env python3
"""
Βελτιωμένο Script για εκκίνηση εφαρμογής με υποστήριξη μεγάλων αρχείων

Αυτό το script βελτιστοποιεί τις ρυθμίσεις για καλύτερη διαχείριση μνήμης
και προσθέτει προηγμένες λειτουργίες για μεγάλα scRNA-seq datasets.

Χρήση:
    python run_large_files_enhanced.py


Ημερομηνία: 2025
"""

import sys
import os
import subprocess
import psutil
import tempfile
from pathlib import Path

def check_system_resources():
    """Έλεγχος πόρων συστήματος και προτάσεις βελτίωσης"""
    
    print("🔍 Έλεγχος πόρων συστήματος...")
    
    # Memory check
    memory = psutil.virtual_memory()
    memory_gb = memory.total / (1024**3)
    available_gb = memory.available / (1024**3)
    
    print(f"💾 Συνολική μνήμη: {memory_gb:.1f} GB")
    print(f"💾 Διαθέσιμη μνήμη: {available_gb:.1f} GB")
    print(f"💾 Χρήση μνήμης: {memory.percent:.1f}%")
    
    # CPU check
    cpu_count = psutil.cpu_count()
    print(f"🖥️  CPU cores: {cpu_count}")
    
    # Disk space check
    disk = psutil.disk_usage('.')
    disk_free_gb = disk.free / (1024**3)
    print(f"💽 Διαθέσιμος χώρος: {disk_free_gb:.1f} GB")
    
    # Προτάσεις βάσει πόρων
    recommendations = []
    
    if available_gb < 4:
        recommendations.append("⚠️  Χαμηλή διαθέσιμη μνήμη (<4GB) - συνιστάται κλείσιμο άλλων εφαρμογών")
    elif available_gb < 8:
        recommendations.append("💡 Μέτρια μνήμη - θα χρησιμοποιηθεί streaming mode για μεγάλα αρχεία")
    else:
        recommendations.append("✅ Επαρκής μνήμη για μεγάλα datasets")
    
    if disk_free_gb < 5:
        recommendations.append("⚠️  Χαμηλός διαθέσιμος χώρος (<5GB) - μπορεί να υπάρξουν προβλήματα με temp files")
    
    if memory.percent > 80:
        recommendations.append("🧹 Υψηλή χρήση μνήμης - συνιστάται restart του συστήματος")
    
    for rec in recommendations:
        print(rec)
    
    return {
        'memory_gb': memory_gb,
        'available_gb': available_gb,
        'cpu_count': cpu_count,
        'disk_free_gb': disk_free_gb,
        'memory_percent': memory.percent
    }

def optimize_streamlit_config():
    """Δημιουργία βελτιστοποιημένου Streamlit config (χωρίς deprecated options)"""
    
    # Δεν δημιουργούμε config στο home directory - χρησιμοποιούμε το τοπικό
    print("⚙️  Χρήση τοπικού Streamlit config (χωρίς deprecated options)")

def set_memory_optimizations():
    """Ρύθμιση environment variables για βελτίωση μνήμης"""
    
    # Python memory optimizations
    os.environ['PYTHONHASHSEED'] = '0'
    os.environ['PYTHONUNBUFFERED'] = '1'
    
    # NumPy optimizations
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
    os.environ['OMP_NUM_THREADS'] = '1'
    
    # Matplotlib backend για memory efficiency
    os.environ['MPLBACKEND'] = 'Agg'
    
    print("🔧 Memory optimizations εφαρμόστηκαν")

def create_temp_directory():
    """Δημιουργία temp directory για μεγάλα αρχεία"""
    
    temp_dir = Path(tempfile.gettempdir()) / 'molecular_biology_app'
    temp_dir.mkdir(exist_ok=True)
    
    os.environ['TMPDIR'] = str(temp_dir)
    
    print(f"📁 Temp directory: {temp_dir}")

def run_streamlit_optimized():
    """Εκκίνηση Streamlit με βελτιστοποιημένες ρυθμίσεις"""
    
    print("\n🚀 Εκκίνηση βελτιστοποιημένης εφαρμογής...")
    print("🌐 Η εφαρμογή θα ανοίξει στη διεύθυνση: http://localhost:8501")
    print("📊 Υποστήριξη αρχείων έως 2GB")
    print("🔄 Automatic streaming mode για μεγάλα datasets")
    print("⏹️  Για τερματισμό πατήστε Ctrl+C")
    print("-" * 60)
    
    try:
        # Εκτέλεση με βελτιστοποιημένες ρυθμίσεις
        subprocess.run([
            sys.executable, '-m', 'streamlit', 'run', 'app.py',
            '--server.port=8501',
            '--server.address=localhost',
            '--server.maxUploadSize=2048',
            '--server.maxMessageSize=2048'
        ])
        
    except KeyboardInterrupt:
        print("\n\n👋 Η εφαρμογή τερματίστηκε από τον χρήστη")
    except Exception as e:
        print(f"\n❌ Σφάλμα κατά την εκκίνηση: {e}")
        print("\n💡 Προτάσεις επίλυσης:")
        print("   • Ελέγξτε ότι όλες οι dependencies είναι εγκατεστημένες")
        print("   • Δοκιμάστε να κλείσετε άλλες εφαρμογές για εξοικονόμηση μνήμης")
        print("   • Χρησιμοποιήστε το run_app.py για βασική εκκίνηση")

def main():
    """Κεντρική συνάρτηση"""
    
    print("🧬 Διαδραστική Εφαρμογή Μοριακής Βιολογίας")
    print("📊 Βελτιστοποιημένη για μεγάλα scRNA-seq datasets")
    print("=" * 60)
    
    # Έλεγχος πόρων συστήματος
    system_info = check_system_resources()
    
    print("\n⚙️  Εφαρμογή βελτιστοποιήσεων...")
    
    # Βελτιστοποιήσεις
    optimize_streamlit_config()
    set_memory_optimizations()
    create_temp_directory()
    
    # Προειδοποιήσεις βάσει πόρων
    if system_info['available_gb'] < 2:
        print("\n⚠️  ΠΡΟΕΙΔΟΠΟΙΗΣΗ: Πολύ χαμηλή διαθέσιμη μνήμη!")
        print("   Συνιστάται να κλείσετε άλλες εφαρμογές πριν συνεχίσετε.")
        response = input("   Θέλετε να συνεχίσετε; (y/N): ")
        if response.lower() != 'y':
            print("👋 Εξαγωγή από την εφαρμογή")
            return
    
    print("\n✅ Όλες οι βελτιστοποιήσεις εφαρμόστηκαν επιτυχώς!")
    
    # Εκκίνηση εφαρμογής
    run_streamlit_optimized()

if __name__ == "__main__":
    main()
