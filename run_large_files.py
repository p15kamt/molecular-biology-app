#!/usr/bin/env python3
"""
Script για εκκίνηση της εφαρμογής με υποστήριξη μεγάλων αρχείων

Αυτό το script ελέγχει τις ρυθμίσεις και εκκινεί την εφαρμογή 
με βελτιστοποιημένες παραμέτρους για μεγάλα αρχεία.

Χρήση:
    python run_large_files.py


Ημερομηνία: 2025
"""

import sys
import os
import subprocess
from pathlib import Path
import psutil

def check_system_requirements():
    """Έλεγχος απαιτήσεων συστήματος για μεγάλα αρχεία"""
    
    print("🔍 Έλεγχος συστήματος για μεγάλα αρχεία...")
    
    # Έλεγχος μνήμης
    memory = psutil.virtual_memory()
    total_memory_gb = memory.total / (1024**3)
    available_memory_gb = memory.available / (1024**3)
    
    print(f"💾 Συνολική μνήμη: {total_memory_gb:.1f} GB")
    print(f"💾 Διαθέσιμη μνήμη: {available_memory_gb:.1f} GB")
    
    if available_memory_gb < 2:
        print("⚠️  ΠΡΟΕΙΔΟΠΟΙΗΣΗ: Χαμηλή διαθέσιμη μνήμη (<2GB)")
        print("   Συνιστάται κλείσιμο άλλων εφαρμογών")
    elif available_memory_gb >= 8:
        print("✅ Άριστη μνήμη για μεγάλα αρχεία")
    else:
        print("✅ Επαρκής μνήμη")
    
    # Έλεγχος χώρου δίσκου
    disk = psutil.disk_usage('.')
    free_space_gb = disk.free / (1024**3)
    
    print(f"💿 Διαθέσιμος χώρος: {free_space_gb:.1f} GB")
    
    if free_space_gb < 1:
        print("⚠️  ΠΡΟΕΙΔΟΠΟΙΗΣΗ: Χαμηλός διαθέσιμος χώρος (<1GB)")
        return False
    
    print("✅ Επαρκής χώρος δίσκου")
    return True

def setup_environment():
    """Ρύθμιση περιβάλλοντος για μεγάλα αρχεία"""
    
    print("\n⚙️  Ρύθμιση περιβάλλοντος...")
    
    # Ρύθμιση environment variables για Streamlit
    os.environ['STREAMLIT_SERVER_MAX_UPLOAD_SIZE'] = '1000'  # 1GB
    os.environ['STREAMLIT_SERVER_ENABLE_CORS'] = 'false'
    os.environ['STREAMLIT_BROWSER_GATHER_USAGE_STATS'] = 'false'
    
    # Ρύθμιση Python memory management
    os.environ['PYTHONHASHSEED'] = '0'
    os.environ['OMP_NUM_THREADS'] = '4'  # Περιορισμός threads για σταθερότητα
    
    print("✅ Περιβάλλον ρυθμίστηκε για μεγάλα αρχεία")

def create_config_file():
    """Δημιουργία βελτιστοποιημένου config file"""
    
    config_dir = Path('.streamlit')
    config_dir.mkdir(exist_ok=True)
    
    config_content = """[server]
headless = true
port = 8501
address = "0.0.0.0"
enableCORS = false
enableXsrfProtection = false
maxUploadSize = 1000
maxMessageSize = 200

[theme]
primaryColor = "#007acc"
backgroundColor = "#ffffff"
secondaryBackgroundColor = "#f0f2f6"
textColor = "#262730"

[browser]
gatherUsageStats = false

[logger]
level = "info"

[runner]
magicEnabled = false
installTracer = false
fixMatplotlib = false
"""
    
    config_path = config_dir / 'config.toml'
    with open(config_path, 'w') as f:
        f.write(config_content)
    
    print(f"✅ Config file δημιουργήθηκε: {config_path}")

def run_streamlit():
    """Εκκίνηση Streamlit με βελτιστοποιημένες παραμέτρους"""
    
    print("\n🚀 Εκκίνηση εφαρμογής...")
    print("🌐 Η εφαρμογή θα ανοίξει στη διεύθυνση: http://localhost:8501")
    print("📁 Υποστήριξη αρχείων έως 1GB")
    print("⏹️  Για τερματισμό πατήστε Ctrl+C")
    print("-" * 60)
    
    try:
        # Εκτέλεση Streamlit με επιπλέον παραμέτρους
        subprocess.run([
            sys.executable, '-m', 'streamlit', 'run', 'app.py',
            '--server.port=8501',
            '--server.address=0.0.0.0',
            '--server.maxUploadSize=1000',
            '--server.enableCORS=false',
            '--browser.gatherUsageStats=false'
        ])
        
    except KeyboardInterrupt:
        print("\n\n👋 Η εφαρμογή τερματίστηκε από τον χρήστη")
    except Exception as e:
        print(f"\n❌ Σφάλμα κατά την εκκίνηση: {e}")

def main():
    """Κεντρική συνάρτηση"""
    
    print("🧬 Διαδραστική Εφαρμογή Μοριακής Βιολογίας")
    print("🗂️  Έκδοση για Μεγάλα Αρχεία")
    print("=" * 60)
    
    # Έλεγχος συστήματος
    if not check_system_requirements():
        print("\n❌ Το σύστημα δεν πληροί τις απαιτήσεις")
        return 1
    
    # Ρύθμιση περιβάλλοντος
    setup_environment()
    
    # Δημιουργία config file
    create_config_file()
    
    # Εκκίνηση εφαρμογής
    run_streamlit()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
