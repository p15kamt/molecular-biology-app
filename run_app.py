#!/usr/bin/env python3
"""
Script για εκκίνηση της εφαρμογής Μοριακής Βιολογίας

Αυτό το script ελέγχει τις dependencies και εκκινεί την εφαρμογή Streamlit.

Χρήση:
    python run_app.py


Ημερομηνία: 2025
"""

import sys
import os
import subprocess
import importlib
from pathlib import Path

def check_python_version():
    """Έλεγχος έκδοσης Python"""
    if sys.version_info < (3, 8):
        print("❌ Απαιτείται Python 3.8 ή νεότερη έκδοση")
        print(f"   Τρέχουσα έκδοση: {sys.version}")
        return False
    print(f"✅ Python έκδοση: {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")
    return True

def check_dependencies():
    """Έλεγχος εγκατεστημένων dependencies"""
    required_packages = [
        'streamlit',
        'pandas', 
        'numpy',
        'plotly',
        'matplotlib'
    ]
    
    missing_packages = []
    
    for package in required_packages:
        try:
            importlib.import_module(package)
            print(f"✅ {package}")
        except ImportError:
            missing_packages.append(package)
            print(f"❌ {package} - Δεν είναι εγκατεστημένο")
    
    if missing_packages:
        print("\n📦 Απαιτείται εγκατάσταση dependencies:")
        print("   pip install -r requirements.txt")
        return False
    
    return True

def check_app_files():
    """Έλεγχος ότι υπάρχουν τα απαραίτητα αρχεία"""
    required_files = [
        'app.py',
        'requirements.txt',
        'modules/data_preprocessing.py',
        'modules/team_info.py'
    ]
    
    missing_files = []
    
    for file_path in required_files:
        if not Path(file_path).exists():
            missing_files.append(file_path)
            print(f"❌ {file_path} - Δεν βρέθηκε")
        else:
            print(f"✅ {file_path}")
    
    if missing_files:
        print("\n🚨 Λείπουν απαραίτητα αρχεία!")
        return False
    
    return True

def run_streamlit_app():
    """Εκκίνηση της Streamlit εφαρμογής"""
    try:
        print("\n🚀 Εκκίνηση της εφαρμογής...")
        print("🌐 Η εφαρμογή θα ανοίξει στη διεύθυνση: http://localhost:8501")
        print("⏹️  Για τερματισμό πατήστε Ctrl+C")
        print("-" * 50)
        
        # Εκτέλεση Streamlit
        subprocess.run([
            sys.executable, '-m', 'streamlit', 'run', 'app.py',
            '--server.port=8501',
            '--server.address=localhost'
        ])
        
    except KeyboardInterrupt:
        print("\n\n👋 Η εφαρμογή τερματίστηκε από τον χρήστη")
    except Exception as e:
        print(f"\n❌ Σφάλμα κατά την εκκίνηση: {e}")

def main():
    """Κεντρική συνάρτηση"""
    print("🧬 Διαδραστική Εφαρμογή Μοριακής Βιολογίας")
    print("=" * 50)
    
    # Έλεγχος συστήματος
    print("\n🔍 Έλεγχος συστήματος...")
    
    if not check_python_version():
        sys.exit(1)
    
    print("\n📦 Έλεγχος dependencies...")
    if not check_dependencies():
        sys.exit(1)
    
    print("\n📁 Έλεγχος αρχείων εφαρμογής...")
    if not check_app_files():
        sys.exit(1)
    
    print("\n✅ Όλοι οι έλεγχοι πέρασαν επιτυχώς!")
    
    # Εκκίνηση εφαρμογής
    run_streamlit_app()

if __name__ == "__main__":
    main()
