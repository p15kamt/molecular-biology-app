"""
Διαδραστική Εφαρμογή Μοριακής Βιολογίας
Ανάλυση δεδομένων scRNA-seq με Streamlit


Ημερομηνία: 2025
"""

import streamlit as st
import pandas as pd
import sys
import os
from pathlib import Path

# Προσθέτουμε το modules directory στο Python path
sys.path.append(str(Path(__file__).parent / "modules"))
sys.path.append(str(Path(__file__).parent / "utils"))

# Imports από τα δικά μας modules
try:
    from data_preprocessing import DataPreprocessingPage
    from data_integration import DataIntegrationPageComplete as DataIntegrationPage
    from deg_analysis import DEGAnalysisPageComplete as DEGAnalysisPage
    from visualization import VisualizationPageComplete as VisualizationPage
    from cell_annotation import CellAnnotationPageComplete as CellAnnotationPage
    from team_info import TeamInfoPage
    from advanced_memory import memory_manager
except ImportError as e:
    st.error(f"Σφάλμα κατά την εισαγωγή modules: {e}")
    st.stop()

# Ρύθμιση της σελίδας
st.set_page_config(
    page_title="Μοριακή Βιολογία - scRNA-seq Ανάλυση",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS για καλύτερη εμφάνιση
st.markdown("""
<style>
    .main > div {
        padding-top: 2rem;
    }
    .stSelectbox > div > div > select {
        background-color: #f0f2f6;
    }
    .upload-section {
        border: 2px dashed #ccc;
        border-radius: 10px;
        padding: 20px;
        text-align: center;
        margin: 20px 0;
    }
    .metric-card {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #007acc;
    }
    h1 {
        color: #007acc;
        text-align: center;
        margin-bottom: 2rem;
    }
    .sidebar-header {
        font-size: 1.2rem;
        font-weight: bold;
        color: #007acc;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

def main():
    """Κεντρική συνάρτηση της εφαρμογής"""
    
    # Τίτλος εφαρμογής
    st.title("🧬 Διαδραστική Εφαρμογή Μοριακής Βιολογίας")
    st.markdown("### Ανάλυση δεδομένων scRNA-seq με σύγχρονα εργαλεία Bioinformatics")
    
    # Sidebar navigation
    st.sidebar.markdown('<div class="sidebar-header">📋 Μενού Πλοήγησης</div>', unsafe_allow_html=True)
    
    # Επιλογές σελίδων
    page_options = {
        "🏠 Αρχική Σελίδα": "home",
        "📊 Προεπεξεργασία Δεδομένων": "preprocessing", 
        "🔗 Ενοποίηση Δεδομένων": "integration",
        "🧬 Ανάλυση Διαφορικής Έκφρασης": "deg_analysis",
        "📈 Οπτικοποιήσεις": "visualization",
        "🏷️ Σχολιασμός Κυττάρων": "cell_annotation",
        "👥 Πληροφορίες Ομάδας": "team_info"
    }
    
    selected_page = st.sidebar.radio(
        "Επιλέξτε σελίδα:",
        list(page_options.keys()),
        index=0
    )
    
    # Memory Monitoring Widget
    st.sidebar.markdown("---")
    st.sidebar.markdown("### 💾 Memory Monitor")
    
    # Real-time memory information
    try:
        memory_info = memory_manager.get_system_memory_info()
        
        # Memory usage bar
        memory_percent = memory_info['percent']
        if memory_percent > 85:
            color = "🔴"
        elif memory_percent > 70:
            color = "🟡"
        else:
            color = "🟢"
        
        st.sidebar.metric(
            f"{color} System RAM", 
            f"{memory_percent:.1f}%",
            f"{memory_info['available_gb']:.1f}GB free"
        )
        
        st.sidebar.metric(
            "Process Memory", 
            f"{memory_info['process_mb']:.0f}MB",
            f"{memory_info['process_percent']:.1f}% of system"
        )
        
        # Memory cleanup button
        if st.sidebar.button("🧹 Cleanup Memory", help="Καθαρισμός μνήμης για καλύτερη απόδοση"):
            memory_manager.cleanup_session_memory()
            st.experimental_rerun()
        
        # Emergency cleanup warning
        if memory_percent > 90:
            st.sidebar.error("⚠️ Κρίσιμα χαμηλή μνήμη!")
            if st.sidebar.button("🚨 Emergency Cleanup"):
                memory_manager.emergency_memory_cleanup()
                st.experimental_rerun()
                
    except Exception as e:
        st.sidebar.warning(f"Memory monitoring unavailable: {str(e)}")
    
    # Πληροφορίες στο sidebar
    st.sidebar.markdown("---")
    st.sidebar.markdown("### 📖 Οδηγίες Χρήσης")
    st.sidebar.markdown("""
    1. **Προεπεξεργασία**: Ανεβάστε και φιλτράρετε δεδομένα
    2. **Ενοποίηση**: Συνδυάστε πολλαπλά datasets
    3. **Ανάλυση**: Εκτελέστε DEG analysis
    4. **Οπτικοποίηση**: Δημιουργήστε plots
    5. **Σχολιασμός**: Αναγνωρίστε τύπους κυττάρων
    """)
    
    # Routing στις διάφορες σελίδες
    page_key = page_options[selected_page]
    
    if page_key == "home":
        show_home_page()
    elif page_key == "preprocessing":
        preprocessing_page = DataPreprocessingPage()
        preprocessing_page.render()
    elif page_key == "integration":
        integration_page = DataIntegrationPage()
        integration_page.render()
    elif page_key == "deg_analysis":
        deg_page = DEGAnalysisPage()
        deg_page.render()
    elif page_key == "visualization":
        viz_page = VisualizationPage()
        viz_page.render()
    elif page_key == "cell_annotation":
        annotation_page = CellAnnotationPage()
        annotation_page.render()
    elif page_key == "team_info":
        team_page = TeamInfoPage()
        team_page.render()

def show_home_page():
    """Εμφάνιση της αρχικής σελίδας"""
    
    st.markdown("## Καλώς ήρθατε στην εφαρμογή ανάλυσης scRNA-seq!")
    
    # Δημιουργία columns για καλύτερη διάταξη
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.markdown("""
        ### 🎯 Στόχος της Εφαρμογής
        
        Αυτή η διαδραστική εφαρμογή σας επιτρέπει να εκτελέσετε ολοκληρωμένη ανάλυση 
        δεδομένων single-cell RNA sequencing (scRNA-seq) χρησιμοποιώντας σύγχρονα 
        εργαλεία bioinformatics και machine learning.
        
        ### 🔬 Βασικές Λειτουργίες
        
        - **Προεπεξεργασία Δεδομένων**: Quality control, filtering, normalization
        - **Ενοποίηση Δεδομένων**: Batch correction με Scanorama
        - **Clustering & Dimensionality Reduction**: UMAP, t-SNE, PCA
        - **Differential Gene Expression**: Στατιστική ανάλυση γονιδιακής έκφρασης
        - **Οπτικοποιήσεις**: Διαδραστικά plots με Plotly
        - **Cell Type Annotation**: Αυτοματοποιημένος σχολιασμός κυττάρων
        
        ### 📋 Υποστηριζόμενα Formats
        
        - **H5AD**: AnnData objects (προτιμώμενο)
        - **CSV**: Comma-separated values
        - **TSV**: Tab-separated values
        - **Excel**: .xlsx και .xls αρχεία
        """)
        
    with col2:
        st.markdown("### 📊 Quick Stats")
        
        # Fake stats για demonstration
        st.markdown('<div class="metric-card">', unsafe_allow_html=True)
        st.metric("Αλγόριθμοι ML", "15+", "")
        st.markdown('</div>', unsafe_allow_html=True)
        
        st.markdown('<div class="metric-card">', unsafe_allow_html=True)
        st.metric("Τύποι Visualization", "10+", "")
        st.markdown('</div>', unsafe_allow_html=True)
        
        st.markdown('<div class="metric-card">', unsafe_allow_html=True)
        st.metric("Supported Formats", "4", "")
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Workflow diagram placeholder
        st.markdown("### 🔄 Workflow")
        st.markdown("""
        ```
        📁 Data Upload
            ⬇️
        🔍 Quality Control
            ⬇️
        🧹 Preprocessing
            ⬇️
        🔗 Integration
            ⬇️
        📊 Analysis
            ⬇️
        📈 Visualization
        ```
        """)
    
    # Tutorial section
    st.markdown("---")
    st.markdown("## 🎓 Οδηγός Χρήσης")
    
    with st.expander("📖 Βήμα προς Βήμα Οδηγίες"):
        st.markdown("""
        ### Βήμα 1: Προετοιμασία Δεδομένων
        1. Προετοιμάστε τα δεδομένα σας σε H5AD ή CSV format
        2. Βεβαιωθείτε ότι τα δεδομένα περιέχουν gene expression matrix
        3. Προσθέστε metadata (αν υπάρχουν)
        
        ### Βήμα 2: Προεπεξεργασία
        1. Πηγαίνετε στη σελίδα "Προεπεξεργασία Δεδομένων"
        2. Ανεβάστε το αρχείο σας
        3. Ρυθμίστε τις παραμέτρους filtering
        4. Εκτελέστε quality control
        
        ### Βήμα 3: Ανάλυση
        1. Επιλέξτε τον τύπο ανάλυσης που θέλετε
        2. Ρυθμίστε τις παραμέτρους
        3. Εκτελέστε την ανάλυση
        4. Εξετάστε τα αποτελέσματα
        
        ### Βήμα 4: Οπτικοποίηση
        1. Δημιουργήστε plots για τα αποτελέσματά σας
        2. Προσαρμόστε τις παραμέτρους visualization
        3. Κατεβάστε τα plots σε υψηλή ανάλυση
        """)
    
    # Contact information
    st.markdown("---")
    st.markdown("### 📞 Επικοινωνία & Υποστήριξη")
    st.markdown("Για ερωτήσεις και υποστήριξη, επικοινωνήστε με την ομάδα ανάπτυξης στη σελίδα 'Πληροφορίες Ομάδας'.")

if __name__ == "__main__":
    main()
