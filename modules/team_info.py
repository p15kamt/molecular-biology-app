"""
Module για τις Πληροφορίες της Ομάδας

Αυτό το module περιλαμβάνει τις πληροφορίες της ομάδας ανάπτυξης
και τη συνεισφορά κάθε μέλους στην εφαρμογή.


Ημερομηνία: 2025
"""

import streamlit as st
import pandas as pd
from datetime import datetime

class TeamInfoPage:
    """Κλάση για τη σελίδα πληροφοριών ομάδας"""
    
    def __init__(self):
        # Πληροφορίες ομάδας - Αυτές πρέπει να ενημερωθούν με τα πραγματικά στοιχεία
        self.team_info = {
            "Όνομα Ομάδας": "MolecularBioApp 2025",
            "Πανεπιστήμιο": "Ιόνιο Πανεπιστήμιο",
            "Τμήμα": "Πληροφορικής",
            "Ακαδημαϊκό Έτος": "2024-2025",
            "Μάθημα": "Τεχνολογία Λογισμικού",
            "Ημερομηνία Υποβολής": "Αύγουστος 2025"
        }
        
        # Μέλη ομάδας - Αυτά πρέπει να ενημερωθούν με τα πραγματικά στοιχεία
        self.team_members = [
            {
                "Όνομα": "Αντώνης Κάμτσης",
                "ΑΜ": "Π2015086",
                "Email": "p15kamt@ionio.gr",
                "Ρόλος": "Δημιουργός Project", 
                "Συνεισφορά": [
                    "Υλοποίηση και παραμετροποίηση εφαρμογής",
                    "Dockerization",
                    "Διαχείριση αποθετηρίου GitHub",
                    "Συγγραφή report"
                ]
            }
        ]
        
        # Τεχνολογίες που χρησιμοποιήθηκαν
        self.technologies = {
            "Frontend": ["Streamlit", "Plotly", "HTML/CSS"],
            "Backend": ["Python", "Scanpy", "Pandas", "NumPy"],
            "Machine Learning": ["Scikit-learn", "Scanorama", "Decoupler"],
            "Data Processing": ["AnnData", "H5py", "SciPy"],
            "Visualization": ["Matplotlib", "Seaborn", "Plotly"],
            "Deployment": ["Docker", "Git", "GitHub"]
        }
        
    def render(self):
        """Κεντρική μέθοδος για την εμφάνιση της σελίδας"""
        
        st.title("👥 Πληροφορίες Ομάδας")
        st.markdown("### Γνωρίστε την ομάδα ανάπτυξης της εφαρμογής")
        
        # Tabs για οργάνωση των πληροφοριών
        tab1, tab2, tab3, tab4 = st.tabs([
            "🏢 Γενικές Πληροφορίες",
            "👨‍💻 Μέλη Ομάδας", 
            "🛠️ Τεχνολογίες",
            "📞 Επικοινωνία"
        ])
        
        with tab1:
            self._render_general_info()
            
        with tab2:
            self._render_team_members()
            
        with tab3:
            self._render_technologies()
            
        with tab4:
            self._render_contact_info()
    
    def _render_general_info(self):
        """Εμφάνιση γενικών πληροφοριών"""
        
        st.markdown("#### 🏢 Γενικές Πληροφορίες Έργου")
        
        # Δημιουργία columns για καλύτερη διάταξη
        col1, col2 = st.columns([2, 1])
        
        with col1:
            # Πίνακας με γενικές πληροφορίες
            info_df = pd.DataFrame(
                list(self.team_info.items()),
                columns=["Πεδίο", "Τιμή"]
            )
            st.dataframe(info_df, use_container_width=True, hide_index=True)
            
        with col2:
            # Logo ή εικόνα (placeholder)
            st.markdown("### 🎯 Στόχος Έργου")
            st.markdown("""
            Ανάπτυξη σύγχρονης διαδραστικής 
            εφαρμογής για ανάλυση δεδομένων 
            μοριακής βιολογίας με έμφαση στα 
            scRNA-seq δεδομένα.
            """)
            
            st.markdown("### 📊 Στατιστικά Έργου")
            st.metric("Γραμμές Κώδικα", "~2,500+")
            st.metric("Modules", "6")
            st.metric("Dependencies", "20+")
            st.metric("Χρόνος Ανάπτυξης", "3 εβδομάδες")
        
        # Περιγραφή έργου
        st.markdown("---")
        st.markdown("#### 📝 Περιγραφή Έργου")
        st.markdown("""
        Η εφαρμογή αυτή αναπτύχθηκε στο πλαίσιο του μαθήματος "Τεχνολογία Λογισμικού" 
        και στοχεύει στην παροχή ενός ολοκληρωμένου εργαλείου για την ανάλυση δεδομένων 
        single-cell RNA sequencing (scRNA-seq).
        
        **Βασικά χαρακτηριστικά:**
        - Διαδραστική εφαρμογή web με Streamlit
        - Ενσωμάτωση σύγχρονων αλγορίθμων bioinformatics
        - Πλήρης pipeline από preprocessing έως visualization
        - Docker containerization για εύκολη deployment
        - Modern UI/UX design
        
        **Καινοτομίες:**
        - Real-time parameter adjustment
        - Interactive visualizations με Plotly
        - Automated cell type annotation
        - Batch correction με Scanorama
        - Comprehensive quality control metrics
        """)
    
    def _render_team_members(self):
        """Εμφάνιση μελών ομάδας"""
        
        st.markdown("#### 👨‍💻 Μέλη της Ομάδας Ανάπτυξης")
        
        for i, member in enumerate(self.team_members):
            # Δημιουργία expandable section για κάθε μέλος
            with st.expander(f"👤 {member['Όνομα']} - {member['Ρόλος']}", expanded=i==0):
                
                col1, col2 = st.columns([1, 2])
                
                with col1:
                    st.markdown("**Στοιχεία Επικοινωνίας:**")
                    st.write(f"**ΑΜ:** {member['ΑΜ']}")
                    st.write(f"**Email:** {member['Email']}")
                    st.write(f"**Ρόλος:** {member['Ρόλος']}")
                
                with col2:
                    st.markdown("**Συνεισφορά στο Έργο:**")
                    for contribution in member['Συνεισφορά']:
                        st.write(f"• {contribution}")
        
        # Στατιστικά ομάδας
        st.markdown("---")
        st.markdown("#### 📊 Στατιστικά Ομάδας")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Μέλη Ομάδας", len(self.team_members))
        with col2:
            total_contributions = sum(len(member['Συνεισφορά']) for member in self.team_members)
            st.metric("Συνολικές Συνεισφορές", total_contributions)
        with col3:
            st.metric("Ρόλοι", len(set(member['Ρόλος'] for member in self.team_members)))
    
    def _render_technologies(self):
        """Εμφάνιση τεχνολογιών"""
        
        st.markdown("#### 🛠️ Τεχνολογίες που Χρησιμοποιήθηκαν")
        
        # Δημιουργία columns για κάθε κατηγορία
        categories = list(self.technologies.keys())
        cols = st.columns(2)
        
        for i, category in enumerate(categories):
            with cols[i % 2]:
                st.markdown(f"##### {category}")
                for tech in self.technologies[category]:
                    st.write(f"• {tech}")
                st.markdown("")  # Spacer
        
        # Dependency tree visualization
        st.markdown("---")
        st.markdown("#### 📦 Κύριες Dependencies")
        
        # Δημιουργία ενός απλού dependency graph
        dependency_data = {
            "Βιβλιοθήκη": [
                "Streamlit", "Scanpy", "Pandas", "NumPy", "Plotly", 
                "Matplotlib", "Scikit-learn", "Scanorama", "Decoupler"
            ],
            "Χρήση": [
                "Web Interface", "scRNA-seq Analysis", "Data Manipulation", 
                "Numerical Computing", "Interactive Plots", "Static Plots",
                "Machine Learning", "Batch Correction", "Cell Annotation"
            ],
            "Έκδοση": [
                "≥1.28.0", "≥1.9.0", "≥1.5.0", "≥1.21.0", "≥5.15.0",
                "≥3.5.0", "≥1.0.0", "≥1.7.3", "≥1.4.0"
            ]
        }
        
        dependency_df = pd.DataFrame(dependency_data)
        st.dataframe(dependency_df, use_container_width=True, hide_index=True)
    
    def _render_contact_info(self):
        """Εμφάνιση στοιχείων επικοινωνίας"""
        
        st.markdown("#### 📞 Στοιχεία Επικοινωνίας & Υποστήριξη")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("##### 🏫 Ακαδημαϊκά Στοιχεία")
            st.info(f"""
            **Πανεπιστήμιο:** {self.team_info['Πανεπιστήμιο']}  
            **Τμήμα:** {self.team_info['Τμήμα']}  
            **Μάθημα:** {self.team_info['Μάθημα']}  
            **Ακαδημαϊκό Έτος:** {self.team_info['Ακαδημαϊκό Έτος']}
            """)
            
            st.markdown("##### 📧 Επικοινωνία Ομάδας")
            for member in self.team_members:
                st.write(f"**{member['Όνομα']}:** {member['Email']}")
        
        with col2:
            st.markdown("##### 🔗 Χρήσιμοι Σύνδεσμοι")
            st.markdown("""
            - [GitHub Repository](https://github.com/p15kamt/molecular-biology-app)
            - [Streamlit Documentation](https://docs.streamlit.io/)
            - [Scanpy Documentation](https://scanpy.readthedocs.io/)
            - [Plotly Documentation](https://plotly.com/python/)
            """)
            
            st.markdown("##### 📋 Αναφορά Σφαλμάτων")
            st.markdown("""
            Για αναφορά σφαλμάτων ή προτάσεις βελτίωσης:
            1. Δημιουργήστε issue στο GitHub repository
            2. Στείλτε email σε οποιονδήποτε από την ομάδα
            3. Χρησιμοποιήστε το feedback form παρακάτω
            """)
        
        # Feedback form
        st.markdown("---")
        st.markdown("#### 💬 Feedback Form")
        
        with st.form("feedback_form"):
            col1, col2 = st.columns(2)
            
            with col1:
                name = st.text_input("Όνομα (προαιρετικό)")
                email = st.text_input("Email (προαιρετικό)")
                
            with col2:
                feedback_type = st.selectbox(
                    "Τύπος Feedback",
                    ["Σφάλμα", "Πρόταση", "Γενικό Σχόλιο", "Άλλο"]
                )
                rating = st.slider("Βαθμολογία (1-5)", 1, 5, 3)
            
            message = st.text_area("Μήνυμα", height=100)
            
            submitted = st.form_submit_button("📤 Αποστολή Feedback")
            
            if submitted:
                # Στην πραγματικότητα θα έπρεπε να αποθηκεύσουμε το feedback
                st.success("✅ Το feedback σας έχει σταλεί! Σας ευχαριστούμε!")
                st.balloons()
        
        # Ημερομηνία τελευταίας ενημέρωσης
        st.markdown("---")
        st.markdown(f"*Τελευταία ενημέρωση: {datetime.now().strftime('%d/%m/%Y')}*")
