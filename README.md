# 🧬 Διαδραστική Εφαρμογή Μοριακής Βιολογίας

Μια σύγχρονη διαδραστική εφαρμογή για ανάλυση δεδομένων single-cell RNA sequencing (scRNA-seq) χρησιμοποιώντας Streamlit και σύγχρονα εργαλεία bioinformatics.

## 🎯 Χαρακτηριστικά

- **Προεπεξεργασία Δεδομένων**: Quality control, filtering, normalization
- **Ενοποίηση Δεδομένων**: Batch correction με Scanorama  
- **Ανάλυση Διαφορικής Έκφρασης**: Στατιστική ανάλυση γονιδιακής έκφρασης
- **Οπτικοποιήσεις**: Διαδραστικά plots με Plotly
- **Σχολιασμός Κυττάρων**: Αυτοματοποιημένη αναγνώριση τύπων κυττάρων
- **Docker Support**: Εύκολη deployment σε οποιοδήποτε περιβάλλον

## 🚀 Γρήγορη Εκκίνηση

### Προαπαιτούμενα
- Python 3.8+ ή Docker
- 4GB+ RAM (συνιστάται)

### Τοπική Εγκατάσταση

1. **Clone το repository**
```bash
git clone https://github.com/[username]/molecular-biology-app.git
cd molecular-biology-app
```

2. **Δημιουργία virtual environment**
```bash
python -m venv venv
source venv/bin/activate  # Linux/Mac
# ή
venv\Scripts\activate     # Windows
```

3. **Εγκατάσταση dependencies**
```bash
pip install -r requirements.txt
```

4. **Εκτέλεση εφαρμογής**

   **Για κανονικά αρχεία (<100MB):**
   ```bash
   streamlit run app.py
   ```

   **Για μεγάλα αρχεία (100MB-1GB):**
   ```bash
   python run_large_files.py
   ```

5. **Άνοιγμα browser** στη διεύθυνση: http://localhost:8501

### Docker Εκτέλεση

1. **Build του image**
```bash
docker build -t molecular-biology-app .
```

2. **Εκτέλεση container**
```bash
docker run -p 8501:8501 molecular-biology-app
```

3. **Άνοιγμα browser** στη διεύθυνση: http://localhost:8501

## 📁 Δομή Project

```
molecular_biology_app/
├── app.py                 # Κεντρική εφαρμογή
├── requirements.txt       # Python dependencies
├── Dockerfile            # Docker configuration
├── README.md             # Αυτό το αρχείο
├── modules/              # Application modules
│   ├── data_preprocessing.py
│   ├── data_integration.py
│   ├── deg_analysis.py
│   ├── visualization.py
│   ├── cell_annotation.py
│   └── team_info.py
├── utils/                # Utility functions
├── data/                 # Sample data
└── assets/               # Static assets
```

## 📊 Υποστηριζόμενα Formats

- **H5AD**: AnnData objects (προτιμώμενο)
- **CSV**: Comma-separated values
- **TSV**: Tab-separated values
- **Excel**: .xlsx και .xls αρχεία

## 🗂️ Υποστήριξη Μεγάλων Αρχείων

### Μεγέθη Αρχείων
- **Μικρά αρχεία** (<100MB): Κανονική λειτουργία
- **Μεσαία αρχεία** (100MB-500MB): Βελτιστοποιημένη φόρτωση
- **Μεγάλα αρχεία** (500MB-1GB): Chunked processing με progress tracking

### Απαιτήσεις Μνήμης
- **4GB RAM**: Ελάχιστη για αρχεία ~200MB
- **8GB RAM**: Συνιστάται για αρχεία ~500MB
- **16GB RAM**: Άνετη χρήση για αρχεία 1GB+

### Βελτιστοποιήσεις
- **Automatic subsample**: Για datasets >10K κύτταρα
- **Memory monitoring**: Real-time παρακολούθηση μνήμης
- **Sparse matrix optimization**: Αυτόματη μετατροπή
- **Progressive loading**: Progress bars για μεγάλα αρχεία
- **Cache management**: Intelligent caching για ταχύτητα

## 🛠️ Τεχνολογίες

- **Frontend**: Streamlit, Plotly
- **Backend**: Python, Scanpy, Pandas
- **ML Libraries**: Scikit-learn, Scanorama, Decoupler
- **Visualization**: Matplotlib, Seaborn, Plotly
- **Deployment**: Docker

## 👥 Ομάδα Ανάπτυξης

- **Αντώνης Κάμτσης**: Δημιουργός Project


## 📖 Οδηγός Χρήσης

1. **Ανέβασμα Δεδομένων**: Επιλέξτε "Προεπεξεργασία Δεδομένων" και ανεβάστε το αρχείο σας
2. **Quality Control**: Εξετάστε τα QC metrics και ρυθμίστε τις παραμέτρους
3. **Preprocessing**: Εκτελέστε filtering και normalization
4. **Ανάλυση**: Χρησιμοποιήστε τα διάφορα modules για περαιτέρω ανάλυση
5. **Αποτελέσματα**: Εξάγετε plots και αποτελέσματα

## 🐛 Αναφορά Σφαλμάτων

Για αναφορά σφαλμάτων ή προτάσεις:
1. Δημιουργήστε issue στο GitHub
2. Χρησιμοποιήστε το feedback form στην εφαρμογή
3. Επικοινωνήστε με την ομάδα

## 📄 Άδεια

Αυτό το project αναπτύχθηκε για εκπαιδευτικούς σκοπούς στο πλαίσιο του μαθήματος "Τεχνολογία Λογισμικού".

## 🔗 Χρήσιμοι Σύνδεσμοι

- [Streamlit Documentation](https://docs.streamlit.io/)
- [Scanpy Documentation](https://scanpy.readthedocs.io/)
- [Plotly Documentation](https://plotly.com/python/)
- [Docker Documentation](https://docs.docker.com/)
