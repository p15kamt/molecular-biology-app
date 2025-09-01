@echo off
echo 🧬 Διαδραστική Εφαρμογή Μοριακής Βιολογίας
echo 🗂️  Έκδοση για Μεγάλα Αρχεία
echo =====================================

echo 🔄 Ενεργοποίηση virtual environment...
call venv\Scripts\activate.bat

echo 🚀 Εκκίνηση εφαρμογής...
echo 🌐 Η εφαρμογή θα ανοίξει στο: http://localhost:8505
echo 📁 Υποστήριξη αρχείων έως 1GB
echo ⏹️  Για τερματισμό πατήστε Ctrl+C
echo.

streamlit run app.py --server.port=8505 --server.address=0.0.0.0 --server.maxUploadSize=1000 --server.maxMessageSize=1000

pause
