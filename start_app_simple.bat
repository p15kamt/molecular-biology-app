@echo off
chcp 65001 > nul
echo.
echo ================================
echo  Molecular Biology App v1.1
echo  Large Files Support
echo ================================
echo.

echo Activating virtual environment...
call venv\Scripts\activate.bat

echo.
echo Starting application...
echo App URL: http://localhost:8505
echo File support: up to 1GB
echo Press Ctrl+C to stop
echo.

streamlit run app.py --server.port=8505 --server.address=0.0.0.0 --server.maxUploadSize=1000 --server.maxMessageSize=1000

echo.
echo Application stopped.
pause
