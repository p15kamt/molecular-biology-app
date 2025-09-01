# PowerShell script for launching the app with advanced memory management

Write-Host "Molecular Biology Interactive App" -ForegroundColor Green
Write-Host "Large Files Edition v1.1" -ForegroundColor Cyan
Write-Host "=======================================" -ForegroundColor Yellow

# Check virtual environment
Write-Host "Activating virtual environment..." -ForegroundColor Blue
try {
    & ".\venv\Scripts\Activate.ps1"
    Write-Host "Virtual environment activated successfully" -ForegroundColor Green
} catch {
    Write-Host "Error activating venv: $($_.Exception.Message)" -ForegroundColor Red
    Write-Host "Try: python -m venv venv; pip install -r requirements.txt" -ForegroundColor Yellow
    exit 1
}

# Check memory
$memory = Get-WmiObject -Class Win32_ComputerSystem
$totalRAM = [math]::Round($memory.TotalPhysicalMemory / 1GB, 1)
Write-Host "Total system memory: $totalRAM GB" -ForegroundColor Cyan

if ($totalRAM -lt 4) {
    Write-Host "Warning: Low memory (<4GB). Use subsample for large files." -ForegroundColor Yellow
} elseif ($totalRAM -ge 8) {
    Write-Host "Excellent memory for large files" -ForegroundColor Green
} else {
    Write-Host "Sufficient memory" -ForegroundColor Green
}

Write-Host ""
Write-Host "Starting app with advanced memory management..." -ForegroundColor Blue
Write-Host "App will open at: http://localhost:8505" -ForegroundColor Magenta
Write-Host "File support up to 1GB" -ForegroundColor Cyan
Write-Host "Advanced memory management: ON" -ForegroundColor Green
Write-Host "Press Ctrl+C to stop" -ForegroundColor Yellow
Write-Host ""

# Start with optimizations
try {
    streamlit run app.py --server.port=8505 --server.address=0.0.0.0 --server.maxUploadSize=1000 --server.maxMessageSize=1000 --server.runOnSave=false --server.fileWatcherType=none
} catch {
    Write-Host "Error starting app: $($_.Exception.Message)" -ForegroundColor Red
    Write-Host "Try: pip install streamlit" -ForegroundColor Yellow
}

Write-Host ""
Write-Host "Application terminated" -ForegroundColor Yellow
Read-Host "Press Enter to exit"
