param([Parameter(Mandatory=$true)][string]$PythonExe)
$ErrorActionPreference = 'Stop'
$projectPath = $PSScriptRoot
$env:TCGA_WORKDIR = Join-Path $projectPath 'data'
$env:MPLCONFIGDIR = Join-Path $projectPath 'data/mplconfig'
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$version = & $PythonExe -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
if ($version -ne '3.12') { throw 'Use Python 3.12 for the frozen environment.' }
& $PythonExe -m venv (Join-Path $projectPath '.venv')
if ($LASTEXITCODE -ne 0) { throw 'Environment creation failed.' }
$analysisPython = Join-Path $projectPath '.venv/Scripts/python.exe'
& $analysisPython -m pip install -r (Join-Path $projectPath 'requirements-lock.txt')
if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed.' }
foreach ($stage in @('00_download_references.py','01_prepare.py','02_expression.py','03_followup.py','04_validate.py','05_report.py')) {
    & $analysisPython (Join-Path $projectPath "scripts/$stage")
    if ($LASTEXITCODE -ne 0) { throw "Analysis stopped at $stage" }
}
Write-Output 'Completed. Open Research_Report.html.'
