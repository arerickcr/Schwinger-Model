$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
try {
    foreach ($pass in 1..2) {
        & pdflatex --disable-installer -interaction=nonstopmode -halt-on-error Schwinger_GPU_guide.tex
        if ($LASTEXITCODE -ne 0) { throw "LaTeX build failed. See Schwinger_GPU_guide.log." }
    }
} finally {
    Pop-Location
}
