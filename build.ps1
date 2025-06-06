# Build and upload script for vgrid
Write-Host "=== Starting vgrid build and upload process ===" -ForegroundColor Green

try {
    # Clean previous builds
    Write-Host "Cleaning previous builds..." -ForegroundColor Yellow
    if (Test-Path "dist") { Remove-Item -Recurse -Force "dist" }
    if (Test-Path "build") { Remove-Item -Recurse -Force "build" }
    if (Test-Path "*.egg-info") { Remove-Item -Recurse -Force "*.egg-info" }

    # Build the package
    Write-Host "Building package..." -ForegroundColor Yellow
    python -m build
    if ($LASTEXITCODE -ne 0) { throw "Build failed" }

    # Upload to TestPyPI (using .pypirc credentials)
    Write-Host "Uploading to TestPyPI..." -ForegroundColor Yellow
    twine upload --config-file .pypirc --repository testpypi dist/*
    if ($LASTEXITCODE -ne 0) { throw "TestPyPI upload failed" }

    # Test installing from TestPyPI
    Write-Host "Testing installation from TestPyPI..." -ForegroundColor Yellow
    pip install --index-url https://test.pypi.org/simple/ vgrid
    if ($LASTEXITCODE -ne 0) { throw "Test installation failed" }

    # Ask for confirmation before uploading to real PyPI
    $confirmation = Read-Host "Do you want to upload to the real PyPI? (y/n)"
    if ($confirmation -eq 'y') {
        Write-Host "Uploading to PyPI..." -ForegroundColor Yellow
        twine upload --config-file .pypirc --repository pypi dist/*
        if ($LASTEXITCODE -ne 0) { throw "PyPI upload failed" }
        Write-Host "Upload complete!" -ForegroundColor Green
    } else {
        Write-Host "Skipping PyPI upload" -ForegroundColor Red
    }

    Write-Host "=== Process complete! ===" -ForegroundColor Green
}
catch {
    Write-Host "Error: $_" -ForegroundColor Red
    exit 1
} 