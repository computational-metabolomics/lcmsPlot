<#
.SYNOPSIS
    Installs the Compounds Chromatograms node (lcmsPlot) into Compound Discoverer.

.DESCRIPTION
    1. Finds R 4.4.0 or later.
    2. Installs the R packages the node needs (lcmsPlot at the version this
       release was built with, and rawrr for Thermo .raw files).
    3. Copies the node into the Scripts folder of each Compound Discoverer
       version found (3.3, 3.4, 3.5, ...), replacing an older copy.
    4. Writes the paths of Rscript and of the node script into node.json.
    Run install.bat (which starts this script as administrator); running it
    again updates the node.

.PARAMETER RscriptPath
    Rscript.exe to use, instead of the newest R 4.4.0 or later found.

.PARAMETER CDScriptsDir
    Compound Discoverer Scripts folder(s) to install into, instead of those
    found under C:\Program Files\Thermo.

.PARAMETER LcmsPlotRef
    Git ref (branch or commit) of lcmsPlot to install instead of the release
    tag, for testing a release before it is tagged.

.PARAMETER SkipRPackages
    Do not install the R packages (for tests).

.PARAMETER NonInteractive
    Install into every Compound Discoverer version found, without asking.
#>
[CmdletBinding()]
param(
    [string]$RscriptPath,
    [string[]]$CDScriptsDir,
    [string]$LcmsPlotRef,
    [switch]$SkipRPackages,
    [switch]$NonInteractive
)

$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$nodeName = 'CompoundsChromatograms'
$minR = [version]'4.4.0'
$logFile = Join-Path $here 'install.log'

function Write-Log([string]$Message) {
    Write-Host $Message
    Add-Content -Path $logFile -Value $Message
}

function Stop-Install([string]$Message) {
    Write-Log ''
    Write-Log "ERROR: $Message"
    Write-Log "Nothing more was changed. See $logFile."
    exit 1
}

# Run a native command. Windows PowerShell 5.1 turns a native command's
# redirected stderr (where R writes its progress) into errors, which stop the
# script under $ErrorActionPreference = 'Stop'.
function Invoke-Native([scriptblock]$Command) {
    $saved = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try { & $Command } finally { $ErrorActionPreference = $saved }
}

# Rscript.exe of every R installed, from the registry and Program Files.
function Find-Rscript {
    $found = @()
    foreach ($key in 'HKLM:\SOFTWARE\R-core\R', 'HKCU:\SOFTWARE\R-core\R') {
        if (Test-Path $key) {
            foreach ($sub in Get-ChildItem $key -ErrorAction SilentlyContinue) {
                $path = (Get-ItemProperty $sub.PSPath -ErrorAction SilentlyContinue).InstallPath
                if ($path) { $found += Join-Path $path 'bin\Rscript.exe' }
            }
        }
    }
    foreach ($root in "$env:ProgramFiles\R", "${env:ProgramFiles(x86)}\R") {
        if ($root -and (Test-Path $root)) {
            foreach ($dir in Get-ChildItem $root -Directory -Filter 'R-*') {
                $found += Join-Path $dir.FullName 'bin\Rscript.exe'
            }
        }
    }
    $found | Where-Object { Test-Path $_ } | Sort-Object -Unique
}

# The version of the R an Rscript.exe belongs to, from its folder (R-4.4.1).
function Get-RVersion([string]$Rscript) {
    $dir = Split-Path -Leaf (Split-Path -Parent (Split-Path -Parent $Rscript))
    if ($dir -match '^R-(\d+\.\d+\.\d+)') { return [version]$Matches[1] }
    $out = Invoke-Native { & $Rscript -e 'cat(as.character(getRversion()))' 2>$null }
    if ($out -match '^\d+\.\d+\.\d+$') { return [version]$out }
    return $null
}

Set-Content -Path $logFile -Value "Compounds Chromatograms node installer, $(Get-Date -Format s)"

# ---------------------------------------------------------------------------
# The node and the lcmsPlot version it was released with
# ---------------------------------------------------------------------------
$source = Join-Path $here $nodeName
$infoPath = Join-Path $source 'node-info.json'
if (-not (Test-Path (Join-Path $source 'node.json'))) {
    Stop-Install "The $nodeName folder was not found next to install.ps1. Unzip the whole download first."
}
$info = Get-Content $infoPath -Raw | ConvertFrom-Json
Write-Log "Node $($info.node_version), lcmsPlot $($info.lcmsplot_version)"

# ---------------------------------------------------------------------------
# 1. R
# ---------------------------------------------------------------------------
Write-Log ''
Write-Log '[1/4] Looking for R 4.4.0 or later'
if (-not $RscriptPath) {
    $candidates = @(Find-Rscript | ForEach-Object {
        [pscustomobject]@{ Path = $_; Version = Get-RVersion $_ }
    } | Where-Object { $_.Version -and $_.Version -ge $minR } |
        Sort-Object Version -Descending)
    if ($candidates.Count -eq 0) {
        Start-Process 'https://cloud.r-project.org/bin/windows/base/'
        Stop-Install 'R 4.4.0 or later is not installed. Install R from the page just opened, then run install.bat again.'
    }
    $RscriptPath = $candidates[0].Path
}
if (-not (Test-Path $RscriptPath)) { Stop-Install "Rscript not found: $RscriptPath" }
Write-Log "      Using $RscriptPath"

# ---------------------------------------------------------------------------
# 2. R packages
# ---------------------------------------------------------------------------
Write-Log ''
Write-Log '[2/4] Installing the R packages (this can take 10-20 minutes the first time)'
if ($SkipRPackages) {
    Write-Log '      Skipped.'
} else {
    Invoke-Native {
        & $RscriptPath (Join-Path $here 'install-packages.R') $info.lcmsplot_version $LcmsPlotRef 2>&1 |
            ForEach-Object {
                # R's stderr lines arrive as error records.
                if ($_ -is [System.Management.Automation.ErrorRecord]) { $_.Exception.Message } else { "$_" }
            } | Tee-Object -FilePath $logFile -Append | Write-Host
    }
    if ($LASTEXITCODE -ne 0) {
        Stop-Install 'The R packages could not be installed (see the messages above).'
    }
}

# ---------------------------------------------------------------------------
# 3. Compound Discoverer
# ---------------------------------------------------------------------------
Write-Log ''
Write-Log '[3/4] Looking for Compound Discoverer'
if (-not $CDScriptsDir) {
    $CDScriptsDir = @(Get-ChildItem "$env:ProgramFiles\Thermo" -Directory `
            -Filter 'Compound Discoverer*' -ErrorAction SilentlyContinue |
        ForEach-Object { Join-Path $_.FullName 'Tools\Scripts' } |
        Where-Object { Test-Path $_ })
    if ($CDScriptsDir.Count -eq 0) {
        Stop-Install "No Compound Discoverer installation was found under $env:ProgramFiles\Thermo."
    }
    if ($CDScriptsDir.Count -gt 1 -and -not $NonInteractive) {
        Write-Host '      Compound Discoverer versions found:'
        for ($i = 0; $i -lt $CDScriptsDir.Count; $i++) {
            Write-Host "        [$($i + 1)] $($CDScriptsDir[$i])"
        }
        $answer = Read-Host '      Install into which? Numbers separated by commas, or Enter for all'
        if ($answer.Trim()) {
            $CDScriptsDir = @($answer.Split(',') | ForEach-Object {
                $CDScriptsDir[[int]$_.Trim() - 1] })
        }
    }
}

# ---------------------------------------------------------------------------
# 4. The node
# ---------------------------------------------------------------------------
Write-Log ''
Write-Log '[4/4] Installing the node'
foreach ($scripts in $CDScriptsDir) {
    New-Item -ItemType Directory -Force -Path $scripts | Out-Null
    $dest = Join-Path $scripts $nodeName
    if (Test-Path $dest) { Remove-Item -Recurse -Force $dest }
    Copy-Item -Recurse $source $dest

    $nodeJson = Join-Path $dest 'node.json'
    $node = Get-Content $nodeJson -Raw | ConvertFrom-Json
    $script = Join-Path $dest "$nodeName.R"
    $node.ScriptProcessorArguments.ExecutablePath = $RscriptPath
    $node.ScriptProcessorArguments.ExecutableCommandLineArguments = "`"$script`" %NODEARGS%"
    # UTF-8 without a byte-order mark.
    [System.IO.File]::WriteAllText($nodeJson, ($node | ConvertTo-Json -Depth 20),
        (New-Object System.Text.UTF8Encoding $false))
    Write-Log "      Installed into $dest"
}

Write-Log ''
Write-Log 'Done. Two steps remain, in Compound Discoverer:'
Write-Log '  1. Help > License Manager > Scan for Missing Features, and confirm.'
Write-Log '  2. Restart Compound Discoverer. The node is under Workflow Nodes >'
Write-Log '     Scripting nodes, as Compounds Chromatograms.'
exit 0
