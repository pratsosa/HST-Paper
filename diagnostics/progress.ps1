# diagnostics/progress.ps1 - AP 2026-09-10
# Progress and ETA for a long pipeline run, from how fast its output files appear.
#
#   rebin:  powershell -File diagnostics\progress.ps1
#   ICA:    powershell -File diagnostics\progress.ps1 -Dir ICA_Plots_Rebin_master_v2 `
#                      -Filter *.png -Log diagnostics\ica_v2_run.log -DonePattern "ICA processing completed"

param(
    [string]$Dir    = "RebinnedSpec_master_v2",
    [string]$Filter = "*.fits",
    [string]$Log    = "diagnostics\rebin_v2_run.log",
    [string]$DonePattern = "Run complete",
    # catalog_bridge yields 615 (object, instrument) pairs, but ~21 have no
    # readable exposures and write no file, so ~594 is the expected final count.
    # The ICA run over those 594 spectra skips a handful it cannot fit.
    [int]   $Total  = 594
)

$files = @(Get-ChildItem $Dir -Filter $Filter -Recurse -ErrorAction SilentlyContinue)
$n = $files.Count

# NB: PowerShell variable names are case-insensitive, so this must NOT be
# called $done -- that would collide with the $DonePattern parameter.
$isDone = $false
if (Test-Path $Log) {
    $isDone = [bool](Select-String -Path $Log -Pattern $DonePattern -Quiet)
}
$alive = @(Get-Process python -ErrorAction SilentlyContinue).Count

if ($n -lt 2) {
    "written : $n / $Total  (too early for an estimate)"
}
else {
    $start   = ($files | Sort-Object LastWriteTime | Select-Object -First 1).LastWriteTime
    $last    = ($files | Sort-Object LastWriteTime | Select-Object -Last 1).LastWriteTime
    $elapsed = (Get-Date) - $start
    $rate    = $n / $elapsed.TotalMinutes            # spectra per minute
    $eta     = [TimeSpan]::FromMinutes(($Total - $n) / $rate)

    "written : {0} / {1}   ({2:N1}%)"          -f $n, $Total, (100 * $n / $Total)
    "elapsed : {0:hh\:mm\:ss}   rate {1:N1}/min" -f $elapsed, $rate
    "remain  : {0:hh\:mm\:ss}   finishes ~{1:HH:mm}" -f $eta, (Get-Date).Add($eta)
    "newest  : {0:HH:mm:ss}  ({1:N1} min ago)" -f $last, ((Get-Date) - $last).TotalMinutes

    # Recent rate reacts faster than the average when the run hits a slow patch.
    # Objects vary from seconds to many minutes (Mrk 817 and NGC 5548 have
    # hundreds of exposures), so treat any ETA as rough.
    $cut    = (Get-Date).AddMinutes(-5)
    $recent = @($files | Where-Object { $_.LastWriteTime -gt $cut }).Count
    if ($recent -gt 0 -and $elapsed.TotalMinutes -gt 5) {
        $rrate = $recent / 5.0
        $reta  = [TimeSpan]::FromMinutes(($Total - $n) / $rrate)
        "last 5m : {0:N1}/min  ->  remain {1:hh\:mm\:ss}, finishes ~{2:HH:mm}" -f $rrate, $reta, (Get-Date).Add($reta)
    }
}

if ($isDone)        { "status  : FINISHED - see $Log" }
elseif ($alive -gt 0) { "status  : running ($alive python processes)" }
else                { "status  : NO python process running - the run has stopped early, check $Log" }
