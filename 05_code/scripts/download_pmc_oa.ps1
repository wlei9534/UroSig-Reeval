# Download PMC Open Access packages from the PMC open-data S3 bucket, in parallel.
# Pass 1: XML, figures, supplementary files. Pass 2: article PDFs. TIFF images are listed, not downloaded.
param([string]$Root = 'D:\article\UroSig-Reeval', [string]$Pass = 'all')
$cand = Import-Csv "$Root\03_extraction\candidates_ai_likely_include_292.csv"
$out = "$Root\03_extraction\fulltext"; New-Item -ItemType Directory -Force $out | Out-Null
$B = 'https://pmc-oa-opendata.s3.amazonaws.com'
$inv = "$Root\03_extraction\fulltext_inventory.csv"
if (-not (Test-Path $inv)) {
  $rows = @()
  foreach ($r in $cand) {
    if (-not $r.pmcid) { continue }
    $x = [xml](curl.exe -s --max-time 60 "$B/?list-type=2&prefix=$($r.pmcid).")
    $keys = @($x.ListBucketResult.Contents)
    if ($keys.Count -eq 0) { $rows += [pscustomobject]@{pmid=$r.pmid; pmcid=$r.pmcid; key=''; size=0; kind='not_in_oa_bucket'}; continue }
    $ver = ($keys | % { ($_.Key -split '/')[0] } | Sort-Object { [int](($_ -split '\.')[1]) } | Select-Object -Last 1)
    foreach ($k in ($keys | ? { $_.Key.StartsWith("$ver/") })) {
      $n = ($k.Key -split '/')[1]
      $kind = if ($n -match '^PMC\d+\.\d+\.pdf$') {'pdf'} elseif ($n -match '\.(tif|tiff)$') {'tiff'} else {'other'}
      $rows += [pscustomobject]@{pmid=$r.pmid; pmcid=$r.pmcid; key=$k.Key; size=[long]$k.Size; kind=$kind}
    }
  }
  $rows | Export-Csv $inv -NoTypeInformation -Encoding UTF8
}
$rows = Import-Csv $inv
$kinds = if ($Pass -eq 'other') {@('other')} elseif ($Pass -eq 'pdf') {@('pdf')} else {@('other','pdf')}
foreach ($kind in $kinds) {
  $cfg = "$Root\03_extraction\curl_$kind.cfg"; $lines = @()
  foreach ($r in ($rows | ? { $_.kind -eq $kind })) {
    $f = Join-Path $out ($r.pmcid + '\' + ($r.key -split '/')[1])
    if ((Test-Path $f) -and ((Get-Item $f).Length -eq [long]$r.size)) { continue }
    New-Item -ItemType Directory -Force (Split-Path $f) | Out-Null
    $lines += "url = `"$B/$($r.key)`""; $lines += "output = `"$($f -replace '\\','/')`""
  }
  if ($lines.Count) { $lines | Set-Content $cfg -Encoding ASCII; curl.exe -s --retry 3 --max-time 600 --parallel --parallel-max 12 -K $cfg }
  "$kind done: " + ($rows | ? { $_.kind -eq $kind }).Count + " files"
}
$rows | Group-Object kind | % { "{0}: {1} files, {2:N0} MB" -f $_.Name, $_.Count, (($_.Group | Measure-Object size -Sum).Sum/1MB) }
"articles in OA bucket: " + (($rows | ? { $_.kind -ne 'not_in_oa_bucket' }).pmcid | Sort-Object -Unique).Count
