## Deterministic test dataset for pairsplit, generated into tests/data/.
##
## Generation is idempotent: a stamp file records the fixture version and
## the data is only rebuilt when it is missing or out of date.

import std/[os, random, strformat, strutils]
import hts/bgzf

const
  FixtureVersion = "4"
  SmallPairs* = 10_000
  UnequalPairs = 1000
  MismatchAt = 500

let dataDir* = currentSourcePath().parentDir / "data"

proc fixture*(name: string): string = dataDir / name

proc randSeq(r: var Rand, n: int): string =
  result = newString(n)
  for i in 0 ..< n:
    result[i] = "ACGT"[r.rand(3)]

proc randQual(r: var Rand, n: int): string =
  result = newString(n)
  for i in 0 ..< n:
    result[i] = char(r.rand(int('!') .. int('J')))

proc fastqText(names, seqs, quals: seq[string]): string =
  for i in 0 ..< names.len:
    result.add &"{names[i]}\n{seqs[i]}\n+\n{quals[i]}\n"

proc pairTexts(prefix: string, n: int, seed: int64): (string, string) =
  ## Generate `n` read pairs named "@<prefix>_NNNNNN/1" and "/2".
  var r = initRand(seed)
  var n1, s1, q1, n2, s2, q2: seq[string]
  for i in 1 .. n:
    let id = &"{prefix}_{i:06}"
    let l1 = r.rand(50 .. 150)
    let l2 = r.rand(50 .. 150)
    n1.add &"@{id}/1"
    s1.add r.randSeq(l1)
    q1.add r.randQual(l1)
    n2.add &"@{id}/2"
    s2.add r.randSeq(l2)
    q2.add r.randQual(l2)
  (fastqText(n1, s1, q1), fastqText(n2, s2, q2))

proc writeCompressed(path, text, mode: string) =
  ## mode "w" writes BGZF, "wg" writes regular gzip.
  var b: BGZ
  b.open(path, mode)
  if text.len > 0:
    doAssert b.write(text) == text.len
  doAssert b.close() == 0

proc writeMultiMember(path, text: string) =
  ## Regular gzip made of two concatenated members (like `cat a.gz b.gz`).
  let lines = text.splitLines(keepEol = true)
  let half = (lines.len div 8) * 4 # split on a record boundary
  let tmpA = path & ".a"
  let tmpB = path & ".b"
  writeCompressed(tmpA, lines[0 ..< half].join, "wg")
  writeCompressed(tmpB, lines[half .. ^1].join, "wg")
  writeFile(path, readFile(tmpA) & readFile(tmpB))
  removeFile(tmpA)
  removeFile(tmpB)

proc generate() =
  createDir(dataDir)

  # Main dataset in every supported input encoding.
  let (r1, r2) = pairTexts("pair", SmallPairs, 20261029)
  writeFile(fixture("small_R1.fastq"), r1)
  writeFile(fixture("small_R2.fastq"), r2)
  writeCompressed(fixture("small_R1.fastq.gz"), r1, "w")
  writeCompressed(fixture("small_R2.fastq.gz"), r2, "w")
  writeCompressed(fixture("small_gzip_R1.fastq.gz"), r1, "wg")
  writeCompressed(fixture("small_gzip_R2.fastq.gz"), r2, "wg")
  writeMultiMember(fixture("small_multi_R1.fastq.gz"), r1)
  writeMultiMember(fixture("small_multi_R2.fastq.gz"), r2)

  let (t1, t2) = pairTexts("tiny", 1, 7)
  writeFile(fixture("tiny_R1.fastq"), t1)
  writeFile(fixture("tiny_R2.fastq"), t2)
  writeFile(fixture("empty_R1.fastq"), "")
  writeFile(fixture("empty_R2.fastq"), "")

  # Format edge cases. Expected parse results live in tests as EdgeRecords.
  writeFile(fixture("edge_R1.fastq"),
    "@edge_1 1:N:0:ACGT\nACGT\n+\nIIII\n" &
    "@edge_2/1\nAACC\n+edge_2/1\nJJJJ\n" &
    "@edge_3/1\nGGTT\n+\n@III\n" &
    "@edge_4/1\nTTAA\n+\n+III") # no trailing newline
  writeFile(fixture("edge_R2.fastq"),
    "@edge_1 2:N:0:ACGT\nTGCA\n+\nHHHH\n" &
    "@edge_2/2\nGGCC\n+edge_2/2\nAAAA\n" &
    "@edge_3/2\nCCAA\n+\n@HHH\n" &
    "@edge_4/2\nAATT\n+\n+HHH\n\n\n") # trailing blank lines
  writeFile(fixture("crlf_R1.fastq"),
    "@crlf_1/1\r\nACGT\r\n+\r\nIIII\r\n@crlf_2/1\r\nGGCC\r\n+\r\nJJJJ\r\n")
  writeFile(fixture("crlf_R2.fastq"),
    "@crlf_1/2\r\nTTGG\r\n+\r\nHHHH\r\n@crlf_2/2\r\nCCAA\r\n+\r\nAAAA\r\n")

  # Error fixtures.
  let (u1, u2) = pairTexts("uneq", UnequalPairs, 11)
  writeFile(fixture("unequal_R1.fastq"), u1)
  writeFile(fixture("unequal_R2.fastq"),
    u2.splitLines(keepEol = true)[0 ..< (UnequalPairs - 1) * 4].join)

  let (m1, m2) = pairTexts("mism", UnequalPairs, 12)
  writeFile(fixture("mismatch_R1.fastq"), m1)
  writeFile(fixture("mismatch_R2.fastq"),
    m2.replace(&"@mism_{MismatchAt:06}/2", "@other_read/2"))

  let (tr1, tr2) = pairTexts("trunc", 20, 13)
  writeFile(fixture("truncated_R1.fastq"),
    tr1.splitLines(keepEol = true)[0 ..< 20 * 4 - 1].join) # drop final qual
  writeFile(fixture("truncated_R2.fastq"), tr2)

  let (b1, b2) = pairTexts("blen", 20, 14)
  var bl = b1.splitLines(keepEol = true)
  bl[10 * 4 + 3] = "I\n" # record 11 quality shorter than sequence
  writeFile(fixture("badlen_R1.fastq"), bl.join)
  writeFile(fixture("badlen_R2.fastq"), b2)

  # Damaged gzip R1 files, each paired with an intact R2: a bad CRC32 in the
  # gzip trailer (all data decodes, the check fails at the end) and a file
  # cut off halfway through.
  writeCompressed(fixture("corrupt_R2.fastq.gz"), r2, "wg")
  copyFile(fixture("corrupt_R2.fastq.gz"), fixture("truncgz_R2.fastq.gz"))
  writeCompressed(fixture("corrupt_R1.fastq.gz"), r1, "wg")
  var c = readFile(fixture("corrupt_R1.fastq.gz"))
  writeFile(fixture("truncgz_R1.fastq.gz"), c[0 ..< c.len div 2])
  c[^8] = char(ord(c[^8]) xor 0xff)
  writeFile(fixture("corrupt_R1.fastq.gz"), c)

  writeFile(fixture("VERSION"), FixtureVersion)

proc ensureFixtures*() =
  let stamp = fixture("VERSION")
  if not (fileExists(stamp) and readFile(stamp) == FixtureVersion):
    generate()

when isMainModule:
  ensureFixtures()
  echo "fixtures in ", dataDir
