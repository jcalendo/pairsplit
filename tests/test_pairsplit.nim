## End-to-end tests: build the pairsplit binary and run it on the fixtures.
## Outputs are checked with the project's reader and independently with
## gzip/awk.

import std/[algorithm, math, os, osproc, re, sequtils, sets, strutils, tables, unittest]
import pairsplit/fastq
import fixtures

ensureFixtures()

let
  root = currentSourcePath().parentDir.parentDir
  outRoot = root / "tests" / "out"
  bin = root / "tests" / "bin" / "pairsplit"

block build:
  let (output, code) = execCmdEx("nim c -d:release --hints:off -o:" &
    quoteShell(bin) & " " & quoteShell(root / "src" / "pairsplit.nim"))
  doAssert code == 0, output

proc sh(cmd: string): (string, int) =
  execCmdEx(cmd, workingDir = root)

proc ps(args: string): (string, int) =
  ## Run pairsplit with `args`; returns combined stdout/stderr and exit code.
  sh(quoteShell(bin) & " " & args)

proc caseDir(name: string): string =
  ## A fresh, empty output directory for one test case.
  result = outRoot / name
  removeDir(result)
  createDir(result)

proc f(name: string): string = quoteShell(fixture(name))

proc readAll(path: string): seq[FastqRecord] =
  let r = openReader(path)
  defer: r.close()
  var rec: FastqRecord
  while r.next(rec):
    result.add rec

proc outputs(dir, prefix: string, splits: int, gz: bool): seq[(string, string)] =
  let ext = if gz: ".fastq.gz" else: ".fastq"
  let width = len($splits)
  for i in 1 .. splits:
    let stem = dir / prefix & "." & align($i, width, '0')
    result.add (stem & "_R1" & ext, stem & "_R2" & ext)

proc filesIn(dir: string): seq[string] =
  for kind, path in walkDir(dir):
    result.add path.extractFilename

proc decompressed(path: string): string =
  ## File content after `gzip -dcf` (independent of htslib).
  let (text, code) = sh("gzip -dcf " & quoteShell(path))
  doAssert code == 0, text
  text

proc awkNames(path: string): seq[string] =
  ## Header lines via gzip + awk, independent of the project's reader.
  let (text, code) = sh("gzip -dcf " & quoteShell(path) & " | awk 'NR%4==1'")
  doAssert code == 0, text
  text.splitLines.filterIt(it.len > 0)

proc reportedSeed(output: string): int64 =
  var m: array[1, string]
  doAssert output.find(re"seed: (-?\d+)", m) >= 0, output
  parseBiggestInt(m[0])

let
  inR1 = readAll(fixture("small_R1.fastq"))
  inR2 = readAll(fixture("small_R2.fastq"))
var inputIndex = initTable[string, int]()
for i, rec in inR1:
  inputIndex[pairId(rec.name)] = i

proc checkSplits(outs: seq[(string, string)]): seq[seq[int]] =
  ## Verify pairing and record integrity against the small dataset and
  ## return, per split, the input indices of its pairs.
  for (p1, p2) in outs:
    let a = readAll(p1)
    let b = readAll(p2)
    check a.len == b.len
    check awkNames(p1) == a.mapIt(it.name)
    check awkNames(p2) == b.mapIt(it.name)
    var idx: seq[int]
    for j in 0 ..< min(a.len, b.len):
      check sameMate(a[j].name, b[j].name)
      let i = inputIndex[pairId(a[j].name)]
      check a[j] == inR1[i]
      check b[j] == inR2[i]
      idx.add i
    result.add idx

proc checkDisjoint(indices: seq[seq[int]]) =
  var seen = initHashSet[int]()
  for idx in indices:
    for i in idx:
      check i notin seen
      seen.incl i

suite "pairing, completeness and evenness":
  for splits in [2, 3, 7]:
    test "all pairs split into " & $splits & " random, disjoint, even splits":
      let d = caseDir("even" & $splits)
      let (output, code) = ps("-s " & $splits & " --seed 11 -o " & quoteShell(d / "x") &
                              " " & f("small_R1.fastq.gz") & " " & f("small_R2.fastq.gz"))
      check code == 0
      let outs = outputs(d, "x", splits, gz = true)
      check filesIn(d).len == 2 * splits
      let indices = checkSplits(outs)
      checkDisjoint(indices)
      check indices.mapIt(it.len).sum == SmallPairs
      let p = 1 / splits
      let sigma = sqrt(SmallPairs.float * p * (1 - p))
      for idx in indices:
        check abs(idx.len.float - SmallPairs.float * p) < 5 * sigma
        # Random, not positional: members span the file in input order and
        # half of them come from each half of the input.
        check idx == sorted(idx)
        check idx[0] < SmallPairs div 20
        check idx[^1] > SmallPairs - SmallPairs div 20
        let firstHalf = idx.countIt(it < SmallPairs div 2)
        check abs(firstHalf.float - idx.len / 2) < 5 * sqrt(idx.len / 4)
      # Independence: pairs `lag` apart share a split with probability
      # 1/splits. Catches periodic or re-seeded assignment.
      var label = newSeq[int](SmallPairs)
      for k, idx in indices:
        for i in idx:
          label[i] = k
      for lag in [1, 2, 7, 1000]:
        let n = SmallPairs - lag
        var same = 0
        for i in 0 ..< n:
          if label[i] == label[i + lag]:
            inc same
        check abs(same.float - n.float * p) < 5 * sqrt(n.float * p * (1 - p))
      check "input pairs: " & $SmallPairs in output
      check "discarded pairs: 0" in output
      for (p1, p2) in outs:
        check sh("gzip -t " & quoteShell(p1) & " " & quoteShell(p2))[1] == 0

  test "-s 1 on plain input reproduces the input byte-for-byte":
    let d = caseDir("single")
    check ps("-s 1 -o " & quoteShell(d / "x") & " " & f("small_R1.fastq") & " " &
             f("small_R2.fastq"))[1] == 0
    let (o1, o2) = outputs(d, "x", 1, gz = false)[0]
    check readFile(o1) == readFile(fixture("small_R1.fastq"))
    check readFile(o2) == readFile(fixture("small_R2.fastq"))

suite "--npairs (exact subsampling)":
  for (splits, n) in [(2, 1000), (3, 3000), (4, 2500), (1, 500)]:
    test "-s " & $splits & " -n " & $n & " gives exactly n pairs per split":
      let d = caseDir("exact" & $splits & "_" & $n)
      let (output, code) = ps("-s " & $splits & " -n " & $n & " --seed 5 -o " &
        quoteShell(d / "x") & " " & f("small_R1.fastq") & " " & f("small_R2.fastq"))
      check code == 0
      let indices = checkSplits(outputs(d, "x", splits, gz = false))
      checkDisjoint(indices)
      for idx in indices:
        check idx.len == n
      check "discarded pairs: " & $(SmallPairs - splits * n) in output
      if splits * n == SmallPairs:
        check indices.concat.toHashSet.len == SmallPairs

  test "n * splits greater than the input is an error and writes nothing":
    let d = caseDir("exact_toomany")
    let (output, code) = ps("-s 3 -n 3334 -o " & quoteShell(d / "x") & " " &
                            f("small_R1.fastq") & " " & f("small_R2.fastq"))
    check code != 0
    check "input has only 10000 pairs" in output
    check filesIn(d).len == 0

  test "-n on empty input is an error":
    let d = caseDir("exact_empty")
    let (output, code) = ps("-n 1 -o " & quoteShell(d / "x") & " " &
                            f("empty_R1.fastq") & " " & f("empty_R2.fastq"))
    check code != 0
    check "no read pairs" in output
    check filesIn(d).len == 0

  test "each pair is selected about equally often across seeds":
    # 40 runs keeping 2 x 500 of 10000 pairs: every first-half/second-half
    # bucket should be hit in proportion.
    var hits = newSeq[int](10)
    for seed in 1 .. 40:
      let d = caseDir("exact_uniform")
      check ps("-s 2 -n 500 --seed " & $seed & " -o " & quoteShell(d / "x") & " " &
               f("small_R1.fastq") & " " & f("small_R2.fastq"))[1] == 0
      for (p1, _) in outputs(d, "x", 2, gz = false):
        for rec in readAll(p1):
          inc hits[inputIndex[pairId(rec.name)] * 10 div SmallPairs]
    # 40 * 1000 = 40000 selections over 10 deciles, 4000 expected each.
    for h in hits:
      check abs(h - 4000).float < 5 * sqrt(4000.0 * 0.9)

suite "seeds":
  test "same seed gives byte-identical outputs; different seeds differ":
    var contents: seq[string]
    for (name, seed) in [("seed_a", 99), ("seed_b", 99), ("seed_c", 100)]:
      let d = caseDir(name)
      check ps("-s 3 --seed " & $seed & " -o " & quoteShell(d / "x") & " " &
               f("small_R1.fastq.gz") & " " & f("small_R2.fastq.gz"))[1] == 0
      contents.add outputs(d, "x", 3, gz = true).mapIt(readFile(it[0]) & readFile(it[1])).join
    check contents[0] == contents[1]
    check contents[0] != contents[2]

  test "the seed reported for a random run reproduces it":
    let d1 = caseDir("seed_random")
    let (output, code) = ps("-s 2 -o " & quoteShell(d1 / "x") & " " &
                            f("small_R1.fastq") & " " & f("small_R2.fastq"))
    check code == 0
    let seed = reportedSeed(output)
    check seed != 0
    let d2 = caseDir("seed_replay")
    check ps("-s 2 --seed " & $seed & " -o " & quoteShell(d2 / "x") & " " &
             f("small_R1.fastq") & " " & f("small_R2.fastq"))[1] == 0
    for (a, b) in zip(outputs(d1, "x", 2, false), outputs(d2, "x", 2, false)):
      check readFile(a[0]) == readFile(b[0])
      check readFile(a[1]) == readFile(b[1])

  test "two random runs use different seeds":
    let d = caseDir("seed_twice")
    let s1 = reportedSeed(ps("-s 2 -f -o " & quoteShell(d / "x") & " " & f("tiny_R1.fastq") &
                             " " & f("tiny_R2.fastq"))[0])
    let s2 = reportedSeed(ps("-s 2 -f -o " & quoteShell(d / "x") & " " & f("tiny_R1.fastq") &
                             " " & f("tiny_R2.fastq"))[0])
    check s1 != s2

suite "compression":
  test "auto output matches input compression":
    for (name, gz) in [("small", true), ("small_gzip", true), ("small_multi", true)]:
      let d = caseDir("auto_" & name)
      check ps("-s 2 --seed 3 -o " & quoteShell(d / "x") & " " & f(name & "_R1.fastq.gz") &
               " " & f(name & "_R2.fastq.gz"))[1] == 0
      for (p1, p2) in outputs(d, "x", 2, gz):
        check readFile(p1)[0 .. 1] == "\x1f\x8b"
        check readFile(p2)[0 .. 1] == "\x1f\x8b"
    let d = caseDir("auto_plain")
    check ps("-s 2 --seed 3 -o " & quoteShell(d / "x") & " " & f("small_R1.fastq") &
             " " & f("small_R2.fastq"))[1] == 0
    for (p1, p2) in outputs(d, "x", 2, gz = false):
      check readFile(p1)[0] == '@'
      check readFile(p2)[0] == '@'

  test "every input encoding and output mode gives the same content":
    var reference = ""
    for (r1, r2, opts) in [
        ("small_R1.fastq", "small_R2.fastq", ""),
        ("small_R1.fastq.gz", "small_R2.fastq.gz", ""),
        ("small_gzip_R1.fastq.gz", "small_gzip_R2.fastq.gz", ""),
        ("small_multi_R1.fastq.gz", "small_multi_R2.fastq.gz", ""),
        ("small_R1.fastq.gz", "small_R2.fastq", ""),          # mixed
        ("small_R1.fastq", "small_R2.fastq.gz", ""),          # mixed
        ("small_R1.fastq", "small_R2.fastq", "-c gzip"),
        ("small_R1.fastq.gz", "small_R2.fastq.gz", "-c none"),
        ("small_R1.fastq.gz", "small_R2.fastq.gz", "-c gzip -t 4")]:
      let d = caseDir("encodings")
      check ps("-s 3 --seed 21 " & opts & " -o " & quoteShell(d / "x") & " " & f(r1) &
               " " & f(r2))[1] == 0
      let gz = if "-c none" in opts: false
               elif "-c gzip" in opts: true
               else: r1.endsWith(".gz")
      let outs = outputs(d, "x", 3, gz)
      check filesIn(d).len == 6
      for (p1, p2) in outs:
        if gz:
          check sh("gzip -t " & quoteShell(p1) & " " & quoteShell(p2))[1] == 0
        else:
          check readFile(p1)[0] == '@'
      let content = outs.mapIt(decompressed(it[0]) & decompressed(it[1])).join
      if reference.len == 0:
        reference = content
      check content == reference

suite "small and unusual inputs":
  test "more splits than pairs leaves valid empty outputs":
    let d = caseDir("tiny")
    let (output, code) = ps("-s 3 -c gzip -o " & quoteShell(d / "x") & " " &
                            f("tiny_R1.fastq") & " " & f("tiny_R2.fastq"))
    check code == 0
    check "input pairs: 1" in output
    let outs = outputs(d, "x", 3, gz = true)
    check outs.mapIt(readAll(it[0]).len).sum == 1
    for (p1, p2) in outs:
      check sh("gzip -t " & quoteShell(p1) & " " & quoteShell(p2))[1] == 0

  test "empty input writes empty outputs":
    let d = caseDir("empty")
    check ps("-s 2 -o " & quoteShell(d / "x") & " " & f("empty_R1.fastq") & " " &
             f("empty_R2.fastq"))[1] == 0
    for (p1, p2) in outputs(d, "x", 2, gz = false):
      check readFile(p1) == ""
      check readFile(p2) == ""

  test "edge-case records keep their pairing and content":
    let d = caseDir("edge")
    check ps("-s 1 -o " & quoteShell(d / "x") & " " & f("edge_R1.fastq") & " " &
             f("edge_R2.fastq"))[1] == 0
    let (p1, p2) = outputs(d, "x", 1, gz = false)[0]
    check readAll(p1) == readAll(fixture("edge_R1.fastq"))
    check readAll(p2) == readAll(fixture("edge_R2.fastq"))
    check "+edge_2/1\n" in readFile(p1)
    check readFile(p1).endsWith("+III\n")

  test "CRLF input is written with LF line endings":
    let d = caseDir("crlf")
    check ps("-s 1 -o " & quoteShell(d / "x") & " " & f("crlf_R1.fastq") & " " &
             f("crlf_R2.fastq"))[1] == 0
    let (p1, _) = outputs(d, "x", 1, gz = false)[0]
    check readFile(p1) == readFile(fixture("crlf_R1.fastq")).replace("\r\n", "\n")

  test "a prefix ending in / is a directory":
    let d = caseDir("prefixdir")
    check ps("-s 2 -o " & quoteShell(d / "sub") & "/ " & f("tiny_R1.fastq") & " " &
             f("tiny_R2.fastq"))[1] == 0
    check filesIn(d / "sub").toHashSet ==
      ["pairsplit.1_R1.fastq", "pairsplit.1_R2.fastq",
       "pairsplit.2_R1.fastq", "pairsplit.2_R2.fastq"].toHashSet

  test "split numbers are zero-padded":
    let d = caseDir("padding")
    check ps("-s 10 -o " & quoteShell(d / "x") & " " & f("tiny_R1.fastq") & " " &
             f("tiny_R2.fastq"))[1] == 0
    check fileExists(d / "x.01_R1.fastq")
    check fileExists(d / "x.10_R2.fastq")

suite "stdin and pipes":
  test "R1 from stdin":
    let d = caseDir("stdin")
    check sh("gzip -dc " & f("small_R1.fastq.gz") & " | " & quoteShell(bin) &
             " -s 2 --seed 8 -o " & quoteShell(d / "x") & " - " & f("small_R2.fastq"))[1] == 0
    let indices = checkSplits(outputs(d, "x", 2, gz = false))
    check indices.mapIt(it.len).sum == SmallPairs

  test "gzipped R2 from stdin":
    let d = caseDir("stdin_gz")
    check sh("cat " & f("small_R2.fastq.gz") & " | " & quoteShell(bin) &
             " -s 2 --seed 8 -o " & quoteShell(d / "x") & " " & f("small_R1.fastq") & " -")[1] == 0
    check checkSplits(outputs(d, "x", 2, gz = false)).mapIt(it.len).sum == SmallPairs

  test "named pipe input":
    let d = caseDir("fifo")
    let fifo = d / "r1.fifo"
    check sh("mkfifo " & quoteShell(fifo))[1] == 0
    # The writer's stderr is detached so a failed run cannot keep the output
    # pipe open, and a non-blocking open afterwards releases a writer that
    # is still waiting for a reader.
    let release = "perl -e 'use Fcntl; sysopen(F, $ARGV[0], O_RDONLY|O_NONBLOCK)' " &
                  quoteShell(fifo)
    let (output, code) = sh("(cat " & f("small_R1.fastq") & " > " & quoteShell(fifo) &
      " 2>/dev/null &) ; perl -e 'alarm 60; exec @ARGV' " & quoteShell(bin) &
      " -s 2 -o " & quoteShell(d / "x") & " " & quoteShell(fifo) & " " &
      f("small_R2.fastq") & "; rc=$?; " & release & "; exit $rc")
    check code == 0
    if code != 0: echo output
    check checkSplits(outputs(d, "x", 2, gz = false)).mapIt(it.len).sum == SmallPairs

  test "--npairs rejects stdin and pipes without reading them":
    let d = caseDir("stdin_exact")
    let fifo = d / "r1.fifo"
    check sh("mkfifo " & quoteShell(fifo))[1] == 0
    # perl alarm guards against a hang if the fifo were ever opened.
    let guard = "perl -e 'alarm 20; exec @ARGV' "
    for input in ["-", quoteShell(fifo)]:
      let (output, code) = sh(guard & quoteShell(bin) & " -n 10 -o " & quoteShell(d / "x") &
                              " " & input & " " & f("small_R2.fastq") & " < /dev/null")
      check code != 0
      check "--npairs requires regular input files" in output
    check filesIn(d) == @["r1.fifo"]

  test "--npairs rejects /dev/stdin redirected from a regular file":
    let d = caseDir("devstdin_exact")
    let (output, code) = ps("-n 10 -o " & quoteShell(d / "x") & " /dev/stdin " &
                            f("small_R2.fastq") & " < " & f("small_R1.fastq"))
    check code != 0
    check "--npairs requires regular input files" in output
    check filesIn(d).len == 0

  test "compression is detected from content, not the file name":
    # gzipped R1 on stdin, and gzipped R1 in a file named .fq
    let d = caseDir("sniff")
    copyFile(fixture("small_gzip_R1.fastq.gz"), d / "in_R1.fq")
    for (name, args) in [("a", "- " & f("small_R2.fastq") & " < " & f("small_R1.fastq.gz")),
                         ("b", quoteShell(d / "in_R1.fq") & " " & f("small_R2.fastq"))]:
      check ps("-s 2 --seed 4 -o " & quoteShell(d / name) & " " & args)[1] == 0
      let outs = outputs(d, name, 2, gz = true)
      for (p1, p2) in outs:
        check readFile(p1)[0 .. 1] == "\x1f\x8b"
        check readFile(p2)[0 .. 1] == "\x1f\x8b"
      check checkSplits(outs).mapIt(it.len).sum == SmallPairs

  test "both inputs from stdin is an error":
    let (output, code) = ps("- - < /dev/null")
    check code != 0
    check "only one input can be read from stdin" in output

suite "errors leave no output behind":
  for (r1, r2, opts, msg) in [
      ("unequal_R1.fastq", "unequal_R2.fastq", "", "has more records than"),
      ("unequal_R2.fastq", "unequal_R1.fastq", "", "has more records than"),
      ("unequal_R1.fastq", "unequal_R2.fastq", "-n 10", "has more records than"),
      ("mismatch_R1.fastq", "mismatch_R2.fastq", "", "read names differ at record 500"),
      ("mismatch_R1.fastq", "mismatch_R2.fastq", "-n 10", "read names differ at record 500"),
      ("truncated_R1.fastq", "truncated_R2.fastq", "", "record 20: truncated record"),
      ("badlen_R1.fastq", "badlen_R2.fastq", "", "record 11: sequence and quality lengths differ"),
      ("corrupt_R1.fastq.gz", "corrupt_R2.fastq.gz", "", "error reading input"),
      ("corrupt_R1.fastq.gz", "corrupt_R2.fastq.gz", "-n 10", "error reading input"),
      ("truncgz_R1.fastq.gz", "truncgz_R2.fastq.gz", "", "error reading input"),
      ("missing_R1.fastq", "small_R2.fastq", "", "input not found"),
      ("small_R1.fastq", "small_R2.fastq", "-s 0", "--splits must be at least 1"),
      ("small_R1.fastq", "small_R2.fastq", "-n -1", "--npairs must be 0"),
      ("small_R1.fastq", "small_R2.fastq", "-t 0", "--threads must be between 1 and 1024"),
      ("small_R1.fastq", "small_R2.fastq", "-t 3000000000", "--threads must be between"),
      ("small_R1.fastq", "small_R2.fastq", "-s 100000000", "open-file limit"),
      ("small_R1.fastq", "small_R2.fastq", "-s 4611686018427387904", "open-file limit"),
      ("small_R1.fastq", "small_R1.fastq", "", "R1 and R2 are the same file"),
      ("small_R1.fastq", "small_R1.fastq.gz", "", "both inputs hold mate 1 reads"),
      ("small_R2.fastq", "small_R2.fastq.gz", "-n 10", "both inputs hold mate 2 reads"),
      ("small_R1.fastq", "small_R2.fastq", "--compress bogus", "bogus"),
      ("small_R1.fastq", "small_R2.fastq", "-s abc", "abc"),
      ("small_R1.fastq", "small_R2.fastq", "--bogus", "bogus")]:
    test r1 & " " & r2 & " " & opts & " -> " & msg:
      let d = caseDir("error")
      let (output, code) = ps(opts & " -o " & quoteShell(d / "x") & " " & f(r1) & " " & f(r2))
      check code != 0
      check msg in output
      check "Traceback" notin output
      check "unhandled exception" notin output
      check filesIn(d).len == 0

  test "an output directory created for a failed run is removed":
    let d = caseDir("error")
    let (output, code) = ps("-o " & quoteShell(d / "new" / "x") & " " &
                            f("mismatch_R1.fastq") & " " & f("mismatch_R2.fastq"))
    check code != 0
    check "read names differ" in output
    check filesIn(d).len == 0

  for sig in ["INT", "TERM"]:
    test "SIG" & sig & " removes temporary outputs":
      let d = caseDir("signal")
      # stdin never delivers data, so pairsplit is blocked mid-run with its
      # temporary outputs open when the signal arrives.
      let cmd = bin & " -o " & d / "x"
      let (output, _) = sh("sleep 5 2>/dev/null | (" & cmd & " - " &
        f("small_R2.fastq") & "; echo rc=$?) & sleep 1; pkill -" & sig & " -f " &
        quoteShell("^" & cmd) & "; wait")
      check "pairsplit: interrupted" in output
      check "rc=130" in output
      check filesIn(d).len == 0

  test "wrong number of input files":
    for files in [f("small_R1.fastq"), f("small_R1.fastq") & " " & f("small_R2.fastq") &
                  " " & f("small_R2.fastq")]:
      let d = caseDir("error")
      let (output, code) = ps("-o " & quoteShell(d / "x") & " " & files)
      check code != 0
      check "expected 2 input files" in output
      check filesIn(d).len == 0

suite "overwrite protection":
  test "existing outputs are kept unless --force":
    let d = caseDir("overwrite")
    let args = " -o " & quoteShell(d / "x") & " " & f("small_R1.fastq") & " " &
               f("small_R2.fastq")
    check ps("-s 2 --seed 1" & args)[1] == 0
    let before = outputs(d, "x", 2, false).mapIt(readFile(it[0]) & readFile(it[1]))
    let (output, code) = ps("-s 2 --seed 2" & args)
    check code != 0
    check "output exists" in output
    check outputs(d, "x", 2, false).mapIt(readFile(it[0]) & readFile(it[1])) == before
    check filesIn(d).len == 4
    check ps("-s 2 --seed 2 --force" & args)[1] == 0
    check outputs(d, "x", 2, false).mapIt(readFile(it[0]) & readFile(it[1])) != before
    check filesIn(d).len == 4

  test "an output that is an input is refused even with --force":
    let d = caseDir("selfoverwrite")
    copyFile(fixture("small_R1.fastq"), d / "x.1_R1.fastq")
    copyFile(fixture("small_R2.fastq"), d / "x.1_R2.fastq")
    let original = readFile(d / "x.1_R1.fastq")
    let args = "-s 1 -o " & quoteShell(d / "x") & " " & quoteShell(d / "x.1_R1.fastq") &
               " " & quoteShell(d / "x.1_R2.fastq")
    let (output, code) = ps(args & " --force")
    check code != 0
    check "output would overwrite input" in output
    check readFile(d / "x.1_R1.fastq") == original
    check "output exists" in ps(args)[0]
    check filesIn(d).len == 2

  test "a failed --force run keeps the previous outputs":
    let d = caseDir("forcefail")
    copyFile(fixture("small_R1.fastq"), d / "x.1_R1.fastq")
    copyFile(fixture("small_R2.fastq"), d / "x.2_R2.fastq")
    let (output, code) = ps("-s 2 -f -o " & quoteShell(d / "x") & " " &
                            f("mismatch_R1.fastq") & " " & f("mismatch_R2.fastq"))
    check code != 0
    check "read names differ" in output
    check readFile(d / "x.1_R1.fastq") == readFile(fixture("small_R1.fastq"))
    check readFile(d / "x.2_R2.fastq") == readFile(fixture("small_R2.fastq"))
    check filesIn(d).toHashSet == ["x.1_R1.fastq", "x.2_R2.fastq"].toHashSet

  test "an output path that is a directory is refused before any work":
    let d = caseDir("dirout")
    createDir(d / "x.2_R2.fastq")
    copyFile(fixture("tiny_R1.fastq"), d / "x.1_R1.fastq")
    let (output, code) = ps("-s 2 -f -o " & quoteShell(d / "x") & " " &
                            f("small_R1.fastq") & " " & f("small_R2.fastq"))
    check code != 0
    check "output path is a directory" in output
    check readFile(d / "x.1_R1.fastq") == readFile(fixture("tiny_R1.fastq"))
    check filesIn(d).toHashSet == ["x.1_R1.fastq", "x.2_R2.fastq"].toHashSet

  test "an input named like a temporary output is never touched":
    let d = caseDir("tmpname")
    copyFile(fixture("small_R1.fastq"), d / "x.1_R1.fastq.tmp")
    copyFile(fixture("small_R1.fastq"), d / "x.1_R1.fastq.tmp-000000000000")
    let (_, code) = ps("-s 1 -o " & quoteShell(d / "x") & " " &
                       quoteShell(d / "x.1_R1.fastq.tmp") & " " & f("small_R2.fastq"))
    check code == 0
    check readFile(d / "x.1_R1.fastq.tmp") == readFile(fixture("small_R1.fastq"))
    check readFile(d / "x.1_R1.fastq") == readFile(fixture("small_R1.fastq"))
    check filesIn(d).toHashSet == ["x.1_R1.fastq.tmp", "x.1_R1.fastq.tmp-000000000000",
                                   "x.1_R1.fastq", "x.1_R2.fastq"].toHashSet

  test "stdin redirected from an output file is refused with --force":
    let d = caseDir("stdinself")
    copyFile(fixture("small_R1.fastq"), d / "x.1_R1.fastq")
    for input in ["-", "/dev/stdin"]:
      let (output, code) = ps("-s 1 -f -o " & quoteShell(d / "x") & " " & input & " " &
                              f("small_R2.fastq") & " < " & quoteShell(d / "x.1_R1.fastq"))
      check code != 0
      check "output would overwrite input" in output
      check readFile(d / "x.1_R1.fastq") == readFile(fixture("small_R1.fastq"))
      check filesIn(d) == @["x.1_R1.fastq"]
