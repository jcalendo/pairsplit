import std/[os, strutils, unittest]
import pairsplit/fastq
import fixtures

ensureFixtures()

proc readAll(path: string): seq[FastqRecord] =
  let r = openReader(path)
  defer: r.close()
  var rec: FastqRecord
  while r.next(rec):
    result.add rec

proc errorOf(path: string): string =
  ## Message raised while reading every record of `path` ("" if none).
  try:
    discard readAll(path)
  except IOError, ValueError:
    return getCurrentExceptionMsg()

let tmpDir = getTempDir() / "pairsplit_test_fastq"
removeDir(tmpDir)
createDir(tmpDir)

suite "FastqReader":
  test "reads every encoding of the main dataset identically":
    let plain = readAll(fixture("small_R1.fastq"))
    check plain.len == SmallPairs
    check plain[0].name == "@pair_000001/1"
    check plain[^1].name == "@pair_010000/1"
    for name in ["small_R1.fastq.gz", "small_gzip_R1.fastq.gz", "small_multi_R1.fastq.gz"]:
      check readAll(fixture(name)) == plain

  test "reports input compression":
    for (name, expected) in [("small_R1.fastq", cmpNone),
                             ("small_R1.fastq.gz", cmpBgzf),
                             ("small_gzip_R1.fastq.gz", cmpGzip),
                             ("small_multi_R1.fastq.gz", cmpGzip)]:
      let r = openReader(fixture(name))
      check r.compression == expected
      r.close()

  test "edge cases parse as expected":
    let r1 = readAll(fixture("edge_R1.fastq"))
    let r2 = readAll(fixture("edge_R2.fastq"))
    check r1.len == 4
    check r2.len == 4
    check r1[0].name == "@edge_1 1:N:0:ACGT"
    check r1[1].plus == "+edge_2/1"
    check r1[2].qual == "@III"
    check r1[3].qual == "+III" # final record without trailing newline
    check r2[3].qual == "+HHH" # trailing blank lines tolerated
    for i in 0 ..< 4:
      check sameMate(r1[i].name, r2[i].name)

  test "CRLF line endings are stripped":
    let recs = readAll(fixture("crlf_R1.fastq"))
    check recs.len == 2
    for rec in recs:
      for field in [rec.name, rec.sequence, rec.plus, rec.qual]:
        check '\r' notin field

  test "empty file has no records":
    check readAll(fixture("empty_R1.fastq")).len == 0

  test "malformed input raises with file and record number":
    let trunc = errorOf(fixture("truncated_R1.fastq"))
    check "truncated_R1.fastq: record 20: truncated record" in trunc
    let badlen = errorOf(fixture("badlen_R1.fastq"))
    check "badlen_R1.fastq: record 11: sequence and quality lengths differ" in badlen
    check "error reading input" in errorOf(fixture("corrupt_R1.fastq.gz"))
    check "error reading input" in errorOf(fixture("truncgz_R1.fastq.gz"))

  test "missing file raises IOError":
    expect IOError:
      discard openReader(fixture("does_not_exist.fastq"))

  test "header must start with @ and blank lines may not precede records":
    let p = tmpDir / "bad.fastq"
    writeFile(p, "pair_1\nACGT\n+\nIIII\n")
    check "does not start with '@'" in errorOf(p)
    writeFile(p, "@a\nACGT\n-\nIIII\n")
    check "does not start with '+'" in errorOf(p)
    writeFile(p, "@a\nACGT\n+\nIIII\n\n@b\nACGT\n+\nIIII\n")
    check "blank line" in errorOf(p)

suite "pair identifiers":
  test "pairId strips mate suffix and comments":
    check pairId("@r1/1") == "r1"
    check pairId("@r1/2 extra") == "r1"
    check pairId("@r1 1:N:0:ACGT") == "r1"
    check pairId("@r1\tx") == "r1"
    check pairId("@r1/3") == "r1/3"

  test "mateNumber":
    check mateNumber("@r1/1") == '1'
    check mateNumber("@r1/2 extra") == '2'
    check mateNumber("@r1 1:N:0:ACGT") == '1'
    check mateNumber("@r1 2:Y:18:ACGT") == '2'
    check mateNumber("@r1") == '\0'
    check mateNumber("@r1 sample=3") == '\0'
    check mateNumber("@r1/3") == '\0'

  test "sameMate":
    check sameMate("@r1/1", "@r1/2")
    check sameMate("@r1 1:N:0:A", "@r1 2:N:0:A")
    check sameMate("@r1/1", "@r1 2:N:0")
    check not sameMate("@r1/1", "@r10/2")
    check not sameMate("@r1/1", "@r2/2")

suite "FastqWriter":
  test "plain round trip is byte-identical":
    let outp = tmpDir / "rt.fastq"
    let w = openWriter(outp, compress = false)
    for rec in readAll(fixture("small_R1.fastq")):
      w.write(rec)
    w.close()
    check readFile(outp) == readFile(fixture("small_R1.fastq"))

  test "compressed round trip is valid gzip with identical content":
    for threads in [1, 4]:
      let outp = tmpDir / ("rt" & $threads & ".fastq.gz")
      let w = openWriter(outp, compress = true, threads = threads)
      for rec in readAll(fixture("small_R2.fastq")):
        w.write(rec)
      w.close()
      check readFile(outp)[0 .. 1] == "\x1f\x8b"
      check execShellCmd("gzip -t " & quoteShell(outp)) == 0
      check readAll(outp) == readAll(fixture("small_R2.fastq"))

  test "empty compressed output is valid gzip":
    let outp = tmpDir / "empty.fastq.gz"
    openWriter(outp, compress = true).close()
    check execShellCmd("gzip -t " & quoteShell(outp)) == 0
    check readAll(outp).len == 0
