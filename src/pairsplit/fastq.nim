## Streaming FASTQ reader and writer built on htslib's BGZF layer.
##
## `bgzf_open` transparently reads plain text, regular (multi-member) gzip
## and BGZF, so one reader handles every supported input format.

import hts/bgzf
import hts/utils
import hts/private/hts_concat

proc bgzf_compression(fp: ptr BGZF): cint {.cdecl, importc: "bgzf_compression".}

type
  Compression* = enum
    ## Compression of an opened input, as reported by htslib.
    cmpNone, cmpGzip, cmpBgzf

  FastqRecord* = object
    name*: string      ## header line, including the leading '@'
    sequence*: string
    plus*: string      ## separator line, kept verbatim (may be "+name")
    qual*: string

  FastqReader* = ref object
    path*: string
    nrec*: int         ## records successfully read so far
    bgz: BGZ
    kstr: kstring_t
    sawBlank: bool

  FastqWriter* = ref object
    path*: string
    bgz: BGZ
    buf: string

proc openReader*(path: string): FastqReader =
  ## Open `path` for reading; "-" reads stdin.
  result = FastqReader(path: path, kstr: kstring_t(l: 0, m: 0, s: nil))
  result.bgz.open(path, "r")

proc compression*(r: FastqReader): Compression =
  case bgzf_compression(r.bgz.cptr)
  of 1: cmpGzip
  of 2: cmpBgzf
  else: cmpNone

proc close*(r: FastqReader) =
  if r.kstr.s != nil:
    free(r.kstr.s)
    r.kstr = kstring_t(l: 0, m: 0, s: nil)
  r.bgz.close()

proc fail(r: FastqReader, msg: string) {.noreturn.} =
  raise newException(ValueError, r.path & ": record " & $(r.nrec + 1) & ": " & msg)

proc readLine(r: FastqReader, line: var string): bool =
  ## Read one line (without '\n' or a trailing '\r'); false at EOF.
  let ret = r.bgz.read_line(r.kstr.addr)
  if ret <= -2:
    raise newException(IOError, r.path & ": error reading input (corrupt or truncated file?)")
  if ret < 0:
    return false
  line.setLen(r.kstr.l)
  if r.kstr.l > 0:
    copyMem(line[0].addr, r.kstr.s, r.kstr.l)
  true

proc next*(r: FastqReader, rec: var FastqRecord): bool =
  ## Read the next record into `rec`. Returns false at end of input.
  ## Blank lines are only tolerated at the end of the file.
  while true:
    if not r.readLine(rec.name):
      return false
    if rec.name.len == 0:
      r.sawBlank = true
      continue
    break
  if r.sawBlank:
    r.fail("unexpected blank line before record")
  if rec.name[0] != '@':
    r.fail("header line does not start with '@'")
  if not r.readLine(rec.sequence):
    r.fail("truncated record (missing sequence line)")
  if not r.readLine(rec.plus):
    r.fail("truncated record (missing '+' line)")
  if rec.plus.len == 0 or rec.plus[0] != '+':
    r.fail("separator line does not start with '+'")
  if not r.readLine(rec.qual):
    r.fail("truncated record (missing quality line)")
  if rec.qual.len != rec.sequence.len:
    r.fail("sequence and quality lengths differ (" & $rec.sequence.len &
           " vs " & $rec.qual.len & ")")
  inc r.nrec
  true

proc idEnd(name: string): int =
  ## End (exclusive) of the pair identifier in a header line: the first
  ## whitespace-delimited token, minus a trailing "/1" or "/2".
  result = 1
  while result < name.len and name[result] notin {' ', '\t'}:
    inc result
  if result >= 3 and name[result - 2] == '/' and name[result - 1] in {'1', '2'}:
    result -= 2

proc mateNumber*(name: string): char =
  ## '1' or '2' from a "/1" "/2" suffix or a Casava "1:N:..." comment;
  ## '\0' when the header does not say.
  var e = 1
  while e < name.len and name[e] notin {' ', '\t'}:
    inc e
  if e >= 3 and name[e - 2] == '/' and name[e - 1] in {'1', '2'}:
    return name[e - 1]
  if e + 2 < name.len and name[e + 1] in {'1', '2'} and name[e + 2] == ':':
    return name[e + 1]
  '\0'

proc pairId*(name: string): string =
  ## Identifier shared by both mates, e.g. "@r1/1 x" and "@r1 1:N:0" -> "r1".
  name[1 ..< idEnd(name)]

proc sameMate*(a, b: string): bool =
  ## True when header lines `a` and `b` name the same read pair.
  let ea = idEnd(a)
  if ea != idEnd(b):
    return false
  for i in 1 ..< ea:
    if a[i] != b[i]:
      return false
  true

proc openWriter*(path: string, compress: bool, threads = 1): FastqWriter =
  ## Open `path` for writing, BGZF-compressed or plain.
  result = FastqWriter(path: path)
  result.bgz.open(path, if compress: "w" else: "wu")
  if compress and threads > 1:
    if bgzf_mt(result.bgz.cptr, cint(threads), 128) != 0:
      raise newException(IOError, path & ": could not start compression threads")

proc write*(w: FastqWriter, rec: FastqRecord) =
  w.buf.setLen(0)
  w.buf.add rec.name
  w.buf.add '\n'
  w.buf.add rec.sequence
  w.buf.add '\n'
  w.buf.add rec.plus
  w.buf.add '\n'
  w.buf.add rec.qual
  w.buf.add '\n'
  if w.bgz.write(w.buf) < 0:
    raise newException(IOError, w.path & ": write failed")

proc close*(w: FastqWriter) =
  if w.bgz.close() != 0:
    raise newException(IOError, w.path & ": error closing file")
