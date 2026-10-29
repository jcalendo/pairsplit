## pairsplit: ingest paired-end FASTQ files and generate even splits
## from subsampled read pairs.

import std/[os, posix, random, strutils, sysrand]
import pairsplit/[assign, fastq]

const
  Version = "0.1.0"
  MaxThreads = 1024
  ReservedFiles = 16   ## descriptors kept free for inputs, stdio, htslib

type
  OutCompression* = enum
    ocAuto = "auto", ocGzip = "gzip", ocNone = "none"

  SplitOutput = object
    r1, r2: string        ## final output paths
    t1, t2: string        ## temporary paths written before the final rename
    w1, w2: FastqWriter
    pairs: int

  FileId = tuple[dev: Dev, ino: Ino]

var tmpRegistry: seq[string]
  ## Temporary outputs that exist right now; removed if we are interrupted.

proc log(msg: string) =
  stderr.writeLine "pairsplit: " & msg

proc isStdinPath(path: string): bool =
  path == "-" or path == "/dev/stdin" or path.startsWith("/dev/fd/")

proc inputStat(path: string, st: var Stat): bool =
  ## stat() of the file an input path actually reads. For stdin and
  ## /dev/fd/N this uses fstat() on the descriptor, because stat() on those
  ## paths reports the device file on macOS, not what is behind it.
  if path == "-" or path == "/dev/stdin":
    return fstat(0, st) == 0
  if path.startsWith("/dev/fd/"):
    try:
      return fstat(cint(parseInt(path[8 .. ^1])), st) == 0
    except ValueError:
      discard
  stat(path.cstring, st) == 0

proc fileId(st: Stat): FileId = (st.st_dev, st.st_ino)

proc openFileLimit(): int =
  var lim: RLimit
  if getrlimit(RLIMIT_NOFILE, lim) != 0:
    return high(int)
  lim.rlim_cur

proc outputPaths(prefix: string, splits: int, gz: bool): seq[(string, string)] =
  ## "<prefix>.<i>_R1.fastq[.gz]" with i zero-padded to the width of
  ## `splits`. A prefix ending in "/" names a directory.
  let base = if prefix.endsWith('/'): prefix / "pairsplit" else: prefix
  let ext = if gz: ".fastq.gz" else: ".fastq"
  let width = len($splits)
  for i in 1 .. splits:
    let stem = base & "." & align($i, width, '0')
    result.add (stem & "_R1" & ext, stem & "_R2" & ext)

proc createTmp(final: string): string =
  ## Create a new, empty, uniquely named temporary file next to `final`.
  ## O_EXCL guarantees an existing file (e.g. an input) is never reused.
  for attempt in 1 .. 100:
    var suffix: array[6, byte]
    if not urandom(suffix):
      raise newException(OSError, "could not read random bytes")
    var hex = ""
    for b in suffix:
      hex.add toHex(b)
    result = final & ".tmp-" & hex.toLowerAscii
    let fd = posix.open(result.cstring, O_WRONLY or O_CREAT or O_EXCL, 0o644)
    if fd >= 0:
      discard posix.close(fd)
      tmpRegistry.add result
      return result
    if errno != EEXIST:
      raiseOSError(osLastError(), result)
  raise newException(IOError, "could not create a temporary file for " & final)

proc removeTmp(path: string) =
  if path.len > 0:
    discard tryRemoveFile(path)
    let i = tmpRegistry.find(path)
    if i >= 0:
      tmpRegistry.del(i)

proc onInterrupt() {.noconv.} =
  for p in tmpRegistry:
    discard tryRemoveFile(p)
  stderr.writeLine "pairsplit: interrupted"
  exitnow(130) # _exit: no exit handlers inside a signal handler

proc onTerminate(sig: cint) {.noconv.} =
  onInterrupt()

proc readPair(a, b: FastqReader, x, y: var FastqRecord): bool =
  ## Read the next mate from each file, checking they stay in lockstep.
  let hasA = a.next(x)
  let hasB = b.next(y)
  if hasA != hasB:
    let (longer, shorter) = if hasA: (a, b) else: (b, a)
    raise newException(ValueError,
      longer.path & " has more records than " & shorter.path & " (" &
      shorter.path & " ended after " & $shorter.nrec & " records)")
  if hasA:
    if not sameMate(x.name, y.name):
      raise newException(ValueError,
        "read names differ at record " & $a.nrec & ": '" & x.name & "' in " &
        a.path & " vs '" & y.name & "' in " & b.path)
    let mate = mateNumber(x.name)
    if mate != '\0' and mate == mateNumber(y.name):
      raise newException(ValueError,
        "both inputs hold mate " & mate & " reads at record " & $a.nrec &
        " ('" & x.name & "'); R1 and R2 must be the two mate files")
  hasA

proc countPairs(path1, path2: string): int =
  ## First pass of -n mode: validate both files and count the pairs.
  let a = openReader(path1)
  defer: a.close()
  let b = openReader(path2)
  defer: b.close()
  var x, y: FastqRecord
  while readPair(a, b, x, y):
    discard
  a.nrec

proc checkOutputs(paths: seq[(string, string)], inputs: openArray[FileId],
                  force: bool) =
  ## Refuse outputs that exist (unless --force), are directories, or are
  ## one of the inputs. Runs before any input is processed.
  for (p1, p2) in paths:
    for p in [p1, p2]:
      var st: Stat
      if lstat(p.cstring, st) != 0:
        continue
      if S_ISDIR(st.st_mode) or dirExists(p):
        raise newException(IOError, "output path is a directory: " & p)
      if not force:
        raise newException(IOError, "output exists (use --force to overwrite): " & p)
      if stat(p.cstring, st) == 0 and fileId(st) in inputs:
        raise newException(IOError, "output would overwrite input: " & p)

proc discardOutputs(outs: var seq[SplitOutput]) =
  ## Close and delete any partially written temporary outputs.
  for o in outs.mitems:
    for w in [o.w1, o.w2]:
      if w != nil:
        try: w.close()
        except IOError: discard
    removeTmp(o.t1)
    removeTmp(o.t2)

proc run(r1Path, r2Path: string, inputs: openArray[FileId], splits, npairs: int,
         prefix: string, seed: int64, compress: OutCompression, threads: int,
         force: bool) =
  let exact = npairs > 0
  let total = if exact: countPairs(r1Path, r2Path) else: 0
  var exactAssigner: ExactAssigner
  if exact:
    if total == 0:
      raise newException(ValueError, "input has no read pairs")
    exactAssigner = initExactAssigner(splits, npairs, total)

  let a = openReader(r1Path)
  defer: a.close()
  let b = openReader(r2Path)
  defer: b.close()
  let gz = case compress
    of ocAuto: a.compression != cmpNone
    of ocGzip: true
    of ocNone: false

  let paths = outputPaths(prefix, splits, gz)
  checkOutputs(paths, inputs, force)
  var outs: seq[SplitOutput]
  for (p1, p2) in paths:
    outs.add SplitOutput(r1: p1, r2: p2)
  let outDir = outs[0].r1.parentDir
  let madeDir = outDir.len > 0 and not dirExists(outDir)
  if madeDir:
    createDir(outDir)

  log "seed: " & $seed
  var rng = initRand(seed)
  let randomAssigner = initRandomAssigner(splits)
  var committed: seq[string]
  try:
    for o in outs.mitems:
      o.t1 = createTmp(o.r1)
      o.t2 = createTmp(o.r2)
      o.w1 = openWriter(o.t1, gz, threads)
      o.w2 = openWriter(o.t2, gz, threads)
    var x, y: FastqRecord
    while readPair(a, b, x, y):
      let k = if exact: exactAssigner.next(rng) else: randomAssigner.next(rng)
      if k < splits:
        outs[k].w1.write(x)
        outs[k].w2.write(y)
        inc outs[k].pairs
      if exact and exactAssigner.filled:
        break
    if exact and not exactAssigner.filled:
      raise newException(ValueError, "the second pass over the inputs found " &
        "fewer pairs than the first (" & $a.nrec & " of " & $total &
        "); inputs must not change while pairsplit runs")
    for o in outs.mitems:
      o.w1.close()
      o.w2.close()
    for o in outs.mitems:
      for (tmp, final) in [(o.t1, o.r1), (o.t2, o.r2)]:
        moveFile(tmp, final)
        committed.add final
      removeTmp(o.t1)
      removeTmp(o.t2)
      o.t1 = ""
      o.t2 = ""
  except CatchableError as e:
    discardOutputs(outs)
    if madeDir and committed.len == 0:
      discard rmdir(outDir.cstring) # only succeeds if still empty
    if committed.len > 0:
      e.msg &= " (outputs already written: " & committed.join(", ") & ")"
    raise

  let inputPairs = if exact: total else: a.nrec
  var written = 0
  log "input pairs: " & $inputPairs
  for i, o in outs:
    log "split " & $(i + 1) & ": " & $o.pairs & " pairs -> " & o.r1 & ", " & o.r2
    written += o.pairs
  log "discarded pairs: " & $(inputPairs - written)

proc pairsplit(files: seq[string], splits = 2, npairs = 0,
               outPrefix = "pairsplit", seed = 0, compress = ocAuto,
               threads = 1, force = false): int =
  ## Randomly split paired-end FASTQ files (R1 R2, gzipped or plain) into
  ## `--splits` disjoint pairs of FASTQ files. Mates always stay together.
  ## By default every pair is assigned to a uniformly random split in one
  ## streaming pass (sizes approximately even; "-" reads stdin). With
  ## `--npairs` N, inputs are counted first and each split receives exactly
  ## N randomly chosen pairs (regular files only).
  try:
    if files.len != 2:
      raise newException(ValueError,
        "expected 2 input files (R1 R2), got " & $files.len)
    if files[0] == "-" and files[1] == "-":
      raise newException(ValueError, "only one input can be read from stdin")
    var inputs: seq[FileId]
    for f in files:
      var st: Stat
      if not inputStat(f, st) or S_ISDIR(st.st_mode):
        raise newException(IOError, "input not found: " & f)
      if npairs > 0 and (isStdinPath(f) or not S_ISREG(st.st_mode)):
        raise newException(ValueError,
          "--npairs requires regular input files (it reads them twice): " & f)
      inputs.add fileId(st)
    if inputs[0] == inputs[1]:
      raise newException(ValueError, "R1 and R2 are the same file")
    if splits < 1:
      raise newException(ValueError, "--splits must be at least 1")
    if npairs < 0:
      raise newException(ValueError, "--npairs must be 0 (all pairs) or positive")
    if threads < 1 or threads > MaxThreads:
      raise newException(ValueError, "--threads must be between 1 and " & $MaxThreads)
    let limit = openFileLimit()
    if splits > (limit - ReservedFiles) div 2:
      raise newException(ValueError, $splits & " splits need two open files " &
        "each, above the open-file limit of " & $limit & " (raise it with `ulimit -n`)")
    let seed = if seed == 0: randomSeed() else: int64(seed)
    setControlCHook(onInterrupt)
    posix.signal(SIGTERM, onTerminate)
    run(files[0], files[1], inputs, splits, npairs, outPrefix, seed, compress,
        threads, force)
  except IOError, ValueError, OSError:
    log "error: " & getCurrentExceptionMsg()
    return 1

when isMainModule:
  import cligen
  clCfg.version = Version
  dispatch pairsplit,
    short = {"splits": 's', "npairs": 'n', "outPrefix": 'o', "seed": '\0',
             "compress": 'c', "threads": 't', "force": 'f'},
    help = {
      "files": "R1 R2",
      "splits": "number of output file pairs (>= 1)",
      "npairs": "exact pairs per split; 0 = split all pairs",
      "outPrefix": "output prefix; writes <prefix>.<i>_R{1,2}.fastq[.gz] " &
                   "(a trailing / means a directory)",
      "seed": "random seed; 0 = pick one (reported on stderr)",
      "compress": "output compression: auto (match R1), gzip or none",
      "threads": "BGZF compression threads per output file (1-1024)",
      "force": "overwrite existing output files"}
