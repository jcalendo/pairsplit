# Link htslib at build time rather than dlopen-ing it at runtime, so the
# binary records an absolute install name and needs no DYLD_LIBRARY_PATH.
# Override the location with HTSLIB_DIR (the htslib install prefix).
import std/strutils

proc htslibLinkFlags(): string =
  let dir = getEnv("HTSLIB_DIR")
  if dir.len > 0:
    return "-L" & dir & "/lib -lhts"
  let (flags, code) = gorgeEx("pkg-config --libs htslib")
  if code == 0 and flags.strip.len > 0:
    return flags.strip
  let (prefix, bcode) = gorgeEx("brew --prefix htslib")
  if bcode == 0 and prefix.strip.len > 0:
    return "-L" & prefix.strip & "/lib -lhts"
  "-lhts"

switch("dynlibOverride", "hts")
switch("passL", htslibLinkFlags())
