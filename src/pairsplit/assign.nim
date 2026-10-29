## Assignment of read pairs to splits.

import std/[random, sysrand]

type
  RandomAssigner* = object
    ## Each pair independently goes to a uniformly random split.
    splits: int

  ExactAssigner* = object
    ## Streams over a known number of pairs and gives each split exactly
    ## `n` of them, discarding the rest. Bucket `splits` is the discard
    ## bucket. Choosing each bucket with probability proportional to its
    ## remaining slots makes every valid labeling equally likely.
    remainingSlots: seq[int]
    remaining: int
    toFill: int

proc initRandomAssigner*(splits: int): RandomAssigner =
  doAssert splits >= 1
  RandomAssigner(splits: splits)

proc next*(a: RandomAssigner, r: var Rand): int =
  r.rand(a.splits - 1)

proc initExactAssigner*(splits, n, total: int): ExactAssigner =
  ## Raises ValueError if `n * splits` exceeds `total`.
  doAssert splits >= 1 and n >= 1 and total >= 0
  if n > total div splits:
    raise newException(ValueError,
      "requested " & $n & " pairs x " & $splits & " splits but input has only " &
      $total & " pairs")
  result.remainingSlots = newSeq[int](splits + 1)
  for i in 0 ..< splits:
    result.remainingSlots[i] = n
  result.remainingSlots[splits] = total - n * splits
  result.remaining = total
  result.toFill = n * splits

proc discardBucket*(a: ExactAssigner): int = a.remainingSlots.len - 1

proc filled*(a: ExactAssigner): bool =
  ## True once every split has its `n` pairs; later pairs are all discarded.
  a.toFill == 0

proc done*(a: ExactAssigner): bool = a.remaining == 0

proc next*(a: var ExactAssigner, r: var Rand): int =
  ## Bucket for the next pair: a split index, or `discardBucket`.
  if a.remaining <= 0:
    raise newException(ValueError, "more pairs than were counted")
  let x = r.rand(a.remaining - 1)
  var cum = 0
  for i, slots in a.remainingSlots:
    cum += slots
    if cum > x:
      dec a.remainingSlots[i]
      dec a.remaining
      if i != a.discardBucket:
        dec a.toFill
      return i
  doAssert false, "unreachable"

proc randomSeed*(): int64 =
  ## A non-zero seed from the OS CSPRNG. `initRand(0)` would silently fall
  ## back to a fixed seed, so 0 is never returned.
  while result == 0:
    var bytes: array[8, byte]
    if not urandom(bytes):
      raise newException(OSError, "could not read random bytes for seed")
    result = cast[int64](bytes) and int64.high
