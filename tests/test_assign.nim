import std/[random, tables, unittest]
import pairsplit/assign

proc chiSquare(counts: openArray[int], expected: float): float =
  for c in counts:
    let d = float(c) - expected
    result += d * d / expected

suite "RandomAssigner":
  test "covers every split roughly uniformly":
    # 6 degrees of freedom: p = 0.001 critical value is 22.46.
    var r = initRand(1)
    let a = initRandomAssigner(7)
    var counts = newSeq[int](7)
    for _ in 1 .. 70_000:
      let k = a.next(r)
      check k in 0 .. 6
      inc counts[k]
    # Bounded below too (p = 0.001 lower quantile 0.381): a round-robin
    # assigner would be suspiciously perfect.
    check chiSquare(counts, 10_000.0) in 0.381 .. 22.46

  test "a single split always gets everything":
    var r = initRand(2)
    let a = initRandomAssigner(1)
    for _ in 1 .. 100:
      check a.next(r) == 0

suite "ExactAssigner":
  test "exactly n per split and the rest discarded":
    for seed in 1 .. 200:
      var r = initRand(seed)
      for (splits, n, total) in [(2, 10, 25), (3, 7, 21), (7, 1, 100), (1, 5, 5)]:
        var a = initExactAssigner(splits, n, total)
        var counts = newSeq[int](splits + 1)
        for _ in 1 .. total:
          inc counts[a.next(r)]
        check a.done
        check a.filled
        for i in 0 ..< splits:
          check counts[i] == n
        check counts[splits] == total - n * splits

  test "filled becomes true exactly when the last split slot is used":
    var r = initRand(3)
    var a = initExactAssigner(2, 3, 50)
    var assigned = 0
    for _ in 1 .. 50:
      check a.filled == (assigned == 6)
      if a.next(r) != a.discardBucket:
        inc assigned
    check a.filled

  test "every labeling is equally likely (T=5, splits=2, n=2)":
    # 5! / (2! 2! 1!) = 30 labelings; 29 dof, p = 0.001 critical value 58.30.
    const Trials = 60_000
    var r = initRand(4)
    var seen = initCountTable[string]()
    var included = newSeq[int](5)
    for _ in 1 .. Trials:
      var a = initExactAssigner(2, 2, 5)
      var key = ""
      for pos in 0 ..< 5:
        let k = a.next(r)
        key.add char(ord('0') + k)
        if k != a.discardBucket:
          inc included[pos]
      seen.inc key
    check seen.len == 30
    var counts: seq[int]
    for _, c in seen:
      counts.add c
    check chiSquare(counts, Trials / 30) in 11.81 .. 58.30 # p = 0.001 both tails
    # Each record is kept with probability n*splits/T = 4/5.
    for c in included:
      check abs(c / Trials - 0.8) < 0.01

  test "rejects n * splits > total":
    expect ValueError:
      discard initExactAssigner(3, 4, 11)
    expect ValueError:
      discard initExactAssigner(1, 1, 0)

  test "raises when called more times than counted":
    var r = initRand(5)
    var a = initExactAssigner(2, 1, 3)
    for _ in 1 .. 3:
      discard a.next(r)
    expect ValueError:
      discard a.next(r)

suite "randomSeed":
  test "never zero and varies":
    let a = randomSeed()
    let b = randomSeed()
    check a != 0
    check b != 0
    check a != b
