# pairsplit

Randomly split paired-end FASTQ files into disjoint, evenly sized pairs of
FASTQ files. Mates always stay together, and every read pair lands in at most
one split.

## Usage

```
pairsplit [options] R1 R2
```

| Option | Default | Meaning |
| --- | --- | --- |
| `-s`, `--splits` | `2` | Number of output file pairs. |
| `-n`, `--npairs` | `0` | Exact number of pairs per split; `0` splits every pair. |
| `-o`, `--outPrefix` | `pairsplit` | Output prefix. A trailing `/` names a directory (`out/` → `out/pairsplit.*`). |
| `--seed` | `0` | Random seed. `0` picks one; the seed used is always printed so a run can be reproduced. |
| `-c`, `--compress` | `auto` | `auto` (match R1), `gzip` or `none`. |
| `-t`, `--threads` | `1` | BGZF compression threads **per output file** (1–1024). |
| `-f`, `--force` | off | Overwrite existing output files. |

Long options are style-insensitive (`--outPrefix`, `--out-prefix` and
`--outprefix` are equivalent).

Example:

```
pairsplit -s 3 --seed 42 -o out/sample sample_R1.fastq.gz sample_R2.fastq.gz
```

writes `out/sample.1_R1.fastq.gz`, `out/sample.1_R2.fastq.gz`, … through
`out/sample.3_R2.fastq.gz`. Split numbers are zero-padded to the width of
`--splits` (`.01` … `.10`). A summary (seed, input pairs, pairs per split,
discarded pairs) is printed on stderr.

## Modes

- **All pairs (default).** One streaming pass. Each pair goes to a uniformly
  random split, so split sizes are approximately equal (binomial). Either
  input may be `-` (stdin) or a named pipe.
- **Exact subsample (`-n N`).** Inputs are read twice. The first pass
  validates and counts the pairs. The second gives each split exactly `N`
  pairs, chosen uniformly at random among all ways to do so, and discards
  the rest. `N * splits` must not exceed the number of input pairs. This mode
  needs regular files, not stdin or pipes.

Within each split, pairs keep their input order.

## Input and output

- Inputs can be plain text, gzip (including multi-member gzip, as made by
  `cat a.gz b.gz`) or BGZF, and R1 and R2 may differ.
- Mates are matched by read name: the first whitespace-delimited token, with
  any trailing `/1` or `/2` removed. Both `@r/1` and Casava-style
  `@r 1:N:0:ACGT` headers work.
- The program fails with a clear message and leaves no output files if any of
  these happen:
  - read names differ
  - one file has more records than the other
  - a record is malformed
  - a gzip stream is corrupt
  - R1 and R2 are the same file, or both hold the same mate (`/1` and `/1`,
    or `1:N:…` and `1:N:…`)
- In `auto` mode, output is BGZF-compressed (`.fastq.gz`, readable by any
  gzip tool) when R1 is compressed, and plain `.fastq` otherwise.
- CRLF input is written with LF line endings.
- Outputs are written to new, uniquely named `<output>.tmp-<random>` files and
  renamed only on success. They are removed on any error and on Ctrl-C or
  SIGTERM (exit 130), so a failed run leaves earlier outputs untouched.
- Existing outputs are never overwritten without `--force`. An output path
  that is a directory, or that is one of the inputs (including an input
  redirected to stdin), is always refused.
- Each compressed output gets its own thread pool, so a run uses up to
  `2 × splits × threads` compression threads.

## Building

Requires Nim ≥ 2.0.4 and htslib.

```
nimble install -y --depsOnly
nimble build
```

htslib is linked at build time instead of being loaded with `dlopen`
(see `htslib.nims`), so the binary needs no `DYLD_LIBRARY_PATH` or
`LD_LIBRARY_PATH`. The library is located through `HTSLIB_DIR` (the install
prefix), then `pkg-config htslib`, then `brew --prefix htslib`.

## Tests

```
nimble test
```

`tests/fixtures.nim` deterministically generates the test dataset in
`tests/data/`: 10,000 read pairs in every input encoding, plus edge-case and
malformed files. The tests check:

- pairing, disjointness and completeness
- split evenness and randomness, including exact uniformity of the `-n`
  sampler
- seed reproducibility
- compression handling
- stdin and named-pipe input
- error handling and overwrite protection

Outputs are verified both with the project's reader and independently with
`gzip` and `awk`.
