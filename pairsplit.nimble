# Package

version       = "0.1.0"
author        = "jcalendo"
description   = "Generate even splits from subsampled read pairs in paired-end FASTQ files"
license       = "MIT"
srcDir        = "src"
bin           = @["pairsplit"]


# Dependencies

requires "nim >= 2.0.4"
requires "hts >= 0.3.25"
requires "cligen >= 1.7.0"
