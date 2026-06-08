# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [2.0.8]

### Added
- FamDB is now a dependency of RepeatModeler and has been updated to support Dfam 4.0
- Support for makeblastdb partitioning of large genomes (e.g. `*.00.nsq`).
- Long format support in 2bit files for large genomes.
- Catch for older versions of RepeatScout.
- Support for RECON 1.10.

### Changed
- Configure now logs without screen clearing making it easier to debug problems.
- RepeatClassifier is now responsible for querying and caching reference 
  libraries from FamDB.

### Fixed
- RepeatModeler/Refiner was not processing consensus extensions properly in 2.0.7
  and in some cases was concatenating forward and reverse strand sequences together,
  adding noise to seed alignments and consensus calls.

## [2.0.7]

### Added
- RepeatAfterMe auto-extension step to Refiner, applied to both RepeatScout and
  RECON generated family models.
- RepeatClassifier now checks for a complete Dfam FamDB database.
- New redundancy filter applied to RepeatScout families.

### Changed
- Lowered masking score threshold to reduce rediscovery.

### Fixed
- Bug with the recoverDir round detection logic.

## [2.0.6]

### Added
- MultAln: new API for calculating Kimura divergence.
- MultAln: new API for slicing operations on MSA.
- MultAln: switch to Smitten V2 for default STK output.
- MultAln: improvements to the Linup reporting format.
- Refiner reorganization for future extension algorithm.

### Changed
- Updated to RepeatScout 1.0.7, which reports extended seed coordinates, reducing
  complexity when searching for exemplars using the consensus sequence.
- Reduced the minimum family instances for RepeatScout families.

### Fixed
- Bug that caused RepeatClassifier to run with only one thread in the last release.
- Bug that caused some RepeatScout-derived exemplars to be reported on the wrong strand.
- Long sequence identifiers from LTR_retriever that altered names longer than 13
  characters are now isolated correctly.
- Fixes for LTR_retriever changes in TRF parameter.

## [2.0.5]

### Added
- `-long` option to `faToTwoBit` to support larger genomes.

### Fixed
- Bug that caused failures in "absolute" reproducibility. Prior to this release,
  use of `-srand` would not guarantee that `consensi.fa.classified` and
  `families-classified.stk` were identical in sequence and sequence order across
  runs. Secondary sorts were added to guarantee a fixed sort order among equally
  scoring results. **Note:** This change applies only to results generated with
  this version and future releases.

## [2.0.4]

### Added
- New `threads` option replacing the `pa` (parallel batches) option, mapping
  directly to the maximum number of threads the program will attempt to use.
- Support for RepeatMasker 4.1.4 and RMBlast 2.13.0 parallel query feature.
- Larger default sample sizes enabled by speed improvements; original sampling
  strategy selectable with the new `quick` option.

### Changed
- Improved gathering of RepeatScout exemplars for building seed alignments.
- Parallelized and improved masking between rounds for faster runs and fewer
  rediscoveries.
- ABBlast is no longer supported.

### Fixed
- Several visual artifacts in the `viewMSA.pl` HTML visualization.

## [2.0.3]

### Added
- Program now generates a logfile in the working directory named
  `<seqfile>-rmod.log`, containing the random seed number and high-level run
  stats for use when reporting problems.

### Fixed
- Problem with orientation/coordinates in Stockholm output format that affected
  a subset of sequences.
- Bug affecting the trim functions of the Linup tool.

## [2.0.2]

### Added
- First release of a set of manual curation tools for use with de-novo generated
  TE libraries.
- `generateSeedAlignments.pl` to generate Dfam-compatible seed alignments given a
  consensus-based TE library and RepeatMasker output.

### Fixed
- Bug in N50 calculation.
- Ruzzo-Tompa maximal scoring subsequences implementation.
- Various small bugfixes.

## [2.0.1]

### Added
- Version of each dependency is now printed in the log.

### Fixed
- Workaround for a bug introduced in NCBI Blast 2.10.0 with version 5 databases.
- RepeatScout in rare cases generates models for very long-period satellites, causing
  Refiner to produce large numbers of off-diagonal alignments; these cases are now
  filtered out.

## [2.0]

### Added
- Container releases (Docker and Singularity).
- Code to remove excess temporary files from the `RM_` directory after they are no
  longer needed.

### Changed
- Seed alignments and consensus sequences previously contained stretches of N's due
  to TRF masking; original genomic sequences are now used instead, preserving
  low-complexity regions in the final consensus.
- Batching system refactored to improve handling of short input sequences common in
  early assemblies or low-coverage genome projects.
- TRF interface refactored to eliminate a major performance bottleneck in the
  simple/tandem repeat masking step.

### Fixed
- Stockholm sequence coordinates did not consistently show strandedness via
  coordinate reversal.
- Temporary file generation in `BuildDatabase` did not use a unique filename,
  risking race conditions when two `BuildDatabase` processes ran in the same
  directory simultaneously.
- Bug that caused some TE families to be labeled with the RepeatMasker
  internal-use label "buffer".
- TRF masking of input sequences caused some low-complexity sequences to
  incorrectly match TE library sequences.
- Conflict resolution often let a lower-scoring match define the class for a family.

## [1.0.11]

### Fixed
- Avoid bailing if Refiner could not produce a consensus from the family sequence set.

## [1.0.10]

### Fixed
- Bug that caused RepeatModeler to exit if no new repeats were detected after
  round 2. This may occur with small sample sizes or repeat-poor sequence sets and
  does not necessarily mean repeats will be absent in later rounds.

## [1.0.9]

### Added
- Random number generator seed is now printed at the start of a run, enabling
  exact reproduction of genome samples via the `-srand ####` flag in future runs.
- Final output files are now placed in the same directory as the input database.
- Additional output file containing the seed alignment for each discovered family,
  stored in a Dfam-compatible Stockholm file. New output files are named
  `<database_name>-families.fa` and `<database_name>-families.stk`.
- Support for Dfam_consensus built into this release. Utilities
  `dfamConsensusTool.pl` and `renameIds.pl` are available in the `util/` directory.

## [1.0.8]

### Added
- New `-recoverDir <dir>` option allowing RepeatModeler to restart where it left
  off after a failure.
- BLAST searches can now be run in parallel on shared-memory multiprocessor
  machines via the new `-pa #` option.

## [1.0.7]

### Added
- Support for RMBlastn version 2.2.27.
- `TRFMask` updated to work with the new RepeatMasker 4.0 release.
- New `viewMSA.pl` utility.
