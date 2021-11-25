# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.1.2]
### Changed
* Use `bamstats` for alignment statistics instead of `stats_from_bam`.
* BED input now optional (default is all chr1-22,X,Y).
* Update clair3 to v0.1.8 for bugfixes.

## [v0.1.1]
# Fixed
* Correct HTML report filename.

## [v0.1.0]
# Added
* First release.
* Port original Bash script to Nextflow to increase parallelism
* Implement a sharded phasing strategy to improve performance further.
