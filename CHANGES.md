# Changes

## v6.3.3

* Singularity specific execution issue corrected

## v6.3.2

* Fix edgecase on RG line handling (#99)
* Ensure attempt to run on completion give exit==0 (#98)

## v6.3.1

* Handle undefined value edge case #95.
  * Thanks to @udvzol for reproducible test data.
* Updated travis build (#89)
* Fixed up dockerfile to reduce build time when testing.

## v6.3.0

* Adds stand alone (supported) docker container

## v6.2.1

* Resolves #85 - command line help
* Resolves #86 - blatSrc code moved

## v6.2.1

* Fixes #80 - allow RGID == 0
* Resolves #84 - use helper threads where possible

## v6.2.0

* Fixed bug where subset of intervals does not have reads -issue#65
* Informative error when .bas file is absent - issue#71
* temp file were now kept after running filter step - issue#73
* Fixed issue#66
* Added script to process centro/telomere data - issue#70
* Fixed bug where intermediate file name matches with one of the sample

## v6.1.2

* Change tabix->query to tabix->query_full

## v6.1.1

* Reduce I/O for small cpu overhead in coverage step
* Fix assembly step to cope with bam.csi and cram.crai
* R lib issues and disable diagnostic plots.

## v6.1.0

* Add travis-ci
* Fix #46
* Correct removal of overlapping reads (small number of bad reads passed to grouping phase).
  * See 508112837271853152ca48d1ecad5a443fadefbb
* Reworked readme and reference dockstore/docker/singularity

## v6.0.5

* Fix `ssearch36` queries attempting to get sequence beyond start of chromosome/contig.

## v6.0.4

* Modify `Rsupport/libInstall.R` to cope with use of R <3.4 for standard builds
  * If you use newer build of R install VGAM as other bioconductor libs
* Change install method in `libInstall.R` to physically exit with error state if library fails to build
* Add a fixed seed (42) to `metropolis_hastings_inversions.R` to ensure reproducible results

## v6.0.2

* Fixes instability in results due to unsorted keys used in outputs.
  * Perls 5.18+ give variability on duplicate runs as hash keys have enforced randomisation. Found all relevant instances, outputs are now stable between runs and perl versions.
* Shortened paths in `*.intermediates.tar.gz`

## v6.0.0

* Modified mormal panel is now handled in more efficient way - makes BRASS GROUP step extremly faster

## v5.4.1

* Fix bug in implementation of cytoband option.

## v5.4.0

* pcf function in copynumber package is now species agnostic
* removed NA rows before applying gam function

## v5.3.4

Merge fix for clang compilation, thanks to @jmarshall

## v5.3.2

Fixed #51

## v5.3.1

Fixed #42, #43, #50

## v5.3.0

* Removed dependency on `bedGraphToBigWig` and use `bg2bw` from cgpBigWig (already dep. for PCAP-core).
