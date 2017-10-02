### v6.0.5
* Fix `ssearch36` queries attempting to get sequence beyond start of chromosome/contig.

### v6.0.4

* Modify `Rsupport/libInstall.R` to cope with use of R <3.4 for standard builds
  * If you use newer build of R install VGAM as other bioconductor libs
* Change install method in `libInstall.R` to physically exit with error state if library fails to build
* Add a fixed seed (42) to `metropolis_hastings_inversions.R` to ensure reproducible results

### v6.0.2

* Fixes instability in results due to unsorted keys used in outputs.
  * Perls 5.18+ give variability on duplicate runs as hash keys have enforced randomisation. Found all relevant instances, outputs are now stable between runs and perl versions.
* Shortened paths in *.intermediates.tar.gz

### v6.0.0
* Modified mormal panel is now handled in more efficient way - makes BRASS GROUP step extremly faster

### v5.4.1
* Fix bug in implementation of cytoband option.

### v5.4.0
* pcf function in copynumber package is now species agnostic
* removed NA rows before applying gam function

### v5.3.4
Merge fix for clang compilation, thanks to @jmarshall

### v5.3.2
Fixed #51

### v5.3.1
Fixed #42, #43, #50

### v5.3.0
* Removed dependency on `bedGraphToBigWig` and use `bg2bw` from cgpBigWig (already dep. for PCAP-core).
