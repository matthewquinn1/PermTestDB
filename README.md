# Permutation Tests for Detecting Differentially Bound Regions with ChIP-seq Data (PermTestDB)

Chromatin immunoprecipitation sequencing (ChIP-seq) is a process that provides a genome-wide profiling of the bindings of different proteins, namely transcription factors (TFs) and histone modifications (HMs), to DNA. ChIP-seq data are commonly used in a wide variety of biological and medical contexts, such as in studying cancer, gene regulation, and genomic organization. Often, investigators may be interested in comparing various cell lines under different experimental conditions to evaluate the presence of differential binding sites for a given TF or HM. The ChIP-seq read counts are usually presented for small "bins" or "windows" along the genome, typically containing on the order of 75-300 base pairs.

Differential binding pipelines offer a way by which an investigator can formally assess the presence and location of differential binding sites. A number of differential binding pipelines exist, but many use relatively ad hoc methods. This is particularly true in the merging of small genomic sites into larger regions, the handling of autocorrelation within samples, and the handling of biological variance across samples. PermTestDB addresses these limitations by using a permutation test procedure to identify differentially bound regions. The technical details are not discussed here, but will be made available in the corresponding paper [[1]](#1), which is currently in progress.

## Current Status

This package is currently in progress. It is **not** suggested that one downloads and installs the version currently available here. This repository primarily exists as a reference until the package is fully prepared for distribution.

## Authors

* **Matthew Quinn** -  <mjq522@g.harvard.edu>
* **Rafael Irizarry** -  <rafa@ds.dfci.harvard.edu>

## References
<a id="1">[1]</a> 
Paper in progress

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


