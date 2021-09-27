# Script archive: Lineage recording in human cerebral organoids
This repository includes scripts to preprocess the iTracer data, as well as scripts to reproduce analysis done on the paper "Lineage recording in human cerebral organoids". In the paper, we describe iTracer, a lineage recorder that combines reporter barcodes with inducible CRISPR/Cas9 scarring, and is compatible with single-cell and spatial transcriptomics. In the paper we also describe how we used the tool to investigate cell lineages in cerebral organoid development.

The paper is published in Nature Methods. Its earlier version is also available in bioRxiv: [Lineage recording reveals dynamics of cerebral organoid regionalization](https://www.biorxiv.org/content/10.1101/2020.06.19.162032v3).

## Folder description
Folder | Description | Related figures
--- | --- | ---
```analysis_lightsheet``` | The scripts to visualize and analyze the lightsheet nuclei-tracked data, given the data frame derived from the Mastodon-exported XML files. | Fig. 3, ExtData Fig. 6
```analysis_microdissection``` | The scripts to analyze the scRNA-seq data of the iTracer-ed microdissected cerebral organoid | Fig. 5, ExtData Fig. 9
```analysis_neuroepithelium``` | The scripts to analyze the scRNA-seq data of the iTracer-ed D15 cerebral organoids | Fig. 4, ExtData Fig. 8
```analysis_perturb``` | The scripts to analyze the iTracer-perturb data | Fig. 6, ExtData Fig. 10
```analysis_visium``` | The notebooks and scripts to analyze the visium data of the iTracer-ed cerebral organoids | Fig. 2, ExtData Fig. 5
```analysis_whole-organoids``` | The scripts to analyze the scRNA-seq data of the iTracer-ed whole cerebral organoids | Fig. 1/4, ExtData Fig. 2/3/4/7
```convert_mastodon_xml2df``` | The python script to convert the Mastodon-exported XML files to a table for further analysis | 
```ext``` | External gene lists used in the analysis |
```extract_barcode_scars``` | The codes to extract barcode and scar information from the post-CellRanger BAM files of the barcode/scar cDNA libraries | 
```r_functions``` | The R functions used in the analysis |

## Related data
The raw data reported in this work can be found in ArrayExpress:
* E-MTAB-10974: scRNA-seq data of iTracer whole-cerebral-organoid samples
* E-MTAB-10973: scRNA-seq data of iTracer early cerebral organoid (d15) samples
* E-MTAB-10972: 10x visium data of iTracer cerebral organoid slices
* E-MTAB-10971: scRNA-seq data of the iTracer-perturb (TSC2) organoids

The processed data is available in Mendeley Data: https://doi.org/10.17632/nj3p3pxv6p. The bulk sequencing raw data of barcode and scar libraries are also available in the same Mendeley Data repository.

## Contact
If you have any question related to the data and the analysis, please contact Dr. Zhisong He (zhisong.he(at)bsse.ethz.ch), Prof. J. Gray Camp (jarrettgrayson.camp(at)unibas.ch), Prof. Barbara Treutlein (barbara.treutlein(at)bsse.ethz.ch)
