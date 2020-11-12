[![Build Status](https://travis-ci.org/pliu55/pram.svg)](https://travis-ci.org/pliu55/pram)
[![bioc](http://www.bioconductor.org/shields/years-in-bioc/pram.svg)](http://bioconductor.org/packages/devel/bioc/html/pram.html)

PRAM: Pooling RNA-seq and Assembling Models
===========================================

Table of Contents
-----------------

* [Introduction](#Introduction)
* [Installation](#Installation)
* [Reference](#Reference)
* [Contact](#Contact)
* [License](#License)

* * *

## <a name='Introduction'></a> Introduction

Pooling RNA-seq and Assembling Models (__PRAM__) is an __Bioconductor__ __R__ 
package that 
utilizes multiple RNA-seq datasets to predict transcript models. The workflow 
of PRAM contains four steps, which is shown in 
the figure below with function names and associated key parameters.  PRAM has a
[vignette](https://bioconductor.org/packages/devel/bioc/vignettes/pram/inst/doc/pram.pdf) that describes each function in details.

<p align='center'>
    <img src="vignettes/workflow_noScreen.jpg" width="400" height="407">
</p>

## <a name='Installation'></a> Installation

### From GitHub

Start __R__ and enter: 

```r
devtools::install_github('pliu55/pram')
```

### From Bioconductor

Start __R__ and enter:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("pram")
```

<!--
- Cufflinks v2.2.1 macOS binary have some issues
- it will report segmentation fault for the same bam file, which Linux 
  Cufflinks runs ok
- Have to use Cufflinks v2.1.1 for macOS instead
-->


<!--
## <a name='Quick-start'></a>Quick start

PRAM provides a function `runPRAM()` to let you run through the whole workflow.

### <a name='predict-only'></a> Predict transcript models only

For a given gene annotation and RNA-seq alignments, you can predict transcript
models in intergenic genomic regions:

```R
runPRAM(in_gtf, in_bamv, out_gtf)
```

- `in_gtf`:  an input GTF file defining genomic coordinates of existing genes. 
             Required to have an attribute of __gene_id__ in the ninth column.
- `in_bamv`:  a vector of input BAM file(s) containing RNA-seq alignments. 
              Currently,
              PRAM only supports strand-specific paired-end RNA-seq with the 
              first mate on the right-most of transcript coordinate, i.e., 
              'fr-firststrand' by Cufflinks definition.
- `out_gtf`:  an output GTF file of predicted transcript models


### <a name='predict-screen'></a> Predict transcript models and screen them by ChIP-seq

If you are interested to predict intergenic transcripts for a particular cell
or tissue type, you can use epigenetic ChIP-seq 
data together with known transcripts and their expression levels to further 
screen intergenic transcript models:
```
runPRAM(in_gtf, in_bamv, out_gtf, in_bedv, training_tpms, training_gtf)
```

- `in_gtf`, `in_bamv`, and `out_gtf` are the same as described above
- `in_bedv`:  A vector of BED file(s) containing ChIP-seq alignments.
- `training_tpms`:  A vector of RSEM quantification results for known
                    transcripts
- `training_gtf`:  A GTF file defining genomic coordinates of known
                   transcripts 

### <a name='Examples'></a> Examples
PRAM has included input examples files in its `extdata/demo/` 
folder.  The table below provides a quick summary of all the example files.

| input argument | file name(s) |
|:--------------:|:------------:|
| `in_gtf`       | [in.gtf](inst/extdata/demo/in.gtf)       |
| `in_bamv`      | [SZP.bam](inst/extdata/demo/SZP.bam), [TLC.bam](inst/extdata/demo/TLC.bam)   |
| `in_bedv`      | H3K79me2.bed.gz, POLR2.bed.gz   |
| `training_tpms`| AED1.isoforms.results, AED2.isoforms.results   |
| `training_gtf` | training.gtf |

You can access example files by `system.file()` in __R__, e.g. for the 
argument `in_gtf`, you can access its example file by

```R
system.file('extdata/demo/in.gtf', package='pram')
```

Below shows usage of `runPRAM()` with example input files: 
##
## Predict transcript models only
##

```R
in_gtf = system.file('extdata/demo/in.gtf', package='pram')

in_bamv = c( system.file('extdata/demo/SZP.bam', package='pram'),
             system.file('extdata/demo/TLC.bam', package='pram') )

pred_out_gtf = tempfile(fileext='.gtf')

runPRAM(in_gtf, in_bamv, pred_out_gtf)
```

##
## Predict transcript models and screen them by ChIP-seq data
##
in_bedv = c( system.file('extdata/demo/H3K79me2.bed.gz', package='pram'),
             system.file('extdata/demo/POLR2.bed.gz',    package='pram') )

training_tpms = c( system.file('extdata/demo/AED1.isoforms.results', package='pram'),
                   system.file('extdata/demo/AED2.isoforms.results', package='pram') )

training_gtf = system.file('extdata/demo/training.gtf', package='pram')

screen_out_gtf = tempfile(fileext='.gtf')

runPRAM(in_gtf, in_bamv, screen_out_gtf, in_bedv, training_tpms, training_gtf)
-->

## <a name="Reference"></a> Reference

[PRAM: a novel pooling approach for discovering intergenic transcripts from large-scale RNA sequencing experiments](https://doi.org/10.1101/gr.252445.119). Peng Liu, Alexandra A. Soukup, Emery H. Bresnick, Colin N. Dewey, and Sündüz Keleş. __Genome Research__ 2020. 30:1655-1666  

For key results reported in the PRAM manuscript and scripts for 
reproducibility, please check out 
[this GitHub repository](https://github.com/pliu55/pram_paper).

## <a name="Contact"></a> Contact

Got a question? Please report it at the [issues tab](https://github.com/pliu55/pram/issues) in this repository.

## <a name="License"></a> License

PRAM is licensed under the [GNU General Public License v3](LICENSE).
