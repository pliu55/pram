# PRAM: Pooling RNA-seq for Assembling Models

## Input

### Transcript model prediction
- RNA-seq data
  - FASTQ files, or
  - BAM files
- A GTF file
  - list genomic ranges where model to be built

### Transcript model screening
- RNA-seq
- ChIP-seq

### Transcript model annotation
- RNA-seq
- ChIP-seq
- Motif
- Enhancer
- Mappability: bigWig
- Conservation: UCSC liftOver chain


## Dependent packages

### STAR precompiled binary
- [v2.4.2a](https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz)
  - bin/Linux_x86_64_static/STAR
  - bin/MacOSX_x86_65/STAR

### Cufflinks precompiled binary
- v2.2.1
  - [Linux](http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz)
  - [Mac OS X](http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.OSX_x86_64.tar.gz)

### StringTie precompiled binary
- v1.3.3
  - [Linux](http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.Linux_x86_64.tar.gz)
  - [Mac OS X](http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.OSX_x86_64.tar.gz)

### RSEM

### EBSeq
- on Bioconductor

### Bowtie

### MOSAiCS
- on Bioconductor

### SPP
- on CRAN

### IDR

### liftOver
