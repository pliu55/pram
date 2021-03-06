% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runPRAM.R
\name{runPRAM}
\alias{runPRAM}
\title{Predict intergenic transcript models from RNA-seq}
\usage{
runPRAM(in_gtf, in_bamv, out_gtf, method, cufflinks = "",
    stringtie = "", taco = "")
}
\arguments{
\item{in_gtf}{An input GTF file for defining genomic coordinates of
existing
genes.  Required to have `gene_id` in the attribute column
(column 9)}

\item{in_bamv}{A character vector of input BAM file(s). If mode 'cf'
or 'st' is used, only one input RNA-seq BAM file is allowed.
Currently, PRAM only supports strand-specific paired-end data
with the first mate on the right-most of transcript
coordinate, i.e., 'fr-firststrand' by Cufflinks's definition.}

\item{out_gtf}{An output GTF file of predicted transcript models}

\item{method}{A character string defining PRAM's model building method.
Current available methods are:
\itemize{
    \item plcf: pooling   + cufflinks
    \item plst: pooling   + stringtie
    \item cfmg: cufflinks + cuffmerge
    \item stmg: stringtie + merging
    \item cftc: cufflinks + taco
    \item cf:   cufflinks
    \item st:   stringtie
}
Default: 'plcf'}

\item{cufflinks}{Cufflinks executable.  Required by mode 'plcf',
'cfmg', and 'cf'.  For mode 'cfmg', executable files of
Cuffmerge, Cuffcompare, and gtf_to_sam from the Cufflinks
suite are assumed to be under the same folder as
Cufflinks.
All the executables are available to download for
Linux \url{http://cole-trapnell-lab.github.io/cufflinks/
assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz}
and MacOS \url{http://cole-trapnell-lab.github.io/
cufflinks/assets/downloads/
cufflinks-2.1.1.OSX_x86_64.tar.gz}. 
Souce code can be obtained from 
\url{http://cole-trapnell-lab.github.io/cufflinks/}. 
Default: ''}

\item{stringtie}{StringTie executable file.  Required by mode
'plst', 'stmg', and 'st'.
Executable can be downloaded for Linux 
\url{http://ccb.jhu.edu/software/stringtie/dl/
stringtie-1.3.3b.Linux_x86_64.tar.gz} and MacOS
\url{http://ccb.jhu.edu/software/stringtie/dl/
stringtie-1.3.3b.OSX_x86_64.tar.gz}.
Souce code can be obtained from 
\url{https://ccb.jhu.edu/software/stringtie/}. 
Default: ''}

\item{taco}{TACO executable file. Required by mode 'cftc'.
Executable can be downloaded for Linux
\url{https://github.com/tacorna/taco/releases/
download/v0.7.0/taco-v0.7.0.Linux_x86_64.tar.gz} and MacOS
\url{https://github.com/tacorna/taco/releases/
download/v0.7.0/taco-v0.7.0.OSX_x86_64.tar.gz}.
Souce code can be obtained from 
\url{https://tacorna.github.io}. 
Default: ''}
}
\value{
None
}
\description{
Predict intergenic transcript models from RNA-seq
}
\examples{

in_gtf = system.file('extdata/demo/in.gtf', package='pram')

in_bamv = c(system.file('extdata/demo/SZP.bam', package='pram'),
            system.file('extdata/demo/TLC.bam', package='pram') )

pred_out_gtf = tempfile(fileext='.gtf')

## assuming the stringtie binary is in folder /usr/local/stringtie-1.3.3/
## you can run runPRAM() by the following example
##
# runPRAM(in_gtf, in_bamv, pred_out_gtf, method='plst',
#         stringtie='/usr/local/stringtie-1.3.3/stringtie')

}
