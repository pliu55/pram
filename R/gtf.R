#' @import data.table

readGTF4Exons <- function(fin) {
  indt <- fread(fin, header=F, sep="\t", colClasses=c('character', 'NULL',
                'character', rep('integer', 2), 'NULL', 'character', 'NULL',
                'NULL'))
  setnames(indt, 1:5, c('chrom', 'feature', 'start', 'end', 'strand'))
  dt <- subset(indt, feature == 'exon')
  dt[, feature:=NULL]
  return(dt)
}

readGTF <- function(fgtf, feature = 'exon',
                    infokeys = c('gene_id', 'transcript_id'), nskip=0) {

  indt <- fread(fgtf, header=F, sep="\t", colClasses=c('character', 'NULL',
                'character', rep('integer', 2), 'NULL', 'character', 'NULL',
                'character'), skip=nskip)
  setnames(indt, 1:6, c('chrom', 'ftr', 'start', 'end', 'strand', 'info'))
  dt <- subset(indt, ftr == 'exon')
  for ( infokey in infokeys ) {
    dt[, eval(infokey) := gsub(paste0('.*', infokey, ' ([^;]+);.*'), '\\1',
                               info, perl=T) ]
    dt[, eval(infokey) := gsub('"', '', get(infokey), fixed=T)]
  }
  dt[, `:=`( ftr=NULL, info=NULL )]
  return(dt)
}


writeDataTable2GTF <- function(dt, fout, append=F,
                                        source='unknown', feature='exon',
                                        infokeys=c('gene_id', 'transcript_id')){
  outdt <- copy(dt)
  outdt[, `:=`( source  = source,
                feature = feature,
                score   = '.',
                frame   = '.'     ) ]

  for ( infokey in infokeys ) {
    outdt[, eval(infokey) := paste0(infokey, ' "', get(infokey), '"')]
  }
  outdt[, irow := .I]
  outdt[, attr:=paste0(paste0(.SD, collapse='; '), ';'),
          by=irow, .SDcols=infokeys]
 #outdt[, `:=`(gene_id = NULL, transcript_id = NULL, irow=NULL)]
  outdt <- outdt[, list(chrom, source, feature, start, end, score, strand,
                        frame, attr)]
 #setcolorder(outdt, c('chrom', 'source', 'feature', 'start', 'end', 'score',
 #                     'strand', 'frame', 'attr'))
  write.table(outdt, fout, quote=F, sep="\t", col.names=F, row.names=F,
              append=append)
  cat('File written:', fout, "\n")
}


readGTF4ExonGeneIDTrID <- function(fin) {
  indt <- fread(fin, header=F, sep="\t", colClasses=c('character', 'NULL',
                'character', rep('integer', 2), 'NULL', 'character', 'NULL',
                'character'))
  setnames(indt, 1:6, c('chrom', 'feature', 'start', 'end', 'strand', 'info'))
  dt <- subset(indt, feature == 'exon')
  dt[, `:=`(geneid = gsub('.*gene_id "([^"]+)";.*', '\\1', info, perl=T),
            trid   = gsub('.*transcript_id "([^"]+)";.*', '\\1', info, perl=T)
           )]

  dt[, `:=`(feature=NULL, info=NULL)]
  return(dt)
}


readGTF4Transcripts <- function(fin) {
  dt <- readGTF4ExonGeneIDTrID(fin)
  trdt <- dt[, list( chrom  = last(chrom),
                     start  = min(start),
                     end    = max(end),
                     strand = last(strand),
                     geneid = last(geneid) ), by=trid ]
  return(trdt)
}


readGTF4Genes <- function(fin) {
  dt <- readGTF4ExonGeneIDTrID(fin)
  genedt <- dt[, list( chrom  = last(chrom),
                       start  = min(start),
                       end    = max(end),
                       strand = last(strand),
                       trids  = paste0(unique(trid), collapse=',')
                     ), by=geneid]
  return(genedt)
}
