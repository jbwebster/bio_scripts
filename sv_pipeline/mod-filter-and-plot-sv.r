# Copied from /gscmnt/gc5111/research/hdang/PCNI/wustl-data/combine/structural-variants/filter-and-plot-sv.r

setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/Pipeline")

annotateSV <- function(z){
  dup = z$svtype == 'DUP'
  del = z$svtype == 'DEL'
  ar.start = 66764464
  ar.stop = 66950461
  ar.exon3.stop = 66863249
  eh.start = 66118920
  eh.stop = 66127943
  t.start = 42836477
  t.stop = 42880085
  e.start = 39739182
  e.stop = 40033704
  
  ar.dup = (z$chr1 == 'chrX' & z$svtype == 'DUP' 
            & z$pos1 < ar.start & z$pos2 > ar.exon3.stop)
  eh.dup = (z$chr1 == 'chrX' & z$svtype == 'DUP' 
            & z$pos1 < eh.start & z$pos2 > eh.stop)
  ar.tail.del = (z$chr1 == 'chrX' & z$svtype == 'DEL' 
                 & z$pos1 > ar.exon3.stop & z$pos1 < ar.stop)
  
  t2e = (z$chr1 == 'chr21' & z$svtype == 'DEL'
         & z$pos1 > e.start & z$pos1 < e.stop 
         & z$pos2 > t.start & z$pos2 < t.stop)
  
  z$annot = ''
  if (any(ar.dup)){z$annot[ar.dup] = paste0(z$annot[ar.dup], 'AR.dup,')}
  if (any(eh.dup)){z$annot[eh.dup] = paste0(z$annot[eh.dup], 'Enhancer.dup,')}
  if (any(ar.tail.del)){z$annot[ar.tail.del] = paste0(z$annot[ar.tail.del],
                                                      'AR.tail.del,')}
  if (any(t2e)){z$annot[t2e] = paste0(z$annot[t2e], 'TMPRSS2-ERG')}
  return(z)
}

swapEnds <- function(bp){
  bp$pos1 = (bp$start1 + bp$stop1)/2
  bp$pos2 = (bp$start2 + bp$stop2)/2
  for (i in 1:nrow(bp)){
    if (bp$chr1[i] != 'chrX' || bp$pos1[i] > maxx || bp$pos1[i] < minx){
      if (bp$chr2[i] == 'chrX' && bp$pos2[i] > minx && bp$pos2[i] < maxx){
        bp[i,c('chr1', 'start1', 'stop1', 'chr2', 'start2', 'stop2',
               'pos1', 'pos2')] = bp[i, c('chr2', 'start2', 'stop2', 'chr1',
                                          'start1', 'stop1', 'pos2', 'pos1')]
        bp[i,c('strand1', 'strand2')] = bp[i, c('strand2', 'strand1')]
      }else{ 
        message('WARN: line ', i, ' out of range')
      } 
    } 
  } 
  return(bp)
} 



sv.files <- list("tulane.rerunHM", "tulane.withHM", "wustl1", 
                 "wustl2", "wustl3")
out.dir <- "custom_ann"
version <- list(T, F)

for(sv.file in sv.files) {
  for(AR in version) {
    # SV calls
    x = read.table(paste0(sv.file, ".tsv"), header=T, stringsAsFactors=F, sep='\t', quote='',
                   comment.char='', na.string='')
    colnames(x) = sub('end', 'stop', colnames(x))
    colnames(x)[colnames(x) == 'name'] = 'sample_svtype'
    colnames(x)[colnames(x) == 'chrom1'] = 'chr1'
    colnames(x)[colnames(x) == 'chrom2'] = 'chr2'
    x$sample = sub('/.*$', '', x$sample_svtype)
    x$svtype = sub('^.*/', '', x$sample_svtype)
    x$pos1 = (x$start1 + x$stop1)/2
    x$pos2 = (x$start2 + x$stop2)/2
    x$loc = paste0(x$chr1, ':', x$start1, '-', x$stop1, ';', x$chr2, ':', x$start2, '-', x$stop2)
    x$bplen1 = (x$stop1 - x$start1)
    x$bplen2 = (x$stop2 - x$start2)
    x$sv.distance = 123456789
    x$sv.distance[x$chr1 == x$chr2] = abs(x$pos1 - x$pos2)[x$chr1 == x$chr2]
    x$id = paste0(x$sample, '/', 1:nrow(x), '/loc=', x$loc, '\nsvlen=',x$sv.distance,
                  '/SR=', x$plasma_split_reads,'/DPE=', x$plasma_pe_reads,
                  '/bplens=', x$bplen1, ';', x$bplen2)
    y = x
    chrs = paste0('chr', c(1:22, 'X', 'Y'))
    y = y[(y$chr1 %in% chrs) & (y$chr2 %in% chrs),]
    
    # hg19
    if (AR){
      # AR & enhancer
      chrom = 'chrX'
      gstart = 66764464
      gstop = 66950461
      estart = 66118920
      estop = 66127943
      prefix = 'AR'
      w = 9
      h = 5
      
      e = read.table('data/tulane/sv/AR.exons.tsv',
                     header=F, stringsAsFactors=F)
      colnames(e) = c('chr', 'src', 'type', 'start', 'stop')
      
    }else{
      
      # TMPRSS2 & ERG
      chrom = 'chr21'
      gstart = 42836477
      gstop = 42880085
      estart = 39739182
      estop = 40033704
      prefix = 'TMPRSS2-ERG'
      w = 8
      h = 6
      e = c()
      
    }
    
    D = 100000
    maxx = gstop + D
    minx = estart - D
    
    exon.pos= c(minx,maxx, e$start, e$stop)
    
    # filtering criteria
    minDistance = 2000
    maxBpLen = 5*100
    maxBpLenBND = 1
    maxDisPairs = 300 # Originally 150
    maxSplitReads = 300 # Originally 150
    minDisPairs = 0 # Handled by pipeline
    minSplitReads = 0 # Handled by pipeline
    svtypes = c('DEL', 'DUP')
    
    z = y[(y$chr1 == chrom & y$pos1 < maxx & y$pos1 > minx) | (y$chr2 == chrom & y$pos2 < maxx & y$pos2 > minx),]
    
    #z = y[y$chr1 == chrom & y$chr2 == chrom,]
    z = z[z$sv.distance > minDistance & !is.na(z$sv.distance),]
    z = z[z$svtype %in% svtypes,]
    z = z[z$bplen1 <= maxBpLen & z$bplen2 <= maxBpLen,]
    z = z[!(z$svtype == 'BND' & (z$bplen1 > maxBpLenBND | z$bplen2 > maxBpLenBND)),]
    z = z[z$plasma_split_reads >= minSplitReads & z$plasma_pe_reads >= minDisPairs,]
    z = z[z$plasma_split_reads <= maxSplitReads & z$plasma_pe_reads <= maxDisPairs,]
    z$plasma_pe_sr_reads <- z$plasma_split_reads + z$plasma_pe_reads
    z = z[z$plasma_pe_sr_reads >= minTotalSupportReads,]
    z = annotateSV(z)
    z$id = paste0(z$id, '---', z$annot)
    #z$patient = sub('(PB\\d+).*$', '\\1', z$sample)
    l = z
    
    res = z[, c('sample', 'annot', 'chr1', 'start1', 'stop1',
                'chr2', 'start2', 'stop2',
                'strand1', 'strand2', 'svtype', 'plasma_pe_reads', 'plasma_split_reads',
                'plasma_pe_sr_reads')]
    
    write.table(res, file=paste0(out.dir, '/', prefix, '.', sv.file, '.tsv'), sep='\t', quote=F, row.names=F)
    
    
  }
}

