# Copied from /gscmnt/gc5111/research/hdang/PCNI/wustl-data/combine/structural-variants/filter-and-plot-sv.r

setwd("/Users/jacewebster/Desktop/GradSchool/Lab/Projects/Pipeline")

#sv.file = 'data/tulane/sv/original_sv.12cols.tsv'
#sv.file = 'data/tulane/sv/hmrerun.v1.12cols.tsv'
#sv.file = 'data/tulane/sv/w93_corrected.12cols.tsv'
#out.dir = 'data/tulane/sv/results'
sv.file = '/Users/jacewebster/Desktop/tmp/w93_hmrerun_share/tulane_hmrerun.12cols.tsv'
out.dir = 'data/tulane/sv/results'

sv.file = '/Users/jacewebster/Desktop/GradSchool/Lab/Projects/Pipeline/data/wustl_validation/outputs/sv/all_batches.12cols.tsv'
out.dir = 'data/wustl_validation/results'

#sv.file = 'tulane.sv.20cols.bedpe'
#out.dir = 'results-tulane'

#sv.file = 'wustl.sv.20cols.bedpe'
#out.dir = 'results-wustl'

#sv.file = 'tumor-and-simulation.sv.20cols.bedpe'
#out.dir = 'results-tumor-simulation'

#dir.create(out.dir)

AR = F

# SV calls
x = read.table(sv.file, header=T, stringsAsFactors=F, sep='\t', quote='',
    comment.char='', na.string='')

# rename some cols (end1, end2 cols to stop1, stop2; name to sample_svtype)
colnames(x) = sub('end', 'stop', colnames(x))
colnames(x)[colnames(x) == 'name'] = 'sample_svtype'
colnames(x)[colnames(x) == 'chrom1'] = 'chr1'
colnames(x)[colnames(x) == 'chrom2'] = 'chr2'


#if ('found.in.healthy' %in% colnames(x)){x=x[!x$found.in.healthy,]}
#colnames(x) = c('chr1', 'start1', 'stop1', 'chr2', 'start2', 'stop2', 'sample.svtype', 'score', 'strand1', 'strand2')
x$sample = sub('/.*$', '', x$sample_svtype)
x$patient = sub('(PB\\d+).*$', '\\1', x$sample)
cat('Number of patients: ', length(unique(x$patient)), '\n')
x$svtype = sub('^.*/', '', x$sample_svtype)
x$pos1 = (x$start1 + x$stop1)/2
x$pos2 = (x$start2 + x$stop2)/2
x$loc = paste0(x$chr1, ':', x$start1, '-', x$stop1, ';', x$chr2, ':', x$start2, '-', x$stop2)
#x$normal.ref.cnt = x$normal.ref.paired.cnt + x$normal.ref.split.cnt
#x$normal.alt.cnt = x$normal.alt.paired.cnt + x$normal.alt.split.cnt 
#x$tumor.ref.cnt = x$tumor.ref.paired.cnt + x$tumor.ref.split.cnt
#x$tumor.alt.cnt = x$tumor.alt.paired.cnt + x$tumor.alt.split.cnt
#x$tumor.vaf = x$tumor.alt.cnt/(x$tumor.ref.cnt + x$tumor.alt.cnt)
x$bplen1 = (x$stop1 - x$start1)
x$bplen2 = (x$stop2 - x$start2)
x$sv.distance = 123456789
x$sv.distance[x$chr1 == x$chr2] = abs(x$pos1 - x$pos2)[x$chr1 == x$chr2]
x$id = paste0(x$sample, '/', 1:nrow(x), '/loc=', x$loc, '\nsvlen=',x$sv.distance,
    '/SR=', x$plasma_split_reads,'/DPE=', x$plasma_pe_reads,
    '/bplens=', x$bplen1, ';', x$bplen2)
y = x

# only keeps events hitting primary chromosomes
chrs = paste0('chr', c(1:22, 'X', 'Y'))
y = y[(y$chr1 %in% chrs) & (y$chr2 %in% chrs),]


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
 


# filtering criteria
minDistance = 2000
maxBpLen = 5*100
maxBpLenBND = 1
maxDisPairs = 150
maxSplitReads = 150
minDisPairs = 0
minSplitReads = 0
minTotalSupportReads = 2#5
svtypes = c('DEL', 'DUP', 'INV', 'BND', 'INS')
#svtypes = c('DEL', 'DUP', 'INV', 'INS')
#svtypes = c('DEL')
svtypes = c('DEL', 'DUP')



# filter
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
z$patient = sub('(PB\\d+).*$', '\\1', z$sample)
l = z

res = z[, c('patient', 'sample', 'annot', 'chr1', 'start1', 'stop1',
    'chr2', 'start2', 'stop2',
    'strand1', 'strand2', 'svtype', 'plasma_pe_reads', 'plasma_split_reads',
    'plasma_pe_sr_reads')]

write.table(res, file=paste0(out.dir, '/', prefix, '.all_batches.tsv'), sep='\t', quote=F, row.names=F)

lcols = c('svtype', 'bplen1', 'bplen2', 'plasma_pe_reads', 'plasma_split_reads', 'id')
mcols = c('svtype', 'bplen1', 'bplen2', 'tumor.alt.paired.cnt', 'tumor.alt.split.cnt', 'id')

#print(l[order(l$svtype), c('svtype', 'bplen1', 'bplen2', 'plasma_pe_reads', 'plasma_split_reads')])
#sel = with(l, plasma_pe_reads > 1 & plasma_split_reads > 1)
#l[sel & l$sample=='PB087C2', c('svtype', 'bplen1', 'bplen2', 'plasma_pe_reads', 'plasma_split_reads')]

#l = l[sel,]
#print(m[order(m$svtype), c('svtype', 'bplen1', 'bplen2', 'tumor.alt.paired.cnt', 'tumor.alt.split.cnt')])
#m[m$sample=='PB087C2', c('svtype', 'bplen1', 'bplen2', 'tumor.alt.paired.cnt', 'tumor.alt.split.cnt')]

prepSV4Plot <- function(x, chrom){
    bnd = x$svtype == 'BND'
    chr1 = x$chr1 == chrom & x$pos1 < maxx & x$pos1 > minx
    chr2 = x$chr2 == chrom & x$pos2 < maxx & x$pos2 > minx
    x$bnd.pos = NA; x$bnd.chr = NA; x$bnd.hjust = NA
    sel = bnd & chr1 & !chr2;
    if (any(sel)){
        x$bnd.pos[sel] = x$pos1[sel]; x$bnd.chr[sel] = x$chr2[sel];
        x$bnd.hjust[sel] = ifelse(x$strand1[sel] == '+', 1, 0)
    }
    sel = bnd & !chr1 & chr2;
    if (any(sel)){
        x$bnd.pos[sel] = x$pos2[sel]; x$bnd.chr[sel] = x$chr1[sel]
        x$bnd.hjust[sel] = ifelse(x$strand2[sel] == '+', 0, 1)
    }
    sel = x$bnd.hjust == 1 & !is.na(x$bnd.hjust);
    if (any(sel)){x$bnd.chr[sel] = paste0(x$bnd.chr[sel], ' -  ')}
    sel = x$bnd.hjust == 0 & !is.na(x$bnd.hjust);
    if (any(sel)){x$bnd.chr[sel] = paste0('  - ', x$bnd.chr[sel])}
    if(any(chr1 & chr2 & bnd)){ #trans within chrom
        cat('WARN: trans within plot region detected. REMOVED\n')
        print(x[chr1 & chr2 & bnd, c('chr1', 'pos1', 'chr2', 'pos2', 'svtype', 'sample')])
        x = x[!(chr1 & chr2 & bnd),]
    }


    return(x)
}


library(ggplot2)
library(reshape2)
library(gridExtra)
colors = c('DEL'='blue', 'DUP'='orange', 'INS'='black', 'INV'='green', 'BND'='cyan')
plotSV <- function(m, title=NULL, left=NULL, right=NULL){
     if (is.null(left)){left = minx}
     if (is.null(right)){right = maxx}
     mp = (ggplot(m)
         + geom_segment(data=m[m$svtype != 'BND',],
            aes(x=pos1, xend=pos2, y=id, yend=id, color=svtype))
         + geom_point(data=m[m$svtype == 'BND',], aes(x=bnd.pos, y=id), shape=4)
         + geom_text(data=m[m$svtype == 'BND',],
            aes(x=bnd.pos, y=id, label=bnd.chr, hjust=bnd.hjust), size=2)
         + geom_vline(xintercept=c(gstart, gstop), color='blue', linetype='dashed', size=0.2)
         + geom_vline(xintercept=c(estart, estop), color='red', linetype='dashed', size=0.2)
         + geom_vline(xintercept=exon.pos, color='darkgray', linetype='dashed',
            size=0.1, alpha=0.7)
         + theme_bw(base_size=9)
         + coord_cartesian(xlim=c(left,right))
         + xlab(chrom)
         + scale_color_manual(values=colors)
         + ggtitle(title)
    )
    return(mp)

}



l = prepSV4Plot(l, chrom)
lp = plotSV(l)
ggsave(lp, file=paste0(out.dir, '/', prefix, '.w93_corrected.pdf'), width=w, height=h)
system(paste0('rclone copy ', out.dir, '/', prefix, '.pdf hdng:tmp/'))

stop()

lp2 = plotSV(l[grepl('PB078', l$sample),])
ggsave(lp2, file=paste0(out.dir, '/', prefix, '-PB078.pdf'), width=w, height=h)
system (paste0('rclone copy ', out.dir, '/TMPRSS2-ERG-PB078.pdf hdng:tmp/'))

lp2 = plotSV(l[grepl('PB202', l$sample),])
ggsave(lp2, file=paste0(out.dir, '/', prefix, '-PB202.pdf'), width=10, height=5)
system (paste0('rclone copy ', out.dir, '/TMPRSS2-ERG-PB202.pdf hdng:tmp/'))

lp3 = plotSV(l[grepl('PB202', l$sample) & l$svtype == 'DUP',], left=64400000, right=67000000)
ggsave(lp3, file=paste0(out.dir, '/', prefix, '-PB202.pdf'), width=10, height=5)
system (paste0('rclone copy ', out.dir, '/', prefix, '-PB202.pdf hdng:tmp/'))



#ggsave(lp, file=paste0('results/', prefix, '.pdf'), width=w, height=40)

l[order(l$svtype), c('plasma_pe_reads', 'plasma_split_reads', 'svtype')]


stop()
 
# merge with CNA results to count AR gain
cn = read.table('/gscmnt/gc5111/research/hdang/PCNI/wustl-data/combine/cna/cnvkit4/results/cna-call-results.tsv', header=T, stringsAsFactors=F,sep='\t')
cn = cn[grepl('^TU|^PB', cn$sample),]
cn = cn[!grepl('UF|HT', cn$sample),]
cn$patient = sub('(PB\\d+).*$', '\\1', cn$sample)

ar.gain.patients = unique(cn$patient[cn$cna.call == 'gain' & cn$gene == 'AR'])
eh.gain.patients = unique(cn$patient[cn$cna.call == 'gain' & cn$gene == 'ARENHCR'])

ar.dup.patients = unique(res$patient[grepl('AR.dup', res$annot)])
eh.dup.patients = unique(res$patient[grepl('Enhancer.dup', res$annot)])


cat('AR dup:', length(ar.dup.patients), '\n')
cat('Enhancer dup:', length(eh.dup.patients), '\n')
cat('AR/Enhancer dup:', length(unique(c(ar.dup.patients, eh.dup.patients))), '\n')
cat('AR gain:', length(ar.gain.patients), '\n')
cat('Enhancer gain:', length(eh.gain.patients), '\n')
cat('AR/Enhancer gain:', length(unique(c(ar.gain.patients, eh.gain.patients))), '\n')
cat('AR/enhancer gain+dup:', length(unique(c(ar.gain.patients, ar.dup.patients,
    eh.gain.patients, eh.dup.patients))), '\n')

cat('# patients w/ TMPRSS2-ERG:', length(unique(res$patient[grepl('TMPRSS2-ERG', res$annot)])), '\n')

pten.loss.patients = unique(cn$patient[cn$cna.call == 'loss' & cn$gene == 'PTEN'])
tp53.loss.patients = unique(cn$patient[cn$cna.call == 'loss' & cn$gene == 'TP53'])
rb1.loss.patients = unique(cn$patient[cn$cna.call == 'loss' & cn$gene == 'RB1'])

cat('PTEN loss:', length(pten.loss.patients), '\n')
cat('TP53 loss:', length(tp53.loss.patients), '\n')
cat('RB1 loss:', length(rb1.loss.patients), '\n')

