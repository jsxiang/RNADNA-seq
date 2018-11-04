% clear

RSIX=RNAseq();
RSIX.load('RSIX_combinedcounts.csv');

%% find control ribozymes
RSIX.findctrls();

%% find validated cdG switches
cdGvalidated={'GCTGTCACCGAGAAACTGTGAAGTATATCTTAAACCTGGGCACTTAAAAGATATATGGAGTTAGTAGTGCAACCTGCTCGGTCTGATGAGTCCCCAATGGGACGAAACAGC',
    'GCTGTCACCGAGAAACTGTGAAGTATATCTTAAACCTGGGCACTTAAAAGATATATGGAGTTAGTAGTGCAACCTGCTCGGTCTGATGAGTCCCCGATGGGACGAAACAGC',
    'GCTGTCACCGAGAAACTGTGAAGTATATCTTAAACCTGGGCACTTAAAAGATATATGGAGTTAGTAGTGCAACCTGCTCGGTCTGATGAGTCCCGGAAGGACGAAACAGC',
    'GCTGTCACCGAGAAACTGTGAAGTATATCTTAAACCTGGGCACTTAAAAGATATATGGAGTTAGTAGTGCAACCTGCTCGGTCTGATGAGTCCCGTAAGGACGAAACAGC',
    'GCTGTCACCGAGAAACTGTGAAGTATATCTTAAACCTGGGCACTTAAAAGATATATGGAGTTAGTAGTGCAACCTGCTCGGTCTGATGAGTCCCGTGAAGGACGAAACAGC',
    'GCTGTCACCGAGAAACTGTGAAGTATATCTTAAACCTGGGCACTTAAAAGATATATGGAGTTAGTAGTGCAACCTGCTCGGTCTGATGAGTCCCGAAAGGACGAAACAGC',
    'GCTGTCACCGAGAAACTGTGAAGTATATCTTAAACCTGGGCACTTAAAAGATATATGGAGTTAGTAGTGCAACCTGCTCGGTCTGATGAGTCCCCAAAGGACGAAACAGC',
    'GCTGTCACCGAGAAACTGTGAAGTATATCTTAAACCTGGGCACTTAAAAGATATATGGAGTTAGTAGTGCAACCTGCTCGGTCTGATGAGTCCCTCAAGGACGAAACAGC'};

val.names={'CCAATG','CCGATG','CGGAA','CGTAA','CGTGAA','CGAAA','CCAAA','CTCAA'};
val.seqs=cdGvalidated;
RSIX.findvalids(val);

%% pair plot RNA and DNA read counts
samplenames={'cdGI.1-','cdGI.1+','cdGI.2-','cdGI.2+','cdGII.1-','cdGII.1+','cdGII.2-','cdGII.2+','cdGI.pDNA','cdGII.pDNA','cdGII3','cdGII3.mut'};
ids=[1 2 3 4 9];

figurename='pair plot cdGI';
% aptaname='cdGII';
RSIX.findpairplots(ids,samplenames,figurename)

%% pair plot RNA and DNA read counts
samplenames={'cdGI.1-','cdGI.1+','cdGI.2-','cdGI.2+','cdGII.1-','cdGII.1+','cdGII.2-','cdGII.2+','cdGI.pDNA','cdGII.pDNA','cdGII3','cdGII3.mut'};
ids=[5:8 10];

figurename='pair plot cdGII';
RSIX.findpairplots(ids,samplenames,figurename)


%% analyze cdGI library
%% how many members in the library were sequenced
cdGI='^GCTGTCACCGCGCACAGGGCAAACCATTCGAAAGAGTGGGACGCAAAGCCTCCGGCCTAAACCAGAAGACATGGTAGGTAGCGGGGTTACCGCGGTCTGATGAGTCC[A|C|G|T][A|C|G|T][A|C|G|T][A|C|G|T][A|C|G|T]GGACGAAACAGC$';
sum(~cellfun('isempty',regexp(RSIX.seqs,cdGI)))


%% find RNA/DNA ratios
minthresh=100;
plotcontrols=true;
plotreplicates=true;
plotgraph=false;
aptamer=cdGI;
clf
cdGI_minus=RSIX.findRNADNAratio(1,3,9,9,'cdGI replicates -lig',aptamer,minthresh,plotcontrols,plotreplicates,plotgraph);

cdGI_plus=RSIX.findRNADNAratio(2,4,9,9,'cdGI replicates +lig',aptamer,minthresh,plotcontrols,plotreplicates,plotgraph);

setfig('cdGI combined replicates');clf
RSIX.plotreplicates(cdGI_minus,cdGI_plus)
title('Cyclic di-GMP-I')
%% find switches
alpha=1e-2;
setfig('combined cdGI switches');clf
filterbyaptamer=true;
aptamer='^GCTGTCACCGCGCACAGGGCAAACCATTCGAAAGAGTGGGACGCAAAGCCTCCGGCCTAAACCAGAAGACATGGTAGGTAGCGGGGTTACCGCGGTCTGATGAGTCC[A|C|G|T][A|C|G|T][A|C|G|T][A|C|G|T][A|C|G|T]GGACGAAACAGC$';

cdGI_all=RSIX.combinesubrun(cdGI_minus,cdGI_plus,filterbyaptamer,aptamer,alpha);
cdGI_all.totalcount=cdGI_minus.totalcount;
cdGI_all.ctrlcounts=cdGI_minus.ctrlcounts;
cdGI_all.ctrlseqs=RSIX.seqs(RSIX.ctrls.ids);
cdGI_all.valcounts=cdGI_minus.valcounts;
cdGI_all.valseqs=RSIX.seqs(RSIX.val.ids);
axis([.5 2.2 .5 2.2])

box on
title('cyclic di-GMP-I')

xlabel('log_1_0(RNA/DNA) (-ligand)')
ylabel('log_1_0(RNA/DNA) (+ligand)')
set(gca,'fontsize',24)

legend('Library sequence','Control ribozymes','Validation candidates','1:1','location','southeast')


%% analyze cdGII library
%% find RNA/DNA ratios
minthresh=100;
plotcontrols=true;
plotreplicates=true;
plotgraph=false;
cdGIIlib='^GCTGTCACCGAGAAACTGTGAAGTATATCTTAAACCTGGGCACTTAAAAGATATATGGAGTTAGTAGTGCAACCTGCTCGGTCTGATGAGTCC[A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G][A|C|T|G]?GGACGAAACAGC$';

aptamer=cdGIIlib;

cdGII_minus=RSIX.findRNADNAratio(5,7,10,10,'cdGII replicates -lig',aptamer,minthresh,plotcontrols,plotreplicates,plotgraph);

cdGII_plus=RSIX.findRNADNAratio(6,8,10,10,'cdGII replicates +lig',aptamer,minthresh,plotcontrols,plotreplicates,plotgraph);

setfig('cdGII combined replicates');clf
RSIX.plotreplicates(cdGII_minus,cdGII_plus)
title('Cyclic di-GMP-II')

%%

setfig('combined cdG switches');clf

filterbyaptamer=true;
aptamer=cdGIIlib;
alpha=1e-2;

cdGII_all=RSIX.combinesubrun(cdGII_minus,cdGII_plus,filterbyaptamer,aptamer,alpha);
cdGII_all.totalcount=cdGII_minus.totalcount;
cdGII_all.ctrlcounts=cdGII_minus.ctrlcounts;
cdGII_all.ctrlseqs=RSIX.seqs(RSIX.ctrls.ids);
cdGII_all.valcounts=cdGII_minus.valcounts;
cdGII_all.valseqs=RSIX.seqs(RSIX.val.ids);

axis([.5 2.2 .5 2.2])

xlabel('log_1_0(RNA/DNA) (-ligand)')
ylabel('log_1_0(RNA/DNA) (+ligand)')

legend('Library sequence','Control ribozymes','Validation candidates','1:1','location','southeast')

%%




%% END

















