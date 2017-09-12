# getGeneInfo
perl script to retrieve various info form a list of HGNC genes, including exon coordinates or domain details.
In particular generates a BED that can be use to design NGS experiences.

## Installation

This progam relies on several data files and one UCSC binary (liftover).

To run it properly you should have inside the root directory of the program:

* a 'data' directory with:
1 one Biomart file with the following columns (tab separated), ideally named 'mart_export.txt':
ENST, ENSP, HGNCid, RefSeq NM, UNIPROT
which can be obtained directly via the ensembl-biomart web [site](http://www.ensembl.org/biomart/):
Select Human genes dataset and attributes: Transcript Stable ID, Protein stable ID in gene attributes, and in external attributes HGNC ID, RefSeq mRNA ID and UniProtKB Gene Name ID .
Or by running the query_biomart.pl file ;before running the script you need to install the biomart API: [http://www.biomart.org/other/install-overview.html]
2 one RefSeq file obtained at [ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene] named LRG_RefSeqGene.txt

* a liftover directory with the liftover binaries and the proper hg19 to hg38 chain file:
1 binaries from [http://hgdownload.soe.ucsc.edu/admin/exe/]
2 chain file from [http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz]

* an active internet connexion as the program needs to connect to togows.org and uniprot.org

* optionnally you woud need [bedtools](http://bedtools.readthedocs.io/en/latest/) installed in your path (default /usr/local/bin, can be modified at the beginning of the script, line $BEDTOOLS = '/usr/local/bin/bedtools')
bedtools is used to merge the bed output with all exons positions.
