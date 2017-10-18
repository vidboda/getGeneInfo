# getGeneInfo
Perl script to retrieve various info from a list of HGNC genes, including exon coordinates or domain details.
In particular generates a BED that can be used to design NGS experiences.

## Installation

This progam relies on several data files and one UCSC binary (liftover).

To run it properly you should have inside the root directory of the program:

* a 'data' directory with:

1. one Biomart file with the following columns (tab separated), ideally named 'mart_export.txt':
ENST, ENSP, HGNCid, RefSeq NM, UNIPROT

which can be obtained directly via the ensembl-biomart web [site](http://www.ensembl.org/biomart/):

Select Human genes dataset and attributes: Transcript Stable ID, Protein stable ID in gene attributes, and in external attributes HGNC ID, RefSeq mRNA ID and UniProtKB Gene Name ID .

Or by running the query_biomart.pl file ;before running the script you need to install the biomart API: <http://www.biomart.org/other/install-overview.html>

2. one RefSeq file obtained at <ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene> named LRG_RefSeqGene.txt

* a liftover directory with the liftover binaries and the proper hg19 to hg38 chain file:


1. binaries from <http://hgdownload.soe.ucsc.edu/admin/exe/>

2. chain file from <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz>

* an active internet connexion as the program needs to connect to togows.org and uniprot.org

* optionnally you would need [bedtools](http://bedtools.readthedocs.io/en/latest/) installed in your path (default /usr/local/bin, can be modified at the beginning of the script, line $BEDTOOLS = '/usr/local/bin/bedtools')
bedtools is used to merge the bed output with all exons positions.

## How to run

Once you have your files in the 'data' directory and the UCSC liftover for your system working in the 'liftover' directory, plus an active internet connection, you can run the script with:

```bash
perl getGeneInfo.pl -l GENE_LIST.txt -g hg19 -o 50 -n
```

with gene_list.txt being a text file with HGNC gene names to process:

USH2A

CLRN1

CFTR

..

the -g option is to fill with hg19 or hg38 human genome assemblies

the -o option defines an offset to be applied to each exon of the genes, if you want a bed with exons +/- 50 intronic base pairs for NGS designs, for example

and the -n option makes the script run even for genes for which the NCBI has no NG_... accesson number.

Another -s option will build sql files ready to be inserted in our in house database system.

If the script encounters critical errors for a gene, this gene will be reported in an gene_list_error.txt file (e.g. when networking issues) which can be used to rerun the script.


## What you will get

Several files in the 'results' directory. GENE_LIST is the name of the file and GENOME the version you provided as input.

* a GENE_LIST_info.txt file which will contain info on accession numbers (NCBI, Ensembl...), exon genomic locations, protein domains for each RefSeq isoform

* a GENE_LIST_LOVD_domains.txt will contain a summary of the Uniprot protein domains, which can be used in LOVD systems to create menus

* a GENE_LIST_exons_GENOME.sorted.merged.bed, is the bed file that you can load in [UCSC genome browser](https://www.genome.ucsc.edu/) which will show the exons. Beware, it is compliant with UCSC [convention for BED](http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/) of being 0-based for start positions and 1-based for end positions. This file can be sent to your NGS probes provider.

* a GENE_LIST_error.txt listing the genes for which the program had serious errors

* optionally a genes_SQL.sql file you probably don't need

## Credits

This script uses togoWS web services, http://togows.org/

Toshiaki Katayama, Mitsuteru Nakao and Toshihisa Takagi: TogoWS: integrated SOAP and REST APIs for interoperable bioinformatics Web services. Nucleic Acids Research 2010, 38:W706-W711. doi:10.1093/nar/gkq386, PMID:20472643 (NAR Web Server Issue 2010)
