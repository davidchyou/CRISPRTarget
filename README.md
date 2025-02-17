# CRISPRTarget

Use case 1: Providing a CRISPRDetect-formatted GFF file of CRISPR arrays, and a set of genomic sequences to be used as DB, the following command will firstly build a BLASTDB and compute an index file for the genomic sequences. Then dinucleotide-shuffle the genomic sequences, build another BLASTDB and compute another index file for the shuffle sequences. Then search for the CRISPR-spacer target and compute P-value using the shuffled BLASTDB. The user can then save the genomic and the shuffled BLASTDBs and index files for later use.

        perl CRISPRTarget.pl -gff sample_crispr_gff/PSA.crispr.gff -user_fasta sample_db/vhdb_selected.fna -dbsize 100000000 -evalue 1 -out test_out -pam_search_all

Use case 2: Providing a CRISPRDetect-formatted GFF of CRISPR arrays, a precomputed BLASTDB and the index file of a set of genomic sequences used in previous call, and those that of the same of the shuffled previously computed, we can search for the CRISPR-spacer target and compute P-value directly with this command.

        perl CRISPRTarget.pl -gff sample_crispr_gff/PSA.crispr.gff -db USER_DB/vhdb_selected.fna -ctrl_db USER_SHUFFLED_DB/vhdb_selected.fna -dbsize 100000000 -evalue 1 -out test_out_2 -pam_search_all

Note that in Use case2, we need BLASTDB and the index file from the genomes and the shuffled sequences.
