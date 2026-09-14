#Make all binaries executable
chmod -R 777 bin

#Download plage and plasmid datbases
mkdir DB
cd DB
wget https://zenodo.org/records/21869551/files/phage.fa.zip
wget https://zenodo.org/records/21869551/files/plasmid.fa.zip
unzip *.zip
cd ../

#Construct BLASTDB and the index file from FASTA and shuffle, then search for targets.
perl CRISPRTarget.pl \
-gff sample_crispr_gff/PSA.crispr.gff \
-user_fasta DB/phage.fa \
-dbsize 100000000 \
-evalue 1 \
-out test_out \
-pam_search_all
-make_user_db

#Search for targets from precomputed BLASTDB and the index file.
perl CRISPRTarget.pl \
-gff sample_crispr_gff/PSA.crispr.gff \
-db USER_DB/vhdb_selected.fna \
-ctrl_db USER_SHUFFLED_DB/vhdb_selected.fna \
-dbsize 100000000 \
-evalue 1 \
-out test_out_2 \
-pam_search_all
