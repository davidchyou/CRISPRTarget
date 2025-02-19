SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
CURDIR=`pwd -P`

INFILE=$1
OUTDIR=$2
BASENAME=$(basename -- $INFILE)

mkdir -p $OUTDIR
perl -pe 's/>(.*)/>\1\t/g; s/\n//g; s/>/\n>/g' $INFILE | grep -v '^$' | \
perl -F'\t' -lane 'BEGIN{$count=0;}{$count++; $header=@F[0]; ($line)=($header=~/^>(.*)/); $id_new="SEQ_".$count."|"."'$BASENAME'"; $header_new=">$id_new $line"; $seq=uc(@F[1]); print "$header_new\n$seq";}' > $OUTDIR/input_genomes.fna

$SCRIPTPATH/bin/esl-shuffle -d -N 4 $OUTDIR/input_genomes.fna > $OUTDIR/$BASENAME 2>/dev/null
rm -f $OUTDIR/input_genomes.fna

$SCRIPTPATH/bin/makeblastdb -in $OUTDIR/$BASENAME -out $OUTDIR/$BASENAME -dbtype nucl -blastdb_version 4 -parse_seqids > /dev/null 2>/dev/null
$SCRIPTPATH/bin/samtools faidx $OUTDIR/$BASENAME
echo $OUTDIR/$BASENAME
