use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $attr_tbl = "";
my $spacer_fasta = ""; #"VHDB/InputData/input_spacers.fna"
my $blastn_out_done = ""; #VHDB/BLASTN/input_array_spacer_seq_blast.out
my $from_tabular = 1;
my $outdir = "BED";
my $genome_seq = "";

my $ind = 0;
foreach(@ARGV) {
	if (@ARGV[$ind] eq '-attr') {
		$attr_tbl = @ARGV[$ind + 1];
		if (! (-e $attr_tbl)) {
			die "cannot open file: " . $attr_tbl . "\n";
		}
		$from_tabular = 1;
	}
	
	if (@ARGV[$ind] eq '-fasta') {
		$spacer_fasta = @ARGV[$ind + 1];
		if (! (-e $spacer_fasta)) {
			die "cannot open file: " . $spacer_fasta . "\n";
		}
		$from_tabular = 0;
	}
	
	if (@ARGV[$ind] eq '-blast_out') {
		$blastn_out_done = @ARGV[$ind + 1];
		if (! (-e $blastn_out_done)) {
			die "cannot open file: " .$blastn_out_done . "\n";
		}
	}
	
	if (@ARGV[$ind] eq '-genome') {
		$genome_seq = @ARGV[$ind + 1];
		if (! (-e $genome_seq)) {
			die "cannot open file: " .$genome_seq . "\n";
		}
	}
	
	if (@ARGV[$ind] eq '-out') {
		$outdir = @ARGV[$ind + 1];
	}
	
	$ind++;
}

if (-d $outdir) {
	system("rm -rf $outdir");
}
mkdir($outdir);

my %spacer_attr = ();

if ($from_tabular > 0) {	
	open (ATTR, $attr_tbl);
	my @contents = <ATTR>;
	close(ATTR);
	
	for (my $i = 0; $i < scalar(@contents); $i++) {
		my $line = $contents[$i];
		$line =~ s/\n//g;
		
		my @toks = split(/[\t]/, $line);
		
		my $spacer_id;
		my $spacer_seq;
		
		if (scalar(@toks) > 2) {
			$spacer_id = $toks[0];
			$spacer_seq = $toks[1];
			
			$spacer_attr{$spacer_id}{CONTIG_ID} = $toks[2];
			$spacer_attr{$spacer_id}{SPACER_LENGTH} = $toks[13];
			$spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} = $toks[14];
			$spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} = $toks[15];						
		} elsif (scalar(@toks) == 2) {
			$spacer_id = $toks[0];
			$spacer_seq = $toks[1];
			
			$spacer_attr{$spacer_id}{CONTIG_ID} = $spacer_id;
			$spacer_attr{$spacer_id}{SPACER_LENGTH} = length($spacer_seq);
			$spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} = 50;
			$spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} = 50;			
		} else {
			next;
		}		
	}
	close(FA);
} else {	
	open (SPACER, $spacer_fasta);
	my @contents = <SPACER>;
	close(SPACER);
	
	for (my $i = 0; $i < scalar(@contents); $i += 2) {
		my $spacer_id = $contents[$i];
		$spacer_id =~ s/\n//g;
		$spacer_id =~ s/\>//g;
		
		my $spacer_seq = $contents[$i + 1];
		$spacer_seq =~ s/\n//g;
			
		$spacer_attr{$spacer_id}{CONTIG_ID} = $spacer_id;
		$spacer_attr{$spacer_id}{SPACER_LENGTH} = length($spacer_seq);
		$spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} = 50;
		$spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} = 50;
	}
}

open(BLAST, $blastn_out_done);
my @contents = <BLAST>;
close(BLAST);

my $bed_out = "$outdir/flank_extract.bed";
my $extract_out = "$outdir/flank_extract.fna";

open(BED, ">$bed_out");
for (my $i = 0; $i < scalar(@contents); $i++) {
	my $line = $contents[$i];
	$line =~ s/\n//g;
		
	my @toks = split(/[\t]/, $line);
	my $genome_id = $toks[0];
	my $spacer_id = $toks[4];
	
	my $partial_protospacer_start = $toks[2];
	my $partial_protospacer_end = $toks[3];
	my $protospacer_strand = $toks[7];
	my $genome_length = $toks[15];
	my $protospacer_evalue = $toks[9]; 
	
	my $spacer_length = $spacer_attr{$spacer_id}{SPACER_LENGTH};
	my $source_flank_length_5p = $spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH};
	my $source_flank_length_3p = $spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH};
	
	my $partial_spacer_start_self = $toks[5];
	my $partial_spacer_end_self = $toks[6];
	
	my $ori = "+";
	if ($protospacer_strand eq "minus") {
		$ori = "-";
		if ($partial_protospacer_start > $partial_protospacer_end) {
			$partial_protospacer_start = $partial_protospacer_start + ($partial_spacer_start_self - 1);
			$partial_protospacer_end = $partial_protospacer_start - ($spacer_length - 1);	
		} else {									
			$partial_protospacer_start = $partial_protospacer_start - ($partial_spacer_start_self - 1);
			$partial_protospacer_end = $partial_protospacer_start + ($spacer_length - 1);	
		}	
	} else {
		$partial_protospacer_start = $partial_protospacer_start - ($partial_spacer_start_self - 1);
		$partial_protospacer_end = $partial_protospacer_start + ($spacer_length - 1);	
	}
	
	my $partial_protospacer_start_with_flank = $partial_protospacer_start;
	my $partial_protospacer_end_with_flank = $partial_protospacer_end;
	
	if($partial_protospacer_start_with_flank > $partial_protospacer_end_with_flank) {
		$partial_protospacer_start_with_flank = $partial_protospacer_start_with_flank + $source_flank_length_3p;
		$partial_protospacer_end_with_flank = $partial_protospacer_end_with_flank - $source_flank_length_5p;	
	} else {
		$partial_protospacer_start_with_flank = $partial_protospacer_start_with_flank - $source_flank_length_3p;
		$partial_protospacer_end_with_flank = $partial_protospacer_end_with_flank + $source_flank_length_5p;
	}
	
	my $blastdbcmd_p_start = $partial_protospacer_start_with_flank;
	my $blastdbcmd_p_end = $partial_protospacer_end_with_flank;
	
	if ($blastdbcmd_p_start > $blastdbcmd_p_end) {
		my $temp = $blastdbcmd_p_start;
		$blastdbcmd_p_start = $blastdbcmd_p_end;
		$blastdbcmd_p_end = $temp;
	}
	
	if ($blastdbcmd_p_start < 1) {
		$blastdbcmd_p_start = 1;
	}
	
	if ($blastdbcmd_p_end > $genome_length) {
		$blastdbcmd_p_end = $genome_length
	}
	
	$bed_name = "HIT_" . $i;
	
	$blastdbcmd_p_start = $blastdbcmd_p_start - 1;
	
	print BED "$genome_id\t$blastdbcmd_p_start\t$blastdbcmd_p_end\t$bed_name\t1\t$ori\n";
}
close(BED);

system("$cd_path/bin/bedtools getfasta -fi $genome_seq -bed $bed_out -nameOnly -s > $extract_out");

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}
