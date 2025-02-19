use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $idnumber = "input_array";
my $attr_tbl = ""; #"out/input_array_spacer_attributes.txt";
my $outdir = ""; #"blast_out";
my $spacer_fasta = "";
my $from_tabular = 1;
my $scoring_flank_length_5p = 8;
my $scoring_flank_length_3p = 8;

my $gapopen = 10;
my $gapextend = 2;
my $reward = 1;
my $penalty = -1;
my $word_size = 7;
my $evalue = 1;
my $cpu = 24;
my $dbsize = 10000;
my $db = ""; #"/DB/PHAGE/phage.fa /DB/REFSEQ_PLASMID/plasmid.fa /DB/IMGVR/test.fna";
my $bed_indexed_fna = "NA";

my $wt_psp_match = 1;
my $wt_psp_mismatch = -1;
my $wt_fp5p_match = 0;
my $wt_fp5p_mismatch = 0;
my $wt_fp3p_match = 0;
my $wt_fp3p_mismatch = 0;
my $show_alignment_string = 0;

my $blastn_out_done = "NA";
my $blastdb_acc_lookup = "NA";
my $no_acc_parsing = 0;
my $parse_acc_from_desc = 0;
my $use_actual_dbsize = 0;

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
	
	if (@ARGV[$ind] eq '-out') {
		$outdir = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-gapopen') {
		$gapopen = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-gapextend') {
		$gapextend = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-reward') {
		$reward = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-penalty') {
		$penalty = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-word_size') {
		$word_size = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-evalue') {
		$evalue = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-cpu') {
		$cpu = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-dbsize') {
		$dbsize = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-db') {
		$db = @ARGV[$ind + 1];
		if ($db eq "") {
			die "BLASTDB not specified.";
		}
		$db =~ s/\'//g;
		$db =~ s/\"//g;
	}
	
	if (@ARGV[$ind] eq '-session') {
		$idnumber = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-show_alignment_string') {
		$show_alignment_string = 1;
	}
	
	if (@ARGV[$ind] eq '-blast_out') {
		$blastn_out_done = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-blastdb_acc_lookup') {
		$blastdb_acc_lookup = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-no_acc_parsing') {
		$no_acc_parsing = 1;
	}
	
	if (@ARGV[$ind] eq '-bed_indexed_fna') {
		$bed_indexed_fna = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-parse_acc_from_desc') {
		$parse_acc_from_desc = 1;
	}
	
	if (@ARGV[$ind] eq '-use_actual_dbsize') {
		$use_actual_dbsize = 1;
	}
	
	$ind++;
}

if (-d $outdir) {
	system("rm -rf $outdir");
}
mkdir($outdir);

my $bed_outdir = "$outdir/BED";

my %new_acc = ();
if (-e $blastdb_acc_lookup) {
	open(ACC, $blastdb_acc_lookup);
	my @acc_contents = <ACC>;
	close(ACC);
	
	for (my $i = 0; $i < scalar(@acc_contents); $i++) {
		my $line = $acc_contents[$i];
		$line =~ s/\n//g;
		my ($dbacc, $orig_acc) = ($line =~ /(\S+)\s(\S+).*/);
		
		if (! exists $new_acc{$orig_acc}) {
			$new_acc{$orig_acc} = $dbacc;
		}
	}
}

my %spacer_attr = ();
my $fasta_in = "$outdir/" . $idnumber . "_spacer_sequences.fasta";

#my $pm = new Parallel::ForkManager($cpu);
if ($from_tabular > 0) {	
	open (ATTR, $attr_tbl);
	my @contents = <ATTR>;
	close(ATTR);
	
	open(FA, ">$fasta_in");
	for (my $i = 0; $i < scalar(@contents); $i++) {
		my $line = $contents[$i];
		$line =~ s/\n//g;
		
		my @toks = split(/[\t]/, $line);
		
		my $spacer_id;
		my $spacer_seq;
		
		if (scalar(@toks) > 2) {
			$spacer_id = $toks[0];
			$spacer_seq = $toks[1];
			
			$spacer_attr{$spacer_id}{SPACER_SEQ} = $spacer_seq;
			$spacer_attr{$spacer_id}{CONTIG_ID} = $toks[2];
			$spacer_attr{$spacer_id}{ARRAY_INDEX} = $toks[3];
			$spacer_attr{$spacer_id}{SEQ_START} = $toks[4];
			$spacer_attr{$spacer_id}{SEQ_END} = $toks[5];
			$spacer_attr{$spacer_id}{ARRAY_START} = $toks[6];
			$spacer_attr{$spacer_id}{ARRAY_END} = $toks[7];
			$spacer_attr{$spacer_id}{CRISPR_TYPE} = $toks[8];
			$spacer_attr{$spacer_id}{DIRECT_REPEAT} = $toks[9];
			$spacer_attr{$spacer_id}{SPACER_INDEX} = $toks[10];
			$spacer_attr{$spacer_id}{SPACER_START} = $toks[11];
			$spacer_attr{$spacer_id}{SPACER_END} = $toks[12];
			$spacer_attr{$spacer_id}{SPACER_LENGTH} = $toks[13];
			$spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} = $toks[14];
			$spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} = $toks[15];
			$spacer_attr{$spacer_id}{SRC_FLANK_5P} = $toks[16];
			$spacer_attr{$spacer_id}{SRC_FLANK_3P} = $toks[17];
			$spacer_attr{$spacer_id}{ORIG_CONTIG_ID} = $toks[18];
						
		} elsif (scalar(@toks) == 2) {
			$spacer_id = $toks[0];
			$spacer_seq = $toks[1];
			
			$spacer_attr{$spacer_id}{SPACER_SEQ} = $spacer_seq;
			$spacer_attr{$spacer_id}{CONTIG_ID} = $spacer_id;
			$spacer_attr{$spacer_id}{ARRAY_INDEX} = 1;
			$spacer_attr{$spacer_id}{SEQ_START} = 1;
			$spacer_attr{$spacer_id}{SEQ_END} = length($spacer_seq);
			$spacer_attr{$spacer_id}{ARRAY_START} = 1;
			$spacer_attr{$spacer_id}{ARRAY_END} = length($spacer_seq);
			$spacer_attr{$spacer_id}{CRISPR_TYPE} = "Unknown";
			$spacer_attr{$spacer_id}{DIRECT_REPEAT} = "NA";
			$spacer_attr{$spacer_id}{SPACER_INDEX} = 1;
			$spacer_attr{$spacer_id}{SPACER_START} = 1;
			$spacer_attr{$spacer_id}{SPACER_END} = length($spacer_seq);
			$spacer_attr{$spacer_id}{SPACER_LENGTH} = length($spacer_seq);
			$spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} = 50;
			$spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} = 50;
			$spacer_attr{$spacer_id}{SRC_FLANK_5P} = stringOfDash(50);
			$spacer_attr{$spacer_id}{SRC_FLANK_3P} = stringOfDash(50);
			$spacer_attr{$spacer_id}{ORIG_CONTIG_ID} = $spacer_id;
			
		} else {
			next;
		}
		
		if (length($spacer_attr{$spacer_id}{SRC_FLANK_5P}) < $spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH}) {
			my $diff = abs($spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} - length($spacer_attr{$spacer_id}{SRC_FLANK_5P}));
			$spacer_attr{$spacer_id}{SRC_FLANK_5P} = fillAlnString($spacer_attr{$spacer_id}{SRC_FLANK_5P}, $diff, 1);
		}
		
		if (length($spacer_attr{$spacer_id}{SRC_FLANK_3P}) < $spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH}) {
			my $diff = abs($spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} - length($spacer_attr{$spacer_id}{SRC_FLANK_3P}));
			$spacer_attr{$spacer_id}{SRC_FLANK_3P} = fillAlnString($spacer_attr{$spacer_id}{SRC_FLANK_3P}, $diff, 0);
		}
		
		if ($spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} > $scoring_flank_length_5p) {
			$spacer_attr{$spacer_id}{SCORING_FLANK_5P} = substr($spacer_attr{$spacer_id}{SRC_FLANK_5P}, -$scoring_flank_length_5p, $scoring_flank_length_5p);
		} elsif ($spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} < $scoring_flank_length_5p) {
			my $diff = abs($scoring_flank_length_5p - $spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH});
			$spacer_attr{$spacer_id}{SCORING_FLANK_5P} = fillAlnString($spacer_attr{$spacer_id}{SRC_FLANK_5P}, $diff, 1);
		} else {
			$spacer_attr{$spacer_id}{SCORING_FLANK_5P} = $spacer_attr{$spacer_id}{SRC_FLANK_5P};
		}

		if ($spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} > $scoring_flank_length_3p) {
			$spacer_attr{$spacer_id}{SCORING_FLANK_3P} = substr($spacer_attr{$spacer_id}{SRC_FLANK_3P}, 0, $scoring_flank_length_3p);
		} elsif ($spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} < $scoring_flank_length_3p) {
			my $diff = abs($scoring_flank_length_3p - $spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH});
			$spacer_attr{$spacer_id}{SCORING_FLANK_3P} = fillAlnString($spacer_attr{$spacer_id}{SRC_FLANK_3P}, $diff, 0);
		} else {
			$spacer_attr{$spacer_id}{SCORING_FLANK_3P} = $spacer_attr{$spacer_id}{SRC_FLANK_3P};
		}
		
		print FA ">$spacer_id\n$spacer_seq\n";		
	}
	close(FA);
} else {
	system("cp $spacer_fasta $fasta_in");
	
	open (SPACER, $spacer_fasta);
	my @contents = <SPACER>;
	close(SPACER);
	
	for (my $i = 0; $i < scalar(@contents); $i += 2) {
		my $spacer_id = $contents[$i];
		$spacer_id =~ s/\n//g;
		$spacer_id =~ s/\>//g;
		
		my $spacer_seq = $contents[$i + 1];
		$spacer_seq =~ s/\n//g;
		
		$spacer_attr{$spacer_id}{SPACER_SEQ} = $spacer_seq;		
		$spacer_attr{$spacer_id}{CONTIG_ID} = $spacer_id;
		$spacer_attr{$spacer_id}{ARRAY_INDEX} = 1;
		$spacer_attr{$spacer_id}{SEQ_START} = 1;
		$spacer_attr{$spacer_id}{SEQ_END} = length($spacer_seq);
		$spacer_attr{$spacer_id}{ARRAY_START} = 1;
		$spacer_attr{$spacer_id}{ARRAY_END} = length($spacer_seq);
		$spacer_attr{$spacer_id}{CRISPR_TYPE} = "Unknown";
		$spacer_attr{$spacer_id}{DIRECT_REPEAT} = "NA";
		$spacer_attr{$spacer_id}{SPACER_INDEX} = 1;
		$spacer_attr{$spacer_id}{SPACER_START} = 1;
		$spacer_attr{$spacer_id}{SPACER_END} = length($spacer_seq);
		$spacer_attr{$spacer_id}{SPACER_LENGTH} = length($spacer_seq);
		$spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} = 50;
		$spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} = 50;
		$spacer_attr{$spacer_id}{SRC_FLANK_5P} = stringOfDash(50);
		$spacer_attr{$spacer_id}{SRC_FLANK_3P} = stringOfDash(50);
		$spacer_attr{$spacer_id}{ORIG_CONTIG_ID} = $spacer_id;
		
		if ($spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} > $scoring_flank_length_5p) {
			$spacer_attr{$spacer_id}{SCORING_FLANK_5P} = substr($spacer_attr{$spacer_id}{SRC_FLANK_5P}, -$scoring_flank_length_5p, $scoring_flank_length_5p);
		} elsif ($spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} < $scoring_flank_length_5p) {
			my $diff = abs($scoring_flank_length_5p - $spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH});
			$spacer_attr{$spacer_id}{SCORING_FLANK_5P} = fillAlnString($spacer_attr{$spacer_id}{SRC_FLANK_5P}, $diff, 1);
		} else {
			$spacer_attr{$spacer_id}{SCORING_FLANK_5P} = $spacer_attr{$spacer_id}{SRC_FLANK_5P};
		}

		if ($spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} > $scoring_flank_length_3p) {
			$spacer_attr{$spacer_id}{SCORING_FLANK_3P} = substr($spacer_attr{$spacer_id}{SRC_FLANK_3P}, 0, $scoring_flank_length_3p);
		} elsif ($spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} < $scoring_flank_length_3p) {
			my $diff = abs($scoring_flank_length_3p - $spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH});
			$spacer_attr{$spacer_id}{SCORING_FLANK_3P} = fillAlnString($spacer_attr{$spacer_id}{SRC_FLANK_3P}, $diff, 0);
		} else {
			$spacer_attr{$spacer_id}{SCORING_FLANK_3P} = $spacer_attr{$spacer_id}{SRC_FLANK_3P};
		}
	}
}

my $blast_out = "$outdir/" . $idnumber . "_spacer_seq_blast.out";
my $report = "$outdir/" . $idnumber . "_report.txt";

#if (scalar(keys(%spacer_attr)) > 0) {
#	my $cmd = "blastn -query $fasta_in -db \'$db\' -num_threads $cpu -word_size $word_size -evalue $evalue -gapopen $gapopen -gapextend $gapextend -penalty $penalty -reward $reward -task blastn-short -outfmt \'6 sseqid stitle sstart send qseqid qstart qend sstrand length evalue mismatch gaps sseq qseq qlen slen\' -dbsize $dbsize > $blast_out";
#	system($cmd);
#} else {
#	exit;
#}

if (scalar(keys(%spacer_attr)) > 0) {
	if (-e $blastn_out_done) {
		system("cp $blastn_out_done $blast_out");
	} else {
		my $cmd = "$cd_path/bin/blastn -query $fasta_in -db \'$db\' -num_threads $cpu -word_size $word_size -evalue $evalue -gapopen $gapopen -gapextend $gapextend -penalty $penalty -reward $reward -task blastn-short -outfmt \'6 sseqid stitle sstart send qseqid qstart qend sstrand length evalue mismatch gaps sseq qseq qlen slen\' ";
		if ($use_actual_dbsize > 0) {
			$cmd .= "-dbsize $dbsize ";
		}
		$cmd .= "> $blast_out";
		system($cmd);
	}
} else {
	exit;
}

if (-e "$bed_indexed_fna" and "$bed_indexed_fna.fai") {
	if ($from_tabular > 0) { 
		system("perl $cd_path/CTFlankExtract.pl -attr $attr_tbl -blast_out $blast_out -genome $bed_indexed_fna -out $bed_outdir");
	} else {
		system("perl $cd_path/CTFlankExtract.pl -fasta $spacer_fasta -blast_out $blast_out -genome $bed_indexed_fna -out $bed_outdir");
	}
}

my @extract_contents = ();
my %extract_lookup = ();
if (-e "$bed_outdir/flank_extract.fna") {
	open (EXT, "$bed_outdir/flank_extract.fna");
	@extract_contents = <EXT>;
	close(EXT);
	
	for (my $i = 0; $i < scalar(@extract_contents); $i += 2) {
		my $seq_id = $extract_contents[$i];
		$seq_id =~ s/\n//g;
		$seq_id =~ s/\>//g;
		
		my ($seq_id_trim) = ($seq_id =~ /(HIT_\d+).*/);
		$seq_id = $seq_id_trim;
		
		my $ext_seq = $extract_contents[$i + 1];
		$ext_seq =~ s/\n//g;
		
		$extract_lookup{$seq_id} = $ext_seq;
	}
}

open(BLAST, $blast_out);
my @contents = <BLAST>;
close(BLAST);

open(REPORT, ">$report");
for (my $i = 0; $i < scalar(@contents); $i++) {
	#$pm->start and next;
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
	
	my $spacer_start = $spacer_attr{$spacer_id}{SPACER_START};
	my $spacer_end = $spacer_attr{$spacer_id}{SPACER_END};
	my $spacer_length = $spacer_attr{$spacer_id}{SPACER_LENGTH};
	my $spacer_seq = $spacer_attr{$spacer_id}{SPACER_SEQ};
	my $source_flank_length_5p = $spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH};
	my $source_flank_length_3p = $spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH};
	my $source_flank_5p = $spacer_attr{$spacer_id}{SRC_FLANK_5P};
	my $source_flank_3p = $spacer_attr{$spacer_id}{SRC_FLANK_3P};
	my $scoring_flank_5p = $spacer_attr{$spacer_id}{SCORING_FLANK_5P};
	my $scoring_flank_3p = $spacer_attr{$spacer_id}{SCORING_FLANK_3P};
	my $contig_id = $spacer_attr{$spacer_id}{CONTIG_ID};
	my $crispr_type = $spacer_attr{$spacer_id}{CRISPR_TYPE};
	my $dr_seq = $spacer_attr{$spacer_id}{DIRECT_REPEAT};
	my $spacer_desc = $spacer_attr{$spacer_id}{ORIG_CONTIG_ID};
	
	my $partial_spacer_start_self = $toks[5];
	my $partial_spacer_end_self = $toks[6];
	
	if ($protospacer_strand eq "minus") {
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
	
	my $blastdbcmd_p_start_orig = $blastdbcmd_p_start;
	my $blastdbcmd_p_end_orig = $blastdbcmd_p_end;
	
	if ($blastdbcmd_p_start < 1) {
		$blastdbcmd_p_start = 1;
	}
	
	if ($blastdbcmd_p_end > $genome_length) {
		$blastdbcmd_p_end = $genome_length
	}
	
	my $db_genome_id = $genome_id;
	if ((-e $blastdb_acc_lookup) && (exists $new_acc{$genome_id})) {
		$db_genome_id = $new_acc{$genome_id};
	}
	
	my $protospacer_seq_with_flanks = "";
	my $bed_name = "HIT_" . $i;
	
	if (exists $extract_lookup{$bed_name}) {
		$protospacer_seq_with_flanks = $extract_lookup{$bed_name};
	} else {	
		my @arr_p_seq = `$cd_path/bin/blastdbcmd -entry '$db_genome_id' -db '$db' -line_length 500 -range '$blastdbcmd_p_start-$blastdbcmd_p_end' -strand '$protospacer_strand' -out - >&1 2>/dev/null`;
		if (scalar(@arr_p_seq) < 2) {
			#$pm->finish;
			next;
		}
	
		$protospacer_seq_with_flanks = $arr_p_seq[1]; 
		chomp $protospacer_seq_with_flanks; 
		$protospacer_seq_with_flanks =~ s/\n//g;
	}
	
	if ($blastdbcmd_p_start_orig < 1) {
		my $diff = abs(1 - $blastdbcmd_p_start_orig);
		$protospacer_seq_with_flanks = fillAlnString($protospacer_seq_with_flanks, $diff, 0);
	}
	
	if ($blastdbcmd_p_end_orig > $genome_length) {
		my $diff = abs($blastdbcmd_p_end_orig - $genome_length);
		$protospacer_seq_with_flanks = fillAlnString($protospacer_seq_with_flanks, $diff, 1);
	}
	
	my $actual_protospacer_seq = substr($protospacer_seq_with_flanks, $source_flank_length_3p, -$source_flank_length_5p);
	my $protospacer_flank_3p_seq = substr($protospacer_seq_with_flanks, 0, $source_flank_length_3p);
	my $protospacer_flank_5p_seq = substr($protospacer_seq_with_flanks, -$source_flank_length_5p, $source_flank_length_5p);
	
	my $actual_protospacer_start = $partial_protospacer_start;
	my $actual_protospacer_end = $partial_protospacer_end;
	if ($actual_protospacer_start > $actual_protospacer_end) {
		my $temp = $actual_protospacer_start;
		$actual_protospacer_start = $actual_protospacer_end;
		$actual_protospacer_end = $temp;
	}
	
	#print "$genome_id\t$spacer_id\t$protospacer_strand\t$source_flank_5p\t$spacer_seq\t$source_flank_3p\t$protospacer_flank_5p_seq\t$actual_protospacer_seq\t$protospacer_flank_3p_seq\n";
	$actual_protospacer_seq = reverse($actual_protospacer_seq); $actual_protospacer_seq =~ tr/ACGT/TGCA/;
	$protospacer_flank_3p_seq = reverse($protospacer_flank_3p_seq); $protospacer_flank_3p_seq =~ tr/ACGT/TGCA/;
	$protospacer_flank_5p_seq = reverse($protospacer_flank_5p_seq); $protospacer_flank_5p_seq =~ tr/ACGT/TGCA/;
	
	my $protospacer_flank_3p_seq_scoring = $protospacer_flank_3p_seq;
	if ($scoring_flank_length_3p < $source_flank_length_3p) {
		$protospacer_flank_3p_seq_scoring = substr($protospacer_flank_3p_seq, 0, $scoring_flank_length_3p);
	}
	
	my $protospacer_flank_5p_seq_scoring = $protospacer_flank_5p_seq;
	if ($scoring_flank_length_5p < $source_flank_length_5p) {
		$protospacer_flank_5p_seq_scoring = substr($protospacer_flank_5p_seq, -$scoring_flank_length_5p, $scoring_flank_length_5p);	
	}
	
	my ($protospacer_contig_id) = ($genome_id =~ /(\S+).*/);
	
	if ($no_acc_parsing < 1) {
		my $acc = "";
		if ($protospacer_contig_id =~ /\S+\|\S+\|/) {
			($acc) = ($protospacer_contig_id =~ /\S+\|(\S+)\|/);
			if ($acc ne "") {
				$protospacer_contig_id = $acc;
			}
		} elsif ($protospacer_contig_id =~ /\S+\|\S+/) {
			($acc) = ($protospacer_contig_id =~ /\S+\|(\S+)/);
			if ($acc ne "") {
				$protospacer_contig_id = $acc;
			}
		}
	}
	
	my $protospacer_desc = $toks[1];
	
	if ($parse_acc_from_desc > 0) {
		($protospacer_contig_id) = ($protospacer_desc =~ /(\S+).*/);
	}
	
	my $actual_protospacer_strand = "-";
	if ($protospacer_strand eq "minus") {
		$actual_protospacer_strand = "+";
	}
	
	my $score = 0;
	my $self_match = 0;
	$spacer_seq =~ tr/ACGT/ACGU/;
	$scoring_flank_5p =~ tr/ACGT/ACGU/;
	$scoring_flank_3p =~ tr/ACGT/ACGU/;
	$source_flank_5p =~ tr/ACGT/ACGU/;
	$source_flank_3p =~ tr/ACGT/ACGU/;
	
	my $actual_protospacer_seq_rev = reverse($actual_protospacer_seq);
	my $protospacer_flank_5p_seq_scoring_rev = reverse($protospacer_flank_5p_seq_scoring);
	my $protospacer_flank_3p_seq_scoring_rev = reverse($protospacer_flank_3p_seq_scoring);
	my ($count_match, $count_mismatch, $count_nn, $count_indel, $cmpstr) = compareAlignedSeqInRC($spacer_seq, $actual_protospacer_seq_rev);
	my ($count_match_5p, $count_mismatch_5p, $count_nn_5p, $count_indel_5p, $cmpstr_5p) = compareAlignedSeqInRC($scoring_flank_5p, $protospacer_flank_5p_seq_scoring_rev);
	my ($count_match_3p, $count_mismatch_3p, $count_nn_3p, $count_indel_3p, $cmpstr_3p) = compareAlignedSeqInRC($scoring_flank_3p, $protospacer_flank_3p_seq_scoring_rev);
	
	$score = ($wt_psp_match * $count_match) + ($wt_psp_mismatch * $count_mismatch) + ($wt_psp_mismatch * $count_indel) + 
	         ($wt_fp5p_match * $count_match_5p) + ($wt_fp5p_mismatch * $count_mismatch_5p) + ($wt_fp5p_mismatch * $count_indel_5p) +
	         ($wt_fp3p_match * $count_match_3p) + ($wt_fp3p_mismatch * $count_mismatch_3p) + ($wt_fp3p_mismatch * $count_indel_3p);
	
	my $num_mismatch = $count_mismatch;
	my $count_mismatch =  $count_mismatch + $count_mismatch_5p + $count_mismatch_3p;
	if ($count_mismatch < 1) {
		$self_match = 1;
	}      
	
	my $comment = "NA";
	if ($show_alignment_string > 0) {
		$comment = $cmpstr;
	} else {
		$comment = $protospacer_evalue;
	}
	
	print REPORT "$contig_id\t$spacer_id\t$protospacer_contig_id\t$actual_protospacer_start\t$actual_protospacer_end\t$protospacer_flank_3p_seq_scoring\t$protospacer_flank_5p_seq_scoring\t$score\t$actual_protospacer_strand\tNA\tNA\tNA\tNA\t$self_match\t$spacer_seq\t$actual_protospacer_seq\t$spacer_desc\t$protospacer_desc\tNA\t$scoring_flank_5p\t$scoring_flank_3p\t$crispr_type\t$dr_seq\tNA\tNA\tNA\tNA\t$comment\t$source_flank_5p\t$source_flank_3p\t$protospacer_flank_5p_seq\t$protospacer_flank_3p_seq\n";
	#$pm->finish;		
}
#$pm->wait_all_children;
close(REPORT);

sub fillAlnString {
	my ($aln, $diff, $prepend) = @_;
	
	my $str = "";
	for (my $i = 0; $i < $diff; $i++) {
		$str .= "-";
	}
	
	my $out = $str . $aln;
	if ($prepend < 1) {
		$out = $aln . $str;
	}
	
	return $out;
}

sub stringOfDash {
	my ($num) = @_;
	
	my $str = "";
	for (my $i = 0; $i < $num; $i++) {
		$str .= "-";
	}
	
	return $str;
}

sub compareAlignedSeqInRC {

	my ($seq1, $seq2) = @_;

	my @arr_seq1=split('',$seq1);
	my @arr_seq2=split('',$seq2);
			
	my %hash_of_complementary_bases;
	$hash_of_complementary_bases{'A'}="T";
	$hash_of_complementary_bases{'C'}="G";
	$hash_of_complementary_bases{'G'}="C";
	$hash_of_complementary_bases{'T'}="A";
	$hash_of_complementary_bases{'U'}="A";
	$hash_of_complementary_bases{'N'}="N";

	my $count_indel = 0;
	my $count_match = 0;
	my $count_mismatch = 0;
	my $count_nn = 0;
	my $cmpstr = "";
	for my $i(0..$#arr_seq1){	
		if(not defined $hash_of_complementary_bases{$arr_seq1[$i]} or not defined $hash_of_complementary_bases{$arr_seq2[$i]}) {
			$count_indel++;
			$cmpstr .= ".";
		} elsif($hash_of_complementary_bases{$arr_seq1[$i]} eq $arr_seq2[$i]) {
			if ($arr_seq2[$i] ne "N") {
				$count_match++;
				$cmpstr .= "+";
			} else {
				$count_nn++;
				$cmpstr .= "*";
			}
		} else {
			$count_mismatch++;
			$cmpstr .= "-";
		}
	}

	return ($count_match, $count_mismatch, $count_nn, $count_indel, $cmpstr);
}

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}