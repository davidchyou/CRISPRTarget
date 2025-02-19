use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $idnumber = "input_array";
my $gff_in = ""; #"NC_002737.crispr.gff";
my $outdir = ""; #"out";
my $known_class = "";
my $spacer_remove_redundancy = 0;
my $user_database_file="";
my $repeat_class_lookup = "NA";
my $repeatType_bin = "NA"; #"/usr/local/bin/anaconda2/envs/cctyper/bin";
my $repeatType_db = "NA"; #/usr/local/bin/anaconda2/envs/cctyper/db;

my $ind = 0;
foreach(@ARGV) {
	if ($ARGV[$ind] eq '-in') {
		$gff_in = $ARGV[$ind + 1];
		if (! (-e $gff_in)) {
			die "cannot open file: " . $gff_in . "\n";
		}
		$is_array = 1;
	}
	
	if ($ARGV[$ind] eq '-out') {
		$outdir = $ARGV[$ind + 1];
	}
	
	if ($ARGV[$ind] eq '-db') {
		$user_database_file = $ARGV[$ind + 1];
	}
	
	if ($ARGV[$ind] eq '-nr') {
		$spacer_remove_redundancy = 1;
	}
	
	if ($ARGV[$ind] eq '-class') {
		$known_class = $ARGV[$ind + 1];
	}
	
	if ($ARGV[$ind] eq '-session') {
		$idnumber =  $ARGV[$ind + 1];
	}
	
	if ($ARGV[$ind] eq '-repeat_lookup') {
		$repeat_class_lookup =  $ARGV[$ind + 1];
	}
	
	if ($ARGV[$ind] eq '-repeatType_bin') {
		$repeatType_bin = $ARGV[$ind + 1];
	}
	
	if ($ARGV[$ind] eq '-repeatType_db') {
		$repeatType_db = $ARGV[$ind + 1];
	}
	
	$ind++;
}

if (-d $outdir) {
	system("rm -rf $outdir");
}
mkdir($outdir);

if (-e $user_database_file) {
	open (FILE, $user_database_file);
	open (WR,">$outdir/$idnumber" . "_" . "$user_database_file");
	
	my $seq_index=1;
	while ( my $line=<FILE> ) {
		chomp $line; $line=~ s/\r//; $line=~ s/\s+//;

		if($line eq ""){next;}

		if($line=~ /^>/) {
			my $o_id=$line;
			chomp $o_id;$o_id=~ s/\r//; $o_id=~ s/^>//;
			my $new_id= $idnumber."_USER_DB_SEQUENCE_".$seq_index;
			$hash_id_lookup_table{$new_id}=$o_id;
			print WR ">$new_id\n";

			$seq_index++;
		} else{
			print WR "$line\n";
		}
	}
	close(WR);
	close(FILE);
	my $user_db = &create_blast_db_from_user_uploaded_sequence("$idnumber" . "_" . "$user_database_file");
}

open(GFF, $gff_in);
my @contents = <GFF>;
close(GFF);

my $proc_spacer_file = "$outdir/" . $idnumber . "_spacer_sequences.fasta";
my $direct_repeat_file = "$outdir/" . $idnumber . "_direct_repeats_list.txt";
my $direct_repeat_file_typed = "$outdir/" . $idnumber . "_direct_repeats_list.typed.txt";
my $spacer_attr_file = "$outdir/" . $idnumber . "_spacer_attributes.txt";
my %spacer_attr = ();
my %nr_repeat_lst = ();
my %array_repeat_lst = ();

open(DR, ">$direct_repeat_file");
for (my $i = 0; $i < scalar(@contents); $i++) {
	my $line = $contents[$i];
	$line =~ s/\n//g;
	
	my @toks = split(/[\t]/, $line);
	my $genome_id = $toks[0];
	my $feature = $toks[2];
	
	my $start = $toks[3];
	my $end = $toks[4];
	my $strand = $toks[6];
	my $comment = $toks[8];
	
	if ($feature eq "binding_site") {
		my $length = 0;
		my ($array_index, $spacer_index, $spacer_start, $spacer_end, $array_start, $array_end, $seq) = ($comment =~ /ID=CRISPR(\d+)_\S+(\d+)_(\d+)_(\d+);Name=\S+;Parent=CRISPR\d+_(\d+)_(\d+);Note=(\S+);Dbxref=\S+;.*/);
		$array_end = $array_end - 1;
		$length = length($seq);
	
		my $spacer_id = $genome_id . "_" . $array_index . "_" . $spacer_index . "|" . $spacer_start . "_" . $spacer_end;
		my $rep = $array_repeat_lst{$genome_id}{$array_index};
		
		if ($strand eq "-") {
			$seq =~ tr/ACGT/TGCA/;
			$seq = reverse($seq);
		}
		
		$spacer_attr{$spacer_id}{SPACER_SEQ} = $seq;
		$spacer_attr{$spacer_id}{CONTIG_ID} = $genome_id;
		$spacer_attr{$spacer_id}{ARRAY_INDEX} = $array_index;
		$spacer_attr{$spacer_id}{SEQ_START} = $array_start;
		$spacer_attr{$spacer_id}{SEQ_END} = $array_end;
		$spacer_attr{$spacer_id}{ARRAY_START} = $array_start;
		$spacer_attr{$spacer_id}{ARRAY_END} = $array_end;
		$spacer_attr{$spacer_id}{CRISPR_TYPE} = "Unknown";
		$spacer_attr{$spacer_id}{DIRECT_REPEAT} = $rep;
		$spacer_attr{$spacer_id}{SPACER_INDEX} = $spacer_index;
		$spacer_attr{$spacer_id}{SPACER_START} = $spacer_end;
		$spacer_attr{$spacer_id}{SPACER_END} = $toks[12];
		$spacer_attr{$spacer_id}{SPACER_LENGTH} = $length;
		$spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} = 50;
		$spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} = 50;
		
		if ($known_class ne "") {
			$spacer_attr{$spacer_id}{CRISPR_TYPE} = $known_class;
		}
		
		my $end_index = -1;
		my $iter_direction = -1;
		if ($strand eq "-") {
			$end_index = scalar(@contents);
			$iter_direction = 1;
		}
		
		my $upstream_joined = "";
		my $upstream_joined_length = 0;
		for (my $j = $i + $iter_direction; $j != $end_index; $j += $iter_direction) {
			my $line = $contents[$j];
			$line =~ s/\n//g;
			
			my @toks = split(/[\t]/, $line);
			my $feature = $toks[2];
			my $comment = $toks[8];
			my ($seq) = ($comment =~ /ID=\S+;Name=\S+;Parent=\S+;Note=(\S+);Dbxref=\S+;.*/);
			
			if ($strand eq "-") {
				$seq =~ tr/ACGT/TGCA/;
				$seq = reverse($seq);
			}
			
			if ($upstream_joined_length < $spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} and $feature ne "repeat_region") {
				$upstream_joined = $seq . $upstream_joined;
				$upstream_joined_length += length($seq);
			} else {
				last;
			}
		}
		
		my $downstream_joined = "";
		my $downstream_joined_length = 0;
		
		$end_index = scalar(@contents);
		$iter_direction = 1;
		if ($strand eq "-") {
			$end_index = -1;
			$iter_direction = -1;
		}
		
		for (my $j = $i + $iter_direction; $j != $end_index; $j += $iter_direction) {
			my $line = $contents[$j];
			$line =~ s/\n//g;
			
			my @toks = split(/[\t]/, $line);
			my $feature = $toks[2];
			my $comment = $toks[8];
			my ($seq) = ($comment =~ /ID=\S+;Name=\S+;Parent=\S+;Note=(\S+);Dbxref=\S+;.*/);
			
			if ($strand eq "-") {
				$seq =~ tr/ACGT/TGCA/;
				$seq = reverse($seq);
			}
			
			if ($downstream_joined_length < $spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} and $feature ne "repeat_region") {
				$downstream_joined .= $seq;
				$downstream_joined_length += length($seq);
			} else {
				last;
			}
		}
		
		if (length($upstream_joined) > $spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH}) {
			$upstream_joined = substr($upstream_joined, -1 * $spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH}, $spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH});
		} elsif (length($upstream_joined) < $spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH}) {
			my $diff = abs($spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} - length($upstream_joined));
			$upstream_joined = fillAlnString($upstream_joined, $diff, 1);
		}
		
		if (length($downstream_joined) > $spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH}) {
			$downstream_joined = substr($downstream_joined, 0, $spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH});
		} elsif (length($downstream_joined) < $spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH}) {
			my $diff = abs($spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} - length($downstream_joined));
			$downstream_joined = fillAlnString($downstream_joined, $diff, 0);
		}
		
		$spacer_attr{$spacer_id}{SRC_FLANK_5P} = $upstream_joined;
		$spacer_attr{$spacer_id}{SRC_FLANK_3P} = $downstream_joined;
		$spacer_attr{$spacer_id}{ORIG_CONTIG_ID} = $genome_id;
		
	} elsif ($feature eq "repeat_region") {
		my ($array_index, $rep) = ($comment =~ /ID=CRISPR(\d+)_\d+_\d+;.*Note=(\S+);Dbxref=\S+;.*/);		
		$array_repeat_lst{$genome_id}{$array_index} = $rep;
		
		if (! exists $nr_repeat_lst{$rep}) {
			$nr_repeat_lst{$rep} = "Unknown";
			if (! ($rep =~/N/)) {
				print DR "$rep\n";
			}
		}
		
		next;
	} else {
		next;
	}
}
close(DR);

if ($known_class eq "") {
	if (-e $repeat_class_lookup) {
		system("cp $repeat_class_lookup $direct_repeat_file_typed");
	} else {
		if ((-d $repeatType_bin) and (-d $repeatType_db)) {
			my $cmd = "$repeatType_bin/repeatType $direct_repeat_file --db $repeatType_db > $direct_repeat_file_typed";
			system($cmd);
		}
	}
	
	if (-e $direct_repeat_file_typed) {
		open(TYPE, $direct_repeat_file_typed);
		@contents = <TYPE>;
		close(TYPE);

		for (my $i = 0; $i < scalar(@contents); $i++) {
			my $line = $contents[$i];
			$line =~ s/\n//g;
	
			my @toks = split(/[\t]/, $line);
			my $rep = $toks[0];
			my $type = $toks[1];
	
			$nr_repeat_lst{$rep} = $type;
		}
	}
}

my %hash_of_spacer_seq = ();

open(FA, ">$proc_spacer_file");
open(ATTR, ">$spacer_attr_file");
foreach my $spacer_id (keys(%spacer_attr)) {
	my $seq = $spacer_attr{$spacer_id}{SPACER_SEQ};
	
	if ($spacer_remove_redundancy > 0) {
		if (exists $hash_of_spacer_seq{$seq}) {
			next;
		} else {
			$hash_of_spacer_seq{$seq} = 1;
		}
	}
	
	if ($known_class eq "") {
		$spacer_attr{$spacer_id}{CRISPR_TYPE} = $nr_repeat_lst{$spacer_attr{$spacer_id}{DIRECT_REPEAT}};
	}

	print ATTR $spacer_id . "\t" . $seq . "\t" .
  	$spacer_attr{$spacer_id}{CONTIG_ID} . "\t" .
  	$spacer_attr{$spacer_id}{ARRAY_INDEX} . "\t" .
  	$spacer_attr{$spacer_id}{SEQ_START} . "\t" .
  	$spacer_attr{$spacer_id}{SEQ_END} . "\t" .
  	$spacer_attr{$spacer_id}{ARRAY_START} . "\t" .
  	$spacer_attr{$spacer_id}{ARRAY_END} . "\t" .
  	$spacer_attr{$spacer_id}{CRISPR_TYPE} . "\t" .
  	$spacer_attr{$spacer_id}{DIRECT_REPEAT} . "\t" .
  	$spacer_attr{$spacer_id}{SPACER_INDEX} . "\t" .
  	$spacer_attr{$spacer_id}{SPACER_START} . "\t".
  	$spacer_attr{$spacer_id}{SPACER_END} . "\t" .
  	$spacer_attr{$spacer_id}{SPACER_LENGTH} . "\t" .
  	$spacer_attr{$spacer_id}{SRC_FLANK_5P_LENGTH} . "\t" .
  	$spacer_attr{$spacer_id}{SRC_FLANK_3P_LENGTH} . "\t" .
  	$spacer_attr{$spacer_id}{SRC_FLANK_5P} . "\t" .
  	$spacer_attr{$spacer_id}{SRC_FLANK_3P} . "\t" .
  	$spacer_attr{$spacer_id}{ORIG_CONTIG_ID} . "\n";
  	
  	print FA ">$spacer_id\n$seq\n";
}
close(ATTR);
close(FA);

#unlink($direct_repeat_file);
#if ($known_class eq "") {	
#	unlink($direct_repeat_file_typed);
#}

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

sub create_blast_db_from_user_uploaded_sequence() {
	my($user_database_file)=@_;
	my $tmp_db_name=$user_database_file . ".db";
	system("$cd_path/bin/makeblastdb -in $outdir/$user_database_file -parse_seqids -dbtype nucl -out $outdir/$tmp_db_name >/dev/null");
	return($tmp_db_name);
}

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}

