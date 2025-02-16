use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $spacer_file = "";
my $spacer_multifasta_file = "";
my $spacer_gff_file = ""; #"NC_002737.crispr.gff";
my $outdir = ""; #"out";
my $known_class = "";
my $spacer_remove_redundancy = 0;
my $user_database_file="";
my $spacer_type="";
my $is_gff = 0;
my $autotype = 0;
my $repeat_class_lookup = "NA";
my $repeatType_bin = "NA"; #"/usr/local/bin/anaconda2/envs/cctyper/bin";
my $repeatType_db = "NA"; #/usr/local/bin/anaconda2/envs/cctyper/db;

my $gapopen = 10;
my $gapextend = 2;
my $reward = 1;
my $penalty = -1;
my $word_size = 7;
my $evalue = 0.1;
my $cpu = 24;
my $dbsize = 1000000000;
my $db = ""; #"USER_DB/phage.fa";
my $algor = "blastn";

my $blastn_out_done = "NA";
my $blastdb_acc_lookup = "NA";
my $bed_indexed_fna = "NA";
my $no_acc_parsing = 0;
my $parse_acc_from_desc = 0;

my $ind = 0;
foreach(@ARGV) {
	if ($ARGV[$ind] eq '-gff') {
		$spacer_gff_file = $ARGV[$ind + 1];
		if (! (-e $spacer_gff_file)) {
			die "cannot open file: " . $spacer_gff_file . "\n";
		}
		$is_gff = 1;
	}

	if ($ARGV[$ind] eq '-array') {
		$spacer_file = $ARGV[$ind + 1];
		if (! (-e $spacer_file)) {
			die "cannot open file: " . $spacer_file . "\n";
		}
		$is_gff = 0;
	}
	
	if ($ARGV[$ind] eq '-fasta') {
		$spacer_multifasta_file = $ARGV[$ind + 1];
		if (! (-e $spacer_multifasta_file)) {
			die "cannot open file: " . $spacer_multifasta_file . "\n";
		}
		$is_gff = 0;
	}
	
	if ($ARGV[$ind] eq '-genome') {
		$user_database_file = $ARGV[$ind + 1];
	}
	
	if ($ARGV[$ind] eq '-nr') {
		$spacer_remove_redundancy = 1;
	}
	
	if ($ARGV[$ind] eq '-out') {
		$outdir = $ARGV[$ind + 1];
	}
	
	if ($ARGV[$ind] eq '-class') {
		$known_class = $ARGV[$ind + 1];
	}
	
	if ($ARGV[$ind] eq '-autotype') {
		$autotype = 1;
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
	
	if (@ARGV[$ind] eq '-algor') {
		$algor = @ARGV[$ind + 1];
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
	
	if ($ARGV[$ind] eq '-repeat_lookup') {
		$repeat_class_lookup =  $ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-bed_indexed_fna') {
		$bed_indexed_fna = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-parse_acc_from_desc') {
		$parse_acc_from_desc = 1;
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

my $preprocess_cmd = "";
if ($is_gff > 0) {
	$preprocess_cmd = "perl $cd_path/CTGFFInput.pl -in $spacer_gff_file -out $outdir/InputData";
} else {
	if ($spacer_file ne "") {
		$preprocess_cmd = "perl $cd_path/CTProcessInput.pl -array $spacer_file -out $outdir/InputData";
	} else {
		$preprocess_cmd = "perl $cd_path/CTProcessInput.pl -fasta $spacer_multifasta_file -out $outdir/InputData";
	}
}

if ($spacer_remove_redundancy > 0) {
	$preprocess_cmd .= " -nr";
}

if ($user_database_file ne "") {
	$preprocess_cmd .= " -db $user_database_file";
}

if ($autotype == 0 and $known_class ne "") {
	$preprocess_cmd .= " -class $known_class";
} elsif ($autotype == 0 and $known_class eq "") {
	$preprocess_cmd .= " -class Unknown";
} elsif ($autotype == 1 and $known_class ne "") {
	$preprocess_cmd .= " -class $known_class";
} else {
	$known_class = "";
}

$preprocess_cmd .= " -repeat_lookup $repeat_class_lookup -repeatType_bin $repeatType_bin -repeatType_db $repeatType_db";

if (! (-e $spacer_multifasta_file)) {
	system("$preprocess_cmd");
} else {
	mkdir("$outdir/InputData");
	system("cp $spacer_multifasta_file $outdir/InputData/input_spacers.fna");
}

my $todo_db = $db;
if ($user_database_file ne "") {
	my @db_files = glob("$outdir/InputData/*.db");
	if (scalar(@db_files) > 0) {
		$todo_db = $db_files[0];
	}
}

my @all_dbs = split(" ", $todo_db);
my $all_db_exist = 1;
for (my $i = 0; $i < scalar(@all_dbs); $i++) {
	if (! (-e $all_dbs[$i])) {
		$all_db_exist *= 0;
	}
}

if ($all_db_exist == 0) {
	print "Some of the BLASTDBs don't exist. Exit.\n";
	exit;
}

my $script_to_use = "$cd_path/TargetPredict.pl";
if (lc($algor) eq "swipe") {
	$script_to_use = "$cd_path/TargetPredictSWIPE.pl";
}

my $action = "-attr $outdir/InputData/input_array_spacer_attributes.txt";
if (-e "$outdir/InputData/input_spacers.fna") {
	$action = "-fasta $outdir/InputData/input_spacers.fna";
}

my $search_cmd = "perl $script_to_use $action -out $outdir/BLASTN -db \'$todo_db\' -gapopen $gapopen -gapextend $gapextend -reward $reward -penalty $penalty -word_size $word_size -evalue $evalue -cpu $cpu -dbsize $dbsize -blast_out $blastn_out_done -blastdb_acc_lookup $blastdb_acc_lookup -bed_indexed_fna $bed_indexed_fna 2> /dev/null";
if ($no_acc_parsing > 0) {
	$search_cmd .= " -no_acc_parsing";
}

if ($parse_acc_from_desc > 0) {
	$search_cmd .= " -parse_acc_from_desc";
}
system("$search_cmd");

my $pam_rescore_cmd = "perl $cd_path/PAMRescore.pl -tab $outdir/BLASTN/input_array_report.txt -out $outdir/BLASTN/input_array_report.pam.txt";
system($pam_rescore_cmd);

system("perl $cd_path/CRISPRTargetGraphics.pl -in $outdir/BLASTN/input_array_report.pam.txt -txt $outdir/report.txt");
system("perl $cd_path/CRISPRTargetGraphics.pl -in $outdir/BLASTN/input_array_report.pam.txt -html $outdir/report.html");
system("cp $outdir/BLASTN/input_array_report.pam.txt $outdir/");

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}
