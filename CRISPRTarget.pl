use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $spacer_gff_file = ""; #"NC_002737.crispr.gff";
my $outdir = ""; #"out";
my $gapopen = 10;
my $gapextend = 2;
my $reward = 1;
my $penalty = -1;
my $word_size = 7;
my $evalue = 01;
my $cpu = 24;
my $dbsize = 1000000000; # not set on web, db is 100Gb this in 10Gb
my $use_actual_dbsize = 0;
my $db = "";
my $ctrl_db = "";
my $user_db = "";
my $use_user_db = 0;
my $min_score = 20;
my $keep_user_db = 0;
my $repeat_class_lookup = "NA";
my $repeatType_bin = "NA"; #"/usr/local/bin/anaconda2/envs/cctyper/bin";
my $repeatType_db = "NA"; #/usr/local/bin/anaconda2/envs/cctyper/db;
my $pam_search_all = 0;
my $pam_search_super = 0;
my $pam_search_exact = 0;

my $ind = 0;
foreach(@ARGV) {
	if ($ARGV[$ind] eq '-gff') {
		$spacer_gff_file = $ARGV[$ind + 1];
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
	
	if (@ARGV[$ind] eq '-use_actual_dbsize') {
		$use_actual_dbsize = 1;
	}
	
	if (@ARGV[$ind] eq '-db') {
		$db = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-ctrl_db') {
		$ctrl_db = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-user_fasta') {
		$user_db = @ARGV[$ind + 1];
		$use_user_db = 1;
	}
	
	if (@ARGV[$ind] eq '-min_score') {
		$min_score = @ARGV[$ind + 1];
	}
	
	if ($ARGV[$ind] eq '-keep_user_db') {
		$keep_user_db = 1;
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
	
	if (@ARGV[$ind] eq '-pam_search_all') {
		$pam_search_all = 1;
		$pam_search_super = 0;
		$pam_search_exact = 0;
	}
	
	if (@ARGV[$ind] eq '-pam_search_super') {
		$pam_search_all = 0;
		$pam_search_super = 1;
		$pam_search_exact = 0;
	}
	
	if (@ARGV[$ind] eq '-pam_search_exact') {
		$pam_search_all = 0;
		$pam_search_super = 0;
		$pam_search_exact = 1;
	}
	
	$ind++;
}

if ($use_user_db == 0) {
	if ($db eq "") {
		die "BLASTDB not specified.";
	}
	if (! -e $db) {
		die "BLASTDB $db dose not exist.";
	}
	if (! -e "$db.fai") {
		die "BLASTDB need to come with a Bedtools-compatible index file (*.fai).";
	}
	if ($ctrl_db eq "") {
		die "Shuffled BLASTDB not specified.";
	}
	if (! -e $ctrl_db) {
		die "Shuffled BLASTDB $ctrl_db dose not exist.";
	}
	if (! -e "$ctrl_db.fai") {
		die "Shuffled BLASTDB need to come with a Bedtools-compatible index file (*.fai).";
	}
} else {
	if ($user_db eq "") {
		die "User FASTA not specified.";
	}
	if (! -e $user_db) {
		die "User FASTA $user_db dose not exist.";
	}

	if ($keep_user_db == 0) {
		system("rm -rf $cd_path/USER_DB");
		system("rm -rf $cd_path/USER_SHUFFLED_DB");
	}

	@rs = `sh $cd_path/MakeUserDB.sh $user_db $cd_path/USER_DB`;
	my $user_db_path = $rs[0]; 
	chomp $user_db_path;
	$db = $user_db_path;
	
	@rs = `sh $cd_path/MakeUserShuffledDB.sh $user_db $cd_path/USER_SHUFFLED_DB`;
	my $user_shuffled_db_path = $rs[0]; 
	chomp $user_shuffled_db_path;
	$ctrl_db = $user_shuffled_db_path;
}

if (-d $outdir) {
	system("rm -rf $outdir");
}
mkdir($outdir);

my $cmd_1 = "perl $cd_path/CRISPRTarget_cmd.pl " . 
            "-gff $spacer_gff_file " . 
            "-db $db -bed_indexed_fna $db " . 
            "-evalue $evalue " . 
            "-gapopen $gapopen " . 
            "-gapextend $gapextend " . 
            "-reward $reward " . 
            "-penalty $penalty " . 
            "-word_size $word_size " .  
            "-cpu $cpu -nr " . 
            "-no_acc_parsing " . 
            "-parse_acc_from_desc " . 
            "-out $outdir/main " .
            "-repeatType_bin $repeatType_bin " .
            "-repeatType_db $repeatType_db " .
            "-repeat_lookup $repeat_class_lookup " .
            "-autotype";

if ($use_actual_dbsize == 0) {
	$cmd_1 .= "-dbsize $dbsize ";
}

system($cmd_1);

my $cmd_2 = "perl $cd_path/CRISPRTarget_cmd.pl " . 
            "-gff $spacer_gff_file " . 
            "-db $ctrl_db -bed_indexed_fna $ctrl_db " . 
            "-evalue $evalue " . 
            "-gapopen $gapopen " . 
            "-gapextend $gapextend " . 
            "-reward $reward " . 
            "-penalty $penalty " .  
            "-word_size $word_size " . 
            "-cpu $cpu -nr " . 
            "-no_acc_parsing " . 
            "-out $outdir/ctrl ";

if ($use_actual_dbsize == 0) {
	$cmd_2 .= "-dbsize $dbsize ";
}

system($cmd_2);

my $cmd_3 = "perl $cd_path/PAMRescore.pl " . 
            "-tab $outdir/main/input_array_report.pam.txt " . 
            "-tab_neg $outdir/ctrl/input_array_report.pam.txt " .
            "-out $outdir/input_array_report.tab.txt " .
            "-min_score $min_score -compute_pval -guess_url";

if ($pam_search_all > 0) {
	$cmd_3 .= " -pam_search_all";
} elsif ($pam_search_super > 0) {
	$cmd_3 .= " -pam_search_super";
} elsif ($pam_search_exact > 0) {
	$cmd_3 .= " -pam_search_exact";
} else {
	$cmd_3 .= "";
}

system($cmd_3);

system("perl $cd_path/CRISPRTargetGraphics.pl -in $outdir/input_array_report.tab.txt -txt $outdir/input_array_report.txt -include_url -include_pvalue");
system("perl $cd_path/CRISPRTargetGraphics.pl -in $outdir/input_array_report.tab.txt -html $outdir/input_array_report.html -include_url -include_pvalue"); 

sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}     
