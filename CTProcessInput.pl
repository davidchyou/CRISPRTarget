use File::Spec;
use Cwd 'getcwd';
use Cwd 'abs_path';

my $cd_path = abs_path($0); 
$cd_path = get_path($cd_path);

my $idnumber = "input_array";

my %hash_id_lookup_table;
my $spacer_file = "";
my $spacer_multifasta_file = "";
my $outdir = "";
my $spacer_remove_redundancy = 0;
my $known_class = "";

my $user_database_file="";
my $spacer_type="";
my $user_spacer_file=$idnumber."_user_spacer_file.out";
my $spacer_source_seq_file=$idnumber."_spacer_source_seq_file.fasta";

my ($new_indentifier,$old_indentifier);
my $is_fasta = 0;
my $is_array = 0;
my $repeat_class_lookup = "NA";

my $ind = 0;
foreach(@ARGV) {
	if ($ARGV[$ind] eq '-array') {
		$spacer_file = $ARGV[$ind + 1];
		if (! (-e $spacer_file)) {
			die "cannot open file: " . $spacer_file . "\n";
		}
		$is_array = 1;
	}
	
	if ($ARGV[$ind] eq '-fasta') {
		$spacer_multifasta_file = $ARGV[$ind + 1];
		if (! (-e $spacer_multifasta_file)) {
			die "cannot open file: " . $spacer_multifasta_file . "\n";
		}
		$is_fasta = 1;
	}
	
	if ($ARGV[$ind] eq '-db') {
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
	
	if ($ARGV[$ind] eq '-session') {
		$idnumber =  $ARGV[$ind + 1];
	}
	
	if ($ARGV[$ind] eq '-repeat_lookup') {
		$repeat_class_lookup =  $ARGV[$ind + 1];
	}
	
	$ind++;
}


if (-d $outdir) {
	system("rm -rf $outdir");
}
mkdir($outdir);

if ($is_array > 0) {
	open (FILE, $spacer_file);
	open (WR,">$outdir/$user_spacer_file");
	
	while ( <FILE> ) {
		if(/\.ORGANISM/) {
			s/\.ORGANISM/\nORGANISM/g;
		}
		print WR;
	}
	close(WR);
	close(FILE);
	
	($new_indentifier,$old_indentifier)=&check_identifiers($user_spacer_file,$outdir,$outdir,\%hash_id_lookup_table);
	
} elsif ($is_fasta > 0) {
	open (FILE, $spacer_multifasta_file);
	open (WR,">$outdir/$user_spacer_file");
	
	while ( my $line=<FILE> ) {
		chomp $line; $line=~ s/\r//; $line=~ s/\s+//;

		if($line eq ""){next;}

		if($line=~ /^>/) {
			my $o_id=$line;
			chomp $o_id;$o_id=~ s/\r//; $o_id=~ s/^>//;
			#$o_id=~ s/[^a-zA-Z0-9_.-]/_/;
			my $time_1=time();
			$time_1=$idnumber.$time_1.int (rand 10000);
			my $new_id="CRISPR_SEQUENCE_".$seq_index.$time_1;
			$hash_id_lookup_table{$new_id}=$o_id;
			print WR ">$new_id\n";

			$seq_index++;
		} else {
			print WR "$line\n";
		}
	}
	close(WR);
	close(FILE);
}

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
				#$o_id=~ s/[^a-zA-Z0-9_.-]/_/;
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

my %array_repeat_lst = ();
if (-e $repeat_class_lookup) {
	open(TYPE, $repeat_class_lookup);
	my @rpt_contents = <TYPE>;
	close(TYPE);

	for (my $i = 0; $i < scalar(@rpt_contents); $i++) {
		my $line = $rpt_contents[$i];
		$line =~ s/\n//g;
	
		my @toks = split(/[\t]/, $line);
		my $rep = $toks[0];
		my $type = $toks[1];
		
		if (! exists $array_repeat_lst{$rep}) {
			$array_repeat_lst{$rep} = $type;
		}
	}
}

my $spacer_first_line=`head -1 $outdir/$user_spacer_file >&1`;
if($spacer_first_line=~/^ORGANISM/) {
	$spacer_type="crt";
} elsif($spacer_first_line=~/^pilercr/) {
	$spacer_type="pilercr";
} elsif($spacer_first_line=~ /########################################/) {
	$spacer_type="CRISPRFinder";
} elsif($spacer_first_line=~ /CRISPRDetect/) {
	$spacer_type="CRISPRDetect";
} elsif($spacer_first_line=~ /^>/) {
	$spacer_type="FASTA";
} elsif($spacer_first_line=~ /^{/) {
	$spacer_type="CRISPRCasFinder";
} else{
	$spacer_type="";
}

my $flag;
my %hash_arrays;
if ($spacer_type ne "") {
	$flag = &create_spacer_sequence_file_from_crispr($user_spacer_file,$spacer_source_seq_file,\%hash_id_lookup_table);
	$flag = &extract_spacers($idnumber, $spacer_type, $user_spacer_file, $spacer_remove_redundancy);
	
	my $id_ref_table = "$outdir/" . $idnumber . "_id_and_name_ref_tab.txt";
	open(WR,">$id_ref_table") or print "$!";
	foreach my $new_id(sort keys %hash_id_lookup_table) {
		print WR "$new_id\t$hash_id_lookup_table{$new_id}\n";
	}
	close(WR);
	
	my $dr_file = "$outdir/" . $idnumber . "_spacer_source_seq_file.fasta.dr.fasta";
	my $dr_file_holer = "$outdir/" . $idnumber . "_spacer_source_seq_file.tmp.txt";
	my $dr_typing = "$outdir/" . $idnumber . "_spacer_source_seq_file.fasta.drtype.txt";
	
	my @contents = ();
	
	if (-e $dr_file) {
		open(DR, "$dr_file");
		@contents = <DR>;
		close(DR);
	}
	
	open(TYPE, ">>$dr_typing");
	for (my $i = 0; $i < scalar(@contents); $i += 2) {
		my $id = $contents[$i];
		my $seq = $contents[$i + 1];
		$id =~ s/\n//g;
		$id =~ s/^>//g;
		$seq =~ s/\n//g;
		
		my $rep_type = "N\A";
		
		if ($known_class ne "") {
			$rep_type = $known_class;
		} elsif ($seq =~ /N/) {
			$rep_type = "Unknown";
		} else {
			if ((-e $repeat_class_lookup) && (exists $array_repeat_lst{$seq})) {
				$rep_type = $array_repeat_lst{$seq};
			} else {
				system("echo \'$seq\' >> $dr_file_holer");
							
				my $cmd_get_type = "/usr/local/bin/anaconda2/envs/cctyper/bin/repeatType $dr_file_holer --db /usr/local/bin/anaconda2/envs/cctyper/db \| /usr/bin/awk '{print \$2}'";
				my @rs = `$cmd_get_type`;
							
				if (scalar(@rs) > 0) {
					$rep_type = $rs[0];
					$rep_type =~ s/\n//g;
				}
				unlink("$dr_file_holer");
			}
		}
		
		#print "$id\t$seq\t$rep_type\n";
		print TYPE "$id\t$seq\t$rep_type\n";
		
		my ($contig_id, $array_index, $array_seq_start, $array_seq_end) = ($id =~ /(\S+)\tCRISPR_(\d+)_(\d+)_(\d+)/);
		$hash_arrays{$contig_id}{$array_index}{SEQ_START} = $array_seq_start;
		$hash_arrays{$contig_id}{$array_index}{SEQ_END} = $array_seq_end;
		$hash_arrays{$contig_id}{$array_index}{DIRECT_REPEAT} = $seq;
		$hash_arrays{$contig_id}{$array_index}{CRISPR_TYPE} = $rep_type;
		$hash_arrays{$contig_id}{$array_index}{ORIG_CONTIG_ID} = $hash_id_lookup_table{$contig_id};
	}
	close(TYPE);
	
	my $proc_spacer_file = "$outdir/" . $idnumber . "_spacer_sequences.fasta";
	my $proc_spacer_scd_file = "$outdir/" . $idnumber . "_spacer_source_seq_file.fasta";
	my $spacer_attr_file = "$outdir/" . $idnumber . "_spacer_attributes.txt";
	
	open(SPACER, $proc_spacer_file);
	@contents = <SPACER>;
	close(SPACER);
	
	open(ATTR, ">$spacer_attr_file");
	for (my $i = 0; $i < scalar(@contents); $i += 2) {
		my $id = $contents[$i];
		my $seq = $contents[$i + 1];
		$id =~ s/\n//g;
		$id =~ s/^>//g;
		$seq =~ s/\n//g;
		
		my ($contig_id, $array_index, $array_start, $array_end, $spacer_index, $spacer_start, $spacer_end, $spacer_len) = ($id =~ /(\S+)\|(\d+)_(\d+)_(\d+)\|\d+_(\d+)_(\d+)_(\d+)\|(\d+)/);
		my $array_seq_start = $hash_arrays{$contig_id}{$array_index}{SEQ_START};
		my $array_seq_end = $hash_arrays{$contig_id}{$array_index}{SEQ_END};
		my $dr_seq = $hash_arrays{$contig_id}{$array_index}{DIRECT_REPEAT};
		my $crispr_type = $hash_arrays{$contig_id}{$array_index}{CRISPR_TYPE};
		my $orig_contig_id = $hash_arrays{$contig_id}{$array_index}{ORIG_CONTIG_ID};
		
		my $ext_start = $spacer_start - 50;
		my $ext_end = $spacer_end + 50;
		if ($spacer_start > $spacer_end) {
			$ext_start = $spacer_start + 50;
			$ext_end = $spacer_end - 50;
		}
		
		my ($spacer_seq_with_flanks) = fetch_crispr_seq($contig_id, $array_index, $proc_spacer_scd_file, $ext_start, $ext_end, 50, 50);
		my $flank_seq_3p = substr($spacer_seq_with_flanks, 0, 50);
		my $flank_seq_5p = substr($spacer_seq_with_flanks, -50, 50);
		my $id_new = $contig_id . "_" . $array_index . "_" . $spacer_index . "|" . $spacer_start . "_" . $spacer_end;
		
		print ATTR "$id_new\t$seq\t$contig_id\t$array_index\t$array_seq_start\t$array_seq_end\t$array_start\t$array_end\t$crispr_type\t$dr_seq\t$spacer_index\t$spacer_start\t$spacer_end\t$spacer_len\t50\t50\t$flank_seq_3p\t$flank_seq_5p\t$orig_contig_id\n";
	}
	close(ATTR);
	
} else {
	print "Error: only CRISPRDetect, CRT, pilercr, CRISPRCasFinder and CRISPRfinder format is supported at the moment.";exit;
}

sub create_blast_db_from_user_uploaded_sequence()
	{
		my($user_database_file)=@_;
		my $tmp_db_name=$user_database_file . ".db";
		system("$cd_path/bin/makeblastdb -in $outdir/$user_database_file -parse_seqids -dbtype nucl -out $outdir/$tmp_db_name >/dev/null");
		return($tmp_db_name);
	}

sub create_spacer_sequence_file_from_crispr()
	{

		my ($spacer_file,$spacer_source_seq_file,$hash_id_lookup_table)=@_;

		open(RD,"$outdir/$spacer_file") or print "$! <br>";
		my @arr_spacer_file=<RD>;
		close(RD);


		my $spacer_type; chomp $arr_spacer_file[0];$arr_spacer_file[0]=~ s/^\s+//;


		if($arr_spacer_file[0]=~ /ORGANISM/)
			{
				$spacer_type="crt";
			}
		elsif($arr_spacer_file[0]=~ /pilercr/)
			{
				$spacer_type="pilercr";
			}
		elsif($arr_spacer_file[0]=~ /########################################/)
			{
				$spacer_type="CRISPRFinder";
			}
		elsif($arr_spacer_file[0]=~ /CRISPRDetect/)
			{
				$spacer_type="CRISPRDetect";
			}
		elsif($arr_spacer_file[0]=~ /^>/)
			{
				$spacer_type="FASTA";
			}
		elsif($arr_spacer_file[0]=~ /^{/)
			{
				$spacer_type="CRISPRCasFinder";
			}
		else{
				print "Error: only CRISPRDetect, CRT, pilercr, CRISPRCasFinder and CRISPRfinder format is supported at the moment.";exit;
			}

		my $warnstr = "";
		my $bad_header = 0;
		if ($spacer_type eq "FASTA") {
			open(OVERWRITE, ">$outdir/$spacer_file");
			for (my $i = 0; $i < scalar(@arr_spacer_file); $i++) {
				my $line = $arr_spacer_file[$i];
				$line =~ s/\n//g;
				if ($line =~ /^>\d+/) {
					$bad_header += 1;
					my ($seq_id) = ($line =~ /^>(\d+)/);
					print OVERWRITE ">SPACER_" . $seq_id . "\n";
				} else {
					print OVERWRITE "$line\n";
				}
			}
			close(OVERWRITE);

			open(RD2,"$outdir/$spacer_file") or print "$! <br>";
			@arr_spacer_file=<RD2>;
			close(RD2);
		}

		if ($bad_header > 0) {
			$warnstr = "Warning: Headers with only digits may produce incomplete output files, and are prepended with \"SPACER_.\"";
		}
		
		open(RWR, ">>$outdir/$spacer_source_seq_file.dr.fasta") or print "$!";
		open(WR,">>$outdir/$spacer_source_seq_file") or print "$!";

		#------ process the crt spacer text -----------------------------------------------------
		if($spacer_type eq "crt")
			{
				my $i;
				for($i=0;$i<=$#arr_spacer_file;$i++)
					{
						if($arr_spacer_file[$i]=~ /^ORGANISM:/)     #------ first level
							{
								#---- check the 4th next line ---------------
								if($arr_spacer_file[$i+4] !~ /^CRISPR/){next;}


								my $line=$arr_spacer_file[$i]; chomp $line;

								my @arr_tmp=split('\:',$line);
								my $org_id=$arr_tmp[$#arr_tmp]; $org_id=~ s/^\s+//;$org_id=~ s/\r//g;$org_id=~ s/\s+$//;

								my $j;
								for($j=$i+1;$j<=$#arr_spacer_file;$j++)
									{
										my $line2=$arr_spacer_file[$j]; chomp $line2;

										if($line2=~ /^CRISPR/)     #------ second level:   'CRISPR 1   Range: 2968265 - 2969027\n'
											{
												$line2=~ s/\s+/ /;$line2=~ s/\r//;


												my @arr_tmp_2=split(' ',$line2);

												my $crispr_index=$arr_tmp_2[1];
												my $sequence_start=$arr_tmp_2[3];
												my $sequence_stop=$arr_tmp_2[5]; chomp $sequence_stop; $sequence_stop=~ s/\r$//;
												my $sequence_strand="";

												my $start_recording=0;
												my $k;
												for ($k=$j+1;$k<=$#arr_spacer_file;$k++)
													{
														my $line3=$arr_spacer_file[$k]; chomp $line3;

														if($line3=~/POSITION/)
															{
																$start_recording=1;
																$k=$k+1; # to skip the first '----' line
																next;

															}
														elsif($line3=~ /-----/ and $start_recording==1)
															{
																$start_recording=0; # to stop recording at the second '----' line
																last;
															}

														elsif($start_recording==1)
															{
																$line3=~ s/\r//;$line3=~ s/\t$//;
																my @arr_tmp_3=split('\t',$line3);

																if($#arr_tmp_3<3)
																	{
																		$sequence_strand=$sequence_strand.$arr_tmp_3[$#arr_tmp_3];
																	}
																else{
																		$sequence_strand=$sequence_strand.$arr_tmp_3[$#arr_tmp_3-2].$arr_tmp_3[$#arr_tmp_3-1]; #---- as $arr_tmp_3[1] contains space [error by crt]
																	}


															}
													}
												print WR ">$org_id\tCRISPR_$crispr_index\_$sequence_start\_$sequence_stop\n$sequence_strand\n";
												$j=$k+1;
											}

										if($line2=~ /^ORGANISM:/){$i=$j;last;}
									}
							}
					}
			}

		#------ process the pilercr spacer text
		elsif($spacer_type eq "pilercr")
			{
				for(my $i=0;$i<=$#arr_spacer_file;$i++)
					{
						my $line=$arr_spacer_file[$i]; chomp $line;
						if($line=~ /^>/)     #------ first level------------------------------------------------------------------------------------------------
							{
								my $org_id=$line;		
								if($org_id=~/^>/){$org_id=~ s/^>//;}

								#----- get CRISPR index from third previous line------------------------------------------------------------------
								my $previous_line=$arr_spacer_file[$i-1];chomp $previous_line;$previous_line=~ s/\r//;

								if($previous_line !~ /Array/)
									{
										next; #----- as thats the summary part
									}

								$previous_line=~ s/Array //;
								my $crispr_index=$previous_line;

								#---- get sequence start from fourth next line -------------------------------------------------------------------

								# check if the line containing the first position is 4th or 5th from the current position
								my $line_containing_first_pos=4;
								if($arr_spacer_file[$i+3]=~ /===/)   # thast the first ======= line
									{
										$line_containing_first_pos=4;
									}
								if($arr_spacer_file[$i+4]=~ /===/)
									{
										$line_containing_first_pos=5;
									}


								my $second_next_line=$arr_spacer_file[$i+$line_containing_first_pos];    # even though it seems like i+4 is the correct number, it works with i+5, don't modify
								chomp $second_next_line; $second_next_line=~ s/\s+/ /g; $second_next_line=~ s/^\s+//;$second_next_line=~ s/\r$//;
								my @arr_second_next_line=split(' ',$second_next_line);

								my $sequence_start=$arr_second_next_line[0];
								my $left_flank=	$arr_second_next_line[4];

								#---- now get the DR and right flank
								my $right_flank="";
								my $direct_repeat="";

								my $line_containing_right_flank;
								for (my $l=$i+$line_containing_first_pos;$l<=$#arr_spacer_file;$l++)
									{
										my $tmp_line=$arr_spacer_file[$l];
										if($tmp_line=~ /===/)
											{
												my $tmp_prev_line=$arr_spacer_file[$l-1]; chomp $tmp_prev_line;$tmp_prev_line=~ s/\s+/ /g;$tmp_prev_line=~ s/\r//g;
												my @arr_tmp_line_prev=split(' ',$tmp_prev_line);
												$right_flank=$arr_tmp_line_prev[$#arr_tmp_line_prev];


												my $tmp_next_line=$arr_spacer_file[$l+1]; chomp $tmp_next_line;$tmp_next_line=~ s/\s+/ /g;
												my @arr_tmp_line_next=split(' ',$tmp_next_line);
												$direct_repeat=$arr_tmp_line_next[$#arr_tmp_line_next];


												$line_containing_right_flank=$l;
												last;
											}
									}
								my @arr_direct_repeat=split('',$direct_repeat);
								#------ now extract the sequence
								my $sequence_strand="";	#$left_flank; #---------- not needed as the spacer itself cotains the needed flanking region

								for (my $l=$i+$line_containing_first_pos;$l<=$line_containing_right_flank;$l++)
									{
										my $tmp_line=$arr_spacer_file[$l]; chomp $tmp_line;$tmp_line=~ s/\r//g;
										if($tmp_line=~ /===/)
											{
												my $sequence_stop=$sequence_start+length($sequence_strand);
												
												if ($direct_repeat ne "") {
													print RWR ">$org_id\tCRISPR_$crispr_index\_$sequence_start\_$sequence_stop\n$direct_repeat\n";
												}
												print WR ">$org_id\tCRISPR_$crispr_index\_$sequence_start\_$sequence_stop\n$sequence_strand\n";
												last;
											}
										else{
												my @arr_tmp_line=split('\s+',$tmp_line);
												my $tmp_repeat=$arr_tmp_line[$#arr_tmp_line-1];
												my $tmp_spacer=$arr_tmp_line[$#arr_tmp_line];

												my @arr_tmp_repeat=split('',$tmp_repeat);
												my $tmp_repeat_dot_replaced="";
												my $dr_index=0;
												foreach my $nt(@arr_tmp_repeat)
													{
														if($nt eq ".")
															{
																$tmp_repeat_dot_replaced=$tmp_repeat_dot_replaced.$arr_direct_repeat[$dr_index];
															}
														else{
																$tmp_repeat_dot_replaced=$tmp_repeat_dot_replaced.$nt;
															}

														$dr_index++
													}
												$sequence_strand=$sequence_strand.$tmp_repeat_dot_replaced.$tmp_spacer;
											}
									}

								next;
							}
					}
			}

		#------ process the CRISPRfinder spacer text
		elsif($spacer_type eq "CRISPRFinder")
			{
				for(my $i=0;$i<=$#arr_spacer_file;$i++)
					{
						my $line=$arr_spacer_file[$i]; chomp $line;

						my($org_id,$crispr_index,$sequence_start,$sequence_stop,$sequence_strand);
						my $direct_repeat="";
						if($line=~ /^# Sequence:/)     #------ first level------------------------------------------------------------------------------------------------
							{
								$line=~ s/^# Sequence://;
								$org_id=$line; $org_id=~ s/\r//g; $org_id=~ s/>//;$org_id=~ s/\s+$//; $org_id=~ s/^\s+//;

								#------ now get the description ------
								my $tmp_description_line=$arr_spacer_file[$i+1]; chomp $tmp_description_line;$tmp_description_line=~ s/\s+$//; $tmp_description_line=~ s/# Description: //;
								my $tmp_id_line=$arr_spacer_file[$i+3]; chomp $tmp_id_line;$tmp_id_line=~ s/\s+$//; $tmp_id_line=~ s/# Id: //;

								my $tmp_total_id_description=$tmp_id_line.$tmp_description_line; $tmp_total_id_description=~ s/\s+/ /g;$tmp_total_id_description=~ s/\s+$//;

								$hash_id_lookup_table->{$org_id}=$tmp_total_id_description;
								$crispr_index=$arr_spacer_file[$i+6]; chomp $crispr_index; $crispr_index=~ s/# Crispr Rank in the sequence: //; $crispr_index=~ s/\r//;

								my $crispr_start_stop_line=$arr_spacer_file[$i+7]; chomp $crispr_start_stop_line;$crispr_start_stop_line=~ s/\r//;
								$crispr_start_stop_line=~ s/^#//;
								$crispr_start_stop_line=~ s/Crispr_begin_position://;
								$crispr_start_stop_line=~ s/Crispr_end_position://;
								$crispr_start_stop_line=~ s/\s+/ /g;
								$crispr_start_stop_line=~ s/^\s+//;
								$crispr_start_stop_line=~ s/\t+/ /;
								my @arr_tmp=split(' ',$crispr_start_stop_line);
								$sequence_start=$arr_tmp[0];
								$sequence_stop=$arr_tmp[1];
								my $dr_line=$arr_spacer_file[$i+8]; chomp $dr_line; $dr_line=~ s/\r//; $dr_line=~ s/# DR: //;
								my($dr,$tmp1)=split('\t',$dr_line);

								$sequence_strand=$dr;
								$i=$i+11;
								while($arr_spacer_file[$i] =~ / +([0-9]+)\t +([0-9]+)\t +([ACGT]+)/)
									{
										$sequence_strand=$sequence_strand.$3.$dr;
										#print "\$1=$1 \$2=$2 \$3=$3<br>";
										if($arr_spacer_file[$i+2] =~ /#====/)
											{
												last;
											}
										$i=$i+2;
									}
								if($sequence_strand ne "")
									{	
										if ($direct_repeat ne "") {
											print RWR ">$org_id\tCRISPR_$crispr_index\_$sequence_start\_$sequence_stop\n$direct_repeat\n";
										}
										print WR ">$org_id\tCRISPR_$crispr_index\_$sequence_start\_$sequence_stop\n$sequence_strand\n";
									}
							}
					}
			}
		#------ process the CRISPRDetect spacer text
		elsif($spacer_type eq "CRISPRDetect")
			{
				for(my $i=0;$i<=$#arr_spacer_file;$i++)
					{
						my $line=$arr_spacer_file[$i]; chomp $line; my $array_seq="";my $left_flank=""; my $right_flank="";
						if($line=~ /^>/)     #------ first level--------------------------------------------------------
							{

								#my $k=$i;
								my $array_direction;my $array_start;my $array_stop;my $model_repeat;
								my ($org_id,$direction)=split('\t+',$line);		#$arr_tmp[$#arr_tmp]; $org_id=~ s/^\s+//;$org_id=~ s/\r//g;$org_id=~ s/\s+$//;
								if($org_id=~/^>/){$org_id=~ s/^>//;}

								#----- get CRISPR index from previous line------------------------------------------------------------------
								my $previous_line=$arr_spacer_file[$i-1];chomp $previous_line;$previous_line=~ s/\r//;
								my $crispr_index=$previous_line;
								if($previous_line =~ /Array (\d+) (\d+)-(\d+)/)
									{
										$crispr_index=$1;
										$array_start=$2;
										$array_stop=$3;
									}
								#---- now get the direction ------------------------------------------------------------------------------

								if($line=~/Array_Orientation: (\S+)/)
									{
										$array_direction=$1; chomp $array_direction; $array_direction=~s/\r+//g;
									}
								#--- first get the model repeat,array start and stop -----------------------------------------------------
								my $j=$i+4;
								while($arr_spacer_file[$j]!~/\/\//)
									{
										if($arr_spacer_file[$j]!~/\S+/){$j++;next;}
										my $current_line=$arr_spacer_file[$j];chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;

										#------- get model repeat-------------------------------------
										if($arr_spacer_file[$j+1]=~/====/)
											{
												my $mr_line=$arr_spacer_file[$j+2];
												$mr_line=~s/^\s+//;
												$mr_line=~s/\s+/\t/g;
												my @arr_t1=split('\t',$mr_line);

												$model_repeat=$arr_t1[4];
											}
										#------- get left and right flanks --------------------------
										if($current_line=~/Left flank :   (\S+)/){$left_flank=$1;if($left_flank=~/^\|/){$left_flank="";}}
										if($current_line=~/Right flank :  (\S+)/){$right_flank=$1; if($right_flank=~/^\|/){$right_flank="";}}
										$j++;
									}
								#---now convert all the dotted repeat, and add to sequence_strand on the go
								my $k=$i+4;
								while($arr_spacer_file[$k]!~/====/)
									{
										if($arr_spacer_file[$k]!~/\S+/){$k++;next;}
										my $current_line=$arr_spacer_file[$k]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;
										my @arr_t1=split('\t',$current_line);
										my $r_start=$arr_t1[0];$r_start=~s/\s+//g;
										my $r_length=$arr_t1[1];$r_length=~s/\s+//g;
										my $s_length=$arr_t1[3];$s_length=~s/\s+//g;
										my $r_seq=$arr_t1[4]; $r_seq=~s/\s+//g; $r_seq=&change_dots_to_bases($r_seq,$model_repeat);
										my $s_seq=$arr_t1[5]; $s_seq=~s/\s+//g;
										my $comment=$arr_t1[6];		$s_length=~s/^\s+//;

										#---- dont forget to check insertion bases in comment ---------------------------
										my %hash_of_insertion_positions;
										my $no_of_insertions=0;
										if($arr_t1[6] and $arr_t1[6]!~/^Del/)
											{
												$comment=$arr_t1[6]; chomp $comment; $comment=~s/^\s+//;
												my @tmp_arr1=split(' ',$comment);
												my $insertion_bases=$tmp_arr1[0];
												my $insertion_positions=$tmp_arr1[1];

												if($insertion_positions)
													{
													   $insertion_positions=~s/\[//g;  $insertion_positions=~s/\]//g;
														my @tmp_arr2=split(',',$insertion_bases);
														my @tmp_arr3=split(',',$insertion_positions);

														for( my $p=0;$p<=$#tmp_arr3;$p++)
															{
																my $pos=$tmp_arr3[$p];
																$hash_of_insertion_positions{$pos}=$tmp_arr2[$p];
															}

														$insertion_bases=~s/,//g;
														$no_of_insertions=length($insertion_bases);
													}
											}
										#--------------------------------------------------------------------------------
										$r_seq=~s/-//g;
										$s_seq=~s/-//g;$s_seq=~s/\|//g;

										#------- split the repeat seq in an array and put the insertion bases back to their positions
										if($no_of_insertions>0)
											{
												my @arr_r_seq=split('',$r_seq);
												my %hash_of_repeat_bases;
												my $index=0;
												foreach my $nt(@arr_r_seq)
													{
														my $current_base_position;

														if($array_direction!~/^Rev/)
															{
																$current_base_position=$r_start+$index;
															}
														else{
																$current_base_position=$r_start-$index;
															}


														if($hash_of_insertion_positions{$current_base_position})
															{
																$hash_of_repeat_bases{$current_base_position}=$hash_of_insertion_positions{$current_base_position}.$nt;
															}
														else{
																$hash_of_repeat_bases{$current_base_position}=$nt;
															}
														$index++;
													}
												#--------- now join the bases to get the corrected repeat --------------------------
												$r_seq="";

												if($array_direction!~/^Rev/)
													{
														foreach my $pos(sort{$a<=>$b} keys %hash_of_repeat_bases)
															{
																$r_seq=$r_seq.$hash_of_repeat_bases{$pos};
															}
													}
												else{
														foreach my $pos(sort{$b<=>$a} keys %hash_of_repeat_bases)
															{
																$r_seq=$r_seq.$hash_of_repeat_bases{$pos};
															}
													}
											}
										#--------------------------------------------------------------------------------------------

										$array_seq=$array_seq.$r_seq.$s_seq;

										$k++;
									}
								#---------------- now add the flanks to the ends----------------------------------------
								my $sequence_strand=$left_flank.$array_seq.$right_flank;

								#------------ at the end write the source sequence line for this array -----------------
								#print "Seq:$array_seq-";exit;
								my $sequence_start;my $sequence_stop;

								if($array_direction=~/^Rev/)
									{
										$sequence_start=$array_start+length($left_flank);
										$sequence_stop=$array_stop-length($right_flank);
									}
								else{
										$sequence_start=$array_start-length($left_flank);
										$sequence_stop=$array_stop+length($right_flank);
									}
								if ($model_repeat ne "") {
									print RWR ">$org_id\tCRISPR_$crispr_index\_$sequence_start\_$sequence_stop\n$model_repeat\n";
								}
								#$sequence_stop=0; #---- $sequence_stop should be set to Zero to be identified that it's a single sequence file, not multi fasta
								print WR ">$org_id\tCRISPR_$crispr_index\_$sequence_start\_$sequence_stop\n$sequence_strand\n";
								$i=$j;
							}
					}
			}
		#------ process the multiFASTA spacer text
		elsif($spacer_type eq "FASTA")
			{
				my($org_id,$crispr_index,$sequence_start,$sequence_stop,$sequence_strand);
				   $crispr_index=0;

				for(my $i=0;$i<=$#arr_spacer_file;$i++)
					{
						my $line=$arr_spacer_file[$i];
						chomp $line; $line=~ s/\r//;	$line=~ s/^\s+//;

						if($line eq ""){next;}

						if($line=~ /^>/)     #------ first level------------------------------------------------------------------------------------------------
							{
								if($i!=0)
									{
										$crispr_index++;
										$sequence_start=1;
										$sequence_stop=length($sequence_strand);
										print WR ">$org_id\tCRISPR_$crispr_index\_$sequence_start\_$sequence_stop\n$sequence_strand\n";
									}

								#----- store the seq id in org_id -----
								#$first_l++;

								$line=~s/\|/_/g;
								$line=~s/\s+/_/g;
								$line=~s/-/_/g;
								$line=~s/\t+/_/g;
								$line=~s/,/_/g;
								$line=~s/_+/_/g;
								$line=~s/:/_/g;
								$line=~s/_+$/_/g;

								$org_id=$line;
								$org_id=~ s/>//; #$org_id=~ s/[^a-zA-Z0-9_.-]/_/; $org_id=~ s/_$//;
								$sequence_strand="";
							}
						else{
								$sequence_strand=$sequence_strand.$line;
							}
					}
				#---- finally write the last record -----
				$crispr_index++;
				$sequence_start=1;
				$sequence_stop=length($sequence_strand);
				print WR ">$org_id\tCRISPR_$crispr_index\_$sequence_start\_$sequence_stop\n$sequence_strand\n";
				$sequence_strand="";

				if (length($warnstr) > 0) {
					print "<font color=red>$warnstr</font><br>";
				}

			}
		elsif($spacer_type eq "CRISPRCasFinder")
			{
				my $path = "$outdir/$spacer_file";
				my $json;

				{
   					local $/;
    				open my $fh, '<', $path or die $!;
   					$json = <$fh>;
				}

				my $perl = decode_json $json;

				my $command_used = $perl->{Command};
				my $gcf <- "NA";

				if (index($command_used, "#GC") >= 0) {
					($gcf) = ($command_used =~ /.*\#(GC\S+)\s+.*/);
				}

				my @seqs = @{$perl->{Sequences}};

				foreach my $seq (@seqs) {
					my $seq_id = $seq->{Id};
					my $org = $seq->{Description};

					$org =~ s/\,//g;
					$org =~ s/\://g;
					$org =~ s/\;//g;
					$org =~ s/\s/_/g;
					$org =~ s/\-/_/g;
					$org =~ s/\./_/g;
					$org =~ s/\#/_/g;
					$org =~ s/\(//g;
					$org =~ s/\)//g;
					$org =~ s/\[//g;
					$org =~ s/\]//g;
					$org =~ s/\%/_/g;

					my @crisprs = @{$seq->{Crisprs}};
					my $count = 0;
					foreach my $crispr (@crisprs) {
						$count++;
						my $org_id = $seq_id . "_" . $org;
						my $arr_start = $crispr->{Start};
						my $arr_end = $crispr->{End};
						my $direct_repeat = $crispr->{DR_Consensus};
						my $crispr_index = $count;

						my @regions = @{$crispr->{Regions}};

						my $region_seq = "";
						foreach my $region (@regions) {

							my $type = $region->{Type};

							if ($type eq "LeftFLANK") {
								next;
							}

							if ($type eq "RightFLANK") {
								next;
							}


							my $region_seq_elem = $region->{Sequence};
							$region_seq .= $region_seq_elem;
						}
						if ($direct_repeat ne "") {
							print RWR ">$org_id\tCRISPR_$crispr_index\_$arr_start\_$arr_end\n$direct_repeat\n";
						}
						print WR ">$org_id\tCRISPR\_$crispr_index\_$arr_start\_$arr_end\n$region_seq\n";
					}
				}
			}

		close(WR);
		close(RWR);

		return 1;
	}

sub change_dots_to_bases()
	{
		my($r_string,$mr_string)=@_;

		if($r_string eq "" or $mr_string eq "" ){return("");}
		my $return_string;

		my @arr_1=split('',$r_string);
		my @arr_2=split('',$mr_string);

		for(my $i=0;$i<=$#arr_2;$i++)
			{
				if($arr_1[$i] eq ".")
					{
						$return_string=$return_string.$arr_2[$i];
					}
				else{
						$return_string=$return_string.$arr_1[$i];
					}
			}

		return($return_string);
	}



sub check_identifiers()
	{
		my($file,$file_path,$dest_folder,$hash_id_lookup_table)=@_;

		open(RD,"$file_path/$file");
		my @array_rd = <RD>;
		close(RD);

		my $old_indentifier="";
		my $new_indentifier="";

		my %hash_of_old_ids;
		#----- check for FASTA seq file -----------------
		my $second_part="";

		if($array_rd[0]=~ />/)
			{
				foreach my $line(@array_rd)
					{
						if($line=~ /^>/)
							{
								#print "\$line=$line<br>";
								my $tmp_old_id;


								if($line=~ /\t/)
									{
										my @arr_tmp_1=split('\t',$line);
										$tmp_old_id=$arr_tmp_1[0];
										$second_part="\t".$arr_tmp_1[1]; chomp $second_part; $second_part=~ s/\r+//g;
									}
								else{
										$tmp_old_id=$line; chomp $tmp_old_id; $tmp_old_id=~ s/\r+//g; $tmp_old_id=~ s/^\s+//; $tmp_old_id=~ s/\s+$//;

									}

								#-------------------------------------------------------------------------------
								$tmp_old_id=~ s/^>//;
								my ($tmp_new_id)=&process_identifier($tmp_old_id);

								$hash_of_old_ids{$tmp_old_id}=$tmp_new_id;
							}
					}

				#------ now create a copy of the file in tmp/ with just the new_id
				open(WR,">$dest_folder/$file");
				foreach my $line(@array_rd)
					{
						if($line=~/>/)
							{
								foreach my $tmp_old_id(sort keys %hash_of_old_ids)
									{
										if($line=~ /$tmp_old_id/)
											{
												my $tmp_new_id=$hash_of_old_ids{$tmp_old_id};

												$line=~ s/\Q$tmp_old_id/$tmp_new_id/i;   #---- the leading \Q quotes all special and reserved charecters like |,'," $ etc in the string

												chomp $line;
												$line=~ s/\r$//; #just to make sure \n wasn't removed
												$line=$line."\n";
												last;
											}
									}
							}

						print WR $line;
					}
				close(WR);
			}
		else{
				#--------------- process the spacer file ---------------------------------------
				if($array_rd[0]=~ /ORGANISM:/)
					{
						#print "---Found CRT spacer file<br>";

						foreach my $line(@array_rd)
							{
								if($line=~ /ORGANISM:/)
									{
										chomp $line; $line=~ s/\r+//; $line=~ s/^\s+//; $line=~ s/\s+$//;

										my($first_part,$tmp_old_id)=split(':  ',$line);

										$tmp_old_id=~ s/^\s+//; $tmp_old_id=~ s/\r+//; chomp $tmp_old_id;
										my $tmp_new_id=&process_identifier($tmp_old_id);

										$tmp_new_id=$tmp_new_id;
										$hash_of_old_ids{$tmp_old_id}=$tmp_new_id;
									}
							}


						open(WR,">$dest_folder/$file");
						foreach my $line(@array_rd)
							{
								if($line=~/ORGANISM:/)
									{
										foreach my $tmp_old_id(sort keys %hash_of_old_ids)
											{
												my $tmp_new_id=$hash_of_old_ids{$tmp_old_id};
												if($line=~ /$tmp_old_id/)
													{
														$line=~ s/\Q$tmp_old_id/$tmp_new_id/i;   #---- the leading \Q quotes all special and reserved charecters like |,'," $ etc in the string

														chomp $line;
														$line=~ s/\r$//; #just to make sure \n wasn't removed
														$line=$line."\n";

														last;
													}

											}

									}

								print WR $line;
							}

						close(WR);
					}
				if($array_rd[0]=~ /pilercr/)
					{
						open(WR,">$dest_folder/$file");
						foreach my $line(@array_rd)
							{
								if($line=~ /^>/)
									{
										chomp $line; $line=~ s/\r+//; $line=~ s/^\s+//; $line=~ s/\s+$//; $line=~ s/^>//;
										my $tmp_old_id=$line;
										$tmp_old_id=~ s/^\s+//; $tmp_old_id=~ s/\r+//; chomp $tmp_old_id;

										my $tmp_new_id=&process_identifier($tmp_old_id);

										$tmp_new_id=$tmp_new_id;
										$hash_of_old_ids{$tmp_old_id}=$tmp_new_id;
										if($tmp_new_id=~/^>/){$tmp_new_id=~ s/^>//;}
										$line=">".$tmp_new_id."\n";
									}
								print WR $line;
							}

						close(WR);
					}
				if($array_rd[0]=~ /CRISPRDetect/)
					{
						open(WR,">$dest_folder/$file");
						foreach my $line(@array_rd)
							{
								if($line=~ /^>/)
									{
										chomp $line; $line=~ s/\r+//; $line=~ s/^\s+//; $line=~ s/\s+$//; $line=~ s/^>//;
										my ($tmp_old_id,$array_orientation)=split('\t+',$line);
										$tmp_old_id=~ s/^\s+//; $tmp_old_id=~ s/\r+//; chomp $tmp_old_id;

										my $tmp_new_id=&process_identifier($tmp_old_id);
										$tmp_new_id=$tmp_new_id;
										$hash_of_old_ids{$tmp_old_id}=$tmp_new_id;
										if($tmp_new_id=~/^>/){$tmp_new_id=~ s/^>//;}
										$line=">".$tmp_new_id."\t\t$array_orientation\n";
									}
								print WR $line;
							}

						close(WR);
					}

				 if($array_rd[0]=~ /{/)
					{
						my $path = "$file_path/$file";
						my $json;

						{
   							local $/;
    						open my $fh, '<', $path or die $!;
   							$json = <$fh>;
						}

						my $perl = decode_json $json;

						my $command_used = $perl->{Command};
						my $gcf <- "NA";

						if (index($command_used, "#GC") >= 0) {
							($gcf) = ($command_used =~ /.*\#(GC\S+)\s+.*/);
						}

						my @seqs = @{$perl->{Sequences}};

						foreach my $seq (@seqs) {
							my $seq_id = $seq->{Id};
							my $org = $seq->{Description};

							$org =~ s/\,//g;
							$org =~ s/\://g;
							$org =~ s/\;//g;
							$org =~ s/\s/_/g;
							$org =~ s/\-/_/g;
							$org =~ s/\./_/g;
							$org =~ s/\#/_/g;
							$org =~ s/\(//g;
							$org =~ s/\)//g;
							$org =~ s/\[//g;
							$org =~ s/\]//g;
							$org =~ s/\%/_/g;

							my @crisprs = @{$seq->{Crisprs}};
							my $count = 0;
							foreach my $crispr (@crisprs) {
								$count++;
								my $org_id = $seq_id . "_" . $org;
								my $tmp_new_id=&process_identifier($org_id);
								$hash_of_old_ids{$org_id}=$tmp_new_id;
							}
						}
					}
			}

		#----- now copy the new and old IDs in lookup table
		foreach my $tmp_old_id(sort keys %hash_of_old_ids)
			{
				my $tmp_new_id=$hash_of_old_ids{$tmp_old_id};
				if(not defined $hash_id_lookup_table->{$tmp_new_id})
					{
						$hash_id_lookup_table->{$tmp_new_id}=$tmp_old_id;
					}
			}

		return($new_indentifier,$old_indentifier);
	}


sub process_identifier()
	{
		my $old_id=shift(@_);

		my $new_id="";


		#--- added by Ambarish on :15/02/2015

		$new_id=$old_id;
		$new_id=~s/\|/_/g;
		$new_id=~s/\s+/_/g;
		$new_id=~s/-/_/g;

		return ($new_id);
		if($old_id=~ /gi\|/)
			{
				if($old_id=~ /\|/)
					{
						my @arr_first_line=split('\|',$old_id);
						$new_id=$arr_first_line[3];
					}
			}
		elsif($old_id=~ /\|/)
			{
				my @arr_first_line=split('\|',$old_id);
				$new_id=$arr_first_line[0];


			}
		elsif($old_id=~ /-/)
			{
				my @arr_first_line=split('-',$old_id);
				$new_id=$arr_first_line[0];


			}
		elsif(length($old_id)>25)
			{
				#------ check if there is any special charecter or spaces, remove with -
				$new_id=substr($old_id,1,25);
				$new_id=~ s/\'//;$new_id=~ s/-//;$new_id=~ s/\"//;$new_id=~ s/\@//;$new_id=~ s/\*//;$new_id=~ s/\///;

			}
		else{
				$new_id=$old_id;
			}


		$new_id=~ s/>//;$new_id=~ s/\r//;
		if($new_id=~ /\./)
			{
				my @arr_tmp=split('\.',$new_id);
				$new_id=$arr_tmp[0];
			}
		return($new_id);
	}
	
sub extract_spacers()
	{
		my($idnumber,$spacer_type,$spacer_file,$spacer_remove_redundancy)=@_;

		my $spacer_seq_fasta_file="$outdir/".$idnumber."_spacer_sequences.fasta";
		open(RD,"$outdir/$spacer_file");
		my @arr_spacer_file=<RD>;
		close(RD);
			
		#------------------------------------------------------------------------------------------
		#------ process the crt spacer text -----------------------------------------------------	
		if($spacer_type eq "crt")
			{		
				open(WR,">$spacer_seq_fasta_file");				
				my $orgName;my $crispr_index;my $rangeStart;my $rangeEnd;
				
				for(my $i=0;$i<scalar(@arr_spacer_file);$i++)
					{						
						if($arr_spacer_file[$i]=~/ORGANISM:  (.*)/)
							{								
								$orgName = $1; chomp $orgName; $orgName=~ s/\r//g;								
							}											
						elsif($arr_spacer_file[$i]=~/CRISPR ([0-9]+).*Range: ([0-9]+) - ([0-9]+)/)
							{								
								$crispr_index=$1; chomp $crispr_index;
								$rangeStart = $2; chomp $rangeStart;
								$rangeEnd = $3; chomp $rangeEnd;
								
								$i+=3;
								my $spacerIndex=0;
								while($arr_spacer_file[$i]=~/([0-9]+)\t\t([ACTG]+)\t([ACTG]+)\t\[ ([0-9]+), ([0-9]+) ]/)
									{									
										if(length($3)<=75)
										{
											$spacerIndex++;
											print WR ">".$orgName."|".$crispr_index."_".$rangeStart."_".$rangeEnd;
											print WR "|".$crispr_index."_".$spacerIndex."_".(int($1)+int($4));
											print WR "_".(int($1)+int($4)+int($5)-1)."|".$5."\n";
											print WR $3."\n\n"; 
										}										
										$i++;																			    
									}
							}		
					}					
				close(WR);
			}
		
		#------ process the pilercr spacer text
		elsif($spacer_type eq "pilercr")
			{		
				open(WR, ">$spacer_seq_fasta_file"); #print "Going to open $pilerCRspacer\n\n";
								
				my $orgName;my $crispr_index;my $rangeStart;my $rangeEnd;

				for(my $i=0;$i<scalar(@arr_spacer_file);$i++)
					{						
						if($arr_spacer_file[$i]=~ />(.*)/)
							{								
								$orgName = $1; chomp $orgName; $orgName=~ s/\r//g;
								
								my $previous_line=$arr_spacer_file[$i-1];chomp $previous_line;$previous_line=~ s/\r//;$previous_line=~ s/\s+$//;
								
								if($previous_line !~ /Array/)
									{										
										next; #----- as thats the summary part
									}	
								$previous_line=~ s/Array //;
								$crispr_index=$previous_line;

								# check if the line containing the first position is 4th or 5th from the current position
								my $line_containing_first_pos=4;
								if($arr_spacer_file[$i+4]=~ /===/)   # thast the first ======= line
									{
										$line_containing_first_pos=4; 
									}
								elsif($arr_spacer_file[$i+5]=~ /===/)
									{
										$line_containing_first_pos=5; 
									}	
																
								my $i=$i+$line_containing_first_pos;
								
								my $spacerIndex=0;
								while($arr_spacer_file[$i]=~/ +([0-9]+) +([0-9]+) +[0-9\.]+ +([0-9]+) +[ACGT]+ +[-ACGT\.]+ + ([ACGT]+)/)
									{										
										$rangeStart = $1; chomp $rangeStart;
										$rangeEnd=$rangeStart+$2+$3-1;
										
										if(length($4)<=75)
										{
											$spacerIndex++;
											print WR ">".$orgName."|".$crispr_index."_".$rangeStart."_".$rangeEnd;
											print WR "|".$crispr_index."_".$spacerIndex."_".(int($1) + int($2));
											print WR "_".(int($1)+int($2)+int($3)-1)."|".$3;
											print WR "\n".$4."\n\n";
										
										if($arr_spacer_file[$i+2]=~ /==/) #---- j+2 as the next line contains the right flank and not a spacer
											{
												last;
											}  # the last line (just before the  ==== line) contains right flank, not spacer 
										}					
										$i++;									
									}
							}
					}
					close(WR);
			}

		#------ process the CRISPRfinder spacer text
		elsif($spacer_type eq "CRISPRFinder")
			{
				open(WR, ">$spacer_seq_fasta_file");
				
				for(my $i=0;$i<=$#arr_spacer_file;$i++)
					{
						my $line=$arr_spacer_file[$i]; chomp $line;						
						if($line=~ /^# Sequence:/)     #------ first level------------------------------------------------------------------------------------------------
							{
								$line=~ s/^# Sequence://;
								
								my($orgName,$crispr_index,$rangeStart,$rangeEnd,$sequence_strand);
								my $direct_repeat="";
								
								$orgName=$line; $orgName=~ s/\r//g; $orgName=~ s/>//;$orgName=~ s/\s+$//; $orgName=~ s/^\s+//;
								
								$crispr_index=$arr_spacer_file[$i+6]; chomp $crispr_index; $crispr_index=~ s/# Crispr Rank in the sequence: //; $crispr_index=~ s/\r//;
								
								my $crispr_start_stop_line=$arr_spacer_file[$i+7]; chomp $crispr_start_stop_line;$crispr_start_stop_line=~ s/\r//;
								$crispr_start_stop_line=~ s/^#//;
								$crispr_start_stop_line=~ s/ Crispr_begin_position: //;
								$crispr_start_stop_line=~ s/ Crispr_end_position: //;
								
								my @arr_tmp=split('\t',$crispr_start_stop_line);
								$rangeStart=$arr_tmp[0];
								$rangeEnd=$arr_tmp[1];
								
								my $dr_line=$arr_spacer_file[$i+8]; chomp $dr_line; $dr_line=~ s/\r//; $dr_line=~ s/# DR: //;
								my($dr,$tmp1)=split('\t',$dr_line); 
								
								$sequence_strand=$dr;
								$i=$i+11;
								my $spacerIndex=0;
								while($arr_spacer_file[$i] =~ / +([0-9]+)\t +([0-9]+)\t +([ACGT]+)/)
									{
										if(length($3)<=75)
											{
												$spacerIndex++;
												print WR ">".$orgName."|".$crispr_index."_".$rangeStart."_".$rangeEnd;
												print WR "|".$crispr_index."_".$spacerIndex."_".$1;
												print WR "_".(int($1)+int($2)-1)."|".$2."\n";
												print WR $3."\n\n"; 
											}
										if($arr_spacer_file[$i+2] =~ /#====/)
											{
												last;
											}
										$i=$i+2;	
									}	
							}
					}	
				close(WR);
			}
		#------ process the CRISPRDetect spacer text
		elsif($spacer_type eq "CRISPRDetect")
			{
				open(WR, ">$spacer_seq_fasta_file");
				
				for(my $i=0;$i<=$#arr_spacer_file;$i++)
					{
						my $line=$arr_spacer_file[$i]; chomp $line; 						
						
						my($orgName,$crispr_index,$rangeStart,$rangeEnd,$sequence_strand);
						if($line=~ /^>/)     #------ first level--------------------------------------------------------
							{
								my $array_seq="";my $crispr_index; my $left_flank=""; my $right_flank="";
								#my $k=$i;
								my $array_direction;my $array_start;my $array_stop;my $model_repeat;
								my ($org_id,$direction)=split('\t+',$line);
								if($org_id=~/^>/)
									{
										$org_id=~ s/^>//;
									}
								my $spacer_index=1;	
								#----- get CRISPR index from previous line------------------------------------------------------------------
								my $previous_line=$arr_spacer_file[$i-1];chomp $previous_line;$previous_line=~ s/\r//;
								if($previous_line =~ /Array (\d+) (\d+)-(\d+)/)
									{
										$crispr_index=$1;
										$array_start=$2;
										$array_stop=$3;
									}
								#---- now get the direction ------------------------------------------------------------------------------
								
								if($line=~/Array_Orientation: (\S+)/)
									{
										$array_direction=$1; chomp $array_direction; $array_direction=~s/\r+//g;
									}
									
								#---- in CRISPRDetect v.2.2 there is a bug creating problem with the array_start position---	
								my $k=$i+4;
								while($arr_spacer_file[$k]!~/====/)
									{
										my $current_line=$arr_spacer_file[$k]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;
										my @arr_t1=split('\t',$current_line);
										my $r_start=$arr_t1[0];$r_start=~s/\s+//g;
										my $r_length=$arr_t1[1];$r_length=~s/\s+//g;
										my $s_length=$arr_t1[3];$s_length=~s/\s+//g;
										my $r_seq=$arr_t1[4]; $r_seq=~s/\s+//g; 
										my $s_seq=$arr_t1[5]; $s_seq=~s/\s+//g;
										my $comment=$arr_t1[6];		$s_length=~s/^\s+//;							
										
										#---- dont forget to check insertion bases in comment ---------------------------
										my %hash_of_insertion_positions;
										my $no_of_insertions=0;
										if($arr_t1[6] and $arr_t1[6]!~/^Del/)
											{
												$comment=$arr_t1[6]; chomp $comment; $comment=~s/^\s+//;											
												my @tmp_arr1=split(' ',$comment);
												my $insertion_bases=$tmp_arr1[0];
												my $insertion_positions=$tmp_arr1[1];
												
												$insertion_bases=~s/,//g;
												$no_of_insertions=length($insertion_bases);													
											}
										
										#---- in CRISPRDetect v.2.2 there is a bug creating problem with the array_start position---	
										if($k==($i+4))
											{
												$array_start=$r_start;
											}	
										if($s_seq=~/\|/) #--- the last repeat --
											{
												$r_seq=~s/-//g;												
												if($array_direction!~/^Rev/)
													{
														$array_stop=$r_start+length($r_seq)+$no_of_insertions;
													}
												else{
														$array_stop=$r_start-length($r_seq)-$no_of_insertions;																										
													}												
											}										
										$k++;
									}			
								
								#---now convert all the dotted repeat, and add to sequence_strand on the go
								$j=$i+4;
								while($arr_spacer_file[$j]!~/====/)
									{
										my $current_line=$arr_spacer_file[$j]; chomp $current_line; $current_line=~s/\r+//g; $current_line=~s/^\s+//;
										my @arr_t1=split('\t',$current_line);
										my $r_start=$arr_t1[0];$r_start=~s/\s+//g;
										my $r_length=$arr_t1[1];$r_length=~s/\s+//g;
										my $s_length=$arr_t1[3];$s_length=~s/\s+//g;
										my $r_seq=$arr_t1[4]; $r_seq=~s/\s+//g; 
										my $s_seq=$arr_t1[5]; $s_seq=~s/\s+//g;
										my $comment=$arr_t1[6];		$s_length=~s/^\s+//;							
										
											
										#---- dont forget to check insertion bases in comment ---------------------------
										my %hash_of_insertion_positions;
										my $no_of_insertions=0;
										if($arr_t1[6] and $arr_t1[6]!~/^Del/)
											{
												$comment=$arr_t1[6]; chomp $comment; $comment=~s/^\s+//;												
												my @tmp_arr1=split(' ',$comment);
												my $insertion_bases=$tmp_arr1[0];
												my $insertion_positions=$tmp_arr1[1];
												
												$insertion_bases=~s/,//g;
												$no_of_insertions=length($insertion_bases);	
												
											}
										#--------------------------------------------------------------------------------
										$r_seq=~s/-//g;
										$s_seq=~s/-//g;
											
										if($s_seq!~/^\|/)
											{	
												my $spacer_start;my $spacer_stop;											
												if($array_direction!~/^Rev/)
													{
														$spacer_start=$r_start+$r_length+$no_of_insertions;
														$spacer_stop=$spacer_start+($s_length-1);
														$array_direction="Forward";
													}
												else{
														$spacer_start=$r_start-$r_length-$no_of_insertions;
														$spacer_stop=$spacer_start-($s_length-1);												
													}
												if(length($s_seq)>12)
													{	
														print WR ">$org_id\|$crispr_index\_$array_start\_$array_stop\|$crispr_index\_$spacer_index\_$spacer_start\_$spacer_stop\|$s_length\n$s_seq\n";	
													}	
												$spacer_index++												
											}	
										#--------------------------------------------------------------------------------------------	

										$j++;
									}
								$i=$j;
							}
					}
					
				close(WR);				
			}
		
		#------ process the CRISPRfinder spacer text
		elsif($spacer_type eq "FASTA")
			{
				open(WR, ">$spacer_seq_fasta_file");				
				my($org_id,$crispr_index,$sequence_start,$sequence_stop,$sequence_strand);
				   $crispr_index=0;
				for(my $i=0;$i<=$#arr_spacer_file;$i++)
					{
						my $line=$arr_spacer_file[$i]; chomp $line; $line=~ s/\r//g;
						
						$line=~ s/^\s+//;
						
						if($line !~/\S/){next;}
						
						$line=~s/\|/_/g;
						$line=~s/\s+/_/g;
						$line=~s/-/_/g;
						$line=~s/\t+/_/g;
						$line=~s/,/_/g;
						$line=~s/_+/_/g;
						$line=~s/:/_/g;
						$line=~s/_+$/_/g;
												
						if($line=~ /^>/)     #------ first level------------------------------------------------------------------------------------------------
							{
								if($i!=0)
									{
										$crispr_index++;
										$sequence_start=1;
										$sequence_stop=length($sequence_strand);
										if($sequence_strand=~/N/i){$sequence_strand=~s/N//gi;}
										if(length($sequence_strand)<=75)
											{
												print WR ">$org_id\|$crispr_index\_$sequence_start\_$sequence_stop\|$crispr_index\_1_$sequence_start\_$sequence_stop\|$sequence_stop\n$sequence_strand\n";
											}	
									}
									
								#----- store the seq id in org_id -----	
								$org_id=$line;
								$org_id=~ s/>//;
								$sequence_strand="";
							}
						else{
								$sequence_strand=$sequence_strand.$line;
							}	
					}
				#---- finally write the last record -----
				$crispr_index++;
				$sequence_start=1;
				$sequence_stop=length($sequence_strand);
				if($sequence_strand=~/N/i){$sequence_strand=~s/N//gi;}
				if(length($sequence_strand)<=75)
					{
						print WR ">$org_id\|$crispr_index\_$sequence_start\_$sequence_stop\|$crispr_index\_1_$sequence_start\_$sequence_stop\|$sequence_stop\n$sequence_strand\n";
					}
				close(WR);	
			}
		elsif($spacer_type eq "CRISPRCasFinder")
			{
				my $path = "$spacer_file";
				my $json;

				{
   					local $/;
    				open my $fh, '<', $path or die $!;
   					$json = <$fh>;
				}

				my $perl = decode_json $json;

				my $command_used = $perl->{Command};
				my $gcf <- "NA";

				if (index($command_used, "#GC") >= 0) {
					($gcf) = ($command_used =~ /.*\#(GC\S+)\s+.*/);
				}
				
				my @seqs = @{$perl->{Sequences}};
				
				open(WR, ">$spacer_seq_fasta_file");
				
				foreach my $seq (@seqs) {
					my $seq_id = $seq->{Id};
					my $org = $seq->{Description};
					
					$org =~ s/\,//g;
					$org =~ s/\://g;
					$org =~ s/\;//g;
					$org =~ s/\s/_/g;
					$org =~ s/\-/_/g;
					$org =~ s/\./_/g;
					$org =~ s/\#/_/g;
					$org =~ s/\(//g;
					$org =~ s/\)//g;
					$org =~ s/\[//g;
					$org =~ s/\]//g;
					$org =~ s/\%/_/g;
					
					my @crisprs = @{$seq->{Crisprs}};
					my $count = 0;
					foreach my $crispr (@crisprs) {
						$count++;
						my $org_id = $seq_id . "_" . $org;
						my $arr_start = $crispr->{Start};
						my $arr_end = $crispr->{End};
						my $crispr_index = $count;
						
						my @regions = @{$crispr->{Regions}};
						
						my $spacer_count = 0;
						foreach my $region (@regions) {
							my $type = $region->{Type};
							
							if ($type ne "Spacer") {
								next;
							}
							
							$spacer_count++;
							my $spacer_start = $region->{Start};
							my $spacer_end = $region->{End};
							my $spacer_seq = $region->{Sequence};
							my $spacer_index = $spacer_count;
							my $spacer_length = length($spacer_seq);
							my $head = ">$org_id\|$crispr_index\_$arr_start\_$arr_end\|$crispr_index\_$spacer_index\_$spacer_start\_$spacer_end\|$spacer_length";
							
							print WR "$head\n";
							print WR "$spacer_seq\n";
						}
					}				
				}
				close(WR);
			}

		#------------- remove redundancy from spacers -----------------------------------------
		if($spacer_remove_redundancy==1)
			{				
				my %hash_of_non_redundant_spacers;
				my %hash_of_redundant_spacers;
				open(RD,"$spacer_seq_fasta_file");
				my @arr_spacers=<RD>;
				close(RD);
				
				for(my $i=0;$i<=$#arr_spacers;$i++)
					{
						my $line=$arr_spacers[$i];
						if($line=~/>/)
							{
								chomp $line;
								$line=~ s/\r//;
								
								my $id_line=$line;
								my $seq=$arr_spacers[$i+1]; chomp $seq; $seq=~ s/\r//;
								
								if(not defined $hash_of_non_redundant_spacers{$seq})
									{
										$hash_of_non_redundant_spacers{$seq}=$id_line;
									}
								else{
										$hash_of_redundant_spacers{$id_line}=$seq;
									}	
							}
					}
				#---- now write the non-redundant spacers to the file again --
				open(WR,">$spacer_seq_fasta_file");
				foreach my $spacer(sort keys %hash_of_non_redundant_spacers)
					{
						print WR "$hash_of_non_redundant_spacers{$spacer}\n$spacer\n";
					}
				close(WR);				
				my $no_of_nonredundant_spacers=keys %hash_of_non_redundant_spacers;
				#----- now write the removed spacers ---------------------------------
				my $removed_spacer_seq_file="$outdir/".$idnumber.".removed_spacer_sequences.fna";
				open(WR,">$removed_spacer_seq_file");
				foreach my $id(sort keys %hash_of_redundant_spacers)
					{
						print WR "$id\n$hash_of_redundant_spacers{$id}\n";
					}
				close(WR);
				my $no_of_redundant_spacers=keys %hash_of_redundant_spacers;
			}
		#-------------------------------------------------------------------------------------------
		return 1;
	}
	
sub fetch_crispr_seq()
	{
		my($id,$sp_array_index,$file,$start,$stop,$flank_3p,$flank_5p)=@_;
		
		open(RD,"$file") or print "$!";
		my @file_content=<RD>;
		close(RD);
		
		my $seq_direction="F";
		my $index=1;
		my $return_seq_strand="";
		
		my $seq_start; my $seq_stop; 
		my $file_type="";
		for( my $i=0;$i<=$#file_content;$i++)
			{
				my $line=$file_content[$i];	chomp $line; $line=~ s/\r//g;
				
				if($line=~ /$id/)
					{
						my ($tmp_id,$pos_det)=split('\t',$line);
						my @arr_tmp=split('_',$pos_det);
						my $array_index = $arr_tmp[1];
						$seq_start = $arr_tmp[2];
						$seq_stop = $arr_tmp[3];
						
						if($arr_tmp[3] == 0){$file_type="SINGLE";last;}
						else{$file_type="MULTI";}
					
						if ("$array_index" ne "$sp_array_index") {next;}

						if($file_type eq "MULTI" )
							{
								my $next_line=$file_content[$i+1]; chomp $next_line; $next_line=~ s/\r//g;
											
								
								
								#----- check if the source sequence is reversed : Arrays predicted by CRISPRDetect may be in Reverse order
								
								if($start>$stop)
									{
										my $t_stop=$stop;
										$stop=$start;
										$start=$t_stop;
									}
								if($seq_start>$seq_stop)
									{
										my $t_stop=$seq_stop;
										$seq_stop=$seq_start;
										$seq_start=$t_stop+1;
										
										#---- now reverse the sequence: not rev_comp
										$next_line=reverse $next_line;
										$seq_direction="R";
									}

								#---------------------------------------------------------------------------------------------------------
								my @tmp_array=split('',$next_line);													
								$index=$seq_start;
						
								if($start>=$seq_start and $stop<=$seq_stop)
									{
										
										foreach my $nt(@tmp_array)
											{									
												if($index>=$start and $index<=$stop)
													{
														$return_seq_strand=$return_seq_strand.$nt;
													}
												if($index>$stop){last;}	
												$index++;
											}		
										last;	
									}									
								elsif(($start+$flank_3p)>=$seq_start and ($stop-$flank_5p)<=$seq_stop)
									{
										
										#---- push the bases in a hash -----
										my %tmp_hash;
										foreach my $nt(@tmp_array)
											{
												$tmp_hash{$index}=$nt;
												$index++	
											}
										for(my $j=$start;$j<=$stop;$j++)
											{
												if(not defined $tmp_hash{$j})
													{
														$return_seq_strand=$return_seq_strand."-";
													}
												else{
														$return_seq_strand=$return_seq_strand.$tmp_hash{$j};
													}	
											}
										last;	
									}
							}								
					}					
									
			}
		
		if($file_type eq "SINGLE")
			{
				shift(@file_content); # remove the first line
					
				foreach my $line(@file_content)
					{
						chomp $line; $line=~ s/\r//g;
						my @tmp_array=split('',$line);
						
						if ($start >= 1 and $stop <= scalar(@tmp_array)) {
							foreach my $nt(@tmp_array)
								{							
									if($index>=$start and $index<=$stop)
										{
											$return_seq_strand=$return_seq_strand.$nt;
										}
									if($index>$stop){last;}	
									$index++;
								}
							last;
							}
						elsif (($start+$flank3p) >= 1 and ($stop-$flank_5p) <= scalar(@tmp_array)) {
							my %tmp_hash;
							foreach my $nt(@tmp_array)
								{
									$tmp_hash{$index}=$nt;
									$index++	
								}
							for(my $j=$start;$j<=$stop;$j++)
								{
									if(not defined $tmp_hash{$j})
										{
											$return_seq_strand=$return_seq_strand."-";
										}
									else{
											$return_seq_strand=$return_seq_strand.$tmp_hash{$j};
										}	
								}
							last;
							}				
					}
			}
		
		if($seq_direction eq "R")
			{
				$return_seq_strand=reverse $return_seq_strand;		### remember, this is spacer sequence, and tested to be 100% correct, so don't mess it up				
			}
		
		return($return_seq_strand);
	}
	
sub get_path() {
	my $dir=shift(@_);
	my @arr_p1=split('\/',$cd_path);
	pop(@arr_p1);
	$dir=join("\/",@arr_p1);
		
	return $dir;
}