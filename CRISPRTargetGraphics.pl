#!/usr/bin/perl

my $file_in = "";
my $file_out = "";
my $html = 1;
my $include_url = 0;
my $include_pvalue = 0;

my $ind = 0;
foreach(@ARGV) {
	if (@ARGV[$ind] eq '-in') {
		$file_in = @ARGV[$ind + 1];
		if (! (-e $file_in)) {
			die "cannot open file: " . $file_in . "\n";
		}
	}
	
	if (@ARGV[$ind] eq '-html') {
		$file_out = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-txt') {
		$file_out = @ARGV[$ind + 1];
		$file_out =~ s/\.html/.txt/g;
		$html = 0;
	}
	
	if (@ARGV[$ind] eq '-include_url') {
		$include_url = 1;
	}
	
	if (@ARGV[$ind] eq '-include_pvalue') {
		$include_pvalue = 1;
	}
	
	$ind++;
}

if (-e $file_out) {
	unlink($file_out);
}

my $flag = 0;
if ($html > 0) {
	$flag = getStyleHTML($file_out);
	$flag = getDisplayHeader($file_out);
}

$flag = createReport($file_in, $file_out);

if ($html > 0) {
	$flag = getDisplayFooter($file_out);
}

sub createReport {
	my ($path, $out) = @_;

	open(TXT, $path);
	my @contents = <TXT>;
	close(TXT);
	
	my $count = 0;
	
	open(OUT, ">>$out");
	
	if ($html < 1) {
		print OUT "# CRISPR target predictions in CRISPRTarget layout (text)\n";
		print OUT "# Prediction data: $file_in\n\n";
	}
	
	for (my $i = 0; $i < scalar(@contents); $i++) {
		my $line = $contents[$i];
		$line =~ s/\n//g;
		
		if ($line =~ /^#.*/) {
			next;
		}
		
		if ($line =~ /^Spacer_accession/) {
			next;
		}
		
		$count++;
		my $hit_id = $count;
		my @toks = split(/[\t]/, $line);
		
		my $strand = $toks[8];
		my $score = $toks[7];
		
		my $fp3p = $toks[5];
		my $fp5p = $toks[6];
		my $psp = $toks[15];
		my $psp_id = $toks[2];
		my $psp_start = $toks[3];
		my $psp_end = $toks[4];
		my $targeted_org = $toks[17];
		
		my $fs3p = $toks[20];
		my $fs5p = $toks[19];
		my $spacer = $toks[14];
		my $spacer_org = $toks[16];
		my $spacer_ind = $toks[1];
		
		my $crispr_type = $toks[21];
		my $dr_seq = $toks[22];
		my $url_link = $toks[18];
		
		my $spacer_flank_length_5p = length($fs5p);
		my $spacer_flank_length_3p = length($fs3p);
		
		my $url_new = "NA";
		my $a_text = "";
		if ($include_url > 0)  {
			($url_new, $a_text) = guessURL($targeted_org, $psp_id, $psp_start, $psp_end, $spacer_flank_length_5p, $spacer_flank_length_3p);
			if ($url_link eq "NA" || $url_link eq "") {
				$url_link = $url_new;
			}
		} else {
			$url_new = "NA";
			$a_text = "";
		}
		
		my $spacer_ind_new = $spacer_ind;
		
		if ($spacer_ind_new =~ /\S+_\d+_\d+\|\d+\-\d+/) {
			my ($id_short, $array_ind, $spacer_ind, $spacer_start, $spacer_end) = ($spacer_ind_new =~ /(\S+)_(\d+)_(\d+)\|(\d+)\-(\d+)/);
			my $id_long = $id_short . "_" . $array_ind . "_" . $spacer_ind;
			$spacer_ind_new = "CRISPR No.$array_ind Spacer No.$spacer_ind \($id_long\) position: $spacer_start-$spacer_end";
		}
		
		my $full_desc = "$hit_id. Match to:  $targeted_org \($psp_id\) position: $psp_start-$psp_end, with $spacer_org $spacer_ind_new, Strand: $strand, Direct Repeat: $dr_seq, Type[cctyper]: $crispr_type";
		my $full_desc_html = "<b>$hit_id</b>. Match to: <b>$targeted_org<\/b> \($psp_id\) position: $psp_start-$psp_end, with <b>$spacer_org<\/b> $spacer_ind_new, Strand: $strand, Direct Repeat: $dr_seq, Type[cctyper]: $crispr_type<br>";
		
		my $spacer_flank_length_5p = length($fs5p);
		my $spacer_flank_length_3p = length($fs3p);
		if ($fs5p eq "-") {$spacer_flank_length_5p = 0;}
		if ($fs3p eq "-") {$spacer_flank_length_3p = 0;}
		
		my $flank_length_5p = length($fp5p);
		my $flank_length_3p = length($fp3p);
		if ($fp3p eq "-") {$flank_length_3p = 0;}
		if ($fp5p eq "-") {$flank_length_5p = 0;}
		
		if ($spacer_flank_length_5p < $flank_length_5p) {
			if ($spacer_flank_length_5p == 0) {
				$fs5p = stringOfDash($flank_length_5p);
			} else {
				$fs5p = fillAlnString($fs5p, abs($flank_length_5p - $spacer_flank_length_5p), 1);
			}
			$spacer_flank_length_5p = $flank_length_5p;
		}
		
		if ($spacer_flank_length_5p > $flank_length_5p) {
			if ($flank_length_5p == 0) {
				$fp5p = stringOfDash($spacer_flank_length_5p);
			} else {
				$fp5p = fillAlnString($fp5p, abs($spacer_flank_length_5p - $flank_length_5p), 0);
			}
			$flank_length_5p = $spacer_flank_length_5p;
		}
		
		if ($spacer_flank_length_3p < $flank_length_3p) {
			if ($spacer_flank_length_3p == 0) {
				$fs3p = stringOfDash($flank_length_3p);
			} else {
				$fs3p = fillAlnString($fs3p, abs($flank_length_3p - $spacer_flank_length_3p), 0);
			}
			$spacer_flank_length_3p = $flank_length_3p;
		}
		
		if ($spacer_flank_length_3p > $flank_length_3p) {
			if ($flank_length_3p == 0) {
				$fp3p = stringOfDash($spacer_flank_length_3p);
			} else {
				$fp3p = fillAlnString($fp3p, abs($spacer_flank_length_3p - $flank_length_3p), 1);
			}
			$flank_length_3p = $spacer_flank_length_3p;
		}
		
		my $fp3p_rc = reverse($fp3p); $fp3p_rc =~ tr/ACGT/TGCA/;
		my $fp5p_rc = reverse($fp5p); $fp5p_rc =~ tr/ACGT/TGCA/;
		my $psp_rc = reverse($psp); $psp_rc =~ tr/ACGT/TGCA/;
		
		my $topline = $fs5p . $spacer . $fs3p;
		my $midline = reverse($fp3p) . reverse($psp) . reverse($fp5p);
		my $btmline = $fp3p_rc . $psp_rc . $fp5p_rc;
		
		my $topline_html = "&nbsp;5<sup>'</sup>&nbsp;" . "<span class=\"flank\">". $fs5p . "</span><span class=\"spacer\">" . $spacer . "</span><span class=\"flank\">" . $fs3p . "</span>&nbsp;3<sup>'</sup>&nbsp;<- CRISPR spacer RNA<br>";
		my $midline_html = "&nbsp;3<sup>'</sup>&nbsp;" . "<span class=\"flank\">". reverse($fp3p) . "</span><span class=\"spacer\">" . $spacer . "</span><span class=\"flank\">" . $fs3p . "</span>&nbsp;5<sup>'</sup>&nbsp;<- Protospacer Sequence<br>";
		my $btmline_html = "";
		
		if ($url_link eq "" || $url_link eq "NA") {
			$btmline_html = "&nbsp;5<sup>'</sup>&nbsp;" . "<span class=\"flank\">". $fp3p_rc . "</span><span class=\"sub\">" . $psp_rc . "</span><span class=\"flank\">" . $fp5p_rc . "</span>&nbsp;3<sup>'</sup>&nbsp;<- []<br>";
		} else {
			$btmline_html = "&nbsp;5<sup>'</sup>&nbsp;" . "<span class=\"flank\">". $fp3p_rc . "</span><span class=\"sub\">" . $psp_rc . "</span><span class=\"flank\">" . $fp5p_rc . "</span>&nbsp;3<sup>'</sup>&nbsp;<- [<a href=\"$url_link\" target='_blank'>$a_text</a>]<br>";
		}
		
		my $pam3pf = $toks[9];
		my $pam5pf = $toks[10];
		my $pam3pr = $toks[11];
		my $pam5pr = $toks[12];
		my $lpam3pf = length($pam3pf);
		my $lpam5pf = length($pam5pf);
		my $lpam3pr = length($pam3pr);
		my $lpam5pr = length($pam5pr);
		if ($pam3pf eq "NA") {$lpam3pf = 0;}
		if ($pam5pf eq "NA") {$lpam5pf = 0;}
		if ($pam3pr eq "NA") {$lpam3pr = 0;}
		if ($pam5pr eq "NA") {$lpam5pr = 0;}
		
		my $pam_type5pf = $toks[23];
		my $pam_type5pr = $toks[24];
		my $pam_type3pf = $toks[25];
		my $pam_type3pr = $toks[26];
		
		my $pam_type5p = $pam_type5pf;
		my $lpam5p = $lpam5pf;
		if ($strand eq "-") {
			$pam_type5p = $pam_type5pr;
			$lpam5p = $lpam5pr;
		}
		
		my $pam_type3p = $pam_type3pf;
		my $lpam3p = $lpam3pf;
		if ($strand eq "-") {
			$pam_type3p = $pam_type3pr;
			$lpam3p = $lpam3pr;
		}
		
		my $midline_pam = $midline;
		my $midline_pam_html = $midline_html;
		
		if ($flank_length_3p > 0 || $flank_length_5p > 0) {
			my $strtmp = "";
			my $strtmp_html = "";
			if ($flank_length_3p > 0) {
				my $p1 = substr($midline, 0, $flank_length_3p - $lpam3p);
				my $p2 = substr($midline, $flank_length_3p - $lpam3p, $lpam3p);
				$strtmp .= $p1 . lc($p2);
				$strtmp_html .= "<span class=\"flank\">" . $p1 . "</span><span class=\"pam\" title='$pam_type3p'>" . $p2 . "</span>";
			}
		
			if ($flank_length_5p > 0) {
				my $p3 = substr($midline, $flank_length_3p, -$flank_length_5p);
				my $p4 = substr($midline, -$flank_length_5p, $lpam5p);
				my $p5 = substr($midline, -$flank_length_5p + $lpam5p, $flank_length_5p - $lpam5p);
				$strtmp .= $p3 . lc($p4) . $p5;
				$strtmp_html .= "<span class=\"spacer\">" . $p3 . "</span><span class=\"pam\" title='$pam_type5p'>" . $p4 . "</span></span><span class=\"flank\">" . $p5 . "</span>";
			} else {
				my $p3 = substr($midline, $flank_length_3p);
				$strtmp_html .= "<span class=\"spacer\">" . $p3 . "</span>";
			}
			$midline_pam = "3' " . $strtmp . " 5' <- Protospacer Sequence";
			$midline_pam_html = "&nbsp;3<sup>'</sup>&nbsp;" . $strtmp_html . "</span>&nbsp;5<sup>'</sup>&nbsp;<- Protospacer Sequence<br>";
		} else {
			$midline_pam = $midline;
			$midline_pam_html = $midline_html;
		}
		
		my $score_line = "Score: $score ";
		my $score_line_html = "<div class=\"features\"><hr></hr>&nbsp;Score: $score";
		
		if ($include_pvalue > 0) {
			my $pvalue = $toks[27];
			$pvalue = sprintf("%.4f", $pvalue);
			$score_line = "Score: $score (P-value: $pvalue)";
			$score_line_html = "<div class=\"features\"><hr></hr>&nbsp;Score: $score (P-value: $pvalue)";
		}
		
		if ($strand eq "+") {
			if ($pam_type5pf ne "NA") {
				$score_line .= ", PAM found in forward strands 5' Tag: $pam_type5pf ($pam5pf)";
				$score_line_html .= "&nbsp;&nbsp;PAM found in forward strands 5' Tag: $pam_type5pf \($pam5pf\) &nbsp;";
			}
		
			if ($pam_type3pf ne "NA") {
				$score_line .= ", PAM found in forward strands 3' Tag: $pam_type3pf ($pam3pf)";
				$score_line_html .= "&nbsp;&nbsp;PAM found in forward strands 3' Tag: $pam_type3pf \($pam3pf\) &nbsp;";
			}
		}
		
		if ($strand eq "-") {
			if ($pam_type5pr ne "NA") {
				$score_line .= ", PAM found in reverse strands 5' Tag: $pam_type5pr ($pam5pr)";
				$score_line_html .= "&nbsp;&nbsp;PAM found in reverse strands 5' Tag: $pam_type5pr \($pam5pr\) &nbsp;";
			}
		
			if ($pam_type3pr ne "NA") {
				$score_line .= ", PAM found in reverse strands 3' Tag: $pam_type3pr ($pam3pr)";
				$score_line_html .= "&nbsp;&nbsp;PAM found in reverse strands 3' Tag: $pam_type3pr \($pam3pr\) &nbsp;";
			}
		}
		$score_line_html .= "</div>";
		
		my $alnstr1 = create_basepairing_string_text($topline, $midline, $flank_length_3p, $flank_length_5p);
		my $alnstr2 = create_basepairing_string_text($midline, $btmline, $flank_length_3p, $flank_length_5p);
		my $alnstr1_html = create_basepairing_string_html($topline, $midline, $flank_length_3p, $flank_length_5p);
		my $alnstr2_html = create_basepairing_string_html($midline, $btmline, $flank_length_3p, $flank_length_5p);
		
		my $bp_part1 = substr($alnstr2, 0, $flank_length_3p);							
		my $bp_part2 = substr($alnstr2,  $flank_length_3p, -$flank_length_5p);
		my $bp_part2_html = $bp_part2;
		my $sp = " "; $bp_part2 =~ s/\|/$sp/g;
		my $sp_html = "&nbsp;"; $bp_part2_html =~ s/\|/$sp_html/g;
									 
		my $bp_part3 = substr($alnstr2, -$flank_length_5p, $flank_length_5p);
		$alnstr2 = $bp_part1 . $bp_part2 . $bp_part3;
		$alnstr2_html = $bp_part1 . $bp_part2_html . $bp_part3;
		
		$alnstr1 = "   " . $alnstr1;
		$alnstr2 = "   " . $alnstr2;
		$alnstr1_html = "&nbsp;&nbsp;&nbsp;&nbsp;" . $alnstr1_html . "<br>";
		$alnstr2_html = "&nbsp;&nbsp;&nbsp;&nbsp;" . $alnstr2_html . "<br>";
		
		$topline = "5' " . $topline . " 3' <- CRISPR spacer RNA";
		$midline = "3' " . $midline . " 5' <- Protospacer Sequence";
		
		if ($url_link eq "" || $url_link eq "NA") {
			$btmline = "5' " . $btmline . " 3' <- []";
		} else {
			$btmline = "5' " . $btmline . " 3' <- [$url_link]";
		}
		
		if ($html < 1) {
			print OUT $full_desc . "\n";
			print OUT $topline . "\n";
			print OUT $alnstr1 . "\n";
			print OUT $midline_pam . "\n";
			print OUT $alnstr2 . "\n";
			print OUT $btmline . "\n";
			print OUT $score_line . "\n";
			print OUT "//\n\n";
		} else {
			print OUT "<div class =\"box\">";
			print OUT $full_desc_html;
			print OUT "<table style='border:0px solid silver;width:100%;'>";
			print OUT "<tr><td style='border-top:1px solid silver;'>";
			print OUT "<div class=\"sequence\"><font face=\"Courier\">";
			print OUT $topline_html;
			print OUT $alnstr1_html;
			print OUT $midline_pam_html;
			print OUT $alnstr2_html;
			print OUT $btmline_html;
			print OUT "</font></div>";
			print OUT "</td></tr>";
			print OUT "</table>";
			print OUT $score_line_html;
			print OUT "</div>\n";
		}
	}
	close(OUT);
	
	return 1;
}

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

sub create_basepairing_string_text {
	my($seq1, $seq2, $flank_3p, $flank_5p) = @_;
			
	chomp $seq1;
	chomp $seq2;
			
	my @arr_seq1=split('',$seq1);
	my @arr_seq2=split('',$seq2);
			
	my %hash_of_complementary_bases;
	$hash_of_complementary_bases{'A'}="T";
	$hash_of_complementary_bases{'C'}="G";
	$hash_of_complementary_bases{'G'}="C";
	$hash_of_complementary_bases{'T'}="A";
	$hash_of_complementary_bases{'U'}="A";
		
	my $basepairing_string="";
	for my $i(0..$#arr_seq1) {
		if(not defined $arr_seq1[$i] or not defined $arr_seq2[$i] or not defined $hash_of_complementary_bases{$arr_seq1[$i]}) {
				$basepairing_string = $basepairing_string . " ";
		} elsif ($hash_of_complementary_bases{$arr_seq1[$i]} eq $arr_seq2[$i]) {
				$basepairing_string = $basepairing_string . "|";
		} else {
			if ($i >= $flank_3p and $i <= ($#arr_seq1 - $flank_5p)) {
				$basepairing_string = $basepairing_string . " ";
			} else {									
				$basepairing_string = $basepairing_string . " ";
			}	
		}	
	}
	return($basepairing_string);
}

sub create_basepairing_string_html {
	my($seq1, $seq2, $flank_3p, $flank_5p) = @_;
			
	chomp $seq1;
	chomp $seq2;
			
	my @arr_seq1=split('',$seq1);
	my @arr_seq2=split('',$seq2);
			
	my %hash_of_complementary_bases;
	$hash_of_complementary_bases{'A'}="T";
	$hash_of_complementary_bases{'C'}="G";
	$hash_of_complementary_bases{'G'}="C";
	$hash_of_complementary_bases{'T'}="A";
	$hash_of_complementary_bases{'U'}="A";
		
	my $basepairing_string="";
	for my $i(0..$#arr_seq1) {
		if(not defined $arr_seq1[$i] or not defined $arr_seq2[$i] or not defined $hash_of_complementary_bases{$arr_seq1[$i]}) {
				$basepairing_string = $basepairing_string . "&nbsp;";
		} elsif ($hash_of_complementary_bases{$arr_seq1[$i]} eq $arr_seq2[$i]) {
				$basepairing_string = $basepairing_string . "|";
		} else {
			if ($i >= $flank_3p and $i <= ($#arr_seq1 - $flank_5p)) {
				$basepairing_string = $basepairing_string . "<span style='background-color:red'>&nbsp;</span>";
			} else {									
				$basepairing_string = $basepairing_string . "&nbsp;";
			}	
		}	
	}
	return($basepairing_string);
}

sub getStyleHTML {
	my ($out) = @_;
	my $html_str = qq~
	<!doctype html>
	<html>
	<STYLE TYPE=\"text/css\">
	 <!-- 
	 body
		{	
			font-family: times;
			font-size: 15px;
			background-color: #ffffff;
			color: #000000;
			text-align: left;
		}
	.flank { background-color: #00FFFF; font-family: courier;} 
	.pam { background-color: green; font-family: courier;} 
	.sub { background-color: white; font-family: courier; }
	.spacer{ background-color: yellow; font-family: courier;} 
	.feature{ clear: left; float: left; width: 30%; min-width: 40px; }
	.sequence{float:left; margin: 5px;margin-left: 10px;}
	.features{clear:left;}	
	.element{width:800px;font-family: courier;font-size:9pt; margin-left: 20px}	
	.box{float:left; width:98%;min-width: 800px; padding: 10px; margin: 5px;-moz-border-radius: 10px; -webkit-border-radius: 10px; border: 1px solid #000;background-color:white;};
	.img{border-style:none;}
	-->
	</STYLE>		
	<body>
	<div id="main" style="position:absolute;top:1px;left:1px;right:1px;bottom:1px;border:1px solid silver;">
	<h1>CRISPR target predictions in CRISPRTarget layout</h1>
	<h2>Prediction data: $file_in</h2>
	~;
	
	open(OUT, ">>$out");
	print OUT $html_str . "\n";
	close(OUT);
	return $html_str;
}

sub getDisplayHeader {
	my ($out) = @_;
	my $html_str = qq~
	<div id="result_div" name="result_div" style="height:'73%';overflow-y:auto;overflow-x:hidden;border:1px solid white;background-color:white;" >
	\n<!-- PRINT_START -->\n
	~;
	open(OUT, ">>$out");
	print OUT $html_str . "\n";
	close(OUT);
	return $html_str;
}

sub getDisplayFooter {
	my ($out) = @_;
	my $html_str = qq~\n<!-- PRINT_STOP -->\n\n\t</div></div></body>\n</html>\n~;
	open(OUT, ">>$out");
	print OUT $html_str . "\n";
	close(OUT);
	return $html_str;
}

sub guessURL {
	my ($targeted_org, $psp_id, $start, $end, $offset_5p, $offset_3p) = @_;
	
	my $link_url = "";
	my $a_text = "";
	if($psp_id =~ /[a-z]+\_[0-9]+\_GID\_[0-9]+/) {
		my @arr_aclame_id = split('_',$psp_id);
		$link_url = "http://aclame.ulb.ac.be/perl/Aclame/Genomes/genes_view.cgi?view=gene&id=gene:$arr_aclame_id[0]:$arr_aclame_id[1]";
		$a_text = "ACLAME Gene";
	} else {
		if (($psp_id=~/Ga\d+\_\d+/) or ($psp_id=~/^IMGVR/) or ($psp_id=~/^DTR/) or ($psp_id=~/^UGV\-GENOME/) or ($psp_id=~/_NODE_/) or ($psp_id=~/Gammaproteobacteria_gi_/)) {
			if ($psp_id =~ /^IMGVR_UViG_\S+\|\S+\|\S+\|\S+-\S+/) {
				my ($uvig) =  ($psp_id =~ /^(IMGVR_UViG_\S+)\|\S+\|\S+\|\S+-\S+/);
				$link_url="https://img.jgi.doe.gov/cgi-bin/vr/main.cgi?section=ViralBrowse&page=uviginfo&uvig_id=$uvig";
				$a_text="JGI Link (UViG=$uvig)";
			} elsif ($psp_id=~/^IMGVR_UViG_\d+_\d+\S+/) {
				my ($uvig) = ($psp_id=~/^(IMGVR_UViG_\d+_\d+)\S+/);
				$link_url="https://img.jgi.doe.gov/cgi-bin/vr/main.cgi?section=ViralBrowse&page=uviginfo&uvig_id=$uvig";
				$a_text="JGI Link (UViG=$uvig)";
			} elsif ($psp_id=~/\S+_gi_\d+\-\d+\-\d+/) {
				my ($giid, $gist) = ($psp_id=~/\S+_gi_(\d+)\-(\d+)\-\d+/);
				my $incr = $gist - 1;
				$link_url="http://www.ncbi.nlm.nih.gov/nuccore/".$giid."?report=graph&v=".(int($start)-$offset_5p).":".(int($end)+$offset_3p)."&c=339966&theme=Default&flip=false&select=null&content=1&color=0&label=1&geneModel=1&decor=0&layout=0&spacing=0&alncolor=on";
				$a_text="JGI mapped onto NCBI (GI=$giid)";
			} elsif ($psp_id=~/\S+_gi_\d+.*/) {
				my ($giid) = ($psp_id=~/\S+_gi_(\d+).*/);
				$link_url="http://www.ncbi.nlm.nih.gov/nuccore/".$giid."?report=graph&v=".(int($start)-$offset_5p).":".(int($end)+$offset_3p)."&c=339966&theme=Default&flip=false&select=null&content=1&color=0&label=1&geneModel=1&decor=0&layout=0&spacing=0&alncolor=on";
				$a_text="JGI mapped onto NCBI (GI=$giid)";
			} elsif ($psp_id=~/Ga\d+\_\d+/) {
				my ($iid) = ($psp_id=~/(Ga\d+)\_\d+/);
				$link_url="https://gold.jgi.doe.gov/analysis_projects?id=$iid";
				$a_text="JGI Link (Genome_assembly=$iid)";	
			} else {
				$link_url = "https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html";
				$a_text="JGI: Click here for further information";
			}	
		} elsif ($targeted_org =~ /accession/i) {
			$link_url = "http://camera.calit2.net/camdata.shtm";
			$a_text = "Link"; 
		} elsif ($psp_id =~ /CAM_/) {
			$link_url = "http://camera.calit2.net/camdata.shtm";
			$a_text = "CAMERA DB: Click here for further information"; 
		} elsif ($psp_id =~ /USER_DB_SEQUENCE/) {
			$link_url = "";
			$a_text = "Data supplied by end-user";
		} elsif ($psp_id =~ /SEQ_\d+/) {
			$link_url = "http://phast.wishartlab.com/";
			$a_text = "PHAST DB: Click here for further information";
		} else {
			$link_url = "http://www.ncbi.nlm.nih.gov/nuccore/".$psp_id."?report=graph&v=".(int($start)-$offset_5p).":".(int($end)+$offset_3p)."&c=339966&theme=Default&flip=false&select=null&content=1&color=0&label=1&geneModel=1&decor=0&layout=0&spacing=0&alncolor=on";
			$a_text = "Entrez Nucleotide";
		}
	}
	return ($link_url, $a_text);
}





