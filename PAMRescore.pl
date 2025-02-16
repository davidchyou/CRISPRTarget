my $ct_tab_in = ""; #"blast_out/input_array_report.txt";
my $out = ""; #"out.txt";

my $scoring_flank_length_5p = 8;
my $scoring_flank_length_3p = 8;
my $wt_psp_match = 1;
my $wt_psp_mismatch = -1;
my $wt_fp5p_match = 0;
my $wt_fp5p_mismatch = 0;
my $wt_fp3p_match = 0;
my $wt_fp3p_mismatch = 0;
my $show_alignment_string = 0;
my $pam_search_all = 0;
my $pam_search_super = 1;
my $pam_search_exact = 0;
my $pam_crispr_class_super = "";
my $pam_crispr_class_exact = "";
my $pam_user_5p = "";
my $pam_user_3p = "";
my $pam_score_wrong_class = 5;
my $pam_score_super_class = 5;
my $pam_score_exact_class = 5;
my $pam_score_user_pam = 5;
my $score_for_pam = 0;
my $chk_range_start = "1-8";
my $chk_range_end = "1-8";
my $chk_range_start_excl = "6";
my $chk_range_end_excl = "";
my $flag_unrealistic_target = 0;

my %hash_of_pams_3p;
my %hash_of_pams_5p;
$hash_of_pams_3p{'I'}{'I-A'}{'I-A,B-NGG'} = "AGG,CGG,GGG,TGG,";
$hash_of_pams_3p{'I'}{'I-B'}{'I-A,B-NGG'} = "AGG,CGG,GGG,TGG,";
$hash_of_pams_3p{'I'}{'I-C'}{'I-C-GAA'} = "GAA,";
$hash_of_pams_3p{'I'}{'I-E'}{'I-E-CAT,CTT,CCT,CTC'} = "CAT,CTT,CCT,CTC,";
$hash_of_pams_3p{'I'}{'I-F'}{'I-F-CC'} = "CC,";
$hash_of_pams_3p{'I'}{'I-F'}{'I-F-GG'} = "GG,";
			
$hash_of_pams_5p{'II'}{'II-A'}{'II-A-WTTCTNN'}="ATTCTAA,ATTCTAC,ATTCTAG,ATTCTAT,ATTCTCA,ATTCTCC,ATTCTCG,ATTCTCT,ATTCTGA,ATTCTGC,ATTCTGG,ATTCTGT,ATTCTTA,ATTCTTC,ATTCTTG,ATTCTTT,TTTCTAA,TTTCTAC,TTTCTAG,TTTCTAT,TTTCTCA,TTTCTCC,TTTCTCG,TTTCTCT,TTTCTGA,TTTCTGC,TTTCTGG,TTTCTGT,TTTCTTA,TTTCTTC,TTTCTTG,TTTCTTT,";
$hash_of_pams_5p{'II'}{'II-A'}{'II-A-TTTYRNNN'}="TTTCAAAA,TTTCAAAC,TTTCAAAG,TTTCAAAT,TTTCAAAN,TTTCAACA,TTTCAACC,TTTCAACG,TTTCAACT,TTTCAACN,TTTCAAGA,TTTCAAGC,TTTCAAGG,TTTCAAGT,TTTCAAGN,TTTCAATA,TTTCAATC,TTTCAATG,TTTCAATT,TTTCAATN,TTTCAANA,TTTCAANC,TTTCAANG,TTTCAANT,TTTCAANN,TTTCACAA,TTTCACAC,TTTCACAG,TTTCACAT,TTTCACAN,TTTCACCA,TTTCACCC,TTTCACCG,TTTCACCT,TTTCACCN,TTTCACGA,TTTCACGC,TTTCACGG,TTTCACGT,TTTCACGN,TTTCACTA,TTTCACTC,TTTCACTG,TTTCACTT,TTTCACTN,TTTCACNA,TTTCACNC,TTTCACNG,TTTCACNT,TTTCACNN,TTTCAGAA,TTTCAGAC,TTTCAGAG,TTTCAGAT,TTTCAGAN,TTTCAGCA,TTTCAGCC,TTTCAGCG,TTTCAGCT,TTTCAGCN,TTTCAGGA,TTTCAGGC,TTTCAGGG,TTTCAGGT,TTTCAGGN,TTTCAGTA,TTTCAGTC,TTTCAGTG,TTTCAGTT,TTTCAGTN,TTTCAGNA,TTTCAGNC,TTTCAGNG,TTTCAGNT,TTTCAGNN,TTTCATAA,TTTCATAC,TTTCATAG,TTTCATAT,TTTCATAN,TTTCATCA,TTTCATCC,TTTCATCG,TTTCATCT,TTTCATCN,TTTCATGA,TTTCATGC,TTTCATGG,TTTCATGT,TTTCATGN,TTTCATTA,TTTCATTC,TTTCATTG,TTTCATTT,TTTCATTN,TTTCATNA,TTTCATNC,TTTCATNG,TTTCATNT,TTTCATNN,TTTCANAA,TTTCANAC,TTTCANAG,TTTCANAT,TTTCANAN,TTTCANCA,TTTCANCC,TTTCANCG,TTTCANCT,TTTCANCN,TTTCANGA,TTTCANGC,TTTCANGG,TTTCANGT,TTTCANGN,TTTCANTA,TTTCANTC,TTTCANTG,TTTCANTT,TTTCANTN,TTTCANNA,TTTCANNC,TTTCANNG,TTTCANNT,TTTCANNN,TTTCGAAA,TTTCGAAC,TTTCGAAG,TTTCGAAT,TTTCGAAN,TTTCGACA,TTTCGACC,TTTCGACG,TTTCGACT,TTTCGACN,TTTCGAGA,TTTCGAGC,TTTCGAGG,TTTCGAGT,TTTCGAGN,TTTCGATA,TTTCGATC,TTTCGATG,TTTCGATT,TTTCGATN,TTTCGANA,TTTCGANC,TTTCGANG,TTTCGANT,TTTCGANN,TTTCGCAA,TTTCGCAC,TTTCGCAG,TTTCGCAT,TTTCGCAN,TTTCGCCA,TTTCGCCC,TTTCGCCG,TTTCGCCT,TTTCGCCN,TTTCGCGA,TTTCGCGC,TTTCGCGG,TTTCGCGT,TTTCGCGN,TTTCGCTA,TTTCGCTC,TTTCGCTG,TTTCGCTT,TTTCGCTN,TTTCGCNA,TTTCGCNC,TTTCGCNG,TTTCGCNT,TTTCGCNN,TTTCGGAA,TTTCGGAC,TTTCGGAG,TTTCGGAT,TTTCGGAN,TTTCGGCA,TTTCGGCC,TTTCGGCG,TTTCGGCT,TTTCGGCN,TTTCGGGA,TTTCGGGC,TTTCGGGG,TTTCGGGT,TTTCGGGN,TTTCGGTA,TTTCGGTC,TTTCGGTG,TTTCGGTT,TTTCGGTN,TTTCGGNA,TTTCGGNC,TTTCGGNG,TTTCGGNT,TTTCGGNN,TTTCGTAA,TTTCGTAC,TTTCGTAG,TTTCGTAT,TTTCGTAN,TTTCGTCA,TTTCGTCC,TTTCGTCG,TTTCGTCT,TTTCGTCN,TTTCGTGA,TTTCGTGC,TTTCGTGG,TTTCGTGT,TTTCGTGN,TTTCGTTA,TTTCGTTC,TTTCGTTG,TTTCGTTT,TTTCGTTN,TTTCGTNA,TTTCGTNC,TTTCGTNG,TTTCGTNT,TTTCGTNN,TTTCGNAA,TTTCGNAC,TTTCGNAG,TTTCGNAT,TTTCGNAN,TTTCGNCA,TTTCGNCC,TTTCGNCG,TTTCGNCT,TTTCGNCN,TTTCGNGA,TTTCGNGC,TTTCGNGG,TTTCGNGT,TTTCGNGN,TTTCGNTA,TTTCGNTC,TTTCGNTG,TTTCGNTT,TTTCGNTN,TTTCGNNA,TTTCGNNC,TTTCGNNG,TTTCGNNT,TTTCGNNN,TTTCRAAA,TTTCRAAC,TTTCRAAG,TTTCRAAT,TTTCRAAN,TTTCRACA,TTTCRACC,TTTCRACG,TTTCRACT,TTTCRACN,TTTCRAGA,TTTCRAGC,TTTCRAGG,TTTCRAGT,TTTCRAGN,TTTCRATA,TTTCRATC,TTTCRATG,TTTCRATT,TTTCRATN,TTTCRANA,TTTCRANC,TTTCRANG,TTTCRANT,TTTCRANN,TTTCRCAA,TTTCRCAC,TTTCRCAG,TTTCRCAT,TTTCRCAN,TTTCRCCA,TTTCRCCC,TTTCRCCG,TTTCRCCT,TTTCRCCN,TTTCRCGA,TTTCRCGC,TTTCRCGG,TTTCRCGT,TTTCRCGN,TTTCRCTA,TTTCRCTC,TTTCRCTG,TTTCRCTT,TTTCRCTN,TTTCRCNA,TTTCRCNC,TTTCRCNG,TTTCRCNT,TTTCRCNN,TTTCRGAA,TTTCRGAC,TTTCRGAG,TTTCRGAT,TTTCRGAN,TTTCRGCA,TTTCRGCC,TTTCRGCG,TTTCRGCT,TTTCRGCN,TTTCRGGA,TTTCRGGC,TTTCRGGG,TTTCRGGT,TTTCRGGN,TTTCRGTA,TTTCRGTC,TTTCRGTG,TTTCRGTT,TTTCRGTN,TTTCRGNA,TTTCRGNC,TTTCRGNG,TTTCRGNT,TTTCRGNN,TTTCRTAA,TTTCRTAC,TTTCRTAG,TTTCRTAT,TTTCRTAN,TTTCRTCA,TTTCRTCC,TTTCRTCG,TTTCRTCT,TTTCRTCN,TTTCRTGA,TTTCRTGC,TTTCRTGG,TTTCRTGT,TTTCRTGN,TTTCRTTA,TTTCRTTC,TTTCRTTG,TTTCRTTT,TTTCRTTN,TTTCRTNA,TTTCRTNC,TTTCRTNG,TTTCRTNT,TTTCRTNN,TTTCRNAA,TTTCRNAC,TTTCRNAG,TTTCRNAT,TTTCRNAN,TTTCRNCA,TTTCRNCC,TTTCRNCG,TTTCRNCT,TTTCRNCN,TTTCRNGA,TTTCRNGC,TTTCRNGG,TTTCRNGT,TTTCRNGN,TTTCRNTA,TTTCRNTC,TTTCRNTG,TTTCRNTT,TTTCRNTN,TTTCRNNA,TTTCRNNC,TTTCRNNG,TTTCRNNT,TTTCRNNN,TTTTAAAA,TTTTAAAC,TTTTAAAG,TTTTAAAT,TTTTAAAN,TTTTAACA,TTTTAACC,TTTTAACG,TTTTAACT,TTTTAACN,TTTTAAGA,TTTTAAGC,TTTTAAGG,TTTTAAGT,TTTTAAGN,TTTTAATA,TTTTAATC,TTTTAATG,TTTTAATT,TTTTAATN,TTTTAANA,TTTTAANC,TTTTAANG,TTTTAANT,TTTTAANN,TTTTACAA,TTTTACAC,TTTTACAG,TTTTACAT,TTTTACAN,TTTTACCA,TTTTACCC,TTTTACCG,TTTTACCT,TTTTACCN,TTTTACGA,TTTTACGC,TTTTACGG,TTTTACGT,TTTTACGN,TTTTACTA,TTTTACTC,TTTTACTG,TTTTACTT,TTTTACTN,TTTTACNA,TTTTACNC,TTTTACNG,TTTTACNT,TTTTACNN,TTTTAGAA,TTTTAGAC,TTTTAGAG,TTTTAGAT,TTTTAGAN,TTTTAGCA,TTTTAGCC,TTTTAGCG,TTTTAGCT,TTTTAGCN,TTTTAGGA,TTTTAGGC,TTTTAGGG,TTTTAGGT,TTTTAGGN,TTTTAGTA,TTTTAGTC,TTTTAGTG,TTTTAGTT,TTTTAGTN,TTTTAGNA,TTTTAGNC,TTTTAGNG,TTTTAGNT,TTTTAGNN,TTTTATAA,TTTTATAC,TTTTATAG,TTTTATAT,TTTTATAN,TTTTATCA,TTTTATCC,TTTTATCG,TTTTATCT,TTTTATCN,TTTTATGA,TTTTATGC,TTTTATGG,TTTTATGT,TTTTATGN,TTTTATTA,TTTTATTC,TTTTATTG,TTTTATTT,TTTTATTN,TTTTATNA,TTTTATNC,TTTTATNG,TTTTATNT,TTTTATNN,TTTTANAA,TTTTANAC,TTTTANAG,TTTTANAT,TTTTANAN,TTTTANCA,TTTTANCC,TTTTANCG,TTTTANCT,TTTTANCN,TTTTANGA,TTTTANGC,TTTTANGG,TTTTANGT,TTTTANGN,TTTTANTA,TTTTANTC,TTTTANTG,TTTTANTT,TTTTANTN,TTTTANNA,TTTTANNC,TTTTANNG,TTTTANNT,TTTTANNN,TTTTGAAA,TTTTGAAC,TTTTGAAG,TTTTGAAT,TTTTGAAN,TTTTGACA,TTTTGACC,TTTTGACG,TTTTGACT,TTTTGACN,TTTTGAGA,TTTTGAGC,TTTTGAGG,TTTTGAGT,TTTTGAGN,TTTTGATA,TTTTGATC,TTTTGATG,TTTTGATT,TTTTGATN,TTTTGANA,TTTTGANC,TTTTGANG,TTTTGANT,TTTTGANN,TTTTGCAA,TTTTGCAC,TTTTGCAG,TTTTGCAT,TTTTGCAN,TTTTGCCA,TTTTGCCC,TTTTGCCG,TTTTGCCT,TTTTGCCN,TTTTGCGA,TTTTGCGC,TTTTGCGG,TTTTGCGT,TTTTGCGN,TTTTGCTA,TTTTGCTC,TTTTGCTG,TTTTGCTT,TTTTGCTN,TTTTGCNA,TTTTGCNC,TTTTGCNG,TTTTGCNT,TTTTGCNN,TTTTGGAA,TTTTGGAC,TTTTGGAG,TTTTGGAT,TTTTGGAN,TTTTGGCA,TTTTGGCC,TTTTGGCG,TTTTGGCT,TTTTGGCN,TTTTGGGA,TTTTGGGC,TTTTGGGG,TTTTGGGT,TTTTGGGN,TTTTGGTA,TTTTGGTC,TTTTGGTG,TTTTGGTT,TTTTGGTN,TTTTGGNA,TTTTGGNC,TTTTGGNG,TTTTGGNT,TTTTGGNN,TTTTGTAA,TTTTGTAC,TTTTGTAG,TTTTGTAT,TTTTGTAN,TTTTGTCA,TTTTGTCC,TTTTGTCG,TTTTGTCT,TTTTGTCN,TTTTGTGA,TTTTGTGC,TTTTGTGG,TTTTGTGT,TTTTGTGN,TTTTGTTA,TTTTGTTC,TTTTGTTG,TTTTGTTT,TTTTGTTN,TTTTGTNA,TTTTGTNC,TTTTGTNG,TTTTGTNT,TTTTGTNN,TTTTGNAA,TTTTGNAC,TTTTGNAG,TTTTGNAT,TTTTGNAN,TTTTGNCA,TTTTGNCC,TTTTGNCG,TTTTGNCT,TTTTGNCN,TTTTGNGA,TTTTGNGC,TTTTGNGG,TTTTGNGT,TTTTGNGN,TTTTGNTA,TTTTGNTC,TTTTGNTG,TTTTGNTT,TTTTGNTN,TTTTGNNA,TTTTGNNC,TTTTGNNG,TTTTGNNT,TTTTGNNN,TTTTRAAA,TTTTRAAC,TTTTRAAG,TTTTRAAT,TTTTRAAN,TTTTRACA,TTTTRACC,TTTTRACG,TTTTRACT,TTTTRACN,TTTTRAGA,TTTTRAGC,TTTTRAGG,TTTTRAGT,TTTTRAGN,TTTTRATA,TTTTRATC,TTTTRATG,TTTTRATT,TTTTRATN,TTTTRANA,TTTTRANC,TTTTRANG,TTTTRANT,TTTTRANN,TTTTRCAA,TTTTRCAC,TTTTRCAG,TTTTRCAT,TTTTRCAN,TTTTRCCA,TTTTRCCC,TTTTRCCG,TTTTRCCT,TTTTRCCN,TTTTRCGA,TTTTRCGC,TTTTRCGG,TTTTRCGT,TTTTRCGN,TTTTRCTA,TTTTRCTC,TTTTRCTG,TTTTRCTT,TTTTRCTN,TTTTRCNA,TTTTRCNC,TTTTRCNG,TTTTRCNT,TTTTRCNN,TTTTRGAA,TTTTRGAC,TTTTRGAG,TTTTRGAT,TTTTRGAN,TTTTRGCA,TTTTRGCC,TTTTRGCG,TTTTRGCT,TTTTRGCN,TTTTRGGA,TTTTRGGC,TTTTRGGG,TTTTRGGT,TTTTRGGN,TTTTRGTA,TTTTRGTC,TTTTRGTG,TTTTRGTT,TTTTRGTN,TTTTRGNA,TTTTRGNC,TTTTRGNG,TTTTRGNT,TTTTRGNN,TTTTRTAA,TTTTRTAC,TTTTRTAG,TTTTRTAT,TTTTRTAN,TTTTRTCA,TTTTRTCC,TTTTRTCG,TTTTRTCT,TTTTRTCN,TTTTRTGA,TTTTRTGC,TTTTRTGG,TTTTRTGT,TTTTRTGN,TTTTRTTA,TTTTRTTC,TTTTRTTG,TTTTRTTT,TTTTRTTN,TTTTRTNA,TTTTRTNC,TTTTRTNG,TTTTRTNT,TTTTRTNN,TTTTRNAA,TTTTRNAC,TTTTRNAG,TTTTRNAT,TTTTRNAN,TTTTRNCA,TTTTRNCC,TTTTRNCG,TTTTRNCT,TTTTRNCN,TTTTRNGA,TTTTRNGC,TTTTRNGG,TTTTRNGT,TTTTRNGN,TTTTRNTA,TTTTRNTC,TTTTRNTG,TTTTRNTT,TTTTRNTN,TTTTRNNA,TTTTRNNC,TTTTRNNG,TTTTRNNT,TTTTRNNN,TTTYAAAA,TTTYAAAC,TTTYAAAG,TTTYAAAT,TTTYAAAN,TTTYAACA,TTTYAACC,TTTYAACG,TTTYAACT,TTTYAACN,TTTYAAGA,TTTYAAGC,TTTYAAGG,TTTYAAGT,TTTYAAGN,TTTYAATA,TTTYAATC,TTTYAATG,TTTYAATT,TTTYAATN,TTTYAANA,TTTYAANC,TTTYAANG,TTTYAANT,TTTYAANN,TTTYACAA,TTTYACAC,TTTYACAG,TTTYACAT,TTTYACAN,TTTYACCA,TTTYACCC,TTTYACCG,TTTYACCT,TTTYACCN,TTTYACGA,TTTYACGC,TTTYACGG,TTTYACGT,TTTYACGN,TTTYACTA,TTTYACTC,TTTYACTG,TTTYACTT,TTTYACTN,TTTYACNA,TTTYACNC,TTTYACNG,TTTYACNT,TTTYACNN,TTTYAGAA,TTTYAGAC,TTTYAGAG,TTTYAGAT,TTTYAGAN,TTTYAGCA,TTTYAGCC,TTTYAGCG,TTTYAGCT,TTTYAGCN,TTTYAGGA,TTTYAGGC,TTTYAGGG,TTTYAGGT,TTTYAGGN,TTTYAGTA,TTTYAGTC,TTTYAGTG,TTTYAGTT,TTTYAGTN,TTTYAGNA,TTTYAGNC,TTTYAGNG,TTTYAGNT,TTTYAGNN,TTTYATAA,TTTYATAC,TTTYATAG,TTTYATAT,TTTYATAN,TTTYATCA,TTTYATCC,TTTYATCG,TTTYATCT,TTTYATCN,TTTYATGA,TTTYATGC,TTTYATGG,TTTYATGT,TTTYATGN,TTTYATTA,TTTYATTC,TTTYATTG,TTTYATTT,TTTYATTN,TTTYATNA,TTTYATNC,TTTYATNG,TTTYATNT,TTTYATNN,TTTYANAA,TTTYANAC,TTTYANAG,TTTYANAT,TTTYANAN,TTTYANCA,TTTYANCC,TTTYANCG,TTTYANCT,TTTYANCN,TTTYANGA,TTTYANGC,TTTYANGG,TTTYANGT,TTTYANGN,TTTYANTA,TTTYANTC,TTTYANTG,TTTYANTT,TTTYANTN,TTTYANNA,TTTYANNC,TTTYANNG,TTTYANNT,TTTYANNN,TTTYGAAA,TTTYGAAC,TTTYGAAG,TTTYGAAT,TTTYGAAN,TTTYGACA,TTTYGACC,TTTYGACG,TTTYGACT,TTTYGACN,TTTYGAGA,TTTYGAGC,TTTYGAGG,TTTYGAGT,TTTYGAGN,TTTYGATA,TTTYGATC,TTTYGATG,TTTYGATT,TTTYGATN,TTTYGANA,TTTYGANC,TTTYGANG,TTTYGANT,TTTYGANN,TTTYGCAA,TTTYGCAC,TTTYGCAG,TTTYGCAT,TTTYGCAN,TTTYGCCA,TTTYGCCC,TTTYGCCG,TTTYGCCT,TTTYGCCN,TTTYGCGA,TTTYGCGC,TTTYGCGG,TTTYGCGT,TTTYGCGN,TTTYGCTA,TTTYGCTC,TTTYGCTG,TTTYGCTT,TTTYGCTN,TTTYGCNA,TTTYGCNC,TTTYGCNG,TTTYGCNT,TTTYGCNN,TTTYGGAA,TTTYGGAC,TTTYGGAG,TTTYGGAT,TTTYGGAN,TTTYGGCA,TTTYGGCC,TTTYGGCG,TTTYGGCT,TTTYGGCN,TTTYGGGA,TTTYGGGC,TTTYGGGG,TTTYGGGT,TTTYGGGN,TTTYGGTA,TTTYGGTC,TTTYGGTG,TTTYGGTT,TTTYGGTN,TTTYGGNA,TTTYGGNC,TTTYGGNG,TTTYGGNT,TTTYGGNN,TTTYGTAA,TTTYGTAC,TTTYGTAG,TTTYGTAT,TTTYGTAN,TTTYGTCA,TTTYGTCC,TTTYGTCG,TTTYGTCT,TTTYGTCN,TTTYGTGA,TTTYGTGC,TTTYGTGG,TTTYGTGT,TTTYGTGN,TTTYGTTA,TTTYGTTC,TTTYGTTG,TTTYGTTT,TTTYGTTN,TTTYGTNA,TTTYGTNC,TTTYGTNG,TTTYGTNT,TTTYGTNN,TTTYGNAA,TTTYGNAC,TTTYGNAG,TTTYGNAT,TTTYGNAN,TTTYGNCA,TTTYGNCC,TTTYGNCG,TTTYGNCT,TTTYGNCN,TTTYGNGA,TTTYGNGC,TTTYGNGG,TTTYGNGT,TTTYGNGN,TTTYGNTA,TTTYGNTC,TTTYGNTG,TTTYGNTT,TTTYGNTN,TTTYGNNA,TTTYGNNC,TTTYGNNG,TTTYGNNT,TTTYGNNN,TTTYRAAA,TTTYRAAC,TTTYRAAG,TTTYRAAT,TTTYRAAN,TTTYRACA,TTTYRACC,TTTYRACG,TTTYRACT,TTTYRACN,TTTYRAGA,TTTYRAGC,TTTYRAGG,TTTYRAGT,TTTYRAGN,TTTYRATA,TTTYRATC,TTTYRATG,TTTYRATT,TTTYRATN,TTTYRANA,TTTYRANC,TTTYRANG,TTTYRANT,TTTYRANN,TTTYRCAA,TTTYRCAC,TTTYRCAG,TTTYRCAT,TTTYRCAN,TTTYRCCA,TTTYRCCC,TTTYRCCG,TTTYRCCT,TTTYRCCN,TTTYRCGA,TTTYRCGC,TTTYRCGG,TTTYRCGT,TTTYRCGN,TTTYRCTA,TTTYRCTC,TTTYRCTG,TTTYRCTT,TTTYRCTN,TTTYRCNA,TTTYRCNC,TTTYRCNG,TTTYRCNT,TTTYRCNN,TTTYRGAA,TTTYRGAC,TTTYRGAG,TTTYRGAT,TTTYRGAN,TTTYRGCA,TTTYRGCC,TTTYRGCG,TTTYRGCT,TTTYRGCN,TTTYRGGA,TTTYRGGC,TTTYRGGG,TTTYRGGT,TTTYRGGN,TTTYRGTA,TTTYRGTC,TTTYRGTG,TTTYRGTT,TTTYRGTN,TTTYRGNA,TTTYRGNC,TTTYRGNG,TTTYRGNT,TTTYRGNN,TTTYRTAA,TTTYRTAC,TTTYRTAG,TTTYRTAT,TTTYRTAN,TTTYRTCA,TTTYRTCC,TTTYRTCG,TTTYRTCT,TTTYRTCN,TTTYRTGA,TTTYRTGC,TTTYRTGG,TTTYRTGT,TTTYRTGN,TTTYRTTA,TTTYRTTC,TTTYRTTG,TTTYRTTT,TTTYRTTN,TTTYRTNA,TTTYRTNC,TTTYRTNG,TTTYRTNT,TTTYRTNN,TTTYRNAA,TTTYRNAC,TTTYRNAG,TTTYRNAT,TTTYRNAN,TTTYRNCA,TTTYRNCC,TTTYRNCG,TTTYRNCT,TTTYRNCN,TTTYRNGA,TTTYRNGC,TTTYRNGG,TTTYRNGT,TTTYRNGN,TTTYRNTA,TTTYRNTC,TTTYRNTG,TTTYRNTT,TTTYRNTN,TTTYRNNA,TTTYRNNC,TTTYRNNG,TTTYRNNT,TTTYRNNN,";
$hash_of_pams_5p{'II'}{'II-B'}{'II-B-CNCCN'}="CACCA,CACCC,CACCG,CACCT,CCCCA,CCCCC,CCCCG,CCCCT,CGCCA,CGCCC,CGCCG,CGCCT,CTCCA,CTCCC,CTCCG,CTCCT,";
$hash_of_pams_5p{'II'}{'II-B'}{'II-B-CCN'}="CCA,CCC,CCG,CCT,";

my $lookup_path = "NA"; #"archaea205.contig.csv";
my $modify_spacer_org = 0;
my $repeat_class_lookup = "NA";
my $update_crispr_class = 0;
my $guess_url = 0;

my $ct_tab_neg = "NA";
my $do_pval = 0;
my $min_score = 0;

my $ind = 0;
foreach(@ARGV) {
	if (@ARGV[$ind] eq '-tab') {
		$ct_tab_in = @ARGV[$ind + 1];
		if (! (-e $ct_tab_in)) {
			die "cannot open file: " . $ct_tab_in . "\n";
		}
	}
	
	if (@ARGV[$ind] eq '-out') {
		$out = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-scoring_flank_length_5p') {
		$scoring_flank_length_5p = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-scoring_flank_length_3p') {
		$scoring_flank_length_3p = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-wt_psp_match') {
		$wt_psp_match = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-wt_psp_mismatch') {
		$wt_psp_mismatch = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-wt_fp5p_match') {
		$wt_fp5p_match = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-wt_fp5p_mismatch') {
		$wt_fp5p_mismatch = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-wt_fp3p_match') {
		$wt_fp5p_match = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-wt_fp3p_mismatch') {
		$wt_fp5p_mismatch = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-show_alignment_string') {
		$show_alignment_string = 1;
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
		$pam_search_super = 1;
		$pam_search_exact = 1;
	}
	
	if (@ARGV[$ind] eq '-pam_crispr_class_super') {
		$pam_crispr_class_super = @ARGV[$ind + 1];
		$pam_search_all = 0;
		$pam_search_super = 1;
		$pam_search_exact = 0;
	}
	
	if (@ARGV[$ind] eq '-pam_crispr_class_exact') {
		$pam_crispr_class_exact = @ARGV[$ind + 1];
		$pam_search_all = 0;
		$pam_search_super = 1;
		$pam_search_exact = 1;
	}
	
	if (@ARGV[$ind] eq '-pam_user_5p') {
		$pam_user_5p = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-pam_user_3p') {
		$pam_user_3p = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-pam_score_wrong_class') {
		$pam_score_wrong_class = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-pam_score_super_class') {
		$pam_score_super_class = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-pam_score_exact_class') {
		$pam_score_exact_class = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-pam_score_user_pam') {
		$pam_score_user_pam = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-score_for_pam') {
		$score_for_pam = 1;
	}
	
	if (@ARGV[$ind] eq '-chk_range_start') {
		$chk_range_start = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-chk_range_end') {
		$chk_range_end = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-chk_range_start_excl') {
		$chk_range_start_excl = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-chk_range_end_excl') {
		$chk_range_end_excl = @ARGV[$ind + 1];
	}
	
	if (@ARGV[$ind] eq '-flag_unrealistic_target') {
		$flag_unrealistic_target = 1;
	}
	
	if (@ARGV[$ind] eq '-spacer_org_lookup') {
		$lookup_path = @ARGV[$ind + 1];
		$modify_spacer_org = 1;
	}
	
	if ($ARGV[$ind] eq '-repeat_lookup') {
		$repeat_class_lookup =  $ARGV[$ind + 1];
		$update_crispr_class = 1;
	}
	
	if ($ARGV[$ind] eq '-guess_url') {
		$guess_url = 1;
	}
	
	if ($ARGV[$ind] eq "-tab_neg") {
		$ct_tab_neg = $ARGV[$ind + 1];
	}
	
	if ($ARGV[$ind] eq "-compute_pval") {
		$do_pval = 1;
	}
	
	if (@ARGV[$ind] eq '-min_score') {
		$min_score = @ARGV[$ind + 1];
	}
	
	$ind++;
}

if (-e $out) {
	unlink($out);
}

if ($pam_user_5p ne "" and $pam_user_5p ne "NA") {
	$hash_of_pams_5p{'USER'}{'USER'}{'USER_PAM'} = $pam_user_5p;
}

if ($pam_user_3p ne "" and $pam_user_3p ne "NA") {
	$hash_of_pams_3p{'USER'}{'USER'}{'USER_PAM'} = $pam_user_3p;
}

my %source_def_lookup = ();
my @lookup_contents = ();

if (-e $lookup_path) {
	open(LOOKUP, $lookup_path);
	@lookup_contents = <LOOKUP>;
	close(LOOKUP);
	
	for (my $i = 1; $i < scalar(@lookup_contents); $i++) {
		my $line = $lookup_contents[$i];
		$line =~ s/\n//g;
		
		my @toks = split(",", $line);
		my $assm = $toks[0];
		my $acc = $toks[1];
		
		my $nn = scalar(@toks);
		my $org = $toks[$nn - 1];
		
		$org =~ s/\s/_/g;
		my $desc = "$acc|$assm|$org";
		
		if (! exists $source_def_lookup{$acc}) {
			$source_def_lookup{$acc} = $desc;
		}
	}
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

my %neg_score_dist = ();
my %neg_count_total = ();
my %neg_score_dist_by_gc_len = ();
my %neg_count_total_by_gc_len = ();

if (-e $ct_tab_neg and $do_pval > 0) {
	open(NEG, $ct_tab_neg);
	my @contents = <NEG>;
	close(NEG);
	
	for (my $i = 0; $i < scalar(@contents); $i++) {
		my $line = $contents[$i];
		$line =~ s/\n//g;
		
		my @toks = split(/[\t]/, $line);
		my $fp_score = $toks[7];
		my $spacer_seq = $toks[14];
		my $spacer_id = $toks[16];
		
		if ($fp_score < -55) {
			$fp_score = -55;
		}
		
		if ($fp_score > 55) {
			$fp_score = 55;
		}
		
		my $spacer_len_cat = int((length($spacer_seq) - 1) / 10.0) + 1;		
		my $gc_perc = gc_percent_int($spacer_seq);
		
		if (exists $neg_count_total_by_gc_len{$spacer_len_cat}{$gc_perc}) {
			$neg_count_total_by_gc_len{$spacer_len_cat}{$gc_perc} += 1;
		} else {
			$neg_count_total_by_gc_len{$spacer_len_cat}{$gc_perc} = 1;
		}
		
		if (exists $neg_count_total{$spacer_id}) {
			$neg_count_total{$spacer_id} += 1;
		} else {
			$neg_count_total{$spacer_id} = 1;
		}
		
		for (my $j = 0; $j < 111; $j++) {
			if ($fp_score >= ($j - 55)) {
				if (exists $neg_score_dist{$spacer_id}{$j}) {
					$neg_score_dist{$spacer_id}{$j} += 1;
				} else {
					$neg_score_dist{$spacer_id}{$j} = 1;
				}
				
				if (exists $neg_score_dist_by_gc_len{$spacer_len_cat}{$gc_perc}{$j}) {
					$neg_score_dist_by_gc_len{$spacer_len_cat}{$gc_perc}{$j} += 1;
				} else {
					$neg_score_dist_by_gc_len{$spacer_len_cat}{$gc_perc}{$j} = 1;
				}
			} else {
				if (exists $neg_score_dist{$spacer_id}{$j}) {
					$neg_score_dist{$spacer_id}{$j} += 0;
				} else {
					$neg_score_dist{$spacer_id}{$j} = 0;
				}
				
				if (exists $neg_score_dist_by_gc_len{$spacer_len_cat}{$gc_perc}{$j}) {
					$neg_score_dist_by_gc_len{$spacer_len_cat}{$gc_perc}{$j} += 0;
				} else {
					$neg_score_dist_by_gc_len{$spacer_len_cat}{$gc_perc}{$j} = 0;
				}
			}
		}
	}
}

open(TXT, $ct_tab_in);
my @contents = <TXT>;
close(TXT);

open(OUT, ">$out");
for (my $i = 0; $i < scalar(@contents); $i++) {
	my $line = $contents[$i];
	$line =~ s/\n//g;
		
	my @toks = split(/[\t]/, $line);
	
	my $spacer_seq = $toks[14]; 
	my $actual_protospacer_seq = $toks[15];
	
	my $spacer_flank_5p_seq = $toks[19];
	my $spacer_flank_3p_seq = $toks[20];
	my $protospacer_flank_5p_seq = $toks[6];
	my $protospacer_flank_3p_seq = $toks[5];
	
	if (scalar(@toks) >= 32) {
		$spacer_flank_5p_seq = $toks[28];
		$spacer_flank_3p_seq = $toks[29];
		$protospacer_flank_5p_seq = $toks[30];
		$protospacer_flank_3p_seq = $toks[31];
	}
	
	my $source_flank_length_3p = length($protospacer_flank_3p_seq);
	my $source_flank_length_5p = length($protospacer_flank_5p_seq);
	
	my $protospacer_flank_3p_seq_scoring = $protospacer_flank_3p_seq;
	if ($scoring_flank_length_3p < $source_flank_length_3p) {
		$protospacer_flank_3p_seq_scoring = substr($protospacer_flank_3p_seq, 0, $scoring_flank_length_3p);
	}
	
	my $protospacer_flank_5p_seq_scoring = $protospacer_flank_5p_seq;
	if ($scoring_flank_length_5p < $source_flank_length_5p) {
		$protospacer_flank_5p_seq_scoring = substr($protospacer_flank_5p_seq, -$scoring_flank_length_5p, $scoring_flank_length_5p);	
	}
	
	my $scoring_flank_3p = $spacer_flank_3p_seq;
	if ($scoring_flank_length_3p < $source_flank_length_3p) {
		$scoring_flank_3p = substr($spacer_flank_3p_seq, 0, $scoring_flank_length_3p);
	}
	
	my $scoring_flank_5p = $spacer_flank_5p_seq;
	if ($scoring_flank_length_5p < $source_flank_length_5p) {
		$scoring_flank_5p = substr($spacer_flank_5p_seq, -$scoring_flank_length_5p, $scoring_flank_length_5p);
	}
	
	my $actual_protospacer_seq_rev = reverse($actual_protospacer_seq);
	my $protospacer_flank_5p_seq_scoring_rev = reverse($protospacer_flank_5p_seq_scoring);
	my $protospacer_flank_3p_seq_scoring_rev = reverse($protospacer_flank_3p_seq_scoring);
	my ($count_match, $count_mismatch, $count_nn, $count_indel, $cmpstr) = compareAlignedSeqInRC($spacer_seq, $actual_protospacer_seq_rev);
	my ($count_match_5p, $count_mismatch_5p, $count_nn_5p, $count_indel_5p, $cmpstr_5p) = compareAlignedSeqInRC($scoring_flank_5p, $protospacer_flank_5p_seq_scoring_rev);
	my ($count_match_3p, $count_mismatch_3p, $count_nn_3p, $count_indel_3p, $cmpstr_3p) = compareAlignedSeqInRC($scoring_flank_3p, $protospacer_flank_3p_seq_scoring_rev);
	
	my $score = ($wt_psp_match * $count_match) + ($wt_psp_mismatch * $count_mismatch) + ($wt_psp_mismatch * $count_indel) + 
	            ($wt_fp5p_match * $count_match_5p) + ($wt_fp5p_mismatch * $count_mismatch_5p) + ($wt_fp5p_mismatch * $count_indel_5p) +
	            ($wt_fp3p_match * $count_match_3p) + ($wt_fp3p_mismatch * $count_mismatch_3p) + ($wt_fp3p_mismatch * $count_indel_3p);
	
	my $self_match = 0;
	my $num_mismatch = $count_mismatch;
	my $count_mismatch = $count_mismatch + $count_mismatch_5p + $count_mismatch_3p;
	if ($count_mismatch < 1) {
		$self_match = 1;
	}
	
	$toks[13] = $self_match;
	if ($show_alignment_string > 0) {
		$toks[27] = $cmpstr;
	}
	
	my $dr_seq = $toks[22];
	if ($update_crispr_class > 0) {
		if (exists $array_repeat_lst{$dr_seq}) {
			$toks[21] = $array_repeat_lst{$dr_seq};
		}
	}
	
	my $pval = 0;
	my $the_spacer = $toks[16];
	my $score_ind = $toks[7] + 55;
	
	if ($score_ind < 0) {
		$score_ind = 0;
	}
	
	if ($score_ind > 110) {
		$score_ind = 110;
	}
	
	my $x = 0;
	my $y = 0;
	if (exists $neg_count_total{$the_spacer} and exists $neg_score_dist{$the_spacer}{$score_ind} and -e $ct_tab_neg and $do_pval > 0) {
		$x = $neg_score_dist{$the_spacer}{$score_ind};
		$y = $neg_count_total{$the_spacer};
		$pval = (($x + 0.0) / ($y + 0.0));
	}
	
	if ($y < 20) {
		my $spacer_len_cat = int((length($spacer_seq) - 1) / 10.0) + 1;		
		my $gc_perc = gc_percent_int($spacer_seq);
		if (exists $neg_count_total_by_gc_len{$spacer_len_cat}{$gc_perc} and exists $neg_score_dist_by_gc_len{$spacer_len_cat}{$gc_perc}{$score_ind} and -e $ct_tab_neg and $do_pval > 0) {
			$x = $neg_score_dist_by_gc_len{$spacer_len_cat}{$gc_perc}{$score_ind};
			$y = $neg_count_total_by_gc_len{$spacer_len_cat}{$gc_perc};
			
			if ($y >= 20) {
				$pval = (($x + 0.0) / ($y + 0.0));
			}
		}
	}
	
	my $offending_spacer = 0;
	if (exists $neg_score_dist{$the_spacer}{75}) {
		if ($neg_score_dist{$the_spacer}{75} > 0) {
			$offending_spacer = 1;
		}
	}
	
	my $proto_ori = $toks[8];
	my $crispr_class = $toks[21];
	my @crispr_class_toks = split("-", $crispr_class);
	my $super_type = $crispr_class_toks[0];
	if ($pam_crispr_class_super ne "all" and $pam_crispr_class_super ne "") {
		$super_type = $pam_crispr_class_super;
	}
	
	if ($pam_crispr_class_exact ne "all" and $pam_crispr_class_exact ne "") {
		$crispr_class = $pam_crispr_class_exact;
	}
	
	my $pam_3p_found = "NA";
	my $pam_3p_seq = "NA";
	my $pam_3p_wrong = 0;
	my $pam_3p_super = 0;
	my $pam_3p_exact = 0;
	my $pam_3p_user = 0;
	foreach my $key (keys(%hash_of_pams_3p)) {
		my %subtypes = %{$hash_of_pams_3p{$key}};
		if ($pam_3p_found ne "NA") {
			last;
		}
		if($pam_search_super == 1 and $key ne $super_type and $pam_search_all == 0 and $key ne "USER") {
			next;
		}
		foreach my $sub (keys(%subtypes)) {
			my %subsubtypes = %{$subtypes{$sub}};
			if ($pam_3p_found ne "NA") {
				last;
			}
			if($pam_search_exact == 1 and $sub ne $crispr_class and $pam_search_all == 0 and $sub ne "USER") {
				next;
			}
			foreach my $tag (keys(%subsubtypes)) {
				my $pams = $subsubtypes{$tag};
				my @all_pams = split(",", $pams);
				if ($pam_3p_found ne "NA") {
					last;
				}
				for (my $j = 0; $j < scalar(@all_pams); $j++) {
					my $the_pam = $all_pams[$j];
					if ($pam_3p_found ne "NA") {
						last;
					}
					if (length($the_pam) > 0) {
						if ($protospacer_flank_3p_seq_scoring =~ /^$the_pam/i) {
							$pam_3p_found = $tag;
							$pam_3p_seq = $the_pam;
							if ($sub eq $crispr_class) {
								$pam_3p_exact += 1;
							} elsif ($key eq $super_type) {
								$pam_3p_super += 1;
							} elsif ($sub eq "USER" or $key eq "USER") {
								$pam_3p_user += 1;
							} else {
								$pam_3p_wrong += 1;
							}
							last;
						}
					}
				}
			}
		}
	}
	
	my $pam_5p_found = "NA";
	my $pam_5p_seq = "NA";
	my $pam_5p_wrong = 0;
	my $pam_5p_super = 0;
	my $pam_5p_exact = 0;
	my $pam_5p_user = 0;
	foreach my $key (keys(%hash_of_pams_5p)) {
		my %subtypes = %{$hash_of_pams_5p{$key}};
		if ($pam_5p_found ne "NA") {
			last;
		}
		if($pam_search_super == 1 and $key ne $super_type and $pam_search_all == 0 and $key ne "USER") {
			next;
		}
		foreach my $sub (keys(%subtypes)) {
			my %subsubtypes = %{$subtypes{$sub}};
			if ($pam_5p_found ne "NA") {
				last;
			}
			if($pam_search_exact == 1 and $sub ne $crispr_class and $pam_search_all == 0 and $sub ne "USER") {
				next;
			}
			foreach my $tag (keys(%subsubtypes)) {
				my $pams = $subsubtypes{$tag};
				my @all_pams = split(",", $pams);
				if ($pam_5p_found ne "NA") {
					last;
				}
				for (my $k = 0; $k < scalar(@all_pams); $k++) {
					my $the_pam = $all_pams[$k];
					if ($pam_5p_found ne "NA") {
						last;
					}
					if (length($the_pam) > 0) {
						if ($protospacer_flank_5p_seq_scoring =~ /$the_pam$/i) {
							$pam_5p_found = $tag;
							$pam_5p_seq = $the_pam;
							if ($sub eq $crispr_class) {
								$pam_5p_exact += 1;
							} elsif ($key eq $super_type) {
								$pam_5p_super += 1;
							} elsif ($sub eq "USER" or $key eq "USER") {
								$pam_5p_user += 1;
							} else {
								$pam_5p_wrong += 1;
							}
							last;
						}
					}
				}
			}
		}
	}
	
	if ($proto_ori eq "+" and $pam_3p_found ne "NA") {
		$toks[9] = $pam_3p_seq;
		$toks[25] = $pam_3p_found;
	}
	
	if ($proto_ori eq "+" and $pam_5p_found ne "NA") {
		$toks[10] = $pam_5p_seq;
		$toks[23] = $pam_5p_found;
	}
	
	if ($proto_ori eq "-" and $pam_3p_found ne "NA") {
		$toks[11] = $pam_3p_seq;
		$toks[26] = $pam_3p_found;
	}
	
	if ($proto_ori eq "-" and $pam_5p_found ne "NA") {
		$toks[12] = $pam_5p_seq;
		$toks[24] = $pam_5p_found;
	}
	
	my $pam_score = ($pam_score_wrong_class * ($pam_5p_wrong + $pam_3p_wrong)) + 
	                ($pam_score_super_class * ($pam_5p_super + $pam_3p_super)) + 
	                ($pam_score_exact_class * ($pam_5p_exact + $pam_3p_exact)) + 
	                ($pam_score_user_pam * ($pam_5p_user + $pam_3p_user));
	
	$score += ($score_for_pam * $pam_score);
	                
	my $b = 1;
	if ($flag_unrealistic_target > 0) {
		my ($start, $end) = (-67676767, -76767676);
		if ($chk_range_start =~ /\d+\-\d+/) {
			($start, $end) = ($chk_range_start =~ /(\d+)\-(\d+)/);
			if ($start > $end) {
				my $temp = $start;
				$start = $end;
				$end = $temp;
			}
			$start = $start - 1;
			$end = $end - 1;
		}
	
		my ($e_start, $e_end) = (-91919191, -19191919);
		if (($chk_range_end =~ /\d+\-\d+/)) {
			($e_start, $e_end) = ($chk_range_end =~ /(\d+)\-(\d+)/);
			$e_start = length($spacer_seq) - $e_start;
			$e_end = length($spacer_seq) - $e_end;
			if ($e_start > $e_end) {
				my $temp = $e_start;
				$e_start = $e_end;
				$e_end = $temp;
			}
		}
	
		my %excl = ();
		my @excl_start = split(",", $chk_range_start_excl);
		for (my $j = 0; $j < scalar(@excl_start); $j++) {
			if ($excl_start[$j] =~ /\d+/) {
				my $ind = $excl_start[$j] - 1;
				$excl{$ind} = 1;
			}
		}
	
		my @excl_end = split(",", $chk_range_end_excl);
		for (my $j = 0; $j < scalar(@excl_end); $j++) {
			if ($excl_start[$j] =~ /\d+/) {
				my $ind = length($spacer_seq) - $excl_end[$j];
				$excl{$ind} = 1;
			}
		}
	
		if ($start >= 0 and $end >= 0 and $e_start >= 0 and $e_end >= 0) {
			for (my $j = 0; $j < length($cmpstr); $j++) {
				my $ch = substr($cmpstr, $j, 1);
		
				if (($j >= $start and $j <= $end) or ($j >= $e_start and $j <= $e_end)) {
					if ($ch eq "+") {
						$b *= 1;
					} else {
						if (exists $excl{$j}) {
							$b *= 1;
						} else {
							$b *= 0;
						}
					}
				} else {
					$b *= 1;
				}
			}
		}
	}
	$b = 1 - $b;
	
	$toks[5] = $protospacer_flank_3p_seq_scoring;
	$toks[6] = $protospacer_flank_5p_seq_scoring;
	$toks[7] = $score;
	
	my $psp_id = $toks[2];
	my $psp_start = $toks[3];
	my $psp_end = $toks[4];
	my $targeted_org = $toks[17];
	
	my $spacer_flank_length_5p = $toks[5];
	my $spacer_flank_length_3p = $toks[6];
	
	if ($guess_url > 0) {
		my ($url_new, $a_text) = guessURL($targeted_org, $psp_id, $psp_start, $psp_end, $spacer_flank_length_5p, $spacer_flank_length_3p);
		$toks[18] = $url_new;
	}
	
	if ($modify_spacer_org > 0) {
		my $spacer_acc = $toks[0];
		$toks[16] = $spacer_acc . "|Unknown|Unknown";
		if (exists $source_def_lookup{$spacer_acc}) {
			$toks[16] = $source_def_lookup{$spacer_acc};
		}
	}
	
	if (-e $ct_tab_neg and $do_pval > 0) {
		$toks[27] = $pval;
	} 
	
	my $strout = ""; #join("\t", @toks);
	if ($flag_unrealistic_target > 0) {
		if (scalar(@toks) >= 33) {
			$toks[32] = $b;
			$strout = join("\t", @toks);
		} elsif (scalar(@toks) == 29) {
			$toks[28] = $b;
			$strout = join("\t", @toks);
		} else {
			$strout = join("\t", @toks);
			$strout .= "\t" . $b;
		}
	} else {
		if (scalar(@toks) >= 33) {
			if (-e $ct_tab_neg and $do_pval > 0) {
				$toks[32] = $num_mismatch;
			}
			$strout = join("\t", @toks);
		} elsif (scalar(@toks) == 29) {
			if (-e $ct_tab_neg and $do_pval > 0) {
				$toks[28] = $num_mismatch;
			}
			$strout = join("\t", @toks);
		} else {
			$strout = join("\t", @toks);
			if (-e $ct_tab_neg and $do_pval > 0) {
				$strout .= "\t" . $num_mismatch;
			}
		}
	}
	#$offending_spacer
	if ($score > $min_score) {
		print OUT "$strout\n";
	}
}
close(OUT);

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

sub gc_percent_int {
	my ($seq1) = @_;
	
	my @arr_seq1=split('',$seq1);
	
	my %hash_of_bases;
	$hash_of_bases{'A'}=0;
	$hash_of_bases{'C'}=0;
	$hash_of_bases{'G'}=0;
	$hash_of_bases{'T'}=0;
	$hash_of_bases{'U'}=0;
	$hash_of_bases{'N'}=0;
	
	my $nbase = scalar(@arr_seq1);
	for (my $i = 0; $i < $nbase; $i++) {
		if (exists $hash_of_bases{$arr_seq1[$i]}) {
			$hash_of_bases{$arr_seq1[$i]} += 1;
		}
	}
	
	my $gc_tot = ($hash_of_bases{'C'} + $hash_of_bases{'G'} + 0.0);
	my $gc_perc = ($gc_tot / $nbase) * 100.0;
	$gc_perc = int($gc_perc);
	
	return($gc_perc);
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
