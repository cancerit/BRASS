diff --git a/perl/bin/filter_with_microbes_and_remapping.pl b/perl/bin/filter_with_microbes_and_remapping.pl
index 4ada730..b4ab412 100644
--- a/perl/bin/filter_with_microbes_and_remapping.pl
+++ b/perl/bin/filter_with_microbes_and_remapping.pl
@@ -74,6 +74,7 @@ my $TMP_DIR = q{};
 my $SCORE_ALG = 'ssearch36';
 my $SEARCH_CORES = 1;
 my $MIN_SUPPORT = 3;
+my $NO_CLEAN = q{};
 GetOptions(
   'blat_mapping_threshold=i' => \$BLAT_MAPPING_THRESHOLD,
   'groups_file=s' => \$GROUPS_FILE,
@@ -90,8 +91,11 @@ GetOptions(
   'score_alg=s' => \$SCORE_ALG,
   'search_cores=i' => \$SEARCH_CORES,
   'min_support=s' => \$MIN_SUPPORT,
+  'noclean' => \$NO_CLEAN
 );
 
+
+
 my %score_method = (emboss => \&pairwise_align_scores_emboss,
                     ssearch36 => \&pairwise_align_scores_ssearch36,
                     );
@@ -103,6 +107,11 @@ my $fa_file = $ARGV[2];
 my $bedpe_out;
 $bedpe_out = $ARGV[3] if(@ARGV == 4);
 
+print Dumper($NO_CLEAN);
+
+exit;
+
+
 if($TMP_DIR eq q{}) {
   undef $TMP_DIR;
 }
@@ -142,7 +151,9 @@ close $IN;
 # viral genomes.
 print STDERR "Performing abnormally paired reads remapping...\n";
 my (@low_end_remap_score, @high_end_remap_score, @low_end_footprint, @high_end_footprint);
-my ($all_reads_fh, $all_reads_file_name) = tempfile('allreads_XXXXXX', DIR => $TMP_DIR, UNLINK=>1);
+
+
+my ($all_reads_fh, $all_reads_file_name) = tempfile('allreads_XXXXXX', DIR => $TMP_DIR, UNLINK=>$NO_CLEAN);
 my %rg_id_of_read;
 
 for my $i(0..$#regions) {
diff --git a/perl/lib/Sanger/CGP/Brass/Implement.pm b/perl/lib/Sanger/CGP/Brass/Implement.pm
index f4d4dd1..70e3c25 100644
--- a/perl/lib/Sanger/CGP/Brass/Implement.pm
+++ b/perl/lib/Sanger/CGP/Brass/Implement.pm
@@ -388,14 +388,16 @@ sub filter {
   my $remap_file = $r_stub.'.r5';
   my $tumour_brm = File::Spec->catfile($tmp, sanitised_sample_from_bam($options->{'tumour'})).'.brm.bam';
   my $remap_micro = $^X.' '._which('filter_with_microbes_and_remapping.pl');
-  $remap_micro .= sprintf ' -virus_db %s -bacterial_db_stub %s -scores_output_file %s -tmpdir %s -score_alg %s -search_cores %d -groups_file %s',
+  $remap_micro .= sprintf ' -virus_db %s -bacterial_db_stub %s -scores_output_file %s -tmpdir %s -score_alg %s -search_cores %d -groups_file %s -noclean %s',
                           $options->{'viral'},
                           $options->{'microbe'},
                           $score_file,
                           File::Spec->catdir($tmp,'remap_micro'),
                           'ssearch36',
                           $options->{'threads'},
-                          $bedpe_no_head;
+                          $bedpe_no_head,
+                          $options->{'noclean'};
+
   $remap_micro .= sprintf ' %s %s %s %s',
                           $abs_bkp_file,
                           $tumour_brm,
