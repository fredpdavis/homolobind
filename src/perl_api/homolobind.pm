=head1 NAME

homolobind.pm

=head1 DESCRIPTION

This module contains routines to annotate binding sites by 
homology transfer. It uses domain assignments from SUPERFAMILY,
alignments from ASTRAL/ASTEROIDS, and template binding sites from
LIGBASE and PIBASE.

=head1 VERSION

1.1

=head1 AUTHOR

Fred P. Davis, HHMI-JFRC (davisf@janelia.hhmi.org)

=head1 LICENCE AND COPYRIGHT

Copyright 2009,2010 Fred P. Davis (davisf@janelia.hhmi.org).
See the file COPYING for copying permission.

This file is part of HOMOLOBIND.

HOMOLOBIND is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

HOMOLOBIND is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with HOMOLOBIND.  If not, see <http://www.gnu.org/licenses/>.

=head1 SUBROUTINES

=cut

package homolobind ;
use strict;
use warnings;

use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw/benchmark_homolobind run_homolobind/ ;

use homolobind::pibase ;
use homolobind::pilig ;
use homolobind::ASTRAL ;
use homolobind::SGE ;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;
use Bit::Vector ;

my $homolobind_specs ;
{

#SET CLUSTER SPECIFICATION: headnode, number of jobs, qstat polling frequency
   $homolobind_specs->{SGE}->{headnode} = 'login-eddy' ;
   $homolobind_specs->{SGE}->{numjobs} = 50 ;
   $homolobind_specs->{SGE}->{qstat_sleep} = 120 ;

#SET DATA DIRECTORIES:
   $homolobind_specs->{data_dir} =
      '/groups/eddy/home/davisf/work/homolobind/data_files' ;

   $homolobind_specs->{benchmark_dir} = $homolobind_specs->{data_dir}.
      '/benchmark' ;
   $homolobind_specs->{benchmark_results} = $homolobind_specs->{benchmark_dir}.
      '/benchmark_homolobind.20100907.out.gz' ;

   $homolobind_specs->{astral}->{dir} = $homolobind_specs->{data_dir}.'/astral';
   $homolobind_specs->{asteroids}->{dir} = $homolobind_specs->{data_dir}.
      '/asteroids' ;
   $homolobind_specs->{pibase}->{dir} = $homolobind_specs->{data_dir}.'/pibase';
   $homolobind_specs->{pilig}->{dir} = $homolobind_specs->{data_dir}.'/pilig';
   $homolobind_specs->{supfam}->{dir} = $homolobind_specs->{data_dir}.
      '/superfamily';
   $homolobind_specs->{scop}->{dir} = $homolobind_specs->{data_dir}.'/scop' ;

#SET DATA SOURCE VERSIONS
   $homolobind_specs->{scop}->{ver} = '1.73' ;
   $homolobind_specs->{astral}->{ver} = '1.73' ;
   $homolobind_specs->{asteroids}->{ver} = '1.73' ;

#SPECIFY URLS FOR REQUIRED EXTERNAL FILES
   $homolobind_specs->{urls}->{pibase} = { #actually also for pilig
      'http://research.janelia.org/davis/homolobind/'.
       'data_files/v1.1/homolobind_data_v1.1.tar.gz' => {
         destination => $homolobind_specs->{data_dir}}
   };


   $homolobind_specs->{urls}->{astral} = {
      'http://astral.berkeley.edu/scopseq-1.73/astral-scopdom-seqres-gd-all-1.73.fa' => {},
      'http://astral.berkeley.edu/scopseq-1.73/astral-scopdom-seqres-gd-sel-gs-bib-verbose-100-1.73.txt' => {},
   } ;

   $homolobind_specs->{urls}->{asteroids} = {
      'http://astral.berkeley.edu/asteroids-1.73/alignments/fam/' => {
         wget_options => '-np -r -l 2 -A aln -nH --cut-dirs=3 ',
         destination => $homolobind_specs->{asteroids}->{dir}.'/alignments/fam'
      },
      'http://astral.berkeley.edu/asteroids-1.73/alignments/sf/' => {
         wget_options => '-np -r -l 2 -A aln -nH --cut-dirs=3 ',
         destination => $homolobind_specs->{asteroids}->{dir}.'/alignments/sf'
      },
      'http://astral.berkeley.edu/asteroids-1.73/sequences/fam/' => {
         wget_options => '-np -r -l 2 -A fa -nH --cut-dirs=3 ',
         destination => $homolobind_specs->{asteroids}->{dir}.'/sequences/fam'
      },
      'http://astral.berkeley.edu/asteroids-1.73/sequences/sf/' => {
         wget_options => '-np -r -l 2 -A fa -nH --cut-dirs=3 ',
         destination => $homolobind_specs->{asteroids}->{dir}.'/sequences/sf'
      },
   } ;

   $homolobind_specs->{urls}->{scop} = {
      'http://scop.mrc-lmb.cam.ac.uk/scop/parse/dir.cla.scop.txt_1.73' => {} } ;


#SET SCOP data
   $homolobind_specs->{scop}->{cla_fn} = $homolobind_specs->{scop}->{dir}.
      '/dir.cla.scop.txt_1.73';

#SET ASTRAL SPECS
   $homolobind_specs->{astral}->{seqcl}->{100} =
      $homolobind_specs->{astral}->{dir}.
      '/astral-scopdom-seqres-gd-sel-gs-bib-verbose-100-'.
      $homolobind_specs->{astral}->{ver}.'.txt';

   $homolobind_specs->{astral}->{gd_seq} =
      $homolobind_specs->{astral}->{dir}.
      '/astral-scopdom-seqres-gd-all-'.
      $homolobind_specs->{astral}->{ver}.'.fa' ;

#SET ASTEROIDS
   $homolobind_specs->{asteroids}->{fam_aln} =
      $homolobind_specs->{asteroids}->{dir}.'/alignments/fam';
   $homolobind_specs->{asteroids}->{sf_aln}=
      $homolobind_specs->{asteroids}->{dir}.'/alignments/sf';
   $homolobind_specs->{asteroids}->{fam_seq} =
      $homolobind_specs->{asteroids}->{dir}.'/sequences/fam' ;
   $homolobind_specs->{asteroids}->{sf_seq} =
      $homolobind_specs->{asteroids}->{dir}.'/sequences/sf' ;

#SET PIBASE table locations
   $homolobind_specs->{pibase}->{tod_dir} =
      $homolobind_specs->{pibase}->{dir}.'/tod' ;

#SET PIBASE/LIGBASE data
   $homolobind_specs->{pilig}->{outfiles}->{liginfo} =
      $homolobind_specs->{pilig}->{dir}."/liginfo.out";
   $homolobind_specs->{pilig}->{outfiles}->{calc_sid_sasa} =
      $homolobind_specs->{pilig}->{dir}.'/calc_sid_sasa.out' ;
   $homolobind_specs->{pilig}->{outfiles}->{assign_pepnuci} =
      $homolobind_specs->{pilig}->{dir}.'/assign_pepnuci.out' ;
   $homolobind_specs->{pilig}->{outfiles}->{assign_exp} =
      $homolobind_specs->{pilig}->{dir}.'/assign_exp.out' ;

   $homolobind_specs->{pilig}->{outfiles}->{assign_pi_clusters} =
      $homolobind_specs->{pilig}->{dir}.'/assign_pi_clusters.out' ;

   $homolobind_specs->{pilig}->{outfiles}->{assign_lig_clusters} =
      $homolobind_specs->{pilig}->{dir}.'/assign_lig_clusters.out' ;
   $homolobind_specs->{pilig}->{outfiles}->{assign_pepnuci_clusters} =
      $homolobind_specs->{pilig}->{dir}.'/assign_pepnuci_clusters.out' ;

#SET SUPFAM data location
   $homolobind_specs->{supfam}->{selfhits_dir} =
      $homolobind_specs->{supfam}->{dir}.'/self_hits' ;
   $homolobind_specs->{supfam}->{fxn} =
      $homolobind_specs->{supfam}->{dir}.'/scop.annotation.1.73.txt' ;
   $homolobind_specs->{supfam}->{fxn_categories} =
      $homolobind_specs->{supfam}->{dir}.'/scop.larger.categories' ;

#SET BENCHMARK SPECS
   $homolobind_specs->{benchmark}->{num_negatives} = 10000 ;

   $homolobind_specs->{fpr_cutoffs} = [0.01, 0.02, 0.05] ;
   $homolobind_specs->{bb_replicates} = 500 ;

}


=head2 run_homolobind()

   Title:       run_homolobind()
   Function:    Maps PIBASE/LIGBASE binding sites onto target sequences
                  annotated with SUPERFAMILY domain assignments 

   Args:        ->{ARGV} = ARGV array reference; parsed to provide:
                ->{ass_fn} = name of SUPERFAMILY domain assignment file
                ->{out_fn} = name of output file
                ->{err_fn} = name of error file
                ->{matrix_fn} = optional substitution matrix, default BLOSUM62
                ->{cluster_fl} = run on an SGE cluster (options in SGE.pm)

   Returns:     NOTHING
   Displays:     1. seq_id
                 2. res_range
                 3. classtype
                 4. class
                 5. bs_type
                 6. bs_template
                 7. partner_descr
                 8. residues
                 9. bs_percseqident
                10. bs_percseqsim
                11. bs_numident
                12. bs_numaln
                13. bs_tmpl_numres
                14. bs_fracaln
                15. wholedom_numident
                16. wholedom_aln_length
                17. wholedom_percseqident
                18. wholedom_percseqsim

=cut

sub run_homolobind {

   my $in = shift;

# Read in command line options
   my $j = 0; my $user_thresh = {};
   while ($j < $#{$in->{ARGV}}) {
      $ARGV[$j] =~ s/^\-// ;
      if ($ARGV[$j] =~ /^thresh/) {
         $user_thresh->{$in->{ARGV}->[$j]} = $in->{ARGV}->[($j+1)] ;
      } else {
         $in->{$in->{ARGV}->[$j]} = $in->{ARGV}->[($j+1)] ;
      }
      $j += 2 ;
   }

# If a (non-run-homolobind) runmode is specified, call appropriate routine
   if (exists $in->{summarize_results}) {
      if (exists $in->{withR} && $in->{withR} == 1) {
         summarize_results_withR({ results_fn => $in->{summarize_results} });
      } else{
         summarize_results_noR({ results_fn => $in->{summarize_results} });
      }
      exit;
   } elsif (exists $in->{make_figures}) {
      make_paper_figures($in) ;
      exit;
   } elsif (exists $in->{plot_annotations}) {
      plot_annotations($in) ;
      exit;
   } elsif (exists $in->{fetch_data}) {
      fetch_data() ;
      exit;
   }


# Set parameters
   my $supfam_specs = set_supfam_specs() ;

# Set usage
   my $usage = "homolobind.pl -ass_fn SUPERFAMILY_assignment_file
 [-cluster_fl 1] to run on an SGE cluster
 [-out_fn output_file] [-err_fn error_file]
 [-matrix_fn substitution_matrix_file] matblas format substitution matrix\n";
# [-thresh_min_L_mw minimum_ligand_MW]\tdefault: 250
# [-thresh_max_L_mw maximum_ligand_MW]\tdefault: 1000
# [-thresh_L_bs_seqid N] ligand binding site sequence identity threshold,
#    default to a bs-specific threshold benchmarked to achieve 1% FPR
# [-thresh_P_bs_seqid N] domain binding site sequence identity threshold,
# [-thresh_p_bs_seqid N] peptide binding site sequence identity threshold;

   my @result_headers = qw/seq_id res_range classtype class bs_type bs_template partner_descr residues bs_percseqident bs_percseqsim bs_numident bs_numaln bs_tmpl_numres bs_fracaln wholedom_numident wholedom_aln_length wholedom_percseqident wholedom_percseqsim/ ;

# Load binding site threholds
   print STDERR "Loading template binding site thresholds: " ;
   my $benchmark_thresholds = readin_benchmark_thresholds({
      fn => $homolobind_specs->{benchmark_results} }) ;
   print STDERR "X\n" ;

# Read in matrix file if user specified
   if (exists $in->{matrix_fn}) {
      $supfam_specs->{substitution_matrix} =
         readin_substitution_matrix({fn => $in->{matrix_fn}}) ;
   }

# Make sure SUPERFAMILY assignment file specified
   if (!exists $in->{ass_fn} || !-s $in->{ass_fn}) {
      die "ERROR: SUPFAM annotation file not found (-ass_fn)\n".
          "\nUSAGE: $usage\n" ; }

# Over-write default thresholds if user specified
   foreach my $type (keys %{$user_thresh}) {
      if ($type eq 'thresh_min_L_mw') {

         $supfam_specs->{thresholds}->{L}->{min_mw} = $user_thresh->{$type} ;

      } elsif ($type eq 'thresh_max_L_mw') {

         $supfam_specs->{thresholds}->{L}->{max_mw} = $user_thresh->{$type} ;

      } elsif ($type eq 'thresh_L_bs_seqid') {

         if ($user_thresh->{type} < 0 || $user_thresh->{type} > 1) {
            die "seqid threshold has to range from 0 to 1\n\nUSAGE: $usage" ; }
         $supfam_specs->{thresholds}->{L}->{bs_seqid} = $user_thresh->{$type} ;

      } elsif ($type eq 'thresh_P_bs_seqid') {

         if ($user_thresh->{type} < 0 || $user_thresh->{type} > 1) {
            die "seqid threshold has to range from 0 to 1\n\nUSAGE: $usage" ; }
         $supfam_specs->{thresholds}->{P}->{bs_seqid} = $user_thresh->{$type} ;

      } elsif ($type eq 'thresh_p_bs_seqid') {

         if ($user_thresh->{type} < 0 || $user_thresh->{type} > 1) {
            die "seqid threshold has to range from 0 to 1\n\nUSAGE: $usage" ; }
         $supfam_specs->{thresholds}->{p}->{bs_seqid} = $user_thresh->{$type} ;

      }
   }


   print STDERR "Predicting binding sites\n";
   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1 ) { # master node:

      print "* run_homolobind() ".localtime()
         if(!exists $in->{quiet_fl});
      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{run_homolobind_in},
       $temp_fn->{run_homolobind_in}) =
       tempfile("splits_run_homolobind_input.XXXXX") ;

      ($temp_fh->{run_homolobind_out},
       $temp_fn->{run_homolobind_out}) =
       tempfile("splits_run_homolobind_SGEout.XXXXX") ;
      ($temp_fh->{run_homolobind_err},
       $temp_fn->{run_homolobind_err}) =
       tempfile("splits_run_homolobind_SGEerr.XXXXX") ;

# Sort ASS file by family identifier - to reduce ASTRAL aln read in time
      my $tcom_sort_input = "cat ".$in->{ass_fn}." | sort -t'	' -k9,9nr >".
                            $temp_fn->{run_homolobind_in} ;
      system($tcom_sort_input) ;

      my $split_dir = tempdir("splits_run_homolobind.XXXXX") ;
      my $numjobs ;
      if (exists $in->{numjobs}) {
         $numjobs = $in->{numjobs} ;
      } else {
         $numjobs = $homolobind_specs->{SGE}->{numjobs} ;
      }

      my $splits = homolobind::SGE::_clust_split_ins({
         fn => $temp_fn->{run_homolobind_in},
         dir => $split_dir,
         numjobs => $numjobs,
      }) ;

      my ($perlscript_fh, $perlscript_fn) =
         tempfile("pb.run_homolobind.XXXXX", SUFFIX =>'.ami.pl') ;

      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase;
use homolobind ;
use Bit::Vector ;

main() ;

sub main {

   homolobind::run_homolobind({
      cluster_fl => 0,
      ARGV => \\\@ARGV,
   }) ;

}
" ;

      my ($sgescript_fh, $sgescript_fn)  =
         tempfile("pb.run_homolobind.XXXXX", SUFFIX => ".SGE.sh") ;
      my $sge_outdir = tempdir("SGEOUT.run_homolobind.XXXXX") ;

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $homolobind_specs->{SGE}->{nodespecs}) {
         print {$sgescript_fh} $homolobind_specs->{SGE}->{nodespecs}."\n" ;}

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1 = \$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn -ass_fn \$input1 

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print STDERR "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = homolobind::SGE::_clust_qsub({
         hostname => $homolobind_specs->{SGE}->{headnode},
         sgescript_fn => $sgescript_fn,
      }) ;
      print STDERR "      job $qsub_job_id\n" ;

      while (1) {
         sleep $homolobind_specs->{SGE}->{qstat_sleep} ;
         my $job_status= homolobind::SGE::_clust_qstat({
            hostname => $homolobind_specs->{SGE}->{headnode},
            job_id => $qsub_job_id});
         if ($job_status) {last;}
      }

      homolobind::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{run_homolobind_out},
         err_fn => $temp_fn->{run_homolobind_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      if (exists $in->{out_fn}) {
# remove internal header lines.
         open(REALOUTF, ">".$in->{out_fn}) ;
         open(FULLOUTF, $temp_fn->{run_homolobind_out}) ;
         my $header_line = '#'.join("\t", @result_headers) ;
         print REALOUTF $header_line."\n" ;
         while (my $line = <FULLOUTF>) {
            chomp $line;
            if ($line eq $header_line) {next;}
            print REALOUTF $line."\n" ;
         }
         close(FULLOUTF) ;
         close(REALOUTF) ;
         unlink $temp_fn->{run_homolobind_out} ;
      }

      if (exists $in->{err_fn}) {
         system("mv ".$temp_fn->{run_homolobind_err}." ".
                $in->{err_fn}) ; }

   } else { #processing node:

      my $out_fh ;
      if (exists $in->{out_fn}) {
         open($out_fh, ">".$in->{out_fn}) ;
      } else {
         open($out_fh, ">-") ;
      }

      if (exists $in->{err_fn}) {
         open(STDERR, ">".$in->{err_fn}) ; }

      my $standardres = $supfam_specs->{standardres};

# Preload ASTRAL data
      my $astral = homolobind::pilig::_pilig_astral_preload({
         specs => $homolobind_specs}) ;

# Load PIBASE and LIGBASE binding site info
      my $pb = homolobind::pilig::_pilig_tod_pibase_preload({
         specs => $homolobind_specs }) ;

      my $liginfo = homolobind::pilig::_pilig_load_liginfo({
         fn => $homolobind_specs->{pilig}->{outfiles}->{liginfo} }) ;
      my $class2alnlength = {};

      my $pepnucibits = homolobind::pilig::readin_pepnuciassignments({
         fn=> $homolobind_specs->{pilig}->{outfiles}->{assign_pepnuci_clusters},
         clustrep_fl => 1,
         standardres => $standardres,
         class2alnlength => $class2alnlength,
      }) ;

      my $ligbits = homolobind::pilig::readin_ligassignments({
         fn => $homolobind_specs->{pilig}->{outfiles}->{assign_lig_clusters},
         clustrep_fl => 1,
         liginfo => $liginfo,
         standardres => $standardres,
         class2alnlength => $class2alnlength,
      }) ;

      my $pibits_both = homolobind::pilig::readin_piassignments({
         fn => $homolobind_specs->{pilig}->{outfiles}->{assign_pi_clusters},
         clustrep_fl => 1,
         pb => $pb,
         dont_read_jdomains_fl => 1,
         standardres => $standardres,
         class2alnlength => $class2alnlength,
      }) ;
      my $pibits = $pibits_both->{pibits} ;
      my $interfaces = $pibits_both->{interfaces} ;

      my $alltogether = homolobind::pilig::combine_piligpepbits({
         pb => $pb,
         liginfo => $liginfo,
         ligbits => $ligbits,
         pibits => $pibits,
         pepbits => $pepnucibits->{pep},
      });
      my $class2bits = $alltogether->{class2bits} ;
      my $class2ligs = $alltogether->{class2ligs} ;
      my $class2allligbits = $alltogether->{class2allligbits} ;

      my $expbits = homolobind::pilig::readin_expassignments({
         fn => $homolobind_specs->{pilig}->{outfiles}->{assign_exp},
         standardres => $standardres,
         class2alnlength => $class2alnlength,
      }) ;

# Create a family indexed list of binding sites
      my $class2sid12_side ;
      my $class2pint;
      foreach my $classtype (@{$supfam_specs->{thresholds}->{class_levels}}) {
         foreach my $sid12 (keys %{$interfaces}) {
            my $sid ;
            ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;
            foreach my $side (1,2) {
               $class2sid12_side->{$pb->{sid2class}->{$classtype}->{$sid->{$side}}}->{$sid12."\t".$side}++ ; }
         }

         foreach my $sid (keys %{$pepnucibits->{pep}->{$classtype}}) {
            $class2pint->{$pb->{sid2class}->{$classtype}->{$sid}}->{$sid}++ ; }
      }


#Read in domain information from scop_cla: fam2sf, px2scopid, scopid2class
      my $scopinfo = readin_scop_cla({
         fn => $homolobind_specs->{scop}->{cla_fn}}) ;


#NOTE: sort ASSF by domain fam/sf so less time spent reading in asteroids alns.
      print STDERR "NOW performing homology transfer\n" ;
      my $cur_aln = {};
      open(ASSF, $in->{ass_fn}) ;

# Print headers
      print {$out_fh} '#'.join("\t", @result_headers)."\n";

# Read in a target-SUPFAM assignment entry
      while (my $line = <ASSF>) {
         chomp $line;
         if ($line =~ /^#/) {next;}
         my @t = split(/\t/, $line) ;
         my $curass ;
         map {$curass->{$supfam_specs->{ass_file_headers}->[$_]} = $t[$_]}
            (0 .. $#t) ;

# 0. Check if this domain class is covered by ASTRAL
         my $tmpl_scopid = $scopinfo->{px2scopid}->{$curass->{px_id}} ;
         my $cur_class;
         $cur_class->{fam} = $scopinfo->{scopid2class}->{$tmpl_scopid} ;
         $cur_class->{sf} = $cur_class->{fam} ;
         $cur_class->{sf} =~ s/\.[0-9]+$// ;

# (ASTRAL alignments only cover classes a-g)
         if ($cur_class->{fam} !~ /[a-g]/) {
            my @outvals = ($curass->{seq_id}, $curass->{res_range},
                           'fam', $cur_class->{fam}, 'all',
                           'domain family not covered by ASTRAL') ;
            print {$out_fh} '#'.join("\t", @outvals)."\n" ;
            next;
         }

# 1. Get seq-SUPFAM template alignment: read line from self_hits/sf_id.hits with ^submodel_id px_id (map from target_seq to tmpl_seq)
         $curass->{sf_id} = $scopinfo->{faid2sfid}->{$curass->{fa_id}} ;

         my $supfam_aln ;
         $supfam_aln->{seq} = $curass->{alnstring} ;

         my $target_seq = $supfam_aln->{seq}; $target_seq =~ s/\-//g ;

         $supfam_aln->{strx} = get_SUPFAM_selfhit({
            selfhits_dir => $homolobind_specs->{supfam}->{selfhits_dir},
            supfam_specs => $supfam_specs,
            sf_id => $curass->{sf_id},
            model_id => $curass->{model_id},
            px_id => $curass->{px_id},
         }) ;

# Exit this domain if the template self-hit string not found
         if (!defined $supfam_aln->{strx}) {
            my @outvals = ($curass->{seq_id}, $curass->{res_range},
                           'fam', $cur_class->{fam}, 'all',
                           'ERROR: self-hit alignment not found for '.
                           "px ".$curass->{px_id}.", ".
                           "sf ".$curass->{sf_id}.", ".
                           "model ".$curass->{model_id}) ;
            print {$out_fh} '#'.join("\t", @outvals)."\n" ;
            next;
         }

# 2. Get SUPFAM template's ASTRAL alignment string from ASTRAL
         my $astral_aln;
         my $target_astral_aln = {} ;
         my $abort_target = 0 ;
         foreach my $classtype (@{$supfam_specs->{thresholds}->{class_levels}}){
            my $class = $cur_class->{$classtype} ;

# Load ASTRAL alignment if necessary (ie, new fam or sf)
            if (!exists $cur_aln->{"cur_".$classtype}  ||
                $cur_aln->{"cur_".$classtype} ne $cur_class->{$classtype}) {
               $cur_aln->{$classtype."_aln"} = 
                homolobind::ASTRAL::load_asteroids_aln({
                  aln_fn => $homolobind_specs->{asteroids}->{$classtype.'_aln'}.
                     '/'.$cur_class->{$classtype}.'.fasta_aln' ,
                  seq_fn => $homolobind_specs->{asteroids}->{$classtype.'_seq'}.
                     '/'.$cur_class->{$classtype}.'.fa' ,
                  seqclcont100 => $astral->{seqcl2cont}->{100},
                  seqcl100 => $astral->{seqcl}->{100},
                  allchains => $pb->{pdbchains}
               }) ;
               $cur_aln->{"cur_".$classtype} = $cur_class->{$classtype} ;
            }

            $astral_aln->{$classtype} = {
               strx => $cur_aln->{$classtype."_aln"}->{aln}->{$tmpl_scopid}} ;

            if (!defined $astral_aln->{$classtype}->{strx}) {
               my @outvals = ($curass->{seq_id}, $curass->{res_range},
                           $classtype, $cur_class->{$classtype}, 'all',
                           "ERROR: no $classtype ASTRAL aln for ".
                           "$tmpl_scopid SUPERFAMILY template domain");
               print {$out_fh} '#'.join("\t", @outvals)."\n";
               next;
            }

# Merge the two alignments.; NOTE: SUPFAM aln is not a true aln - includes
#  lower case characters that are not actually part of the alignment.
            $target_astral_aln->{$classtype} = merge_SUPFAM_ASTRAL_alignments({
               superfam_aln => $supfam_aln,
               astral_aln => $astral_aln->{$classtype},
               common => 'strx',
               target => 'seq',
            }) ;

            if (exists $target_astral_aln->{$classtype}->{error_fl}) {

               my @outvals = ($curass->{seq_id}, $curass->{res_range},
                           $classtype, $cur_class->{$classtype}, 'all',
                           "ERROR in merging SUPFAM/ASTRAL alignments, ".
                           "px ".$curass->{px_id}.", ".
                           "sf ".$curass->{sf_id}.", ".
                           "model ".$curass->{model_id}.": ".
                           $target_astral_aln->{$classtype}->{error_fl});
               print {$out_fh} '#'.join("\t", @outvals)."\n";

               $abort_target = 1; 
               last;
            }
         }

         if ($abort_target) {next;}


# Transfer exposed aln positions onto the target sequence domain
         my $warning_noexp = 0;
         foreach my $classtype ('fam') {
            if (!exists $expbits->{$classtype}->{$cur_class->{$classtype}}) {
               print STDERR "WARNING: no assign-exp entries for $classtype ".
                            $cur_class->{$classtype}."\n" ;
               $warning_noexp++ ;
               next;
            }
            my $curexpbits = $expbits->{$classtype}->{$cur_class->{$classtype}};
            my @exp_alnpos = $curexpbits->Index_List_Read() ;
            my $target_exp = {residues => []} ;
            foreach my $exp_alnpos (@exp_alnpos) {
               if (!exists
               $target_astral_aln->{$classtype}->{alnpos2resno}->{$exp_alnpos}){
                  next;}

               my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$exp_alnpos} ;
               push @{$target_exp->{residues}}, $target_res ;
            }
            my $target_exp_reslist = join(", ", @{$target_exp->{residues}}) ;

            my @outvals = ($curass->{seq_id}, $curass->{res_range},
                           $classtype, $cur_class->{$classtype},
                           'exp', '', '', $target_exp_reslist) ;
            print {$out_fh} join("\t", @outvals)."\n" ;
         }

# Iterate over the class level for homology transfer, from (fam, sf)
         foreach my $classtype (@{$supfam_specs->{thresholds}->{class_levels}}){
            my $class = $cur_class->{$classtype} ;

# Initialize warnings to display if no annotations are made
         my $fail_type = {
            L => 'no templates',
            p => 'no templates',
            P => 'no templates',
         } ;

# Iterate over all ligand binding sites in this class;
#   homology transfer to target seq and calc bs/dom-seqid
#  (code adapated from pibase::pilig::collate_per_instance)
         my $target_bs ;
         if (exists $class2allligbits->{$classtype}->{$class}) {
            my $curligs = $class2allligbits->{$classtype}->{$class};

# Iterate over all ligand binding sites in this class;
            foreach my $j (0 .. $#{$curligs}) {

               my ($sid_origdom_l) =
                  ($curligs->[$j]->[2] =~ /SCOP\.(.+)/) ;
               my $curligsig = $curligs->[$j]->[0] ;
               my (undef, $ligcod, undef) = split(/\t/, $curligsig) ;
               my $outcurligsig = $curligsig; $outcurligsig =~ s/\t/:/g ;
               my $curligbits = $curligs->[$j]->[1] ;
               my $curligsid = $curligs->[$j]->[2] ; #fpd090505_1049 


# Check to make sure minimum number of residues met.
               if (exists $supfam_specs->{thresholds}->{min_bs_size} &&
                   $curligbits->Norm() <
                   $supfam_specs->{thresholds}->{min_bs_size}) {next;}

# Check if ligand is within molecular weight range, if specified
               if (exists $supfam_specs->{thresholds}->{L}->{min_mw} &&
                    $liginfo->{mw}->{$ligcod} <
                       $supfam_specs->{thresholds}->{L}->{min_mw}) {next;}

               if (exists $supfam_specs->{thresholds}->{L}->{max_mw} &&
                    $liginfo->{mw}->{$ligcod} >
                       $supfam_specs->{thresholds}->{L}->{max_mw}) {next;}

               $fail_type->{L} = 'sub-threshold templates';

               my $bs_sig = $sid_origdom_l.":".$outcurligsig ;
               $target_bs->{L}->{$bs_sig}->{num_aln} = 0 ;
               if (exists $liginfo->{mw}->{$ligcod} &&
                   $liginfo->{mw}->{$ligcod} ne '') {
                  $target_bs->{L}->{$bs_sig}->{lig_mw} =
                     sprintf("%.3f",$liginfo->{mw}->{$ligcod});
               } else {
                  $target_bs->{L}->{$bs_sig}->{lig_mw} = '' ;
               }
               $target_bs->{L}->{$bs_sig}->{descr_field} =
                  'MW='.$target_bs->{L}->{$bs_sig}->{lig_mw} ;

# Iterate over all alignment positions to calc whole dom seqid
               $target_bs->{L}->{$bs_sig}->{wholedom_aln_length} = 
                 length($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_l});
               $target_bs->{L}->{$bs_sig}->{wholedom_sum_simscore} = 0 ;
               $target_bs->{L}->{$bs_sig}->{wholedom_num_ident} = 0 ;

               foreach my $cur_alnpos (0 ..
                  length($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_l})
                         - 1) {

                  if (!exists 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos}) { next; }

                  my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos} ;

                  my $cur_tmpl_char = substr(
                     $cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_l},
                     $cur_alnpos,1) ;

                  my $cur_targ_char = substr($target_seq, ($target_res - 1),1);
                  if ($cur_tmpl_char eq $cur_targ_char) {
                     $target_bs->{L}->{$bs_sig}->{wholedom_num_ident}++ ; }

                  if (exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char} &&
                      exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char}) {
                     $target_bs->{L}->{$bs_sig}->{wholedom_sum_simscore} +=
$supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char};
                  }
               }

# Iterate over the binding alignment positions and transfer them onto target
               my @bs_alnpos = $curligbits->Index_List_Read() ;
               $target_bs->{L}->{$bs_sig}->{residues} = [] ;
               foreach my $bs_alnpos (@bs_alnpos) {
                  if (!exists
                   $target_astral_aln->{$classtype}->{alnpos2resno}->{$bs_alnpos}){
#                     print STDERR "   no target res aligned to bs alnpos $bs_alnpos ".__LINE__."\n";
                     next;}

                  my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$bs_alnpos} ;
                  push @{$target_bs->{L}->{$bs_sig}->{residues}}, $target_res ;
                  $target_bs->{L}->{$bs_sig}->{num_aln}++ ;

                  my $cur_tmpl_char = substr(
                     $cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_l},
                     $bs_alnpos,1) ;
                  my $cur_targ_char = substr($target_seq, ($target_res - 1),1);

# Check if identical
                  if ($cur_tmpl_char eq $cur_targ_char) {
                     $target_bs->{L}->{$bs_sig}->{num_ident}++ ; }

                  if (exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char} &&
                      exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char}) {
                     $target_bs->{L}->{$bs_sig}->{sum_simscore} +=
$supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char};
                  }
               }

               if ($#{$target_bs->{L}->{$bs_sig}->{residues}} >= 0) {
                  if (!exists $target_bs->{L}->{$bs_sig}->{num_ident} ) {
                     $target_bs->{L}->{$bs_sig}->{num_ident} = 0 ;}
                  if (!exists $target_bs->{L}->{$bs_sig}->{sum_simscore} ) {
                     $target_bs->{L}->{$bs_sig}->{sum_simscore} = 0 ;}
                  $target_bs->{L}->{$bs_sig}->{tmpl_bs_size} = $#bs_alnpos + 1 ;
                  $target_bs->{L}->{$bs_sig}->{targ_bs_size} =
                     $#{$target_bs->{L}->{$bs_sig}->{residues}} + 1 ;
                  $target_bs->{L}->{$bs_sig}->{reslist} = join(", ",
                     @{$target_bs->{L}->{$bs_sig}->{residues}}) ;

               } else {
                  delete $target_bs->{L}->{$bs_sig} ;
               }
            }
         }

# Iterate over all protein binding sites in this class;
#   homology transfer to target seq and calc bs/dom-seqid
         if (exists $class2sid12_side->{$class}) {
            foreach my $sid12_side (keys %{$class2sid12_side->{$class}}) {
               my ($sid, $side) ;
               ($sid->{1}, $sid->{2}, $side) = split(/\t/, $sid12_side);
               my $sid12 = $sid->{1}."\t".$sid->{2} ;
               my $classes ;
               $classes->{1} = $pb->{sid2class}->{'fam'}->{$sid->{1}} ;
               $classes->{2} = $pb->{sid2class}->{'fam'}->{$sid->{2}} ;
               my ($sid_origdom_p) = ($sid->{$side} =~ /SCOP\.(.+)/) ;
               my $bs_sig = $sid12_side; $bs_sig =~ s/\t/:/g ;

               $target_bs->{P}->{$bs_sig}->{partner_class} = $classes->{2} ;
               if ($side == 2) {
                  $target_bs->{P}->{$bs_sig}->{partner_class} = $classes->{1}; }

               $target_bs->{P}->{$bs_sig}->{num_aln} = 0 ;
               $target_bs->{P}->{$bs_sig}->{residues} = [] ;

# Get alignment positions for this binding site.
               if (!exists $interfaces->{$sid12}->{$side} ||
                !exists $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}||
                ($interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Norm()
                    == 0)) {next;}

# Check to make sure minimum number of residues met.
               if (exists $supfam_specs->{thresholds}->{min_bs_size} &&
   $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Norm() <
                   $supfam_specs->{thresholds}->{min_bs_size}) {next;}


               $fail_type->{P} = 'sub-threshold templates';

               $target_bs->{P}->{$bs_sig}->{chains} =
                  $interfaces->{$sid12}->{chains} ;
               $target_bs->{P}->{$bs_sig}->{descr_field} =
                  'chains='.$target_bs->{P}->{$bs_sig}->{chains}.';'.
                  'class='.$target_bs->{P}->{$bs_sig}->{partner_class} ;

               my @bs_alnpos =
   $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Index_List_Read() ;

# Iterate over all alignment positions to calc whole dom seqid
               $target_bs->{P}->{$bs_sig}->{wholedom_aln_length} =
                 length($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p});
               $target_bs->{P}->{$bs_sig}->{wholedom_num_ident} = 0 ;
               $target_bs->{P}->{$bs_sig}->{wholedom_sum_simscore} = 0 ;

               foreach my $cur_alnpos (0 ..
                  length($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p})
                         - 1) {

                  if (!exists 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos}) { next; }
                  my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos} ;

                  my $cur_tmpl_char = substr(
                     $cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p},
                     $cur_alnpos,1) ;
                  my $cur_targ_char = substr($target_seq, ($target_res - 1),1);

                  if ($cur_tmpl_char eq $cur_targ_char) {
                     $target_bs->{P}->{$bs_sig}->{wholedom_num_ident}++ ; }

                  if (exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char} &&
                      exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char}) {
                     $target_bs->{P}->{$bs_sig}->{wholedom_sum_simscore} +=
$supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char};
                  }
               }


# Iterate over alnpos2targetseqresno to homology transfer the bs
               foreach my $bs_alnpos (@bs_alnpos) {
                  if (!exists
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$bs_alnpos}){ next;}

                  my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$bs_alnpos} ;
                  push @{$target_bs->{P}->{$bs_sig}->{residues}}, $target_res ;
                  $target_bs->{P}->{$bs_sig}->{num_aln}++ ;

                  my $cur_tmpl_char = substr(
                     $cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p},
                     $bs_alnpos,1) ;
                  my $cur_targ_char = substr($target_seq, ($target_res - 1),1);

# Check if identical
                  if ($cur_tmpl_char eq $cur_targ_char) {
                     $target_bs->{P}->{$bs_sig}->{num_ident}++ ; }

                  if (exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char} &&
                      exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char}) {
                     $target_bs->{P}->{$bs_sig}->{sum_simscore} +=
$supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char};
                  }
               }

               if ($#{$target_bs->{P}->{$bs_sig}->{residues}} >= 0) {
                  if (!exists $target_bs->{P}->{$bs_sig}->{num_ident} ) {
                     $target_bs->{P}->{$bs_sig}->{num_ident} = 0 ;}
                  if (!exists $target_bs->{P}->{$bs_sig}->{sum_simscore} ) {
                     $target_bs->{P}->{$bs_sig}->{sum_simscore} = 0 ;}
                  $target_bs->{P}->{$bs_sig}->{tmpl_bs_size} = $#bs_alnpos + 1 ;
                  $target_bs->{P}->{$bs_sig}->{targ_bs_size} =
                     $#{$target_bs->{P}->{$bs_sig}->{residues}} + 1 ;
                  $target_bs->{P}->{$bs_sig}->{reslist} = join(", ",
                     @{$target_bs->{P}->{$bs_sig}->{residues}}) ;
               } else {
                  delete $target_bs->{P}->{$bs_sig} ;
               }
            }
         }


# Iterate over all peptide binding sites in this class;
#   homology transfer to target seq and calc bs/dom-seqid
         if (exists $class2pint->{$class}) {
            foreach my $sid (keys %{$class2pint->{$class}}) {
               my ($sid_origdom_p) = ($sid =~ /SCOP\.(.+)/) ;
            foreach my $targetch (
                keys %{$pepnucibits->{pep}->{$classtype}->{$sid}} ){
                if ($targetch eq 'cumulative') {next;}

#Check to make sure minimum number of residues met.
               if (exists $supfam_specs->{thresholds}->{min_bs_size} &&
   $pepnucibits->{pep}->{$classtype}->{$sid}->{$targetch}->Norm() <
                   $supfam_specs->{thresholds}->{min_bs_size}) {next;}

               my $bs_sig = $sid.":".$targetch ;
               $target_bs->{p}->{$bs_sig}->{num_aln} = 0 ;
               $target_bs->{p}->{$bs_sig}->{residues} = [] ;

# Get alignment positions for this binding site.
               my @bs_alnpos = 
   $pepnucibits->{pep}->{$classtype}->{$sid}->{$targetch}->Index_List_Read();

               $fail_type->{p} = 'sub-threshold templates';

# Iterate over all alignment positions to calc whole dom seqid
               $target_bs->{p}->{$bs_sig}->{wholedom_aln_length} = 
                 length($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p});
               $target_bs->{p}->{$bs_sig}->{wholedom_num_ident} = 0 ;
               $target_bs->{p}->{$bs_sig}->{wholedom_sum_simscore} = 0 ;

               $target_bs->{p}->{$bs_sig}->{chain_length} =
                  $pepnucibits->{chain_info}->{$targetch}->{chain_length};
               $target_bs->{p}->{$bs_sig}->{descr_field} =
                  'peptide_length='.$target_bs->{p}->{$bs_sig}->{chain_length};

               foreach my $cur_alnpos (0 ..
                  length($cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p})
                         - 1) {

                  if (!exists 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos}) { next; }
                  my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$cur_alnpos} ;

                  my $cur_tmpl_char = substr(
                     $cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p},
                     $cur_alnpos,1) ;
                  my $cur_targ_char = substr($target_seq, ($target_res - 1),1);

                  if ($cur_tmpl_char eq $cur_targ_char) {
                     $target_bs->{p}->{$bs_sig}->{wholedom_num_ident}++ ; }

                  if (exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char} &&
                      exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char}) {
                     $target_bs->{p}->{$bs_sig}->{wholedom_sum_simscore} +=
$supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char};
                  }
               }

# Iterate over alnpos2targetseqresno to homology transfer the bs
               foreach my $bs_alnpos (@bs_alnpos) {
                  if (!exists
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$bs_alnpos}){ next;}

                  my $target_res = 
   $target_astral_aln->{$classtype}->{alnpos2resno}->{$bs_alnpos} ;

                  push @{$target_bs->{p}->{$bs_sig}->{residues}}, $target_res ;
                  $target_bs->{p}->{$bs_sig}->{num_aln}++ ;

                  my $cur_tmpl_char = substr(
                     $cur_aln->{$classtype."_aln"}->{aln}->{$sid_origdom_p},
                     $bs_alnpos,1) ;
                  my $cur_targ_char = substr($target_seq, ($target_res - 1),1);

# Check if identical
                  if ($cur_tmpl_char eq $cur_targ_char) {
                     $target_bs->{p}->{$bs_sig}->{num_ident}++ ; }

                  if (exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char} &&
                      exists $supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char}) {
                     $target_bs->{p}->{$bs_sig}->{sum_simscore} +=
$supfam_specs->{substitution_matrix}->{nl}->{$cur_tmpl_char}->{$cur_targ_char};
                  }
               }

               if ($#{$target_bs->{p}->{$bs_sig}->{residues}} >= 0) {
                  if (!exists $target_bs->{p}->{$bs_sig}->{num_ident} ) {
                     $target_bs->{p}->{$bs_sig}->{num_ident} = 0 ;}
                  if (!exists $target_bs->{p}->{$bs_sig}->{sum_simscore} ) {
                     $target_bs->{p}->{$bs_sig}->{sum_simscore} = 0 ;}
                  $target_bs->{p}->{$bs_sig}->{tmpl_bs_size} = $#bs_alnpos + 1 ;
                  $target_bs->{p}->{$bs_sig}->{targ_bs_size} =
                     $#{$target_bs->{p}->{$bs_sig}->{residues}} + 1 ;
                  $target_bs->{p}->{$bs_sig}->{reslist} = join(", ",
                     @{$target_bs->{p}->{$bs_sig}->{residues}}) ;
               } else {
                  delete $target_bs->{p}->{$bs_sig} ;
               }
            }
            }
         }

# Display transferred binding sites
         foreach my $bs_type (qw/L p P/) {
            if (! exists $target_bs->{$bs_type}) {
                  my @outvals = ($curass->{seq_id}, $curass->{res_range},
                                 $classtype, $cur_class->{$classtype},
                                 $bs_type, $fail_type->{$bs_type}) ;
               print {$out_fh} '#'.join("\t", @outvals)."\n" ;
               next;
            }

            my $num_transf = 0 ;
            foreach my $bs_sig (keys %{$target_bs->{$bs_type}}) {

               if (! exists $target_bs->{$bs_type}->{$bs_sig}->{num_aln} ||
                   $target_bs->{$bs_type}->{$bs_sig}->{num_aln} == 0) {
                  next; }
               my $wholedom_percident = sprintf("%.3f",
                  ($target_bs->{$bs_type}->{$bs_sig}->{wholedom_num_ident} / 
                   $target_bs->{$bs_type}->{$bs_sig}->{wholedom_aln_length})) ;
               my $bs_percident = sprintf("%.3f",
                  ($target_bs->{$bs_type}->{$bs_sig}->{num_ident} / 
                   $target_bs->{$bs_type}->{$bs_sig}->{tmpl_bs_size})) ;

               my $wholedom_percsim = sprintf("%.3f",
                  ($target_bs->{$bs_type}->{$bs_sig}->{wholedom_sum_simscore} / 
                   $target_bs->{$bs_type}->{$bs_sig}->{wholedom_aln_length})) ;
               my $bs_percsim = sprintf("%.3f",
                  ($target_bs->{$bs_type}->{$bs_sig}->{sum_simscore} / 
                   $target_bs->{$bs_type}->{$bs_sig}->{tmpl_bs_size})) ;


# If user thresholds defined, will supercede benchmark thresholds
               if (exists $supfam_specs->{thresholds}->{$bs_type}->{bs_seqid}){
                  if ($bs_percident < 
                     $supfam_specs->{thresholds}->{$bs_type}->{bs_seqid}){next;}

# Skip if threshold not determined for template OR seqid < threshold
               } elsif (!exists $benchmark_thresholds->{$bs_sig} ||
                   $benchmark_thresholds->{$bs_sig} eq 'UNDEF' ||
                   $bs_percident < $benchmark_thresholds->{$bs_sig}) {
                  next;
               }

# If a minimum fraction aligned binding site residues defined, check it
               if (exists $supfam_specs->{thresholds}->{frac_bs_aligned} &&
                  $target_bs->{$bs_type}->{$bs_sig}->{targ_bs_size} <
                  ($target_bs->{$bs_type}->{$bs_sig}->{tmpl_bs_size} *
                   $supfam_specs->{thresholds}->{frac_bs_aligned})) {
                  next;
               }

               $num_transf++ ;
               my $bs_fracaln = sprintf("%.3f",
                  ($target_bs->{$bs_type}->{$bs_sig}->{num_aln} / 
                   $target_bs->{$bs_type}->{$bs_sig}->{tmpl_bs_size})) ;
               my @outvals = ($curass->{seq_id}, $curass->{res_range},
                              $classtype, $cur_class->{$classtype},
                              $bs_type, $bs_sig,
                              $target_bs->{$bs_type}->{$bs_sig}->{descr_field},
                              $target_bs->{$bs_type}->{$bs_sig}->{reslist},
                              $bs_percident,
                              $bs_percsim,
                              $target_bs->{$bs_type}->{$bs_sig}->{num_ident},
                              $target_bs->{$bs_type}->{$bs_sig}->{num_aln},
                              $target_bs->{$bs_type}->{$bs_sig}->{tmpl_bs_size},
                              $bs_fracaln,
   $target_bs->{$bs_type}->{$bs_sig}->{wholedom_num_ident},
   $target_bs->{$bs_type}->{$bs_sig}->{wholedom_aln_length},
                              $wholedom_percident,
                              $wholedom_percsim
                             );
               print {$out_fh} join("\t", @outvals)."\n" ;
            }
            if ($num_transf == 0) {
                  my @outvals = ($curass->{seq_id}, $curass->{res_range},
                                 $classtype, $cur_class->{$classtype},
                                 $bs_type, $fail_type->{$bs_type}) ;
               print {$out_fh} '#'.join("\t", @outvals)."\n" ;
            }
         }
         }
      }
   }

}



=head2 set_supfam_specs

   Title:       set_supfam_specs()
   Function:    Sets configuration parameters
   Args:        None
   Returns:     $_->{option} = value; hash of parameters

=cut

sub set_supfam_specs {

   my $params ;

   my $pilig_specs = homolobind::pilig::set_pilig_specs() ;
   $params->{pilig_specs} = $pilig_specs ;

   $params->{standardres} = homolobind::pilig::_list_standardres() ;

   $params->{thresholds} = {
      min_bs_size => 5,       # at least 5 residues
      frac_bs_aligned => 0.5, # at least 50% of template residues aligned
      class_levels => ['fam'],
      'L' => { min_mw => 250, max_mw => 1000 },
      'P' => { },
      'p' => { },
   } ;

   $params->{scop_cla_fn} =
      $homolobind_specs->{scop}->{cla_fn};

   $params->{ass_file_headers} = [
      'sp',
      'seq_id',
      'model_id',
      'res_range',
      'sf_eval',
      'alnstring',
      'fa_eval',
      'px_id',
      'fa_id'
   ] ;
   map {$params->{ass_file_f2i}->{$params->{ass_file_headers}->[$_]} = $_;}
      (0 .. $#{$params->{ass_file_headers}}) ;

   $params->{selfhits_headers} = [
      'model_id',
      'px_id',
      'alnstring',
   ] ;
   map {$params->{selfhits_f2i}->{$params->{selfhits_headers}->[$_]} = $_;}
      (0 .. $#{$params->{selfhits_headers}}) ;

# BLOSUM62 matrix: http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
   my $blosum62_matrix = "#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
" ;

   $params->{substitution_matrix} = readin_substitution_matrix({
      string => $blosum62_matrix}) ;

   return $params ;

}


=head2 readin_substitution_matrix

   Title:       readin_substitution_matrix()
   Function:    Parses a substitution matrix in matblas format
   Args:        ->{matrix_fn} = filename of substitution matrix
                ->{string} = string containing contents of matrix file

   Returns:     ->{raw}->{aa1}->{aa2} = substitution matrix score for aa1,aa2
                ->{nl}->{aa1}->{aa2} = normalized aa similarity score
                  (see karling_normalize_matrix() for normalization scheme)

=cut

sub readin_substitution_matrix {
   my $in = shift ;

   my @matrix_string ;
   if (exists $in->{string}) {
      @matrix_string = split(/\n/, $in->{string});
   } elsif (exists $in->{matrix_fn}) {
      open(MATF, $in->{matrix_fn});
      @matrix_string = <MATF> ;
      close(MATF);
   } else {
      die "substitution matrix not specified" ;
   }

   my @aa_order ;
   my $raw_matrix ;
   foreach my $mat_line (@matrix_string) {
      $mat_line =~ s/\s$// ;
      if ($mat_line =~ /^#/) {
         next; }
      my @t = split(' ', $mat_line) ;
      if ($mat_line !~ /[0-9]/) {
         @aa_order = @t ;
      } else {
         foreach my $j ( 1 .. $#t) {
            $raw_matrix->{$t[0]}->{$aa_order[($j -1)]} =
               $t[$j] ;
         }
      }
   }

   my $nl_matrix = karlin_normalize_matrix({
      matrix => $raw_matrix}); 

   return {
      raw => $raw_matrix,
      nl => $nl_matrix
   } ;

}


=head2 karlin_normalize_matrix

   Title:       karlin_normalize_matrix()
   Function:    Normalizes a substitution matrix per Karlin and Brocchieri,
                  J Bacteriol 1996; and rescaled to range from 0-1:

               0.5 * (1 + ( mat(i,j) / sqrt(|mat(i,i) * mat(j,j)|)))

   Args:        ->{matrix}->{aa1}->{aa2} = raw substitution matrix
   Returns:     ->{aa1}->{aa2} = normalized similarity score

=cut

sub karlin_normalize_matrix {

   my $in = shift ;
   my $raw_matrix = $in->{matrix};

   my $nl_matrix ;
   foreach my $aa1 ( keys %{$raw_matrix} ) {
      foreach my $aa2 ( keys %{$raw_matrix->{$aa1}} ) {
         $nl_matrix->{$aa1}->{$aa2} = ($raw_matrix->{$aa1}->{$aa2}) /
            sqrt( abs($raw_matrix->{$aa1}->{$aa1} *
                      $raw_matrix->{$aa2}->{$aa2})) ;

         $nl_matrix->{$aa1}->{$aa2} = ($nl_matrix->{$aa1}->{$aa2} + 1) / 2 ;
      }
   }

   return $nl_matrix ;
}


=head2 readin_benchmark_thresholds

   Title:       readin_benchmark_thresholds()
   Function:    Reads in sequence identity cutoffs determined by binding site
                 library benchmarking to yield 1% FPR
   Args:        ->{fn} = benchmark file location
   Returns:     ->{binding site identifier} = seqid cutoff for 1% FPR

=cut

sub readin_benchmark_thresholds {

#BS      Pintra  dom	a.1.2.1 BDP67107-0_SCOP.d2b76n1 BDP67107-0_SCOP.d2b76n2 1       20      10028   0.250   0.002   0.200   0.016   0.200   0.016
#BS      Pinter  dom	a.1.1.1 BDP32771-0_SCOP.d1ngkd_ BDP32771-0_SCOP.d1ngkg_ 2       13      10050   0.308   0.003   0.308   0.003   0.231   0.027
#BS      p       dom	a.1.1.2 BDP59974-0_SCOP.d1yiha1 BDP59974-4_CHAIN-D      13      11853   0.308   0.005   0.308   0.005   0.231   0.032
#BS      L       dom	a.1.1.1 1ux8    HEM     001     26      10050   0.231   0.003   0.192   0.012   0.154   0.050

   my $in = shift;

   my $thresh; 
   if ($in->{fn} =~ /\.gz$/) {
      open(BMF, "zcat ".$in->{fn}." |") ;
   } else {
      open(BMF, $in->{fn}) ;
   }
   while (my $line = <BMF>) {
      chomp $line;
      if ($line !~ /^BS/) {next;}
      my @t =split(/\t/, $line) ;
      my $btype = $t[1] ;
      my $sid = $t[2] ;
      my $fam = $t[3] ;

      if ($btype eq 'Pintra' || $btype eq 'Pinter') {
         my $bssig = $t[4].":".$t[5].":".$t[6] ;
         $thresh->{$bssig} = $t[9] ;
      } elsif ($btype eq 'p') {
         my $bssig = $t[4].":".$t[5] ;
         $thresh->{$bssig} = $t[8] ;
      } elsif ($btype eq 'L') {
         my ($osid) = ($t[2] =~ /SCOP\.(.+)/) ;
         my $bssig = $osid.":".$t[4].":".$t[5].":".$t[6] ;
         $thresh->{$bssig} = $t[9] ;
      }
   }
   close(BMF) ;

   return $thresh ;

}


=head2 fetch_data

   Title:       fetch_data()
   Function:    Fetch data from ASTRAL, ASTEROIDS, HOMOLOBIND, and SCOP servers.

   Args:        NOTHING
   Returns:     NOTHING

=cut

sub fetch_data {

   my $orig_dir = `pwd` ;
   foreach my $data_source (keys %{$homolobind_specs->{urls}}) {
      foreach my $url (keys %{$homolobind_specs->{urls}->{$data_source}}) {
         print STDERR "$data_source\t$url\n" ;
         my $info = $homolobind_specs->{urls}->{$data_source}->{$url};
         my $dest = $homolobind_specs->{$data_source}->{dir} ;
         if (exists $info->{destination}) {
            $dest = $info->{destination} ; }
         if (!-s $dest) {
            mkpath($dest) ; }
         chdir $dest ;

         my $wget_com = "wget ";
         if (exists $info->{wget_options}) {
            $wget_com .= $info->{wget_options}.' ' ; }
         $wget_com .= $url ;

         system($wget_com) ;

         my $fn = $url ; $fn =~ s/^.*\///g ;

         if ($fn =~ /\.tgz$/ || $fn =~ /\.tar\.gz$/) {
            my $uncomp_com = "tar xvfz ".$fn ;
            system($uncomp_com) ;
         }
      }
   }
   chdir $orig_dir ;

   if (!-s $homolobind_specs->{supfam}->{dir}) {
      mkpath($homolobind_specs->{supfam}->{dir}) ; }

   return ;

}


=head2 merge_SUPFAM_ASTRAL_alignments

   Title:       merge_SUPFAM_ASTRAL_alignments()
   Function:    Merges SUPERFAMILY alignment string with ASTRAL alignment to
                 get SUPERFAMILY annotated target sequence in the ASTRAL
                 alignment frame (where it can receive binding site annotations)

   Args:        ->{superfam_aln} = SUPERFAMILY alignment string
                ->{astral_aln} = ASTRAL alignment string
                ->{target} = target sequence name
                ->{common} = template sequence name present in both alignments

   Returns:     ->{resno2alnpos}->{target resno} = alignment position
                ->{alnpos2resno}->{alignment position} = target resno

=cut

sub merge_SUPFAM_ASTRAL_alignments {

   my $in = shift ;
   my $aln1= $in->{superfam_aln} ;
   my $aln2 = $in->{astral_aln} ;
   my $target = $in->{target} ;
   my $common = $in->{common} ;

# CHECK 1. MAKE SURE THE COMMON SEQUENCES ARE ACTUALLY IDENTICAL!
   my $seq1 = $aln1->{$common} ; $seq1 =~ s/\-//g ; $seq1 = uc($seq1) ;
   my $seq2 = $aln2->{$common} ; $seq2 =~ s/\-//g ; $seq2 = uc($seq2) ;
   if ($seq1 ne $seq2) {
      return {error_fl => "the \"common\" sequence differs between the two ".
              "alignments; CANNOT MERGE!\tseq1: $seq1\tseq2: $seq2"} ;
   }

# ALIGNMENT 1: SUPFAM alignment string
# NOTE NOT A TRUE ALIGNMENT: skip residues in lower case
   my $aln1_resno2alnpos = {}; my $aln1_alnpos2resno = {};
   foreach my $seq ($target, $common) {
      my $cur_resno = 0 ; my $cur_alnpos = 0 ;
      foreach my $j (0 .. (length($aln1->{$seq}) - 1)) {
         my $alnchar = substr($aln1->{$seq}, $j, 1) ;
         if ($alnchar ne '-') {$cur_resno++;}
         if ($alnchar =~ /[a-z]/) {next;}
         $cur_alnpos++ ;

         if ($alnchar eq '-') {next;}

         $aln1_resno2alnpos->{$seq}->{$cur_resno} = $cur_alnpos ;
         $aln1_alnpos2resno->{$seq}->{$cur_alnpos} = $cur_resno ;
      }
   }

# Parse aln2 for seqresno_2_alnpos mapping for common
# ALIGNMENT 2: ASTRAL alignment string
   my $aln2_resno2alnpos = {}; my $aln2_alnpos2resno = {};
   my $cur_resno = 0 ;
   foreach my $alnpos (0 .. (length($aln2->{$common}) - 1)) {
      my $alnchar = substr($aln2->{$common}, $alnpos, 1) ;
      if ($alnchar eq '-') {next; }
      $cur_resno++ ;
      $aln2_resno2alnpos->{$common}->{$cur_resno} = $alnpos ;
      $aln2_alnpos2resno->{$common}->{$alnpos} = $cur_resno ;
   }

# Combine seqresno_2_alnpos mappings to get target_2_alnpos2 mapping
   foreach my $target_resno (sort {$a <=> $b}
                             keys %{$aln1_resno2alnpos->{$target}}) {

      my $aln1_pos = $aln1_resno2alnpos->{$target}->{$target_resno} ;
      if (!exists $aln1_alnpos2resno->{$common}->{$aln1_pos}) {next;}
      my $common_resno = $aln1_alnpos2resno->{$common}->{$aln1_pos} ;

      if (!exists $aln2_resno2alnpos->{$common}->{$common_resno}) {next;}
      my $aln2_pos = $aln2_resno2alnpos->{$common}->{$common_resno} ;

      $aln2_resno2alnpos->{$target}->{$target_resno} = $aln2_pos ;
      $aln2_alnpos2resno->{$target}->{$aln2_pos} = $target_resno ;
   }

   return {
      resno2alnpos => $aln2_resno2alnpos->{$target},
      alnpos2resno => $aln2_alnpos2resno->{$target},
   } ;

}


=head2 readin_scop_cla

   Title:       readin_scop_cla()
   Function:    Reads in SCOP cla file
   Args:        ->{fn} = SCOP cla filename
   Returns:     ->{px2scopid}->{px_id} = scopid
                ->{scopid2class}->{scopid} = class
                ->{faid2sfid}->{fa_id} = sf_id

=cut

sub readin_scop_cla {

   my $in = shift ;
   my $scopinfo ;

   open(SCOPCLA, $in->{fn}) ;
   while (my $line = <SCOPCLA>) { #parsing logic from pibase::SCOP
      if ($line =~ /^\#/) {next;}
      chomp $line;
      my ($sid, $pdb_id, $raw_domaindef, $sccs, $raw_px_id, $raw_sunid ) =
         split(/\t/, $line) ;

      my ($cl_id, $cf_id, $sf_id, $fa_id, $dm_id, $sp_id, $px_id) =
( $raw_sunid =~ /cl=(.+),cf=(.+),sf=(.+),fa=(.+),dm=(.+),sp=(.+),px=(.+)/ ) ;

      $scopinfo->{px2scopid}->{$px_id} = $sid ;
      $scopinfo->{scopid2class}->{$sid} = $sccs ;
      $scopinfo->{faid2sfid}->{$fa_id} = $sf_id ;
   }
   close(SCOPCLA) ;

   return $scopinfo ;
}


=head2 get_SUPFAM_selfhit

   Title:       get_SUPFAM_selfhit()
   Function:    Retrieves alignment string from SUPFAM self-hit files for
                  a particular SCOP template domain

   Args:        ->{supfam_specs} = configuration parameters
                ->{px_id} = SCOP px id
                ->{model_id} = SUPERFAMILY model id

   Returns:     self-hit alignment string

=cut

sub get_SUPFAM_selfhit {

   my $in = shift ;
   my $supfam_specs = $in->{supfam_specs} ;
   my $selfhits_dir = $in->{selfhits_dir} ;

   open(SELFHITF, $selfhits_dir.'/'.$in->{sf_id}.'.tab') ;
   my @selfhit_headers ;
   my $alnstring ;
   while (my $line = <SELFHITF>) {
      chomp $line;
      my @t = split(/\t/, $line) ;
      if ($t[$supfam_specs->{selfhits_f2i}->{px_id}] eq $in->{px_id} &&
          $t[$supfam_specs->{selfhits_f2i}->{model_id}] eq $in->{model_id}) {
         $alnstring = $t[$supfam_specs->{selfhits_f2i}->{alnstring}] ;
         last;
      }
   }

   return $alnstring ;

}


=head2 summarize_results_withR()

   Title:       summarize_results_withR()
   Function:    Parses run_homolobind() output and reports prediction summary.

                Uses R through RSPerl interface to calculate significance
                of overlap between predicted ligand and protein binding sites

   Args:        ->{results_fn} = run_homolobind() output file
   Returns:     nothing
   Displays:    table describing numbers of proteins,domains,families,residues
                 with each kind of annotation

=cut

sub summarize_results_withR {

   require R;
   require RReferences ;
   R::initR("--silent") ;

   my $in = shift ;
   open(RESF, $in->{results_fn}) ;
   my $f2i ;
   my $header_line ;
   {
      my $headers = <RESF> ; $header_line = $headers ;
      $headers =~ s/^\#// ;
      my @headers = split(/\t/, $headers) ;
      map {$f2i->{$headers[$_]} = $_} (0 .. $#headers) ;
   }

   my $protein2dom ;
   my $protein2domfam ;
   my $protein2dom2bs_res ;
   my $protein2dom2exp_res ;
   while (my $line = <RESF>) {
      chomp $line;
      if ($line eq $header_line) { #skip if an (internal) header line
         next;
      } elsif ($line =~ /^#/) { #domain was skipped

         $line =~ s/^\#// ;
         my @t = split(/\t/, $line) ;
         $protein2dom->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$t[$f2i->{'classtype'}]} = $t[$f2i->{'class'}] ;
         $protein2domfam->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'class'}]}++ ;

      } elsif ($line =~ /	exp	/) { #exposed residue list

         my @t = split(/\t/, $line) ;
         if ($#t >= $f2i->{'residues'} && $t[$f2i->{'residues'}] ne '') {
            my @res = split(/\, /, $t[$f2i->{'residues'}]) ;
            map {$protein2dom2exp_res->{$t[$f2i->{'classtype'}]}->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$_}++ } @res ;
         }

      } else { #annotation line

         my @t = split(/\t/, $line) ;
         $protein2dom->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$t[$f2i->{'classtype'}]} = $t[$f2i->{'class'}] ;
         $protein2domfam->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'class'}]}++ ;

         my @res = split(/\, /, $t[$f2i->{'residues'}]) ;

         map {$protein2dom2bs_res->{$t[$f2i->{'classtype'}]}->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$t[$f2i->{'bs_type'}]}->{$_}++ } @res ;
      }
   }
   close(RESF) ;

   my $results ;
   my @classtypes = keys %{$protein2dom2bs_res} ;
   $results->{all}->{num_proteins} = keys %{$protein2dom} ;
   $results->{all}->{num_domains} = 0 ;
   $results->{all}->{num_residues} = 0 ;


   {
      my @summary_headers= ('OVERLAP_SUMMARY', 'sequence_id',
                            'num_exposed_res',
                            'num_ligand_binding_res',
                            'num_protein_binding_res',
                            'num_bifunctional_res',
                            'Pvalue_right', 'Pvalue_left',
                            'overlap_score',
                            'domain_families') ;
      print '#'.join("\t", @summary_headers)."\n" ;
   }

# Iterate over all proteins
   foreach my $seq_id (sort keys %{$protein2dom}) {

      $results->{all}->{num_domains} += keys %{$protein2dom->{$seq_id}} ;

      my $res_list ; #tally of binding site annotated residues
      my $curseq_nums = {} ; # Keep track of E, L, P, and LP positions

# Count number of bi-functional positions
      $curseq_nums->{fam} = {
            'L_and_protein' => 0,
            'E' => 0,
            'L' => 0,
            'protein' => 0,
            'P' => 0,
            'p' => 0,
      } ;

# Iterate over all domains
      my $noexpres_fl = 0 ;
      foreach my $res_range (keys %{$protein2dom->{$seq_id}}) {
         my $curdom_res_list ; #res_list for the current domain

# Count all residues in domain
         foreach my $subrange (split(',', $res_range)) {
            if ($subrange =~ /\-/) {
               my ($start, $stop) = ($subrange =~ /^([0-9]+)\-([0-9]+)$/);
               map {$res_list->{all}->{$_}++} ($start .. $stop) ;
            } else {
               $res_list->{all}->{$subrange}++ ;
            }
         }

# Count family/superfamily type
         foreach my $classtype (keys %{$protein2dom->{$seq_id}->{$res_range}}){
            $results->{all}->{'list_'.$classtype}->{$protein2dom->{$seq_id}->{$res_range}->{$classtype}}++ ;
         }

# Count up annotated binding sites: protein, domains, fam/sf, residues
         foreach my $classtype (keys %{$protein2dom2bs_res}) {
            if (!exists $protein2dom2bs_res->{$classtype}->{$seq_id} ||
             !exists $protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}){
               next;}
            foreach my $bs_type (keys %{$protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}}) {

               my @cur_bs_types = ('annotated', $bs_type);
               if ($bs_type eq 'P' || $bs_type eq 'p') {
                  push @cur_bs_types, 'protein' ; }

               foreach my $cur_bs_type (@cur_bs_types) {
                  map {$curdom_res_list->{bs}->{$classtype}->{$cur_bs_type}->{$_}++; }
                 (keys
   %{$protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}->{$bs_type}}) ;

                  $results->{bs}->{$classtype}->{$cur_bs_type}->{'list_proteins'}->{$seq_id}++;
                  $results->{bs}->{$classtype}->{$cur_bs_type}->{'list_domains'}->{$seq_id."\t".$res_range}++ ;
                  $results->{bs}->{$classtype}->{$cur_bs_type}->{"list_$classtype"}->{$protein2dom->{$seq_id}->{$res_range}->{$classtype}}++ ;
               }
            }

# count up bi-functional = L and (P or p) residues
            if (exists $curdom_res_list->{bs}->{$classtype}->{L} &&
                exists $curdom_res_list->{bs}->{$classtype}->{protein}){
            foreach my $res (keys %{$curdom_res_list->{bs}->{$classtype}->{L}}) {
               if (exists $curdom_res_list->{bs}->{$classtype}->{protein}->{$res}){
                  $curdom_res_list->{bs}->{$classtype}->{"L_and_protein"}->{$res}++ ;
               }
            }
            }

# If any bi-functional residues, count protein/domain/family
            if (exists $curdom_res_list->{bs}->{$classtype}->{"L_and_protein"}){
               $results->{bs}->{$classtype}->{"L_and_protein"}->{'list_proteins'}->{$seq_id}++;
               $results->{bs}->{$classtype}->{"L_and_protein"}->{'list_domains'}->{$seq_id."\t".$res_range}++ ;
               $results->{bs}->{$classtype}->{"L_and_protein"}->{"list_$classtype"}->{$protein2dom->{$seq_id}->{$res_range}->{$classtype}}++ ;
            }

#Transfer residues onto master res-list
            foreach my $bs_type (keys %{$curdom_res_list->{bs}->{$classtype}}){
               map {$res_list->{bs}->{$classtype}->{$bs_type}->{$_}++}
                  keys %{$curdom_res_list->{bs}->{$classtype}->{$bs_type}} ;}

         }

#100608_1101 - Count number of exposed and predicted binding site residues
         foreach my $classtype (keys %{$protein2dom2exp_res}) {
            if (exists $protein2dom2exp_res->{$classtype}->{$seq_id} &&
            exists $protein2dom2exp_res->{$classtype}->{$seq_id}->{$res_range}){
               $curseq_nums->{$classtype}->{E} += keys
                 %{$protein2dom2exp_res->{$classtype}->{$seq_id}->{$res_range}};
            } else {
               print "#WARNING: exposed residues not reported for ".
                     "$seq_id ($res_range)\n" ;
               $noexpres_fl = 1 ;
            }
         }
      }


      foreach my $classtype (keys %{$res_list->{bs}}) {
# count up residues per binding site type
         foreach my $bs_type (keys %{$res_list->{bs}->{$classtype}}) {
            $curseq_nums->{$classtype}->{$bs_type} =
               keys %{$res_list->{bs}->{$classtype}->{$bs_type}} ;

            $results->{bs}->{$classtype}->{$bs_type}->{num_residues} +=
               keys %{$res_list->{bs}->{$classtype}->{$bs_type}} ;
         }
      }

# Calculate p-value of observed overlap or non-overlap for this protein

      my $overlap_score = ''; 
      my $overlap_pvals = {
         "L_P" => '',
         "L_P_less" => '',
      } ;
      if ($noexpres_fl == 0) {
      { # code from pibase::pilig.pm snapshot
         my $tnorm ;
         $tnorm->{lp} = $curseq_nums->{fam}->{"E"} -
                        $curseq_nums->{fam}->{"L"} -
                        $curseq_nums->{fam}->{"protein"} +
                        $curseq_nums->{fam}->{"L_and_protein"} ;

         $tnorm->{Lp} = $curseq_nums->{fam}->{"L"} -
                        $curseq_nums->{fam}->{"L_and_protein"} ;

         $tnorm->{lP} = $curseq_nums->{fam}->{"protein"} -
                        $curseq_nums->{fam}->{"L_and_protein"} ;

         $tnorm->{LP} = $curseq_nums->{fam}->{"L_and_protein"} ;

# Check to ensure any negative values set to 0; occurs in cases where
#  a residue not considered exposed is part of a binding site;
#  in the n=30K human proteins, this happens n=3 times: a single residue
#  part of the binding site is not considered exposed.
# (ENSP00000356121, ENSP00000356122, ENSP00000375801)
         foreach my $type (keys %{$tnorm}) {
            if ($tnorm->{$type} < 0) {
               print STDERR "Note: $seq_id |$type| = ".
                               $tnorm->{$type}."; setting to 0.\n";
               $tnorm->{$type} = 0 ;
            }
         }

         my @x = &R::eval("capture.output(fisher.test(matrix(c($tnorm->{lp}, $tnorm->{Lp}, $tnorm->{lP}, $tnorm->{LP}), 2, 2),workspace=2e7,alternative='g'))") ;
         my $x = join("\n", @x) ;
         my ($pval) = ($x =~ /p-value [\<=] (.*)\nalternative/) ;
         $overlap_pvals->{"L_P"} = $pval ;

         my @y = &R::eval("capture.output(fisher.test(matrix(c($tnorm->{lp}, $tnorm->{Lp}, $tnorm->{lP}, $tnorm->{LP}), 2, 2),workspace=2e7,alternative='l'))") ;
         my $y = join("\n", @y) ;
         my ($pval_less) = ($y =~ /p-value [\<=] (.*)\nalternative/) ;
         $overlap_pvals->{"L_P_less"} = $pval_less ;
      }

      $overlap_score = 0 ;
      if ($curseq_nums->{fam}->{L} > 0 && $curseq_nums->{fam}->{protein} > 0) {
         my $overlap_score_raw = (($curseq_nums->{fam}->{"L_and_protein"} *
                                   $curseq_nums->{fam}->{E}) /
                                  ($curseq_nums->{fam}->{L} *
                                   $curseq_nums->{fam}->{protein})) ;
         $overlap_score = sprintf("%.3f", $overlap_score_raw) ;
      }
      }

      my @curseq_outvals = ("OVERLAP_SUMMARY", $seq_id,
                            $curseq_nums->{fam}->{E},
                            $curseq_nums->{fam}->{L},
                            $curseq_nums->{fam}->{protein},
                            $curseq_nums->{fam}->{"L_and_protein"},
                            $overlap_pvals->{"L_P"},
                            $overlap_pvals->{"L_P_less"},
                            $overlap_score,
                            join(";",sort keys
                                 %{$protein2domfam->{$seq_id}})) ;
      print '#'.join("\t", @curseq_outvals)."\n" ;

      $results->{all}->{num_residues} += keys %{$res_list->{all}} ;
   }

# Count final summary statistics
   foreach my $classtype (keys %{$results->{bs}}) {

# count up number of families/superfamilies/domains/proteins
      foreach my $bstype (keys %{$results->{bs}->{$classtype}}) {
         $results->{bs}->{$classtype}->{$bstype}->{num_proteins} =
            keys %{$results->{bs}->{$classtype}->{$bstype}->{list_proteins}} ;

         $results->{bs}->{$classtype}->{$bstype}->{num_domains} =
            keys %{$results->{bs}->{$classtype}->{$bstype}->{list_domains}} ;

         $results->{bs}->{$classtype}->{$bstype}->{'num_'.$classtype} =
           keys %{$results->{bs}->{$classtype}->{$bstype}->{'list_'.$classtype}};
      }
   }

# Display overall counts
   print "\n\n#OVERALL COUNTS\t(note: p=peptide, P=domain-domain, L=small molecule ligand):\n" ;
   foreach my $field  (keys %{$results->{all}}) {
      if ($field =~ /list/) {
         my $newfield = $field ;$newfield =~ s/^list_/num_/ ;
         my $val = keys %{$results->{all}->{$field}} ;
         $results->{all}->{$newfield} = $val ;
         print "#ALL $newfield: $val\n" ;
      } else {
         print "#ALL $field: ".$results->{all}->{$field}."\n";
      }
   }

# Display counts for each kind of binding site
   foreach my $class_type (keys %{$results->{bs}}) {
      foreach my $bs_type (keys %{$results->{bs}->{$class_type}}) {
         foreach my $field (keys %{$results->{bs}->{$class_type}->{$bs_type}}){
            if ($field !~ /^num/) {next;}
         print "#BS: $class_type $bs_type $field ".
            $results->{bs}->{$class_type}->{$bs_type}->{$field}."\n"; }}}

# Format summary numbers as a LaTeX table
   print "\n\n" ;
   print '#'.join(" & ", 'CLASSTYPE', 'group', 'num_proteins',
               'num_domains',
               'num_fam',
               'num_residues')."\\\\\n" ;
   print '#'.join(" & ", '', 'annotated_domains',
               $results->{all}->{num_proteins},
               $results->{all}->{num_domains},
               $results->{all}->{num_fam},
               $results->{all}->{num_residues})."\\\\\n";
   foreach my $class_type (keys %{$results->{bs}}) {
      foreach my $bs_type (keys %{$results->{bs}->{$class_type}}) {
         print '#'.join(" & ", $class_type, $bs_type,
         $results->{bs}->{$class_type}->{$bs_type}->{num_proteins},
         $results->{bs}->{$class_type}->{$bs_type}->{num_domains},
         $results->{bs}->{$class_type}->{$bs_type}->{num_fam},
         $results->{bs}->{$class_type}->{$bs_type}->{num_residues},
         )."\\\\\n" ;
      }
   }

}


=head2 summarize_results_noR()

   Title:       summarize_results_noR()
   Function:    Parses run_homolobind() output and reports prediction summary.
   Args:        ->{results_fn} = run_homolobind() output file
   Returns:     nothing
   Displays:    table describing numbers of proteins,domains,families,residues
                 with each kind of annotation

=cut

sub summarize_results_noR {

   my $in = shift ;
   open(RESF, $in->{results_fn}) ;
   my $f2i ;
   my $header_line ;
   {
      my $headers = <RESF> ; $header_line = $headers ;
      $headers =~ s/^\#// ;
      my @headers = split(/\t/, $headers) ;
      map {$f2i->{$headers[$_]} = $_} (0 .. $#headers) ;
   }

   my $protein2dom ;
   my $protein2domfam ;
   my $protein2dom2bs_res ;
   my $protein2dom2exp_res ;
   while (my $line = <RESF>) {
      chomp $line;
      if ($line eq $header_line) { #skip if an (internal) header line
         next;
      } elsif ($line =~ /^#/) { #domain was skipped

         $line =~ s/^\#// ;
         my @t = split(/\t/, $line) ;
         $protein2dom->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$t[$f2i->{'classtype'}]} = $t[$f2i->{'class'}] ;
         $protein2domfam->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'class'}]}++ ;

      } elsif ($line =~ /	exp	/) { #exposed residue list

         my @t = split(/\t/, $line) ;
         if ($#t >= $f2i->{'residues'} && $t[$f2i->{'residues'}] ne '') {
            my @res = split(/\, /, $t[$f2i->{'residues'}]) ;
            map {$protein2dom2exp_res->{$t[$f2i->{'classtype'}]}->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$_}++ } @res ;
         }

      } else { #annotation line

         my @t = split(/\t/, $line) ;
         $protein2dom->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$t[$f2i->{'classtype'}]} = $t[$f2i->{'class'}] ;
         $protein2domfam->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'class'}]}++ ;

         my @res = split(/\, /, $t[$f2i->{'residues'}]) ;

         map {$protein2dom2bs_res->{$t[$f2i->{'classtype'}]}->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$t[$f2i->{'bs_type'}]}->{$_}++ } @res ;
      }
   }
   close(RESF) ;

   my $results ;
   my @classtypes = keys %{$protein2dom2bs_res} ;
   $results->{all}->{num_proteins} = keys %{$protein2dom} ;
   $results->{all}->{num_domains} = 0 ;
   $results->{all}->{num_residues} = 0 ;


   {
      my @summary_headers= ('OVERLAP_SUMMARY', 'sequence_id',
                            'num_exposed_res',
                            'num_ligand_binding_res',
                            'num_protein_binding_res',
                            'num_bifunctional_res',
                            'Pvalue_right', 'Pvalue_left',
                            'overlap_score',
                            'domain_families') ;
      print '#'.join("\t", @summary_headers)."\n" ;
   }

# Iterate over all proteins
   foreach my $seq_id (sort keys %{$protein2dom}) {

      $results->{all}->{num_domains} += keys %{$protein2dom->{$seq_id}} ;

      my $res_list ; #tally of binding site annotated residues
      my $curseq_nums = {} ; # Keep track of E, L, P, and LP positions

# Count number of bi-functional positions
      $curseq_nums->{fam} = {
            'L_and_protein' => 0,
            'E' => 0,
            'L' => 0,
            'protein' => 0,
            'P' => 0,
            'p' => 0,
      } ;

# Iterate over all domains
      my $noexpres_fl = 0 ;
      foreach my $res_range (keys %{$protein2dom->{$seq_id}}) {
         my $curdom_res_list ; #res_list for the current domain

# Count all residues in domain
         foreach my $subrange (split(',', $res_range)) {
            if ($subrange =~ /\-/) {
               my ($start, $stop) = ($subrange =~ /^([0-9]+)\-([0-9]+)$/);
               map {$res_list->{all}->{$_}++} ($start .. $stop) ;
            } else {
               $res_list->{all}->{$subrange}++ ;
            }
         }

# Count family/superfamily type
         foreach my $classtype (keys %{$protein2dom->{$seq_id}->{$res_range}}){
            $results->{all}->{'list_'.$classtype}->{$protein2dom->{$seq_id}->{$res_range}->{$classtype}}++ ;
         }

# Count up annotated binding sites: protein, domains, fam/sf, residues
         foreach my $classtype (keys %{$protein2dom2bs_res}) {
            if (!exists $protein2dom2bs_res->{$classtype}->{$seq_id} ||
             !exists $protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}){
               next;}
            foreach my $bs_type (keys %{$protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}}) {

               my @cur_bs_types = ('annotated', $bs_type);
               if ($bs_type eq 'P' || $bs_type eq 'p') {
                  push @cur_bs_types, 'protein' ; }

               foreach my $cur_bs_type (@cur_bs_types) {
                  map {$curdom_res_list->{bs}->{$classtype}->{$cur_bs_type}->{$_}++; }
                 (keys
   %{$protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}->{$bs_type}}) ;

                  $results->{bs}->{$classtype}->{$cur_bs_type}->{'list_proteins'}->{$seq_id}++;
                  $results->{bs}->{$classtype}->{$cur_bs_type}->{'list_domains'}->{$seq_id."\t".$res_range}++ ;
                  $results->{bs}->{$classtype}->{$cur_bs_type}->{"list_$classtype"}->{$protein2dom->{$seq_id}->{$res_range}->{$classtype}}++ ;
               }
            }

# count up bi-functional = L and (P or p) residues
            if (exists $curdom_res_list->{bs}->{$classtype}->{L} &&
                exists $curdom_res_list->{bs}->{$classtype}->{protein}){
            foreach my $res (keys %{$curdom_res_list->{bs}->{$classtype}->{L}}) {
               if (exists $curdom_res_list->{bs}->{$classtype}->{protein}->{$res}){
                  $curdom_res_list->{bs}->{$classtype}->{"L_and_protein"}->{$res}++ ;
               }
            }
            }

# If any bi-functional residues, count protein/domain/family
            if (exists $curdom_res_list->{bs}->{$classtype}->{"L_and_protein"}){
               $results->{bs}->{$classtype}->{"L_and_protein"}->{'list_proteins'}->{$seq_id}++;
               $results->{bs}->{$classtype}->{"L_and_protein"}->{'list_domains'}->{$seq_id."\t".$res_range}++ ;
               $results->{bs}->{$classtype}->{"L_and_protein"}->{"list_$classtype"}->{$protein2dom->{$seq_id}->{$res_range}->{$classtype}}++ ;
            }

#Transfer residues onto master res-list
            foreach my $bs_type (keys %{$curdom_res_list->{bs}->{$classtype}}){
               map {$res_list->{bs}->{$classtype}->{$bs_type}->{$_}++}
                  keys %{$curdom_res_list->{bs}->{$classtype}->{$bs_type}} ;}

         }

#100608_1101 - Count number of exposed and predicted binding site residues
         foreach my $classtype (keys %{$protein2dom2exp_res}) {
            if (exists $protein2dom2exp_res->{$classtype}->{$seq_id} &&
            exists $protein2dom2exp_res->{$classtype}->{$seq_id}->{$res_range}){
               $curseq_nums->{$classtype}->{E} += keys
                 %{$protein2dom2exp_res->{$classtype}->{$seq_id}->{$res_range}};
            } else {
               print "#WARNING: exposed residues not reported for ".
                     "$seq_id ($res_range)\n" ;
               $noexpres_fl = 1 ;
            }
         }
      }


      foreach my $classtype (keys %{$res_list->{bs}}) {
# count up residues per binding site type
         foreach my $bs_type (keys %{$res_list->{bs}->{$classtype}}) {
            $curseq_nums->{$classtype}->{$bs_type} =
               keys %{$res_list->{bs}->{$classtype}->{$bs_type}} ;

            $results->{bs}->{$classtype}->{$bs_type}->{num_residues} +=
               keys %{$res_list->{bs}->{$classtype}->{$bs_type}} ;
         }
      }

# Calculate p-value of observed overlap or non-overlap for this protein

      my $overlap_score = ''; 
      my $overlap_pvals = {
         "L_P" => '',
         "L_P_less" => '',
      } ;
      if ($noexpres_fl == 0) {
      $overlap_score = 0 ;
      if ($curseq_nums->{fam}->{L} > 0 && $curseq_nums->{fam}->{protein} > 0) {
         my $overlap_score_raw = (($curseq_nums->{fam}->{"L_and_protein"} *
                                   $curseq_nums->{fam}->{E}) /
                                  ($curseq_nums->{fam}->{L} *
                                   $curseq_nums->{fam}->{protein})) ;
         $overlap_score = sprintf("%.3f", $overlap_score_raw) ;
      }
      }

      my @curseq_outvals = ("OVERLAP_SUMMARY", $seq_id,
                            $curseq_nums->{fam}->{E},
                            $curseq_nums->{fam}->{L},
                            $curseq_nums->{fam}->{protein},
                            $curseq_nums->{fam}->{"L_and_protein"},
                            $overlap_pvals->{"L_P"},
                            $overlap_pvals->{"L_P_less"},
                            $overlap_score,
                            join(";",sort keys
                                 %{$protein2domfam->{$seq_id}})) ;
      print '#'.join("\t", @curseq_outvals)."\n" ;

      $results->{all}->{num_residues} += keys %{$res_list->{all}} ;
   }

# Count final summary statistics
   foreach my $classtype (keys %{$results->{bs}}) {

# count up number of families/superfamilies/domains/proteins
      foreach my $bstype (keys %{$results->{bs}->{$classtype}}) {
         $results->{bs}->{$classtype}->{$bstype}->{num_proteins} =
            keys %{$results->{bs}->{$classtype}->{$bstype}->{list_proteins}} ;

         $results->{bs}->{$classtype}->{$bstype}->{num_domains} =
            keys %{$results->{bs}->{$classtype}->{$bstype}->{list_domains}} ;

         $results->{bs}->{$classtype}->{$bstype}->{'num_'.$classtype} =
           keys %{$results->{bs}->{$classtype}->{$bstype}->{'list_'.$classtype}};
      }
   }

# Display overall counts
   print "\n\n#OVERALL COUNTS\t(note: p=peptide, P=domain-domain, L=small molecule ligand):\n" ;
   foreach my $field  (keys %{$results->{all}}) {
      if ($field =~ /list/) {
         my $newfield = $field ;$newfield =~ s/^list_/num_/ ;
         my $val = keys %{$results->{all}->{$field}} ;
         $results->{all}->{$newfield} = $val ;
         print "#ALL $newfield: $val\n" ;
      } else {
         print "#ALL $field: ".$results->{all}->{$field}."\n";
      }
   }

# Display counts for each kind of binding site
   foreach my $class_type (keys %{$results->{bs}}) {
      foreach my $bs_type (keys %{$results->{bs}->{$class_type}}) {
         foreach my $field (keys %{$results->{bs}->{$class_type}->{$bs_type}}){
            if ($field !~ /^num/) {next;}
         print "#BS: $class_type $bs_type $field ".
            $results->{bs}->{$class_type}->{$bs_type}->{$field}."\n"; }}}

# Format summary numbers as a LaTeX table
   print "\n\n" ;
   print '#'.join(" & ", 'CLASSTYPE', 'group', 'num_proteins',
               'num_domains',
               'num_fam',
               'num_residues')."\\\\\n" ;
   print '#'.join(" & ", '', 'annotated_domains',
               $results->{all}->{num_proteins},
               $results->{all}->{num_domains},
               $results->{all}->{num_fam},
               $results->{all}->{num_residues})."\\\\\n";
   foreach my $class_type (keys %{$results->{bs}}) {
      foreach my $bs_type (keys %{$results->{bs}->{$class_type}}) {
         print '#'.join(" & ", $class_type, $bs_type,
         $results->{bs}->{$class_type}->{$bs_type}->{num_proteins},
         $results->{bs}->{$class_type}->{$bs_type}->{num_domains},
         $results->{bs}->{$class_type}->{$bs_type}->{num_fam},
         $results->{bs}->{$class_type}->{$bs_type}->{num_residues},
         )."\\\\\n" ;
      }
   }

}


sub _list_standardres_20aa {

   return {
   'A' => 1,
   'R' => 1,
   'N' => 1,
   'D' => 1,
   'C' => 1,
   'Q' => 1,
   'E' => 1,
   'G' => 1,
   'H' => 1,
   'I' => 1,
   'L' => 1,
   'K' => 1,
   'M' => 1,
   'F' => 1,
   'P' => 1,
   'S' => 1,
   'T' => 1,
   'W' => 1,
   'Y' => 1,
   'V'=> 1,
   } ;

}


sub make_paper_figures {
# fpd100614_1150: Routine to generate MS text numbers, figures, and tables using
#                 homolobind (1) benchmark # results and (2) predictions

# -------------------------------------------------------
# |                    Goal data                        |
# -------------------------------------------------------
#  Abs-N1: Avg TPR for ligand bs @ 1%FPR
#  Abs-N2: Avg TPR for protein bs @ 1%FPR
#  Abs-N3: # human proteins with sig high # bi-func
#  Res-N1: # lig BS templates
#  Res-N2: # pep BS templates
#  Res-N3: # inter-dom BS templates
#  Res-N4: # intra-dom BS templates
#  Res-N5: avg coverage per family
#  Res-N6: avg seqid thresh for lig BS @ 1%FPR
#  Res-N7: avg seqid thresh for pep BS @ 1%FPR
#  Res-N8: avg seqid thresh for inter-dom BS @ 1%FPR
#  Res-N9: avg seqid thresh for intra-dom BS @ 1%FPR
#  Res-N10: Table 1 #s
#  Res-N11: Abs-N1 (# human seq with sig high # bi-func)
#  Res-N12: # human proteins with sig low # bi-func
#  Dis-N0: Thresholds for template/transfer
#  Dis-N1: Res-N1
#  Dis-N2: Res-N2
#  Dis-N3: Res-N3
#  Dis-N4: Res-N4
#  Dis-N5: # lig BS without 1%FPR threshold
#  Dis-N6: # pep BS without 1%FPR threshold
#  Dis-N7: # inter-dom BS without 1%FPR threshold
#  Dis-N8: # intra-dom BS without 1%FPR threshold
#  Fig-1D: density plot of TPR confidence interval width
# -------------------------------------------------------

   my $in = shift ;
   my $standardres_20aa =  _list_standardres_20aa();

# R pipe for making plots
   my $r_pipe ;
#   open($r_pipe, "|R --no-save > /dev/null");
   open($r_pipe, "|R --no-save ") ;

   my $gnuplot_pipe ;
   open($gnuplot_pipe, "|gnuplot - > /dev/null");

# R pipe for calculating p-values
   require R;
   require RReferences ;
   R::initR("--silent") ;

# SET NECESSARY FILE LOCATIONS
   my $pilig_specs = homolobind::pilig::set_pilig_specs() ;
   my $data_fn ;
   $data_fn->{benchmark_fn} = $homolobind_specs->{benchmark_results} ;

# MAKE FIGURE 2 - BENCHMARKING RESULTS
# Read in benchmark results (need for Fig 2A-D)
# Also: Call R to figure out mean of actual TPR mean/lower/upper

   if (exists $in->{fig2}) {

#    Prep data files
      my ($t_fn, $t_fh) ;
      foreach my $fig (qw/fam_coverage tpr seqid_thresh/) {
         foreach my $btype (qw/L p Pinter Pintra/) {
            ($t_fh->{$fig."-".$btype}, $t_fn->{$fig."-".$btype}) =
               tempfile("$fig-$btype.XXXXX", SUFFIX =>".temp"); }}

      if ($data_fn->{benchmark_fn} =~ /.gz$/) {
         open(BMF, "zcat ".$data_fn->{benchmark_fn}." |") ;
      } else {
         open(BMF, $data_fn->{benchmark_fn}) ;
      }

      my $num_fpr1_undef = {};
      while (my $line = <BMF>) {
         chomp $line;
         my @t = split(/\t/, $line) ;
         my $btype = $t[1] ;

         if ($line =~ /^FAM/) {
            print {$t_fh->{"fam_coverage-$btype"}} $t[5]."\n" ;

            if ($#t > 5) { # ie, if more than 1 BS in fam so TPR calc possible
               print {$t_fh->{"tpr-$btype"}} join("\t",
                  $t[6], $t[7], $t[8])."\n";
            }

         } elsif ($line =~ /^BS/) {
            if ($t[-6] ne 'UNDEF') {
               print {$t_fh->{"seqid_thresh-$btype"}} $t[-6]."\n" ;
            } else {
               $num_fpr1_undef->{$t[1]}++ ;
            }
         }
      }

# Dis-N5: # lig BS without 1%FPR threshold
# Dis-N6: # pep BS without 1%FPR threshold
# Dis-N7: # inter-dom BS without 1%FPR threshold
# Dis-N8: # intra-dom BS without 1%FPR threshold
      map { print "Number of $_ BS without a seqid thresh yielding 1% FPR:".
                  $num_fpr1_undef->{$_}."\n";} (keys %{$num_fpr1_undef}) ;
      map { close($t_fh->{$_}) ; } (keys %{$t_fh}) ;

# Plot Fig 2A: Coverage histogram
# Abs-N1: Avg TPR for ligand bs @ 1%FPR
# Abs-N2: Avg TPR for protein bs @ 1%FPR
# Res-N5: avg coverage per family
# Res-N1: # lig BS templates
# Res-N2: # pep BS templates
# Res-N3: # inter-dom BS templates
# Res-N4: # intra-dom BS templates
print {$r_pipe} 'coverage_lig<-read.table("'.$t_fn->{"fam_coverage-L"}.'");' ;
print {$r_pipe} 'coverage_pinter<-read.table("'.$t_fn->{"fam_coverage-Pinter"}.'");' ;
print {$r_pipe} 'coverage_pintra<-read.table("'.$t_fn->{"fam_coverage-Pintra"}.'");' ;
print {$r_pipe} 'coverage_pep<-read.table("'.$t_fn->{"fam_coverage-p"}.'");' ;

print {$r_pipe} 'cat("Mean coverage of tmpl ligand BS", mean(coverage_lig$V1),"\n");';
print {$r_pipe} 'cat("Mean coverage of tmpl peptide BS", mean(coverage_pep$V1),"\n");';
print {$r_pipe} 'cat("Mean coverage of tmpl inter-dom BS", mean(coverage_pinter$V1),"\n");';
print {$r_pipe} 'cat("Mean coverage of tmpl intra-dom BS", mean(coverage_pintra$V1),"\n");';

print {$r_pipe} 'cat("Number of tmpl ligand BS with covg data", length(coverage_lig$V1),"\n");';
print {$r_pipe} 'cat("Number of tmpl peptide BS with covg data", length(coverage_pep$V1),"\n");';
print {$r_pipe} 'cat("Number of tmpl inter-dom BS with covg data", length(coverage_pinter$V1),"\n");';
print {$r_pipe} 'cat("Number of tmpl intra-dom BS with covg data", length(coverage_pintra$V1),"\n");';

print {$r_pipe} 'pdf(file="fam_coverage_histo.pdf" ,height=3.5, width=5);'  ;
print {$r_pipe} 'par(mar=c(3,3,0.5,0.1), mgp=c(2,1,0));'  ;
print {$r_pipe} 'plot(density(coverage_lig$V1),lty=1,xlab="covered fraction of binding residues", main="",cex.lab=1.2,cex.axis=1.3, col="black", lwd=2);' ;
print {$r_pipe} 'lines(density(coverage_pintra$V1), lty=1,col="orange", lwd=2);' ;
print {$r_pipe} 'lines(density(coverage_pinter$V1), lty=1,col="purple", lwd=2);'  ;
print {$r_pipe} 'lines(density(coverage_pep$V1), lty=1,col="cyan", lwd=2);'  ;
print {$r_pipe} 'legend("left", legend = c("Small molecule", "Inter-molecular domain", "Intra-molecular domain", "Domain-peptide"), lty = c(1, 1, 1, 1), lwd = c(3 ,3, 3, 3), pch = c(NA,NA, NA, NA), col = c("black","purple","orange","cyan"), pt.bg = c(NA,NA, NA, NA), pt.cex = 1.5, pt.lwd = 1, inset = .01, cex = 1.2, adj = 0,box.lty=0);'  ;
print {$r_pipe} 'rug(coverage_lig$V1,col="black", pos=7.5);'  ;
print {$r_pipe} 'rug(coverage_pinter$V1,col="purple", pos=7.1);'  ;
print {$r_pipe} 'rug(coverage_pintra$V1,col="orange",pos=6.7);'  ;
print {$r_pipe} 'rug(coverage_pep$V1,col="cyan",pos=6.3);'  ;
print {$r_pipe} 'dev.off();'  ;


# Plot Fig 2B: Sequence identity threshold
print {$r_pipe} 'seqid_lig<-read.table("'.$t_fn->{"seqid_thresh-L"}.'");' ;
print {$r_pipe} 'seqid_pinter<-read.table("'.$t_fn->{"seqid_thresh-Pinter"}.'");' ;
print {$r_pipe} 'seqid_pintra<-read.table("'.$t_fn->{"seqid_thresh-Pintra"}.'");' ;
print {$r_pipe} 'seqid_pep<-read.table("'.$t_fn->{"seqid_thresh-p"}.'");' ;

# Res-N6: avg seqid thresh for lig BS @ 1%FPR
# Res-N7: avg seqid thresh for pep BS @ 1%FPR
# Res-N8: avg seqid thresh for inter-dom BS @ 1%FPR
# Res-N9: avg seqid thresh for intra-dom BS @ 1%FPR
print {$r_pipe} 'cat("Mean seqid thresh for lig BS", mean(seqid_lig$V1),"\n");';
print {$r_pipe} 'cat("Mean seqid thresh for pep BS", mean(seqid_pep$V1),"\n");';
print {$r_pipe} 'cat("Mean seqid thresh for inter-dom BS", mean(seqid_pinter$V1),"\n");';
print {$r_pipe} 'cat("Mean seqid thresh for intra-dom BS", mean(seqid_pintra$V1),"\n");';

print {$r_pipe} 'cat("Number of tmpl ligand BS with 1%FPR: ", length(seqid_lig$V1),"\n");';
print {$r_pipe} 'cat("Number of tmpl peptide BS with 1%FPR: ", length(seqid_pep$V1),"\n");';
print {$r_pipe} 'cat("Number of tmpl inter-dom BS with 1%FPR: ", length(seqid_pinter$V1),"\n");';
print {$r_pipe} 'cat("Number of tmpl intra-dom BS with 1%FPR: ", length(seqid_pintra$V1),"\n");';


print {$r_pipe} 'pdf(file="seqid_thresh_histo.pdf" ,height=3.5, width=5);'  ;
print {$r_pipe} 'par(mar=c(3,3,0.5,0.1), mgp=c(2,1,0));'  ;
print {$r_pipe} 'plot(density(seqid_pinter$V1),lty=1,xlab="sequence identity threshold", main="",cex.lab=1.2,cex.axis=1.3, lwd=2, col="purple",xlim=c(0,1));'  ;
print {$r_pipe} 'lines(density(seqid_lig$V1), lty=1,col="black", lwd=2);' ;
print {$r_pipe} 'lines(density(seqid_pintra$V1), lty=1,col="orange", lwd=2);' ;
print {$r_pipe} 'lines(density(seqid_pep$V1), lty=1,col="cyan", lwd=2);' ;
print {$r_pipe} 'legend("right", legend = c("Small molecule", "Inter-molecular domain", "Intra-molecular domain", "Domain-peptide"), lty = c(1, 1, 1, 1), lwd = c(3, 3, 3, 3), pch = c(NA,NA, NA, NA), col = c("black","purple","orange","cyan"), pt.bg = c(NA,NA, NA, NA), pt.cex = 1.5, pt.lwd = 1, inset = .01, cex = 1.2, adj = 0,box.lty=0);'  ;
print {$r_pipe} 'rug(seqid_lig$V1,col="black", pos=5.5);'  ;
print {$r_pipe} 'rug(seqid_pinter$V1,col="purple", pos=5.2);'  ;
print {$r_pipe} 'rug(seqid_pintra$V1,col="orange",pos=4.9);'  ;
print {$r_pipe} 'rug(seqid_pep$V1,col="cyan",pos=4.6);'  ;
print {$r_pipe} 'dev.off();'."\n"  ;

# Plot Fig 2C: TPR
print {$r_pipe} 'tpr_lig<-read.table("'.$t_fn->{"tpr-L"}.'");
tpr_pinter<-read.table("'.$t_fn->{"tpr-Pinter"}.'");
tpr_pintra<-read.table("'.$t_fn->{"tpr-Pintra"}.'");
tpr_pep<-read.table("'.$t_fn->{"tpr-p"}.'");

pdf(file="tpr_histo.pdf" ,height=3.5, width=5);
par(mar=c(3,3,0.5,0.1), mgp=c(2,1,0));
plot(density(tpr_pep$V1),lty=1,xlab="true positive rate", main="",cex.lab=1.2,cex.axis=1.3, lwd=2, col="cyan");
lines(density(tpr_pinter$V1), lty=1,col="purple", lwd=2);
lines(density(tpr_pintra$V1), lty=1,col="orange", lwd=2);
lines(density(tpr_lig$V1), lty=1,col="black", lwd=2);
legend("left", legend = c("Small molecule", "Inter-molecular domain", "Intra-molecular domain", "Domain-peptide"), lty = c(1, 1, 1, 1), lwd = c(3 ,3, 3, 3), pch = c(NA,NA, NA, NA), col = c("black","purple","orange","cyan"), pt.bg = c(NA,NA, NA, NA), pt.cex = 1.5, pt.lwd = 1, inset = .01, cex = 1.2, adj = 0,box.lty=0);
rug(tpr_lig$V1,col="black", pos=26);
rug(tpr_pinter$V1,col="purple", pos=24.5);
rug(tpr_pintra$V1,col="orange",pos=23);
rug(tpr_pep$V1,col="cyan",pos=21.5);
dev.off();'."\n" ;

# Plot Fig 2D: TPR CI interval width density
print {$r_pipe} 'pdf(file="tpr_95CI_width_histo.pdf" ,height=3.5, width=5);
par(mar=c(3,3,0.5,0.1), mgp=c(2,1,0));

tpr_pinter_ciwidth<-(tpr_pinter$V3 - tpr_pinter$V2);
tpr_pintra_ciwidth<-(tpr_pintra$V3 - tpr_pintra$V2);
tpr_pep_ciwidth<-(tpr_pep$V3 - tpr_pep$V2);
tpr_lig_ciwidth<-(tpr_lig$V3 - tpr_lig$V2);

plot(density(tpr_pep_ciwidth),lty=1,xlab="width of true positive rate 95% confidence interval", main="",cex.lab=1.2,cex.axis=1.3, lwd=2, col="cyan");
lines(density(tpr_pinter_ciwidth), lty=1,col="purple", lwd=2);
lines(density(tpr_pintra_ciwidth), lty=1,col="orange", lwd=2);
lines(density(tpr_lig_ciwidth), lty=1,col="black", lwd=2);
legend("right", legend = c("Small molecule", "Inter-molecular domain", "Intra-molecular domain", "Domain-peptide"), lty = c(1, 1, 1, 1), lwd = c(3 ,3, 3, 3), pch = c(NA,NA, NA, NA), col = c("black","purple","orange","cyan"), pt.bg = c(NA,NA, NA, NA), pt.cex = 1.5, pt.lwd = 1, inset = .01, cex = 1.2, adj = 0,box.lty=0);
rug(tpr_lig_ciwidth,col="black", pos=30);
rug(tpr_pinter_ciwidth ,col="purple", pos=28.5);
rug(tpr_pintra_ciwidth,col="orange",pos=27);
rug(tpr_pep_ciwidth,col="cyan",pos=25.5);
dev.off();'."\n" ;

   }


   if (exists $in->{fig3} ||
       exists $in->{fig3a} || exists $in->{fig3b} || exists $in->{fig3c} ||
       exists $in->{tab1} || exists $in->{tab2} || exists $in->{tab3} ||
       exists $in->{tab4}) {

# Abs-N3: # human proteins with sig high # bi-func
# Res-N10: Table 1 #s
# Res-N11: # human proteins with sig high # bi-func
# Res-N12: # human proteins with sig low # bi-func


      my $histograms ;
# 1. Read in SUPERFAMILY function assignments
      my $superfam_fxn = readin_superfam_fxn() ;

# Takes as input: (1) SUPERFAMILY assignment file, (2) HOMOLOBIND output file
      my $ass_fn = $in->{ass_fn} ;
      my $res_fn = $in->{results_fn} ;

      my $seqid_list = {};

# Step 1. Read in predictions
# Sort ASS file by seqid so that sequence can be read in sequentially
      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{sorted_ass}, $temp_fn->{sorted_ass}) = tempfile(
        "sorted_ass.XXXXX", SUFFIX => ".temp") ; close($temp_fh->{sorted_ass});
      system("sort -k2,2 $ass_fn > ".$temp_fn->{sorted_ass}) ;

# Read in header line for results file
      my $f2i ;
      my $header_line ;
      open(RESF, $in->{results_fn}) ;
      {
         my $headers = <RESF> ; $header_line = $headers ;
         $headers =~ s/^\#// ;
         my @headers = split(/\t/, $headers) ;
         map {$f2i->{$headers[$_]} = $_} (0 .. $#headers) ;
      }

# 1a. Read and store all residue sets for each target domain
      my $protein2dom ;
      my $protein2domfam ;
      my $protein2dom2bs_res ;
      my $protein2dom2exp_res ;
      while (my $line = <RESF>) {
         chomp $line;
         if ($line eq $header_line) { #skip if an (internal) header line
            next;
         } elsif ($line =~ /^#/) { #domain was skipped

            $line =~ s/^\#// ;
            my @t = split(/\t/, $line) ;
            $protein2dom->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$t[$f2i->{'classtype'}]} = $t[$f2i->{'class'}] ;
            $protein2domfam->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'class'}]}++ ;

         } elsif ($line =~ /	exp	/) { #exposed residue list

            my @t = split(/\t/, $line) ;
            if ($#t >= $f2i->{'residues'} && $t[$f2i->{'residues'}] ne '') {
               my @res = split(/\, /, $t[$f2i->{'residues'}]) ;
               map {$protein2dom2exp_res->{$t[$f2i->{'classtype'}]}->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$_}++ } @res ;
            }

         } else { #annotation line

            my @t = split(/\t/, $line) ;
            $protein2dom->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$t[$f2i->{'classtype'}]} = $t[$f2i->{'class'}] ;
            $protein2domfam->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'class'}]}++ ;

            my @res = split(/\, /, $t[$f2i->{'residues'}]) ;

            map {$protein2dom2bs_res->{$t[$f2i->{'classtype'}]}->{$t[$f2i->{'seq_id'}]}->{$t[$f2i->{'res_range'}]}->{$t[$f2i->{'bs_type'}]}->{$_}++ } @res ;
         }
      }
      close(RESF) ;

# Step 2. Read in ASS file for sequence info; if any annotations for that
# Numbers summary code from summarize-results
# Should handle: Table 1/S1, S2, Fig 2B
# Have to adjust to handle residue (F2A, TS3) and function (F2C, TS4) propensity
      my $results ;
      my @classtypes = keys %{$protein2dom2bs_res} ;
      $results->{all}->{num_proteins} = keys %{$protein2dom} ;
      $results->{all}->{num_domains} = 0 ;
      $results->{all}->{num_residues} = 0 ;

# Open sorted ASS file to read in sequences sequentially
      open(ASSF, $temp_fn->{sorted_ass}) ;

# Set up array to hold pvalues for top-ten table (S2)
      my $overlap_pvals_histos ;
# Set up files to hold overlap-scores for all/sig over-/sig under-lap
      my ($t_fn, $t_fh) ;
      foreach my $fig (qw/overlap_score overlap_score_sig_over overlap_score_sig_under/) {
         ($t_fh->{$fig}, $t_fn->{$fig}) =
            tempfile("$fig.XXXXX", SUFFIX =>".temp"); }


# Iterate over all proteins
      {

#must set use locale so that the LC_COLLATE is respected
# and the hash sort is in the same order as the unix sorted ASS file

      use locale;
      foreach my $seq_id (sort keys %{$protein2dom}) {
# Load protein sequence from ASS file
         my $curseq_info = {};
#         print STDERR "looking up seq for $seq_id\n" ;
         while (1) {
            my $ass_line = <ASSF>; chomp $ass_line;
            my @t = split(/\t/, $ass_line) ;
            if ($t[1] lt $seq_id) {
               next;
            } elsif ($t[1] gt $seq_id) {
               die "ERROR: no sequence found for $seq_id; now on $ass_line\n" ;
            }
            $curseq_info->{aaseq} = $t[5] ;
            $curseq_info->{aaseq} =~ s/\-//g ;
            $curseq_info->{aaseq} = uc($curseq_info->{aaseq}) ;
            last;
         }

         $results->{all}->{num_domains} += keys %{$protein2dom->{$seq_id}} ;
         my $res_list ; #tally of binding site annotated residues

# Keep track of E, L, P, and LP positions
         my $curseq_nums = {} ;
# Count number of bi-functional positions
         $curseq_nums->{fam} = {
            'L_and_protein' => 0,
            'E' => 0,
            'L' => 0,
            'protein' => 0,
            'P' => 0,
            'p' => 0,
         } ;

# Iterate over all domains
         my $noexpres_fl = 0 ;
         foreach my $res_range (sort keys %{$protein2dom->{$seq_id}}) {
            my $curdom_res_list; #res_list for the current_domain

# Count all residues in domain
            foreach my $subrange (split(',', $res_range)) {
               if ($subrange =~ /\-/) {
                  my ($start, $stop) = ($subrange =~ /^([0-9]+)\-([0-9]+)$/);
                  map {$res_list->{all}->{$_}++} ($start .. $stop) ;
               } else {
                  $res_list->{all}->{$subrange}++ ;
               }
            }

# Count family/superfamily type
            foreach my $classtype (keys %{$protein2dom->{$seq_id}->{$res_range}}){
               $results->{all}->{'list_'.$classtype}->{$protein2dom->{$seq_id}->{$res_range}->{$classtype}}++ ;
            }

# Count up annotated binding sites: protein, domains, fam/sf, residues
            foreach my $classtype (keys %{$protein2dom2bs_res}) {
               if (!exists $protein2dom2bs_res->{$classtype}->{$seq_id} ||
                !exists $protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}){
                  next;}
               foreach my $bs_type (keys %{$protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}}) {

                  my @cur_bs_types = ('annotated', $bs_type);
                  if ($bs_type eq 'P' || $bs_type eq 'p') {
                     push @cur_bs_types, 'protein' ; }

                  foreach my $cur_bs_type (@cur_bs_types) {
                     map {$curdom_res_list->{bs}->{$classtype}->{$cur_bs_type}->{$_}++}
                        (keys %{$protein2dom2bs_res->{$classtype}->{$seq_id}->{$res_range}->{$bs_type}});

                     $results->{bs}->{$classtype}->{$cur_bs_type}->{'list_proteins'}->{$seq_id}++;
                     $results->{bs}->{$classtype}->{$cur_bs_type}->{'list_domains'}->{$seq_id."\t".$res_range}++ ;
                     $results->{bs}->{$classtype}->{$cur_bs_type}->{"list_$classtype"}->{$protein2dom->{$seq_id}->{$res_range}->{$classtype}}++ ;
                  }
               }

# count up L and (P or p) residues (can do protein/domain/family check post-hoc)
            if (exists $curdom_res_list->{bs}->{$classtype}->{L} &&
                exists $curdom_res_list->{bs}->{$classtype}->{protein}) {
            foreach my $res (keys %{$curdom_res_list->{bs}->{$classtype}->{L}}) {
               if (exists $curdom_res_list->{bs}->{$classtype}->{protein}->{$res}) {
                  $curdom_res_list->{bs}->{$classtype}->{"L_and_protein"}->{$res}++ ; }
            }
            }

# If any bi-functional residues, count protein/domain/family
            if (exists $curdom_res_list->{bs}->{$classtype}->{"L_and_protein"}){
               $results->{bs}->{$classtype}->{"L_and_protein"}->{'list_proteins'}->{$seq_id}++;
               $results->{bs}->{$classtype}->{"L_and_protein"}->{'list_domains'}->{$seq_id."\t".$res_range}++ ;
               $results->{bs}->{$classtype}->{"L_and_protein"}->{"list_$classtype"}->{$protein2dom->{$seq_id}->{$res_range}->{$classtype}}++ ;
            }

#Transfer residues onto master res-list
            foreach my $bs_type (keys %{$curdom_res_list->{bs}->{$classtype}}){
               map {$res_list->{bs}->{$classtype}->{$bs_type}->{$_}++}
                  keys %{$curdom_res_list->{bs}->{$classtype}->{$bs_type}} ;}

            }

#100608_1101 - Count number of exposed and predicted binding site residues
            foreach my $classtype (keys %{$protein2dom2exp_res}) {
               if (exists $protein2dom2exp_res->{$classtype}->{$seq_id} &&
               exists $protein2dom2exp_res->{$classtype}->{$seq_id}->{$res_range}){
                  $curseq_nums->{$classtype}->{E} += keys
                    %{$protein2dom2exp_res->{$classtype}->{$seq_id}->{$res_range}};
               } else {
                  print "#WARNING: exposed residues not reported for ".
                        "$seq_id ($res_range)\n" ;
                  $noexpres_fl = 1 ;
               }
            }
         }

# Tally residue type count for exposed residues
         foreach my $classtype (keys %{$protein2dom2exp_res}) {
            foreach my $res_range (keys %{$protein2dom2exp_res->{$classtype}->{$seq_id}}) {
            map {$histograms->{aa_freq}->{$classtype}->{"E"}->{substr($curseq_info->{aaseq}, $_, 1)}++ ;} (keys %{$protein2dom2exp_res->{$classtype}->{$seq_id}->{$res_range}}); }}

         foreach my $classtype (keys %{$res_list->{bs}}) {

# using current L/protein/L_and_protein list also make L-only and P-only lists
#  for residue propensity calcs.

# L-only
            if (exists $res_list->{bs}->{$classtype}->{L}) {
               if (!exists $res_list->{bs}->{$classtype}->{"L_and_protein"}){
                  map {$res_list->{bs}->{$classtype}->{"L-only"}->{$_}++;}
                     (keys %{$res_list->{bs}->{$classtype}->{L}}) ;
               } else {
                  foreach my $res (keys %{$res_list->{bs}->{$classtype}->{L}}) {
                     if (!exists
                     $res_list->{bs}->{$classtype}->{"L_and_protein"}->{$res}) {
         $res_list->{bs}->{$classtype}->{"L-only"}->{$res}++;}
                  }
               }
            }

# P-only
            if (exists $res_list->{bs}->{$classtype}->{protein}) {
               if (!exists $res_list->{bs}->{$classtype}->{"L_and_protein"}){
                  map {$res_list->{bs}->{$classtype}->{"protein-only"}->{$_}++;}
                     (keys %{$res_list->{bs}->{$classtype}->{protein}}) ;
               } else {
                  foreach my $res (keys %{$res_list->{bs}->{$classtype}->{protein}}) {
                     if (!exists
                     $res_list->{bs}->{$classtype}->{"L_and_protein"}->{$res}) {
         $res_list->{bs}->{$classtype}->{"protein-only"}->{$res}++;} }

               }
            }


# count up residues per binding site type
            foreach my $bs_type (keys %{$res_list->{bs}->{$classtype}}) {
# Tally amino acid frequency per class
#               print STDERR "trying to tally residue type from positions: ".join(", ", sort {$a <=> $b} keys %{$res_list->{bs}->{$classtype}->{$bs_type}})."\n" ;
#               print STDERR " FROM seq: ".$curseq_info->{aaseq}."\n\n";
               map {$histograms->{aa_freq}->{$classtype}->{$bs_type}->{substr($curseq_info->{aaseq}, ($_ - 1), 1)}++ ;} (keys %{$res_list->{bs}->{$classtype}->{$bs_type}}) ;

               $curseq_nums->{$classtype}->{$bs_type} =
                  keys %{$res_list->{bs}->{$classtype}->{$bs_type}} ;

               $results->{bs}->{$classtype}->{$bs_type}->{num_residues} +=
                  keys %{$res_list->{bs}->{$classtype}->{$bs_type}} ;
            }
         }

# Calculate p-value of observed overlap or non-overlap for this protein
         my $overlap_pvals = {
            "L_P" => '',
            "L_P_less" => '',
         };
         my $overlap_score = '' ;

         if ($curseq_nums->{fam}->{L} > 0 && $curseq_nums->{fam}->{protein} > 0){
            $seqid_list->{"L_and_protein_binding"}->{$seq_id} = {} ;}

         if ($noexpres_fl == 0) {
         { # code from pibase::pilig.pm snapshot
            my $tnorm ;
            $tnorm->{lp} = $curseq_nums->{fam}->{"E"} -
                           $curseq_nums->{fam}->{"L"} -
                           $curseq_nums->{fam}->{"protein"} +
                           $curseq_nums->{fam}->{"L_and_protein"} ;

            $tnorm->{Lp} = $curseq_nums->{fam}->{"L"} -
                           $curseq_nums->{fam}->{"L_and_protein"} ;

            $tnorm->{lP} = $curseq_nums->{fam}->{"protein"} -
                           $curseq_nums->{fam}->{"L_and_protein"} ;

            $tnorm->{LP} = $curseq_nums->{fam}->{"L_and_protein"} ;

# Check to ensure any negative values set to 0; occurs in cases where
#  a residue not considered exposed is part of a binding site;
#  in the n=30K human proteins, this happens n=3 times: a single residue
#  part of the binding site is not considered exposed.
# (ENSP00000356121, ENSP00000356122, ENSP00000375801)
            foreach my $type (keys %{$tnorm}) {
               if ($tnorm->{$type} < 0) {
                  print STDERR "Note: $seq_id |$type| = ".
                               $tnorm->{$type}."; setting to 0.\n";
                  $tnorm->{$type} = 0 ;
               }
            }

#         print STDERR "E=".$curseq_nums->{fam}->{"E"}." ".
#                      "lp=".$tnorm->{lp}." ".
#                      "Lp=".$tnorm->{Lp}." ".
#                      "lP=".$tnorm->{lP}." ".
#                      "LP=".$tnorm->{LP}."\n" ;

            my @x = &R::eval("capture.output(fisher.test(matrix(c($tnorm->{lp}, $tnorm->{Lp}, $tnorm->{lP}, $tnorm->{LP}), 2, 2),workspace=2e7,alternative='g'))") ;
            my $x = join("\n", @x) ;
            my ($pval) = ($x =~ /p-value [\<=] (.*)\nalternative/) ;
            $overlap_pvals->{"L_P"} = $pval ;

            my @y = &R::eval("capture.output(fisher.test(matrix(c($tnorm->{lp}, $tnorm->{Lp}, $tnorm->{lP}, $tnorm->{LP}), 2, 2),workspace=2e7,alternative='l'))") ;
            my $y = join("\n", @y) ;
            my ($pval_less) = ($y =~ /p-value [\<=] (.*)\nalternative/) ;
            $overlap_pvals->{"L_P"."_less"} = $pval_less ;
         }

         $overlap_score = 0 ;
         if ($curseq_nums->{fam}->{L} > 0 && $curseq_nums->{fam}->{protein} > 0){

            my $overlap_score_raw =
               (($curseq_nums->{fam}->{"L_and_protein"} *
                 $curseq_nums->{fam}->{E}) /
                ($curseq_nums->{fam}->{L} *
                 $curseq_nums->{fam}->{protein})) ;
            $overlap_score = sprintf("%.3f", $overlap_score_raw) ;
            push @{$overlap_pvals_histos->{label}}, $seq_id ;
            push @{$overlap_pvals_histos->{pval}}, $overlap_pvals->{"L_P"};
            push @{$overlap_pvals_histos->{pval_less}},
               $overlap_pvals->{"L_P_less"};

#ORIGINAL POSITION -100916_1415 
            $seqid_list->{"L_and_protein_binding"}->{$seq_id} = {
               overlap => $overlap_score,
               pval => $overlap_pvals->{"L_P"},
               pval_less => $overlap_pvals->{"L_P_less"}
            } ;

            print {$t_fh->{overlap_score}} $overlap_score_raw."\n";
            if ($overlap_pvals->{"L_P_less"} < 0.01) {
               $seqid_list->{sig_under}->{$seq_id} = {
                  overlap => $overlap_score,
                  pval => $overlap_pvals->{"L_P"},
                  pval_less => $overlap_pvals->{"L_P_less"}
               } ;
               print {$t_fh->{overlap_score_sig_under}}
                  $overlap_score_raw."\n";
            }

            if ($overlap_pvals->{"L_P"} < 0.01) {
               $seqid_list->{sig_over}->{$seq_id} = {
                  overlap => $overlap_score,
                  pval => $overlap_pvals->{"L_P"},
                  pval_less => $overlap_pvals->{"L_P_less"}
               } ;
               print {$t_fh->{overlap_score_sig_over}}
                  $overlap_score_raw."\n";
            }
         }
         }

         my @curseq_outvals = ("OVERLAP_SUMMARY", $seq_id,
                               $curseq_nums->{fam}->{E},
                               $curseq_nums->{fam}->{L},
                               $curseq_nums->{fam}->{protein},
                               $curseq_nums->{fam}->{"L_and_protein"},
                               $overlap_pvals->{"L_P"},
                               $overlap_pvals->{"L_P_less"},
                               $overlap_score,
                               join(";",sort keys
                                    %{$protein2domfam->{$seq_id}})) ;
         print '#'.join("\t", @curseq_outvals)."\n" ;

         $results->{all}->{num_residues} += keys %{$res_list->{all}} ;
      }
      map {close($t_fh->{$_});} (keys %{$t_fh}) ;
      close(ASSF) ;
      }

# Make Table 3
      if (exists $in->{tab3}) {

      my @pval_right_sorted = sort {
         $seqid_list->{sig_over}->{$a}->{pval} <=>
         $seqid_list->{sig_over}->{$b}->{pval} ||
         $seqid_list->{sig_over}->{$b}->{overlap} <=>
         $seqid_list->{sig_over}->{$a}->{overlap} }
                              (keys %{$seqid_list->{sig_over}}) ;

      my @pval_left_sorted = sort {
         $seqid_list->{sig_under}->{$a}->{pval_less} <=>
         $seqid_list->{sig_under}->{$b}->{pval_less} ||
         $seqid_list->{sig_under}->{$a}->{overlap} <=>
         $seqid_list->{sig_under}->{$b}->{overlap} }
                              (keys %{$seqid_list->{sig_under}}) ;


      print "Table 3.\n" ;
      print '\begin{tabular}{|llrr|}'."\n";
      print '\hline'."\n";
      print 'Protein & overlap & significance\\\\'."\n";
      print '\hline'."\n";
      print '\hline'."\n";
      print '\multicolumn{3}{|l|}{\it Most significant overlapping proteins} & P-value (right)\\\\'."\n";
      print '\hline'."\n";
      foreach my $j ( 0 .. 49) {
         my $curseq = $pval_right_sorted[$j] ;
         print join(" & ", $curseq, ' ',
            $seqid_list->{sig_over}->{$curseq}->{overlap},
            $seqid_list->{sig_over}->{$curseq}->{pval})."\\\\\n" ; }

      print '\hline'."\n";
      print '\hline'."\n";
      print '\multicolumn{3}{|l|}{\it Most significant non-overlapping proteins} & P-value (left)'."\n";
      print '\hline'."\n";
      foreach my $j ( 0 .. 49) {
         my $curseq = $pval_left_sorted[$j] ;
         print join(" & ", $curseq, ' ',
            $seqid_list->{sig_under}->{$curseq}->{overlap},
            $seqid_list->{sig_under}->{$curseq}->{pval_less})."\\\\\n" ; }
      print '\hline'."\n";
      print '\end{tabular}'."\n";
      }

# Count final summary statistics
      foreach my $classtype (keys %{$results->{bs}}) {

# count up number of families/superfamilies/domains/proteins
         foreach my $bstype (keys %{$results->{bs}->{$classtype}}) {
            $results->{bs}->{$classtype}->{$bstype}->{num_proteins} =
               keys %{$results->{bs}->{$classtype}->{$bstype}->{list_proteins}} ;

            $results->{bs}->{$classtype}->{$bstype}->{num_domains} =
               keys %{$results->{bs}->{$classtype}->{$bstype}->{list_domains}} ;

            $results->{bs}->{$classtype}->{$bstype}->{'num_'.$classtype} =
              keys %{$results->{bs}->{$classtype}->{$bstype}->{'list_'.$classtype}};
         }
      }

# Display overall counts
      foreach my $field  (keys %{$results->{all}}) {
         if ($field =~ /list/) {
            my $newfield = $field ;$newfield =~ s/^list_/num_/ ;
            my $val = keys %{$results->{all}->{$field}} ;
            $results->{all}->{$newfield} = $val ;
            print "ALL $newfield: $val\n" ;
         } else {
            print "ALL $field: ".$results->{all}->{$field}."\n";
         }
      }

# Display counts for each kind of binding site
      foreach my $class_type (keys %{$results->{bs}}) {
         foreach my $bs_type (keys %{$results->{bs}->{$class_type}}) {
            foreach my $field (keys %{$results->{bs}->{$class_type}->{$bs_type}}){
               if ($field !~ /^num/) {next;}
            print "BS: $class_type $bs_type $field ".
               $results->{bs}->{$class_type}->{$bs_type}->{$field}."\n"; }}}

# Format summary numbers as a LaTeX table
      print "\n\n" ;
      print join(" & ", 'CLASSTYPE', 'group', 'num_proteins',
                  'num_domains',
                  'num_fam',
                  'num_residues')."\\\\\n" ;
      print join(" & ", '', 'annotated_domains',
                  $results->{all}->{num_proteins},
                  $results->{all}->{num_domains},
                  $results->{all}->{num_fam},
                  $results->{all}->{num_residues})."\\\\\n";
      foreach my $class_type (keys %{$results->{bs}}) {
         foreach my $bs_type (keys %{$results->{bs}->{$class_type}}) {
            print join(" & ", $class_type, $bs_type,
            $results->{bs}->{$class_type}->{$bs_type}->{num_proteins},
            $results->{bs}->{$class_type}->{$bs_type}->{num_domains},
            $results->{bs}->{$class_type}->{$bs_type}->{num_fam},
            $results->{bs}->{$class_type}->{$bs_type}->{num_residues},
            )."\\\\\n" ;
         }
      }

# Plot Fig 3A: Residue type propensity
      if (exists $in->{fig3} || exists $in->{fig3a}) {
      foreach my $btype (keys %{$histograms->{aa_freq}->{fam}}) {
      $histograms->{aa_freq}->{fam}->{$btype}->{total} = 0 ;
      foreach my $aatype (keys %{$histograms->{aa_freq}->{fam}->{$btype}}) {
         if ($aatype eq 'total') {next;}
         if (!exists $standardres_20aa->{$aatype}) {next;}
         $histograms->{aa_freq}->{fam}->{$btype}->{total} +=
            $histograms->{aa_freq}->{fam}->{$btype}->{$aatype} ; }}

      foreach my $btype (keys %{$histograms->{aa_freq}->{fam}}) {
         ($t_fh->{aa_freq}->{$btype},
          $t_fn->{aa_freq}->{$btype}) = tempfile("aahisto_nongap.$btype.XXXXX",
                                                  SUFFIX => ".temp") ;
         print STDERR "now on aa propensities for $btype\n" ;
         foreach my $aatype (sort
            keys %{$histograms->{aa_freq}->{fam}->{$btype}}) {
            if (!exists $standardres_20aa->{$aatype}) {next;}

            if ($aatype eq 'total') {next;}
            my $p_aa_given_btype =
                $histograms->{aa_freq}->{fam}->{$btype}->{$aatype} /
                $histograms->{aa_freq}->{fam}->{$btype}->{total} ;

            my $p_aa_given_exp =
                ($histograms->{aa_freq}->{fam}->{'E'}->{$aatype} /
                $histograms->{aa_freq}->{fam}->{'E'}->{total}) ;

            my $prop_aa_btype = $p_aa_given_btype / $p_aa_given_exp ;
            my @outvals = ($aatype,
                           sprintf("%.4f",$p_aa_given_btype),
                           sprintf("%.4f",$p_aa_given_exp),
                           sprintf("%.4f",$prop_aa_btype)) ;
            print {$t_fh->{aa_freq}->{$btype}} join("\t", @outvals)."\n" ;
         }
         close($t_fh->{aa_freq}->{$btype}) ;
      }

      my $gcom_resprop = '' ;

      $gcom_resprop .='set term postscript enhanced solid "Helvetica, 22"'."\n";
      $gcom_resprop .='set output "residue_propensity.'.$$.'.eps'."\n";

      $gcom_resprop .='set style data histogram'."\n";
      $gcom_resprop .='set style histogram cluster gap 1'."\n";
      $gcom_resprop .='set style fill solid 1 noborder'."\n";
      $gcom_resprop .='set boxwidth 0.75'."\n";
      $gcom_resprop .='set xlabel "amino acid residue type"'."\n";
      $gcom_resprop .='set ylabel "propensity vs all exposed residues"'."\n";
      $gcom_resprop .='show label'."\n";
      $gcom_resprop .='set ylabel offset 2'."\n";
      $gcom_resprop .='set yrange [0.6:1.55]'."\n";
      $gcom_resprop .='set y2range [-0.4:0.55]'."\n";

      $gcom_resprop .='set grid ytics'."\n";
      $gcom_resprop .='set tics scale 0'."\n";
      $gcom_resprop .='set arrow 1 from -0.5,0.6 to -0.5, 1.5 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 2 from 0.5,0.6 to 0.5, 1.5 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 3 from 1.5,0.6 to 1.5, 1.5 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 4 from 2.5,0.6 to 2.5, 1.5 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 5 from 3.5,0.6 to 3.5, 1.5 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 6 from 4.5,0.6 to 4.5, 1.5 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 7 from 5.5,0.6 to 5.5, 1.5 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 8 from 6.5,0.6 to 6.5, 1.5 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 9 from 7.5,0.6 to 7.5, 1.5 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 10 from 8.5,0.6 to 8.5, 1.35 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 11 from 9.5,0.6 to 9.5, 1.35 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 12 from 10.5,0.6 to 10.5, 1.35 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 13 from 11.5,0.6 to 11.5, 1.35 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 14 from 12.5,0.6 to 12.5, 1.35 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 15 from 13.5,0.6 to 13.5, 1.35 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 16 from 14.5,0.6 to 14.5, 1.35 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 17 from 15.5,0.6 to 15.5, 1.35 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 18 from 16.5,0.6 to 16.5, 1.35 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 19 from 17.5,0.6 to 17.5, 1.5 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 20 from 18.5,0.6 to 18.5, 1.5 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 21 from 19.5,0.6 to 19.5, 1.5 nohead lt 1'."\n";
      $gcom_resprop .='set arrow 22 from -0.5,1.4 to 7.5, 1.4 nohead ls 0'."\n";
      $gcom_resprop .='set arrow 23 from 17.5,1.4 to 19.5, 1.4 nohead ls 0'."\n";
      $gcom_resprop .='set key at 17.7,1.505 font "Helvetica, 18"'."\n";
      $gcom_resprop .='plot [-0.5: 19.5] 1 t "" lw 1 lt rgb "black" axis x1y1, "'.$t_fn->{aa_freq}->{"L-only"}.'" using ($4 - 1):xtic(1) t "ligand binding only" lt rgb "cyan" axis x1y2, "'.$t_fn->{aa_freq}->{"protein-only"}.'" using ($4 - 1):xtic(1) t "protein binding only" lt rgb "orange" axis x1y2, "'.$t_fn->{aa_freq}->{L_and_protein}.'" using ($4 - 1):xtic(1) lt rgb "black" axis x1y2 t "ligand and protein binding"'."\n";
      $gcom_resprop .='set terminal x11'."\n";

print "GNUPLOT COMMAND FOR Fig 3A:\n$gcom_resprop\n\n";
print {$gnuplot_pipe} $gcom_resprop ;

      my $num_restype = {} ;
      foreach my $restype (keys %{$histograms->{aa_freq}->{fam}}) {
         $num_restype->{$restype} =
            $histograms->{aa_freq}->{fam}->{$restype}->{total} ; }

# Table 2. Bootstrap the residue propensities to get confidence intervals.
# Collect all commands into a single string, then send to the
# R.pm interface to capture.output(). then parse into latex table format.
      my $rcom_resprop_bootstrap = '
prob_exp <- read.table("'.$t_fn->{aa_freq}->{"E"}.'")
prob_lp <- read.table("'.$t_fn->{aa_freq}->{"L_and_protein"}.'")
prob_l <- read.table("'.$t_fn->{aa_freq}->{"L-only"}.'")
prob_p <- read.table("'.$t_fn->{aa_freq}->{"protein-only"}.'")

aa<-levels(prob_exp$V1)

orig_propensity_l <- prob_l$V4
orig_propensity_p <- prob_p$V4
orig_propensity_lp <- prob_lp$V4

propensity_distr_l <- vector(mode="list",length=length(aa))
propensity_distr_p <- vector(mode="list",length=length(aa))
propensity_distr_lp <- vector(mode="list",length=length(aa))

for (trial in 1:1000) {
   lp_sample<-sample(prob_lp$V1,'.$num_restype->{"L_and_protein"}.', replace=TRUE,prob_lp$V2)
   l_sample<-sample(prob_l$V1,'.$num_restype->{"L-only"}.',replace=TRUE,prob_l$V2)
   p_sample<-sample(prob_p$V1,'.$num_restype->{"protein-only"}.',replace=TRUE,prob_p$V2)
   exp_sample<-sample(prob_exp$V1,'.$num_restype->{"E"}.',replace=TRUE,prob_exp$V2)

   ft_lp <- table(lp_sample)
   ft_p <- table(p_sample)
   ft_l <- table(l_sample)
   ft_exp <- table(exp_sample)

   for (aatype in 1:length(aa)) {
      count_lp <- 0
      if (! is.na(ft_lp[aa[aatype]])){
         count_lp <- ft_lp[aa[aatype]] }

      count_l <- 0
      if (! is.na(ft_l[aa[aatype]])){
         count_l <- ft_l[aa[aatype]] }

      count_p <- 0
      if (! is.na(ft_p[aa[aatype]])){
         count_p <- ft_p[aa[aatype]] }

      count_exp <- 0
      if (! is.na(ft_exp[aa[aatype]])){
         count_exp <- ft_exp[aa[aatype]] }

      if (count_exp == 0) {
         if (count_l == 0) {
            curprop_l <- 1
         } else {
            curprop_l <- Inf
         }

         if (count_p == 0) {
            curprop_p <- 1
         } else {
            curprop_p <- Inf
         }

         if (count_lp == 0) {
            curprop_lp <- 1
         } else {
            curprop_lp <- Inf
         }
      } else {
         curprop_l <- (count_l /length(l_sample))/ (count_exp / length(exp_sample))
         curprop_lp <- (count_lp /length(lp_sample))/ (count_exp / length(exp_sample))
         curprop_p <- (count_p /length(p_sample))/ (count_exp / length(exp_sample))
      }

      propensity_distr_l[[aatype]]<-c(propensity_distr_l[[aatype]], curprop_l)
      propensity_distr_lp[[aatype]]<-c(propensity_distr_lp[[aatype]], curprop_lp)
      propensity_distr_p[[aatype]]<-c(propensity_distr_p[[aatype]], curprop_p)
   }
}

for (aatype in 1:length(aa)) {
   q_l_lower <- quantile(propensity_distr_l[[aatype]],0.025)
   q_l_upper <- quantile(propensity_distr_l[[aatype]],0.975)

   q_p_lower <- quantile(propensity_distr_p[[aatype]],0.025)
   q_p_upper <- quantile(propensity_distr_p[[aatype]],0.975)

   q_lp_lower <- quantile(propensity_distr_lp[[aatype]],0.025)
   q_lp_upper <- quantile(propensity_distr_lp[[aatype]],0.975)

   sig_l <- ""
   if ((q_l_lower < 1 && q_l_upper < 1) ||
       (q_l_lower > 1 && q_l_upper > 1)) {sig_l <- "*"}

   sig_p <- ""
   if ((q_p_lower < 1 && q_p_upper < 1) ||
       (q_p_lower > 1 && q_p_upper > 1)) {sig_p <- "*"}

   sig_lp <- ""
   if ((q_lp_lower < 1 && q_lp_upper < 1) ||
       (q_lp_lower > 1 && q_lp_upper > 1)) {sig_lp <- "*"}

   cat(aa[[aatype]], "&",
       round(orig_propensity_l[aatype],3), "&", sig_l, "(", round(q_l_lower, 3), ",", round(q_l_upper, 3), ")", "&",
       round(orig_propensity_p[aatype],3), "&", sig_p, "(", round(q_p_lower, 3), ",", round(q_p_upper, 3), ")", "&",
       round(orig_propensity_lp[aatype],3), "&", sig_lp, "(", round(q_lp_lower, 3), ",", round(q_lp_upper, 3), ")\\\\",
   "\n")

}
' ;
      print "R COMMANDS FOR Table 2:\n\n$rcom_resprop_bootstrap\n";

      print {$r_pipe} $rcom_resprop_bootstrap ;
# should print table to STDOUT

      }

      if (exists $in->{fig3} || exists $in->{fig3b}) {
# Plot Fig 3B: Overlap score distributions
         my $rcom_fig3b = '';
 $rcom_fig3b .= 'overlap_all<-read.table("'.$t_fn->{"overlap_score"}.'");'."\n" ;
 $rcom_fig3b .= 'overlap_over<-read.table("'.$t_fn->{"overlap_score_sig_over"}.'");'."\n" ;
 $rcom_fig3b .= 'overlap_under<-read.table("'.$t_fn->{"overlap_score_sig_under"}.'");'."\n" ;

 $rcom_fig3b .= 'pdf(file="overlap_score_histo.pdf" ,height=3.5, width=5);'."\n" ;
 $rcom_fig3b .= 'par(mar=c(3.1,3.1,0.1,0.1), mgp=c(2,1,0));'."\n"  ;
 $rcom_fig3b .= 'hist(overlap_all$V1,freq=FALSE, col="darkgrey",main="",xlab="protein and ligand binding site overlap score", ylim=c(0,2.7),xlim=c(0,10),border="darkgrey",cex.lab=1.3,cex.axis=1.2, breaks=c(100));'."\n" ;
 $rcom_fig3b .= 'rug(overlap_all$V1);'."\n" ;

 $rcom_fig3b .= 'lines(density(overlap_over$V1),lty=1,lwd=2);'."\n" ;
 $rcom_fig3b .= 'lines(density(overlap_under$V1),lty=2,lwd=2);'."\n" ;
 $rcom_fig3b .= 'legend("right", legend = c("all proteins that bind ligands and proteins", "overlapping proteins (p < 0.01)", "non-overlapping proteins (p < 0.01)"), lty = c(0, 1, 2), lwd = c(0, 2, 2), pch = c(22, NA, NA), pt.bg = c("grey", NA, NA), pt.cex = 1.5, pt.lwd = 1, inset = .01, cex = 1.0, adj = 0,box.lty=0);'."\n" ;
 $rcom_fig3b .= "dev.off();\n" ;
      print "R COMMANDS FOR Fig 3B:\n\n$rcom_fig3b\n";
      print {$r_pipe} $rcom_fig3b ;
      }


      if (exists $in->{fig3} || exists $in->{fig3c}) {
# Calc Function propensities
         my $fam_fxncounts ;
         my $fam_fxnbroadcounts ;
         foreach my $seq_type (keys %{$seqid_list}) {
            foreach my $seq_id (keys %{$seqid_list->{$seq_type}}) {
               foreach my $resrange (keys %{$protein2dom->{$seq_id}}) {
                  my $curfam = $protein2dom->{$seq_id}->{$resrange}->{fam} ;
                  my ($cursf) = ($curfam =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
         if (!exists $superfam_fxn->{sf2fxn_broad}->{$cursf}) {
print STDERR "Warning: no function found for $seq_id $resrange $curfam ($cursf)\n";
         } else {
            $fam_fxnbroadcounts->{$seq_type}->{$superfam_fxn->{sf2fxn_broad}->{$cursf}}++;
            $fam_fxnbroadcounts->{"totals_$seq_type"}++ ;

            $fam_fxncounts->{$seq_type}->{$superfam_fxn->{sf2fxn}->{$cursf}}++;
            $fam_fxncounts->{"totals_$seq_type"}++ ;
         }
               }
            }
         }

      ($t_fh->{fxnprop_detail}, $t_fn->{fxnprop_detail}) =
         tempfile("fxn_propensities_detailed.sigunder_sigover.XXXXX",
                  SUFFIX => ".temp") ;
   my $fxn_propensities ;
   foreach my $t_fxn ( sort keys %{$fam_fxncounts->{"L_and_protein_binding"}}) {
      foreach my $o_type (qw/sig_under sig_over/) {
         if (!exists $fam_fxncounts->{$o_type}->{$t_fxn}) {
            $fam_fxncounts->{$o_type}->{$t_fxn} = 0 ;
            $fxn_propensities->{$o_type}->{$t_fxn} = 0 ;
         } else {
            $fxn_propensities->{$o_type}->{$t_fxn} =
            ($fam_fxncounts->{$o_type}->{$t_fxn} / 
             $fam_fxncounts->{"totals_".$o_type}) /
            ($fam_fxncounts->{"L_and_protein_binding"}->{$t_fxn} / 
             $fam_fxncounts->{"totals_L_and_protein_binding"}) ;
         }
      }
      print {$t_fh->{fxnprop_detail}} join("\t", $t_fxn,
                           $fam_fxncounts->{sig_under}->{$t_fxn},
                           $fam_fxncounts->{totals_sig_under},
                           $fam_fxncounts->{sig_over}->{$t_fxn},
                           $fam_fxncounts->{totals_sig_over},
                           $fam_fxncounts->{"L_and_protein_binding"}->{$t_fxn},
                           $fam_fxncounts->{"totals_L_and_protein_binding"},
                           $fxn_propensities->{sig_under}->{$t_fxn},
                           $fxn_propensities->{sig_over}->{$t_fxn})."\n";
   }
   close($t_fh->{fxnprop_detail}) ;

      ($t_fh->{fxnprop_broad}, $t_fn->{fxnprop_broad}) =
         tempfile("fxn_propensities_broad.sigunder_sigover.XXXXX",
                  SUFFIX => ".temp") ;
   my $fxnbroad_propensities ;
   foreach my $t_fxn ( sort keys %{$fam_fxnbroadcounts->{"L_and_protein_binding"}}) {
      foreach my $o_type (qw/sig_under sig_over/) {
         if (!exists $fam_fxnbroadcounts->{$o_type}->{$t_fxn}) {
            $fam_fxnbroadcounts->{$o_type}->{$t_fxn} = 0 ;
            $fxnbroad_propensities->{$o_type}->{$t_fxn} = 0 ;
         } else {
            $fxnbroad_propensities->{$o_type}->{$t_fxn} =
            ($fam_fxnbroadcounts->{$o_type}->{$t_fxn} / 
             $fam_fxnbroadcounts->{"totals_".$o_type}) /
            ($fam_fxnbroadcounts->{"L_and_protein_binding"}->{$t_fxn} / 
             $fam_fxnbroadcounts->{"totals_L_and_protein_binding"}) ;
         }
      }
      print {$t_fh->{fxnprop_broad}} join("\t", $t_fxn,
                           $fam_fxnbroadcounts->{sig_under}->{$t_fxn},
                           $fam_fxnbroadcounts->{totals_sig_under},
                           $fam_fxnbroadcounts->{sig_over}->{$t_fxn},
                           $fam_fxnbroadcounts->{totals_sig_over},
                     $fam_fxnbroadcounts->{"L_and_protein_binding"}->{$t_fxn},
                     $fam_fxnbroadcounts->{"totals_L_and_protein_binding"},
                           $fxnbroad_propensities->{sig_under}->{$t_fxn},
                           $fxnbroad_propensities->{sig_over}->{$t_fxn})."\n";
   }
   close($t_fh->{fxnprop_broad}) ;

my $gcom_fig3c = '
set term postscript enhanced solid "Helvetica, 18"
set output "fxn_propensity.'.$$.'.eps"

set style data histogram
set style histogram cluster gap 1
set style fill solid 1 border
set boxwidth=0.75
set xlabel "function"
set ylabel "propensity"
show label
set yrange [0.2:1.8]
set y2range [-0.8:0.8]

set grid ytics
set xtics nomirror rotate by -45
set tics scale 0
set arrow 1 from -0.5,0.2 to -0.5, 1.8 nohead lt 1
set arrow 2 from 0.5,0.2 to 0.5, 1.8 nohead lt 1
set arrow 3 from 1.5,0.2 to 1.5, 1.8 nohead lt 1
set arrow 4 from 2.5,0.2 to 2.5, 1.8 nohead lt 1
set arrow 5 from 3.5,0.2 to 3.5, 1.6 nohead lt 1
set arrow 6 from 4.5,0.2 to 4.5, 1.6 nohead lt 1
set arrow 7 from 5.5,0.2 to 5.5, 1.6 nohead lt 1
set key at 6.6,1.805 font "Helvetica, 18"

plot [-0.5:6.5] 1 t "" lw 1 lt rgb "black" axis x1y1, "'.$t_fn->{fxnprop_broad}.'" using ($9 - 1):xtic(1) t "significantly overlapping" lt rgb "grey" axis x1y2, "'.$t_fn->{fxnprop_broad}.'" using ($8 - 1):xtic(1) t "significantly non-overlapping" axis x1y2 
set terminal x11
' ;
      print "GNUPLOT commands for Fig 3C: $gcom_fig3c\n\n" ;
      print {$gnuplot_pipe} $gcom_fig3c ;

      my $rcom_tab4 = '
a <- read.table("'.$t_fn->{fxnprop_broad}.'")
p_under<-(a$V2/a$V3)
p_over<-(a$V4/a$V5)
p_all<-(a$V6/a$V7)

orig_propensity_under <- ((a$V2/a$V3) / (a$V6 / a$V7))
orig_propensity_over <- ((a$V4/a$V5) / (a$V6 / a$V7))

fxns<-levels(a$V1)
propensity_distr_under <- vector(mode="list",length=length(fxns))
propensity_distr_over <- vector(mode="list",length=length(fxns))
for (trial in 1:1000) {
   underfam_sample<-sample(fxns,a$V3[[1]], replace=TRUE,p_under)
   overfam_sample<-sample(fxns,a$V5[[1]], replace=TRUE,p_over)
   allfam_sample<-sample(fxns,a$V7[[1]], replace=TRUE,p_all)

   ft_under <- table(underfam_sample)
   ft_over <- table(overfam_sample)
   ft_all <- table(allfam_sample)

   for (fxntype in 1:length(fxns)) {
      count_allfam <- 0
      if (! is.na(ft_all[fxns[fxntype]])){
         count_allfam <- ft_all[fxns[fxntype]] }

      count_overfam <- 0
      if (! is.na(ft_over[fxns[fxntype]])){
         count_overfam <- ft_over[fxns[fxntype]] }

      count_underfam <- 0
      if (! is.na(ft_under[fxns[fxntype]])){
         count_underfam <- ft_under[fxns[fxntype]] }

      if (count_allfam == 0) {
         if (count_underfam == 0) {
            curprop_under <- 1
         } else {
            curprop_under <- Inf
         }

         if (count_overfam == 0) {
            curprop_over <- 1
         } else {
            curprop_over <- Inf
         }
      } else {
         curprop_under <- (count_underfam /length(underfam_sample))/ (count_allfam / length(allfam_sample))
         curprop_over <- (count_overfam /length(overfam_sample))/ (count_allfam / length(allfam_sample))
      }

      propensity_distr_under[[fxntype]]<-c(propensity_distr_under[[fxntype]], curprop_under)
      propensity_distr_over[[fxntype]]<-c(propensity_distr_over[[fxntype]], curprop_over)
   }
}

for (fxntype in 1:length(fxns)) {
   q_under_lower <- quantile(propensity_distr_under[[fxntype]],0.025)
   q_under_upper <- quantile(propensity_distr_under[[fxntype]],0.975)

   q_over_lower <- quantile(propensity_distr_over[[fxntype]],0.025)
   q_over_upper <- quantile(propensity_distr_over[[fxntype]],0.975)

   sig_under <- ""
   if ((q_under_lower < 1 && q_under_upper < 1) ||
       (q_under_lower > 1 && q_under_upper > 1)) {sig_under <- "*"}

   sig_over <- ""
   if ((q_over_lower < 1 && q_over_upper < 1) ||
       (q_over_lower > 1 && q_over_upper > 1)) {sig_over <- "*"}

   cat(fxns[[fxntype]], "&",
       round(orig_propensity_under[fxntype],3), " & ", sig_under, "(", round(q_under_lower,3), ",", round(q_under_upper,3), ")", "&",
       round(orig_propensity_over[fxntype],3), "&", sig_over, "(", round(q_over_lower,3), ",", round(q_over_upper,3), ")\\\\\n"
   )
}
' ;
      print "R COMMANDS FOR Table 4:\n\n$rcom_tab4\n\n";
      print {$r_pipe} $rcom_tab4 ;
      }

   }

# 4. Read in sequence info from SUPERFAMILY assignment file

# Table  1. Overview of predicted binding sites
# Table  2. Residue type propensities
# Table  3. Top ten most/least overlap proteins
# Table  4. Function propensities
#   Fig 2A. Coverage histogram
#   Fig 2B. Sequence identity cutoff histogram
#   Fig 2C. TPR histogram
#   Fig 2D. TPR 95% CI width histogram
#   Fig 3A. Residue type propensity histogram
#   Fig 3B. Overlap score distribution
#   Fig 3C. Function propensity

}


=head2 readin_superfam_fxn()

   Title:       readin_superfam_fxn()
   Function:    Read in SUPERFAMILY function annotation file
   Args:        ->{results_fn} = run_homolobind() output file
   Returns:     ->{sf2fxn}->{sf} = function
                ->{fxn2broad}->{detailed fxn} = broad fxn
                ->{sffxn_broad}->{sf} = broad fxn

=cut

sub readin_superfam_fxn {
# Purpose: read in SUPERFAMILY function annotation

   my $in ;

   my $sf2fxn; my $fxn2broad ; my $sf2fxn_broad;
   open(SUPERFAM_FXN, $homolobind_specs->{supfam}->{fxn}) ;
   while (my $line = <SUPERFAM_FXN>) {
      chomp $line ;
      my ($fxn, $sunid, $classtype, $class, undef,$name) = split(/\t/, $line) ;
      if ($fxn eq "S'") {$fxn = 'S';}
      $sf2fxn->{$class} = $fxn ;
   }
   close(SUPERFAM_FXN) ;

   open(SUPERFAM_FXN_CATEGORIES, $homolobind_specs->{supfam}->{fxn_categories});
   while (my $line = <SUPERFAM_FXN_CATEGORIES>) {
      chomp $line ;
      my ($broad, $fxn, $fxn_name, $fxn_descr) = split(/\t/, $line) ;
      if ($fxn eq "S'") {$fxn = 'S';}
      $fxn2broad->{$fxn} = $broad;
   }
   close(SUPERFAM_FXN_CATEGORIES) ;

   foreach my $sf (keys %{$sf2fxn}) {
      if (!exists $fxn2broad->{$sf2fxn->{$sf}}) {
         $sf2fxn_broad->{$sf} = 'N_A' ;
         print STDERR "WARNING: Fxn class for $sf ".
                      "(".$sf2fxn->{$sf}.") not found\n";
      } else {
         $sf2fxn_broad->{$sf} = $fxn2broad->{$sf2fxn->{$sf}} ;
      }
   }

   return {
      sf2fxn => $sf2fxn,
      fxn2broad => $fxn2broad,
      sf2fxn_broad => $sf2fxn_broad,
   }

}


=head2 plot_annotations()

   Title:       plot_annotations()

   Function:    Generates a postscript diagram depicting domain annotations
                and transferred binding sites; given (1) assignment file,
                (2) homolobind results file, and (3) sequence identifier

   Args:        ->{results_fn} = run_homolobind() output file
   Returns:     ->{sf2fxn}->{sf} = function
                ->{fxn2broad}->{detailed fxn} = broad fxn
                ->{sffxn_broad}->{sf} = broad fxn

=cut

sub plot_annotations {

   my $diagram_specs = {
      linewidth => 20,
      residue_scaling => 10,      # pixels per residue
      tic_interval => 10, # place tic marks every nth residue
      tic_labels => 5,    # label every nth tic mark
      tic_height => 50,
      tic_width => 10,
      tic_labypos => -120,
      domain_box_height => 600,
      color_ligand => "0.000 0.000 0.000",
      color_peptide => "0.000 0.933 0.933",
      color_Pinter => "0.569 0.173 0.933",
      color_Pintra => "0.933 0.604 0.000",
      color_bifunc => "0.000 0.933 0.463",
      bstype2layout => [ 'bifunc', 'Pintra', 'Pinter', 'peptide', 'ligand']
   } ;
   $diagram_specs->{bsrow_height} = $diagram_specs->{domain_box_height} /
                                    ($#{$diagram_specs->{bstype2layout}} + 1);
   $diagram_specs->{circle_radius} = $diagram_specs->{residue_scaling} / 2;

# 0. figure out which sequences to analyze
   my $in = shift ;
   my $goal_seqids ;
   if (exists $in->{seq_id})  {
      $goal_seqids->{$in->{seq_id}}++ ;
   } elsif (exists $in->{seq_id_fn}) {
      open(SEQIDF, $in->{seq_id_fn}) ;
      while (my $line = <SEQIDF>) {
         chomp $line;
         $goal_seqids->{$line}++ ;
      }
      close(SEQIDF) ;
   } else {
      die "Must specify desired sequences with either -seq_id or -seq_id_fn" ;
   }
      my $supfam_specs = set_supfam_specs() ;

   my $fig_range = {}; #residue range if specified on command line
   if (exists $in->{range}) {
      ($fig_range->{start}, $fig_range->{end}) =
         ($in->{range} =~ /([0-9]+)\-([0-9]+)/); }


# 1. Parse domain assignments from SUPERFAMILY ass file
      if (!exists $in->{ass_fn}) {
         die "Must specify -ass_fn"; }
      my $seqid_info ;
      open(ASSF, $in->{ass_fn}) ;
      while (my $line = <ASSF>) {
         chomp $line;
         if ($line =~ /^#/) {next;}
         my @t = split(/\t/, $line) ;
         my $curass ;
         map {$curass->{$supfam_specs->{ass_file_headers}->[$_]} = $t[$_]}
            (0 .. $#t) ;
         if (!exists $goal_seqids->{$curass->{seq_id}}) {next;}


         my $cur_aaseq = $curass->{alnstring} ;
         $cur_aaseq =~ s/\-//g ;
         my $seqlength = length($cur_aaseq) ;

         $seqid_info->{$curass->{seq_id}}->{seqlength} = $seqlength ;
         $seqid_info->{$curass->{seq_id}}->{domains}->{$curass->{res_range}}++ ;
      }
      close(ASSF) ;

# 2. Parse predicted binding sites from HOMOLOBIND output file
      if (!exists $in->{results_fn}) {
         die "Must specify -results_fn" ; }
      {
         my $f2i ;
         my $header_line ;
         open(RESF, $in->{results_fn}) ;
         {
            my $headers = <RESF> ; $header_line = $headers ;
            $headers =~ s/^\#// ;
            my @headers = split(/\t/, $headers) ;
            map {$f2i->{$headers[$_]} = $_} (0 .. $#headers) ;
         }

         my $protein2dom2bs_res ;
         while (my $line = <RESF>) {
            chomp $line;
            if ($line eq $header_line) { #skip if an (internal) header line
               next;
            } elsif ($line =~ /^#/) { #domain was skipped
               next;
            } elsif ($line =~ /	exp	/) { #exposed residue list
               next;
            } else { #annotation line
               my @t = split(/\t/, $line) ;
               if (!exists $goal_seqids->{$t[$f2i->{'seq_id'}]}) { next; }
               my @res = split(/\, /, $t[$f2i->{'residues'}]) ;

               my $bstypes = [ $t[$f2i->{'bs_type'}]];

               if ($bstypes->[0]  eq 'p') {
                  $bstypes->[0] = 'peptide' ;
                  $bstypes->[1] = 'protein' ;
               } elsif ($bstypes->[0] eq 'L') {
                  $bstypes->[0] = 'ligand' ;
               } elsif ($bstypes->[0] eq 'P' &&
                        $t[$f2i->{'partner_descr'}] =~ /^chains=diff/) {
                  $bstypes->[0] = 'Pinter' ;
                  $bstypes->[1] = 'protein' ;
               } elsif ($bstypes->[0] eq 'P' &&
                        $t[$f2i->{'partner_descr'}] =~ /^chains=same/) {
                  $bstypes->[0] = 'Pintra' ;
                  $bstypes->[1] = 'protein' ;
               } else {
                  die "DONT KNOW WHAT KIND OF BS THIS IS!!: $line\n" ;
               }
               foreach my $bstype (@{$bstypes}) {
                  map {$seqid_info->{$t[$f2i->{'seq_id'}]}->{bsres}->{$bstype}->{$_}++ } @res ; }
            }
         }
         close(RESF) ;
      }

# 3. Draw plots for each sequence.
      foreach my $seq_id (keys %{$seqid_info}) {
         my ($t_fh, $t_fn) = tempfile("$seq_id.XXXXX",
                                       SUFFIX => ".homolobind_annotation.eps") ;

# set up canvas
         my $fig_res_length = $seqid_info->{$seq_id}->{seqlength} ;
         my $canvas_start = 0 ;
         if (exists $fig_range->{start}) {
            $canvas_start = ($fig_range->{start} - 1) *
                            $diagram_specs->{residue_scaling} ;
            $fig_res_length = $fig_range->{end} - $fig_range->{start} + 1;
         }

         my $canvas_width = $diagram_specs->{residue_scaling} * $fig_res_length;
         my $canvas_height = $diagram_specs->{domain_box_height} +
                             2 * (abs($diagram_specs->{tic_labypos}));
         my $bbox ;
         $bbox->{xmin} = $canvas_start - ($diagram_specs->{linewidth} * 2) ;
         $bbox->{xmax} = $canvas_width + ($diagram_specs->{linewidth} * 2) +
                         800;
#         $bbox->{ymin} = 0 - ($diagram_specs->{linewidth} * 2) ;
         $bbox->{ymin} = 0 - (abs($diagram_specs->{tic_labypos}) * 2) ;
         $bbox->{ymax} = $canvas_height - ($diagram_specs->{linewidth} * 2) ;

         my $date = `date` ; chomp $date ;
         print {$t_fh} '%!PS-Adobe-3.0 EPSF-3.0'."\n";
         print {$t_fh} '%%Creator: graph2ps.pl'."\n";
         print {$t_fh} '%%Title: graph layout'."\n";
         print {$t_fh} '%%CreationDate: '.$date."\n";
         print {$t_fh} '%%DocumentData: Clean7Bit'."\n";
         print {$t_fh} '%%Origin: [0 0]'."\n";
         print {$t_fh} '%%BoundingBox: '.join(' ',($bbox->{xmin},
               $bbox->{ymin}, $bbox->{xmax}, $bbox->{ymax}))."\n";
         print {$t_fh} '%%LanguageLevel: 2'."\n";
         print {$t_fh} '%%Pages: 1'."\n";
         print {$t_fh} '%%Page: 1 1'."\n";
         print {$t_fh} '%%EOF'."\n";

# i. Draw line for sequence
         print {$t_fh} $diagram_specs->{linewidth}." setlinewidth\n" ;
         print {$t_fh} "0 0 0 setrgbcolor\n" ;
         print {$t_fh} "newpath\n" ;

         my $vertmid = $diagram_specs->{domain_box_height} / 2 ;

         print {$t_fh} "0 $vertmid moveto\n" ;
         print {$t_fh} "$canvas_width $vertmid lineto\n" ;
         print {$t_fh} "stroke\n" ;

         print {$t_fh} "/Helvetica findfont 96 scalefont setfont\n" ;

# Provide legend for each bstype
         foreach my $j ( 0 .. $#{$diagram_specs->{bstype2layout}}) {
            my $key_text = $diagram_specs->{bstype2layout}->[$j] ;
            if ($key_text eq 'bifunc') {
               $key_text = 'bi-functional' ;
            } elsif ($key_text eq 'Pintra') {
               $key_text = 'intra-mol domain' ;
            } elsif ($key_text eq 'Pinter') {
               $key_text = 'inter-mol domain' ;
            }
            my $key_x = $canvas_width + 20 ;
            my $key_y = ($j + 0.5) * $diagram_specs->{bsrow_height} ;
            print {$t_fh} "$key_x $key_y moveto ($key_text) show\n" ;
         }

# make a hash list of all tic marks necessary, then figure
# Mark those within domain box ranges, then finish the remaning ones
         my $all_tic_marks ;
         map {$all_tic_marks->{$_} = 1;} (
            0 .. POSIX::floor($seqid_info->{$seq_id}->{seqlength} /
                               $diagram_specs->{tic_interval})) ;

# ii. Draw boxes for each annotated domain
         my @res_ranges = split(/\,/, join(',',
                           keys %{$seqid_info->{$seq_id}->{domains}})) ;
         foreach my $res_range (@res_ranges) {
            my ($start_res, $end_res) = split('-', $res_range)  ;

            if (exists $fig_range->{start}) {
               if ($end_res < $fig_range->{start} ||
                   $start_res > $fig_range->{end}) {next;}
               if ($start_res < $fig_range->{start}) {
                  $start_res = $fig_range->{start}; }
               if ($end_res > $fig_range->{end}) {
                  $end_res = $fig_range->{end}; }
            }

            my $points = [
             [($start_res * $diagram_specs->{residue_scaling}), 0],
             [($start_res * $diagram_specs->{residue_scaling}),
               $diagram_specs->{domain_box_height}],
             [($end_res * $diagram_specs->{residue_scaling}),
               $diagram_specs->{domain_box_height}],
             [($end_res * $diagram_specs->{residue_scaling}), 0],
            ] ;

            print {$t_fh} $diagram_specs->{linewidth}." setlinewidth\n" ;
            print {$t_fh} "0 0 0 setrgbcolor\n" ;
            print {$t_fh} "newpath\n" ;
            print {$t_fh} $points->[0]->[0]." ".$points->[0]->[1]." moveto\n" ;
            print {$t_fh} $points->[1]->[0]." ".$points->[1]->[1]." lineto\n" ;
            print {$t_fh} $points->[2]->[0]." ".$points->[2]->[1]." lineto\n" ;
            print {$t_fh} $points->[3]->[0]." ".$points->[3]->[1]." lineto\n" ;
            print {$t_fh} "closepath\n" ;
            print {$t_fh} "gsave\n" ;
            print {$t_fh} "stroke\n" ;
            print {$t_fh} "grestore\n" ;
            print {$t_fh} "1 1 1 setrgbcolor\n" ;
            print {$t_fh} "fill\n" ;

# Draw tic marks
            my ($start_tic, $end_tic) ;
            if ($start_res % $diagram_specs->{tic_interval} > 0) {
               $start_tic = POSIX::floor($start_res /
                           $diagram_specs->{tic_interval}) + 1 ;
            } else {
               $start_tic = POSIX::floor($start_res /
                           $diagram_specs->{tic_interval}) ;
            }
            $end_tic = POSIX::floor($end_res /
                           $diagram_specs->{tic_interval}) ;

            foreach my $tic ($start_tic .. $end_tic) {
               my $tic_x = $tic * $diagram_specs->{tic_interval} *
                           $diagram_specs->{residue_scaling} ;
               my $tic_y1 = 0 - $diagram_specs->{tic_height} ;
               my $tic_y2 = 0;

               if ($tic % $diagram_specs->{tic_labels} == 0) {
                  print {$t_fh} ($diagram_specs->{tic_width} * 2).
                                 " setlinewidth\n" ;
                  my $resnum = $tic * $diagram_specs->{residue_scaling} ;
                  print {$t_fh} ($tic_x - 80)." ".$diagram_specs->{tic_labypos}.
                                " moveto ($resnum) show\n" ;
                  print {$t_fh} ($tic_x - 80)." ".$diagram_specs->{tic_labypos}.
                                " moveto ($resnum) show\n" ;
               } else {
                  print {$t_fh} $diagram_specs->{tic_width}." setlinewidth\n" ;
               }

               print {$t_fh} "0 0 0 setrgbcolor\n" ;
               print {$t_fh} "newpath\n" ;
               print {$t_fh} $tic_x." ".$tic_y1." moveto\n" ;
               print {$t_fh} $tic_x." ".$tic_y2." lineto\n" ;
               print {$t_fh} "stroke\n" ;

# Code to place thin lines across the domain box for every tic
#                  print {$t_fh} $diagram_specs->{tic_width}." setlinewidth\n" ;
#                  print {$t_fh} "newpath\n" ;
#                  print {$t_fh} $tic_x." ".$tic_y1." moveto\n" ;
#                  print {$t_fh} $tic_x." ".$diagram_specs->{domain_box_height}.
#                                " lineto\n" ;
#                  print {$t_fh} "stroke\n" ;
               $all_tic_marks->{$tic} = 2 ;
            }

# draw separating lines for the bsres
            print {$t_fh} ($diagram_specs->{linewidth} / 10)." setlinewidth\n" ;
            print {$t_fh} "0 0 0 setrgbcolor\n" ;
            foreach my $j ( 1 .. $#{$diagram_specs->{bstype2layout}}) {
               my $line_y = $j * $diagram_specs->{bsrow_height} ;
               print {$t_fh} "newpath\n" ;
               print {$t_fh} $points->[0]->[0]." $line_y moveto\n" ;
               print {$t_fh} $points->[2]->[0]." $line_y lineto\n" ;
               print {$t_fh} "stroke\n" ;
            }
         }

         foreach my $tic (sort {$a <=> $b} keys %{$all_tic_marks}) {
            if ($all_tic_marks->{$tic} == 2) {next;}
            my $tic_res = $tic * $diagram_specs->{tic_interval} ;

            if (exists $fig_range->{start}) {
               if ($tic_res < $fig_range->{start}) {next;}
               if ($tic_res > $fig_range->{end})   {last;} }

            my $tic_x = $tic_res * $diagram_specs->{residue_scaling} ;
            my $tic_y1 = $vertmid - $diagram_specs->{tic_height} ;
            my $tic_y2 = $vertmid ;
            if ($tic % $diagram_specs->{tic_labels} == 0) {
                  print {$t_fh} ($diagram_specs->{tic_width} * 2).
                                 " setlinewidth\n" ;
#                  print {$t_fh} ($tic_x - 80)." ".
#                                ($tic_y2 + $diagram_specs->{tic_labypos}).
#                                " moveto ($tic_res) show\n" ;
                  if ($tic > 0) {
                  print {$t_fh} ($tic_x - 80)." ".
                                $diagram_specs->{tic_labypos}.
                                " moveto ($tic_res) show\n" ;
                  }
            } else {
               print {$t_fh} $diagram_specs->{tic_width}." setlinewidth\n" ;
            }

            print {$t_fh} "0 0 0 setrgbcolor\n" ;
            print {$t_fh} "newpath\n" ;
            print {$t_fh} $tic_x." ".$tic_y1." moveto\n" ;
            print {$t_fh} $tic_x." ".$tic_y2." lineto\n" ;
            print {$t_fh} "stroke\n" ;
            $all_tic_marks->{$tic} = 2 ;
         }

# Make list of bi-functional residues
         if (exists $seqid_info->{$seq_id}->{bsres}->{ligand} &&
             exists $seqid_info->{$seq_id}->{bsres}->{protein}) {
            foreach my $bs_res (keys %{$seqid_info->{$seq_id}->{bsres}->{ligand}}){
               if (exists $seqid_info->{$seq_id}->{bsres}->{protein}->{$bs_res}) {
                  $seqid_info->{$seq_id}->{bsres}->{bifunc}->{$bs_res}++; }}}

# Plot filled squares for each predicted binding residue
         foreach my $j (0 .. $#{$diagram_specs->{bstype2layout}}) {
            my $bstype = $diagram_specs->{bstype2layout}->[$j] ;
            if (!exists $seqid_info->{$seq_id}->{bsres}->{$bstype}) {next;}

         if ($bstype eq 'ass') {
            my $mid_y = ($j + 0.5) * $diagram_specs->{bsrow_height} ;
            print {$t_fh} $diagram_specs->{"color_bifunc"}." setrgbcolor\n" ;
            foreach my $bsres (keys %{$seqid_info->{$seq_id}->{bsres}->{bifunc}}){
               if (exists $fig_range->{start} &&
                   ($bsres < $fig_range->{start} ||
                    $bsres > $fig_range->{end})) {next;}

               my $x = ($bsres - 0.5) * $diagram_specs->{residue_scaling} ;
               print {$t_fh} "newpath $x $mid_y ".$diagram_specs->{circle_radius}.
                             " 0 360 arc closepath fill\n" ;
            }

         } else {

            my $min_y = $j * $diagram_specs->{bsrow_height} ;
            my $max_y = $min_y + $diagram_specs->{bsrow_height} ;

            print {$t_fh} $diagram_specs->{"color_$bstype"}." setrgbcolor\n" ;
            foreach my $bsres (keys %{$seqid_info->{$seq_id}->{bsres}->{$bstype}}){
               if (exists $fig_range->{start} &&
                   ($bsres < $fig_range->{start} ||
                    $bsres > $fig_range->{end})) {next;}
               my $min_x = ($bsres - 1 ) * $diagram_specs->{residue_scaling} ;
               my $max_x = $min_x + $diagram_specs->{residue_scaling} ;
               print {$t_fh} "newpath\n" ;
               print {$t_fh} $min_x." ".$min_y." moveto\n" ;
               print {$t_fh} $min_x." ".$max_y." lineto\n" ;
               print {$t_fh} $max_x." ".$max_y." lineto\n" ;
               print {$t_fh} $max_x." ".$min_y." lineto\n" ;
               print {$t_fh} "closepath\n" ;
               print {$t_fh} "fill\n" ;
            }
         }
      }


      print {$t_fh} "showpage\n" ;
      close($t_fh) ;
   }

}


=head2 benchmark_homolobind()

   Title:       benchmark_homolobind()
   Function:    Estimate the family-specific binding site coverage and 
                template BS-specific accuracy at varying seqid thresholds

   Args:        ->{ARGV} = ARGV array reference; parsed to provide:
                ->{ass_fn} = name of SUPERFAMILY domain assignment file
                ->{out_fn} = name of output file
                ->{err_fn} = name of error file
                ->{matrix_fn} = optional substitution matrix, default BLOSUM62
                ->{cluster_fl} = run on an SGE cluster (options in SGE.pm)

   Returns:     NOTHING
   Displays:     1. seq_id

=cut

sub benchmark_homolobind {

   require Math::Random::MT ;

   my $in = shift;

# Read in command line options
   my $j = 0;
   while ($j < $#{$in->{ARGV}}) {
      $ARGV[$j] =~ s/^\-// ;
      $in->{$in->{ARGV}->[$j]} = $in->{ARGV}->[($j+1)] ;
      $j += 2 ;
   }

# Set parameters
   my $supfam_specs = set_supfam_specs() ;
   my $standardres = $supfam_specs->{standardres};
   
# Preload ASTRAL data
   print STDERR "loading ASTRAL\n" ;
   my $astral = homolobind::pilig::_pilig_astral_preload({
      specs => $homolobind_specs}) ;

# Load PIBASE and LIGBASE binding site info
   my $pb = homolobind::pilig::_pilig_tod_pibase_preload({
      specs => $homolobind_specs }) ;
   my $liginfo = homolobind::pilig::_pilig_load_liginfo({
      fn => $homolobind_specs->{pilig}->{outfiles}->{liginfo} }) ;
   my $class2alnlength = {};

   print STDERR "reading in pepnucibits\n" ;
   my $pepnucibits = homolobind::pilig::readin_pepnuciassignments({
      fn => $homolobind_specs->{pilig}->{outfiles}->{assign_pepnuci_clusters},
      clustrep_fl => 1,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   print STDERR "reading in ligbits\n" ;
   my $ligbits = homolobind::pilig::readin_ligassignments({
      fn => $homolobind_specs->{pilig}->{outfiles}->{assign_lig_clusters},
      clustrep_fl => 1,
      liginfo => $liginfo,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;

   print STDERR "reading in protein interfaces\n" ;
   my $pibits_both = homolobind::pilig::readin_piassignments({
      fn => $homolobind_specs->{pilig}->{outfiles}->{assign_pi_clusters},
      clustrep_fl => 1,
      pb => $pb,
      dont_read_jdomains_fl => 1,
      standardres => $standardres,
      class2alnlength => $class2alnlength,
   }) ;
   my $pibits = $pibits_both->{pibits} ;
   my $interfaces = $pibits_both->{interfaces} ;

   my $alltogether = homolobind::pilig::combine_piligpepbits({
      pb => $pb,
      liginfo => $liginfo,
      ligbits => $ligbits,
      pibits => $pibits,
      pepbits => $pepnucibits->{pep},
   });
   my $class2bits = $alltogether->{class2bits} ;
   my $class2ligs = $alltogether->{class2ligs} ;
   my $class2allligbits = $alltogether->{class2allligbits} ;

# Create a family indexed list of binding sites
   my $class2sid12_side ;
   my $class2pint;
   foreach my $classtype (@{$supfam_specs->{thresholds}->{class_levels}}) {
      foreach my $sid12 (keys %{$interfaces}) {
         my $sid ;
         ($sid->{1}, $sid->{2}) = split(/\t/, $sid12) ;
         foreach my $side (1,2) {
            $class2sid12_side->{$pb->{sid2class}->{$classtype}->{$sid->{$side}}}->{$sid12."\t".$side}++ ; }
      }

      foreach my $sid (keys %{$pepnucibits->{pep}->{$classtype}}) {
         $class2pint->{$pb->{sid2class}->{$classtype}->{$sid}}->{$sid}++ ; }
   }

# Read in domain information from scop_cla: fam2sf, px2scopid, scopid2class
   my $scopinfo = readin_scop_cla({
      fn => $homolobind_specs->{scop}->{cla_fn}}) ;

# Make a list of families with at least one binding site
   my $fam_with_bs = {} ;
   map { $fam_with_bs->{$_}++; } (keys %{$class2allligbits->{fam}}) ; # L
   map { $fam_with_bs->{$_}++; } (keys %{$class2sid12_side}) ;        # P
   map { $fam_with_bs->{$_}++; } (keys %{$class2pint}) ;              # pep

# remove any families outside of a-g
   {
      my $delete_fam = {} ;
      foreach my $fam (keys %{$fam_with_bs}) {
         if ($fam !~ /[a-g]/) {$delete_fam->{$fam}++;} }
      map { delete $fam_with_bs->{$_};} (keys %{$delete_fam}) ;
   }


# If specified, load list of families to process
   my $fam_list = {};
   if (exists $in->{fam_list_fn}) {
      open(FAMLISTF, $in->{fam_list_fn}) ;
      while (my $line = <FAMLISTF>) {
         chomp $line;
         $fam_list->{$line}++ ;
      }
      close(FAMLISTF) ;
   } else {
      map {$fam_list->{$_}++;} (keys %{$fam_with_bs}) ;
   }


   if (exists $in->{cluster_fl} && $in->{cluster_fl} == 1 ) { # master node:

      print "* benchmark_homolobind() ".localtime()
         if(!exists $in->{quiet_fl});
      my ($temp_fh, $temp_fn) ;
      ($temp_fh->{bm_homolobind_in},
       $temp_fn->{bm_homolobind_in}) =
       tempfile("splits_bm_homolobind_input.XXXXX") ;

      ($temp_fh->{bm_homolobind_out},
       $temp_fn->{bm_homolobind_out}) =
       tempfile("splits_bm_homolobind_SGEout.XXXXX") ;
      ($temp_fh->{bm_homolobind_err},
       $temp_fn->{bm_homolobind_err}) =
       tempfile("splits_bm_homolobind_SGEerr.XXXXX") ;

# Set up split files
      my $split_dir = tempdir("splits_bm_homolobind.XXXXX") ;
      my $numjobs ;
      if (exists $in->{numjobs}) {
         $numjobs = $in->{numjobs} ;
      } else {
         $numjobs = $homolobind_specs->{SGE}->{numjobs} ;
      }

      if (exists $in->{fam_list_fn}) {
         system("cp ".$in->{fam_list_fn}." ".$temp_fn->{bm_homolobind_in}) ;
      } else {
         open(TEMPFH, ">".$temp_fn->{bm_homolobind_in}) ;
         foreach my $fam (sort keys %{$fam_list}) {
            print TEMPFH $fam."\n" ;
         }
         close(TEMPFH) ;
      }

      my $splits = homolobind::SGE::_clust_split_ins({
         fn => $temp_fn->{bm_homolobind_in},
         dir => $split_dir,
         numjobs => $numjobs,
      }) ;

      my ($perlscript_fh, $perlscript_fn) =
         tempfile("pb.bm_homolobind.XXXXX", SUFFIX =>'.ami.pl') ;


      print {$perlscript_fh} '#!/usr/local/bin/perl'."
use strict;
use warnings;
use pibase;
use homolobind ;
use Bit::Vector ;

main() ;

sub main {

   homolobind::benchmark_homolobind({
      cluster_fl => 0,
      ARGV => \\\@ARGV,
   }) ;

}
" ;

      my ($sgescript_fh, $sgescript_fn)  =
         tempfile("pb.bm_homolobind.XXXXX", SUFFIX => ".SGE.sh") ;
      my $sge_outdir = tempdir("SGEOUT.bm_homolobind.XXXXX") ;

      print {$sgescript_fh} "#!/bin/csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -o $sge_outdir
#\$ -e $sge_outdir
#\$ -r y\n" ;

      if (exists $homolobind_specs->{SGE}->{nodespecs}) {
         print {$sgescript_fh} $homolobind_specs->{SGE}->{nodespecs}."\n" ;}

      print {$sgescript_fh} "#\$ -t 1-$splits->{numjobs}

set tasks1=( $splits->{tasklist} )
set input1 = \$tasks1[\$SGE_TASK_ID\]

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
set scratchdir=/tmp/fred/\$input1.\$\$

rm -rf \$scratchdir
mkdir -p \$scratchdir

cp $perlscript_fn \$scratchdir
cp $split_dir/\$input1 \$scratchdir

cd \$scratchdir

echo \"#sgejob run started on \$curhost at \$curtime\"
perl $perlscript_fn -fam_list_fn \$input1 

set curtime=`date`
echo \"#sgejob run finished on \$curhost at \$curtime\"

rm -f \$scratchdir/\$input1 \$scratchdir/$perlscript_fn
cd \$curdir
rmdir \$scratchdir\n" ;
      close($sgescript_fh) ;

      print STDERR "   submitted $sgescript_fn ".localtime() if
         (!exists $in->{quiet_fl});
      my $qsub_job_id = homolobind::SGE::_clust_qsub({
         hostname => $homolobind_specs->{SGE}->{headnode},
         sgescript_fn => $sgescript_fn,
      }) ;
      print STDERR "      job $qsub_job_id\n" ;

      while (1) {
         sleep $homolobind_specs->{SGE}->{qstat_sleep} ;
         my $job_status= homolobind::SGE::_clust_qstat({
            hostname => $homolobind_specs->{SGE}->{headnode},
            job_id => $qsub_job_id});
         if ($job_status) {last;}
      }

      homolobind::SGE::_clust_merge_outs({
         script_fn => $sgescript_fn,
         out_fn => $temp_fn->{bm_homolobind_out},
         err_fn => $temp_fn->{bm_homolobind_err},
         job_id => $qsub_job_id,
         outdir => $sge_outdir,
         numjobs => $splits->{numjobs}
      }) ;

      if (exists $in->{out_fn}) {
         system("cp ".$temp_fn->{bm_homolobind_out}." ".
                $in->{out_fn}) ;
         unlink $temp_fn->{bm_homolobind_out} ;
      }

      if (exists $in->{err_fn}) {
         system("mv ".$temp_fn->{bm_homolobind_err}." ".
                $in->{err_fn}) ; }

   } else { #processing node:

# Set up Mersene-Twister random number generator for bayesian bootstrapping
      my $mtgen_seed = time ;
      my $mtgen = Math::Random::MT->new($mtgen_seed) ;

# Setup output file
      my $out_fh ;
      if (exists $in->{out_fn}) {
         open($out_fh, ">".$in->{out_fn}) ;
      } else {
         open($out_fh, ">-") ;
      }
   
# Iterate over each family with at least one binding site
      my $fam_coverage ;
      foreach my $class (sort keys %{$fam_list}) {
         if (!exists $fam_with_bs->{$class}) {next;}
   
         print STDERR "NOW ON $class\n" ;
         my $classtype = 'fam' ;
   
# Load whole family ASTRAL alignment
         my $fam_aln = homolobind::ASTRAL::load_asteroids_aln({
            aln_fn => $homolobind_specs->{asteroids}->{$classtype.'_aln'}.
                        '/'.$class.'.fasta_aln' ,
            seq_fn => $homolobind_specs->{asteroids}->{$classtype.'_seq'}.
                        '/'.$class.'.fa' ,
            seqclcont100 => $astral->{seqcl2cont}->{100},
            seqcl100 => $astral->{seqcl}->{100},
            allchains => $pb->{pdbchains}
         }) ;
   
# Make a generic structure to unify code for benchmarking L/P/p site transfer
# {L|P|p}->{ligsig|sid12_side|sid_chain}->{source_ind} = sid entry in cur_aln
#                                       ->{bits} = bit vector
# (eg only load ligbs that match MW thresholds
         my $curfam_bs = {} ;
         {
         if (exists $class2allligbits->{fam}->{$class}) {
            my $curligs = $class2allligbits->{fam}->{$class};
            foreach my $j ( 0 .. $#{$curligs}) {
               my (undef, $ligcod, undef) = split(/\t/, $curligs->[$j]->[0]) ;
               if (exists $supfam_specs->{thresholds}->{min_bs_size} &&
                   $curligs->[$j]->[1]->Norm() < 
                   $supfam_specs->{thresholds}->{min_bs_size}) {
                  next; }
   
               if (exists $supfam_specs->{thresholds}->{L}->{min_mw} &&
                   $liginfo->{mw}->{$ligcod} <
                   $supfam_specs->{thresholds}->{L}->{min_mw}) {next;}
   
               if (exists $supfam_specs->{thresholds}->{L}->{max_mw} &&
                   $liginfo->{mw}->{$ligcod} >
                   $supfam_specs->{thresholds}->{L}->{max_mw}) {next;}
               
               my $sid = $curligs->[$j]->[2] ;
               my ($osid) = ($sid =~ /SCOP\.(.+)/) ;
               push @{$curfam_bs->{L}}, {
                  name => $curligs->[$j]->[0],  #identifier for the BS
                  bits => $curligs->[$j]->[1],
                  sid => $sid,
                  osid => $osid,
               } ;
            }
         }
   
         if (exists $class2sid12_side->{$class}) {
            foreach my $sid12_side (keys %{$class2sid12_side->{$class}}) {
               my ($sid1, $sid2, $side) = split(/\t/, $sid12_side) ;
               my $sid12 = $sid1."\t".$sid2 ;
               my $sid = $sid1 ;
               if ($side == 2) {$sid = $sid2 ;}
               my ($osid) = ($sid =~ /SCOP\.(.+)/) ;
               my $btype = 'P';
               if ($interfaces->{$sid12}->{chains} eq 'diff') {
                  $btype .= 'inter' ;
               } else {
                  $btype .= 'intra' ;
               }
   
#100601_1514 - if cluster_rep only have been read in, is possible that this
#  interface is a clusterrep for superfamily level, but not at family level
#  so that the fam bits are undefined; must skip these.
# eg BDP23296-0_SCOP.d1jjua4-BDP23296-0_SCOP.d1jjuc_ is sf, but not fam defined
               if (!defined $interfaces->{$sid12}->{$side}->{pibits}->{fam} ||
                   (exists $supfam_specs->{thresholds}->{min_bs_size} &&
                    $interfaces->{$sid12}->{$side}->{pibits}->{fam}->Norm() <
                    $supfam_specs->{thresholds}->{min_bs_size})) {
                  next; }
   
               push @{$curfam_bs->{$btype}}, {
                  name => $sid12_side,
                  bits => $interfaces->{$sid12}->{$side}->{pibits}->{fam},
                  sid => $sid,
                  osid => $osid,
               } ;
            }
         }
   
         if (exists $class2pint->{$class}) {
            foreach my $sid (keys %{$class2pint->{$class}}) {
               my ($osid) = ($sid =~ /SCOP\.(.+)/) ;
               foreach my $target_ch (keys %{$pepnucibits->{pep}->{fam}->{$sid}}) {
                  if ($target_ch eq 'cumulative') {next;}
                  if (exists $supfam_specs->{thresholds}->{min_bs_size} &&
                      $pepnucibits->{pep}->{fam}->{$sid}->{$target_ch}->Norm() <
                      $supfam_specs->{thresholds}->{min_bs_size}) {
                     next; }
                  push @{$curfam_bs->{p}}, {
                     name => $sid."\t".$target_ch,
                     bits => $pepnucibits->{pep}->{fam}->{$sid}->{$target_ch},
                     sid => $sid,
                     osid => $osid,
                  } ;
               }
            }
         }
         }
   
   
# Iterate over binding site types: L, P, p
         foreach my $btype (keys %{$curfam_bs}) {
# If this family only has 1 binding site, can't estimate coverage or acc
            if ($#{$curfam_bs->{$btype}} == 0) {
               print STDERR "SKIPPING $class $btype: only has 1 binding site\n";
               next;
            }
   
# Iterate over all binding sites in this family
# - homology transfer from other binding sites in this family w diff PDB
# 1. covered bs residue = one aligned to a bs in a different PDB entry
   
            my $bs_coverage = [] ;
            my $fam_roc_scores = {}; #family-wide ROC curve
   
            my $uniq_pdb = {};
            print STDERR " total: ".($#{$curfam_bs->{$btype}} + 1).
                         " binding sites\n" ;
   
# Generate 10,000 negative sequences - shuffle each seq x (10K/numfammembers)
            my $num_neg = $homolobind_specs->{benchmark}->{num_negatives} ;
            my $shuffled_aln_seq = []; #shuffle alignment string when needed
            {
               my $numfam_members = keys %{$fam_aln->{aln}} ;
               my $num_shuffle_each = POSIX::ceil($num_neg / $numfam_members) ;
               $num_neg = $num_shuffle_each * $numfam_members ;
   
               foreach my $osid (keys %{$fam_aln->{aln}}) {
                  my $cur_shuffle = shuffle_aln_string_keepgapstrx({
                     num_shuffle => $num_shuffle_each,
                     orig => $fam_aln->{aln}->{$osid},
                  }) ;
                  push @{$shuffled_aln_seq}, @{$cur_shuffle} ;
               }
            }
   
# Estimate family-wide coverage and binding site-specific accuracy
# Iterate over all binding sites in this class;
   
# Two-pass scheme:
# 1. BS as template: establish (i) seqid cutoffs and (ii) FPR/CI at FPR cutoffs
# 2. BS as target: establish family-wide TPR (and BB CI)
   
# Pass1. Iterate over each ligand binding site as template
            my $skip_templates ;   # don't use a template BS if any of the
                                   # strictest goal FPR is not achievable.
            my $cur_cutoffs ;      # sequence identity thresholds for goal FPR
            foreach my $j (0 .. $#{$curfam_bs->{$btype}}) {
               my $bs = $curfam_bs->{$btype}->[$j] ;
   
               my @bs_alnpos = $bs->{bits}->Index_List_Read() ;
               my $transf_bs = {} ;
# Transfer onto all shuffled sequences:
               foreach my $k ( 0 .. $#{$shuffled_aln_seq}) {
                  $transf_bs->{$k}->{residues} = [] ;
                  $transf_bs->{$k}->{num_ident} = 0 ;
                  $transf_bs->{$k}->{num_aln} = 0 ;

                  foreach my $bs_alnpos (@bs_alnpos) {
                     my $cur_tmpl_char = substr(
                        $fam_aln->{aln}->{$bs->{osid}},$bs_alnpos,1);
                     my $cur_targ_char = substr( $shuffled_aln_seq->[$k],
                                                 $bs_alnpos,1 );
                     if ($cur_targ_char eq '-') {next;}
                     push @{$transf_bs->{$k}->{residues}}, $bs_alnpos ;
                     $transf_bs->{$k}->{num_aln}++ ;
   
                     if ($cur_tmpl_char eq $cur_targ_char) {
                        $transf_bs->{$k}->{num_ident}++ ; }
                  }
                  $transf_bs->{$k}->{bs_seqid} = $transf_bs->{$k}->{num_ident} /
                                                 ($#bs_alnpos + 1) ;
# fpd100824_1428 - changed to correctly compute frac aln
#                  if ($#{$transf_bs->{$k}->{residues}} >= 
#                      $supfam_specs->{thresholds}->{frac_bs_aligned}) {
#                     $transf_bs->{$k}->{tmpl_bs_size} = $#bs_alnpos + 1 ;
#                  } else {
#                     delete $transf_bs->{$k} ;
#                  }

# IMPOSE MINIMUM NUMBER OF aligned residues to template bs
                  if ($transf_bs->{$k}->{num_aln} >= 
                      ($supfam_specs->{thresholds}->{frac_bs_aligned} * 
                       ($#bs_alnpos + 1)) ) {
                     $transf_bs->{$k}->{tmpl_bs_size} = $#bs_alnpos + 1 ;
                  } else {
                     delete $transf_bs->{$k} ;
                  }
               }
   
# Sort negative scores to determine seqid cutoff for 1/2/5 FPR
# False positive = binding site? or 
   
# Build FP-curves to determine seqid cutoffs for 1/2/5a FPR
# 1. get seqid range that gives each FPR; orig/mean/CI
# 2. get FPR mean/CI at the identified seqid cutoffs.
   
               my @neg_scores ;
               map {push @neg_scores, $transf_bs->{$_}->{bs_seqid}}
                  (keys %{$transf_bs}) ;
   
# Estimate seqid cutoffs for goal FPRs
               my $bs_fpc = build_fpcurve({
                  goal_fpr      => $homolobind_specs->{fpr_cutoffs},
                  scores        => \@neg_scores,
                  sort_order    => 'desc'
               }) ;
   
# Store seqid thresholds at each fpr for this BS
               $cur_cutoffs->[$j] = $bs_fpc->{cutoffs} ;
   
                {
                  my @outvals ;
                  push @outvals, ("BS", $btype, $bs->{sid}, $class, $bs->{name},
              $bs->{bits}->Norm(),
              ($#neg_scores + 1)) ;
                  foreach my $fpr (@{$homolobind_specs->{fpr_cutoffs}}) {
                     if (!defined $bs_fpc->{cutoffs}->{"fpr$fpr"}) {
                        push @outvals, "UNDEF", "UNDEF" ;
                        $skip_templates->{$j}++ ;
                     } else {
                        push @outvals, 
              (sprintf("%.3f", $bs_fpc->{cutoffs}->{"fpr$fpr"}->{score}),
               sprintf("%.3f", $bs_fpc->{cutoffs}->{"fpr$fpr"}->{fpr}));
                     }
                  }
                  print {$out_fh} join("\t", @outvals)."\n" ;
               }
            }
   
# Pass2. Iterate over each ligand binding site as target to cross-annotate
#        from all other binding sites in the family using just-established
#        thresholds
   
   
# goal: hash list of {target seq}->{BS residue} = [(tmpl, seqid)]
   
# make a family-wide array of "true" scores across all family
# post-hoc build a TP analog of the fpcurve (call same routine actually)
   
# keep a running tally of covered positives that pass or don't pass
# seqid cutoffs establisehd above: $cur_cutoffs->[$k]->{fprFPR}->{score}
   
            my $truepos = {}; my $num_covered = 0;
            foreach my $fpr (@{$homolobind_specs->{fpr_cutoffs}}) {
               $truepos->{$fpr} = 0;}
   
            foreach my $j (0 .. $#{$curfam_bs->{$btype}}) {
               my $bs1 = $curfam_bs->{$btype}->[$j] ;
   
# Keep track of union of all other BS - for coverage calc
               my $allother_union = $bs1->{bits}->Shadow() ;
   
               my $transf_bs ; #index is over templates - all except self
   
# Iterate over all other lig bs in this class
               foreach my $k (0 .. $#{$curfam_bs->{$btype}}) {
                  if (exists $skip_templates->{$k}) { next; }
   
                  if ($j == $k) {next;} #ignore if same
                  my $bs2 = $curfam_bs->{$btype}->[$k] ;
                  if (!exists $fam_aln->{aln}->{$bs2->{osid}}) {
                     print STDERR "WARNING: COULDNT FIND ALN STRING for ".
                                  $bs2->{osid}." ($btype binding site) \n" ;
                     next;
                  }
   
# Skip this comparison if they are from the same PDB file.
                  $uniq_pdb->{$pb->{sid2pdb}->{$bs1->{sid}}}++ ;
                  $uniq_pdb->{$pb->{sid2pdb}->{$bs2->{sid}}}++ ;
                  if ($pb->{sid2pdb}->{$bs1->{sid}} eq
                      $pb->{sid2pdb}->{$bs2->{sid}}) { next;}
   
                  $allother_union->Union($allother_union,$bs2->{bits}) ;
   
# Set transferred BS specs
                  my $bs2_sig = $bs2->{name} ;
   
# Homology transfer from bs2 onto bs1
# Calculate seqid over template bs positions - for family-wide ROC curves
# * Get TPR at the seqid cutoffs identified above
                  $transf_bs->{$k}->{source_ind} = $k;
                  $transf_bs->{$k}->{residues} = [] ;
                  $transf_bs->{$k}->{num_ident}= 0;
                  $transf_bs->{$k}->{num_aln}= 0;
                  my @bs2_alnpos = $bs2->{bits}->Index_List_Read() ;
                  foreach my $bs2_alnpos (@bs2_alnpos) {
   
                     my $cur_tmpl_char = substr($fam_aln->{aln}->{$bs2->{osid}},
                                                $bs2_alnpos,1);
   
                     my $cur_targ_char = substr($fam_aln->{aln}->{$bs1->{osid}},
                                                $bs2_alnpos,1);
   
                     if ($cur_targ_char eq '-') {next;}
   
                     push @{$transf_bs->{$k}->{residues}}, $bs2_alnpos ;
                     $transf_bs->{$k}->{num_aln}++ ;
   
#check if identical
                     if ($cur_tmpl_char eq $cur_targ_char) {
                       $transf_bs->{$k}->{num_ident}++;}
                  }
   
                  $transf_bs->{$k}->{bs_seqid} = $transf_bs->{$k}->{num_ident} /
                                                 ($#bs2_alnpos + 1) ;
# fpd100824_1430 - changed to correctly compute frac aln
#                  if ($#{$transf_bs->{$k}->{residues}} >= 
#                      $supfam_specs->{thresholds}->{frac_bs_aligned}) {
#                     $transf_bs->{$k}->{tmpl_bs_size} = $#bs2_alnpos + 1;
#                  } else {
#                     delete $transf_bs->{$k} ;
#                  }

# IMPOSE MINIMUM NUMBER OF aligned residues to template bs
                  if ($transf_bs->{$k}->{num_aln} >= 
                      ($supfam_specs->{thresholds}->{frac_bs_aligned} * 
                       ($#bs2_alnpos + 1)) ) {
                     $transf_bs->{$k}->{tmpl_bs_size} = $#bs2_alnpos + 1 ;
                  } else {
                     delete $transf_bs->{$k} ;
                  }

               }
   
#compute coverage for this targbs: how well is bs1 covered by other templates?
               my $covered_bits = $bs1->{bits}->Clone() ;
               $covered_bits->Intersection($bs1->{bits}, $allother_union);
               $bs_coverage->[$j] = $covered_bits->Norm() /
                                       $bs1->{bits}->Norm() ;
# If coverage is 0, no point in addressing accuracy (no P examples)
               if ($bs_coverage->[$j] == 0) {next;}
               $num_covered += $covered_bits->Norm() ;
   
# Determine accuracy for bs1 as target: (fam-wide ROC)
# * seqid over bs2 pos (P = covered bs1 res)
               foreach my $fpr (@{$homolobind_specs->{fpr_cutoffs}}) {
               my $cur_tp_res = {} ;
               foreach my $k (keys %{$transf_bs}) {
                  my $cur_bsseqid = $transf_bs->{$k}->{bs_seqid};
   
# If seqid cutoff not met, skip to next template
                  if ($cur_bsseqid < $cur_cutoffs->[$k]->{"fpr$fpr"}->{score}) {
                     next;}
   
                  foreach my $res ( @{$transf_bs->{$k}->{residues}}) {
   
# Only consider BS residues transfered onto actual target BS
                     if ($bs1->{bits}->bit_test($res) == 0) {next;}
                     $cur_tp_res->{$res}++ ;
                  }
               }
               $truepos->{$fpr} += keys %{$cur_tp_res} ;
               }
            }
   
# average over coverage values for each binding site to get family-wide value;
            {
               $fam_coverage->{$class}->{$btype} = 0 ;
               my $n ;
               foreach my $j ( 0 .. $#{$curfam_bs->{$btype}}) {
                  $fam_coverage->{$class}->{$btype} += $bs_coverage->[$j] ; }
               $fam_coverage->{$class}->{$btype} /= ($#{$curfam_bs->{$btype}} + 1);
            }
   
            my $num_uniq_pdb = keys %{$uniq_pdb} ;
            if ($num_uniq_pdb == 1) {
               print STDERR "SKIPPING: $class only has $btype binding sites ".
                            "from only 1 PDB code\n" ;
               next;}
            my @outvals = ("FAM", $btype, $class,
                           ($#{$curfam_bs->{$btype}} + 1), $num_uniq_pdb,
                           sprintf("%.3f",$fam_coverage->{$class}->{$btype})) ;
   
            if ($fam_coverage->{$class}->{$btype} > 0) {
   
# After all "TP" BS are collected, figure out what fraction are let through
# Bayesian bootstrap the TPR also, by pulling from a dirichlet.
   
#Simple calc:
# 1. count number of TP (= P - FN)
   
# 2. generate numtrials x dirichlet distributions of size P
            my $bb_tprs = {} ;
            foreach my $trial (1 .. $homolobind_specs->{bb_replicates}) {
               my $weights = generaterand_dirichlet({
                  mtgen => $mtgen, n => $num_covered }) ;
               foreach my $fpr ( @{$homolobind_specs->{fpr_cutoffs}}) {
   
# 3. calc TPR estimates by adding the first TP cells in the dirichlet distr.
                  my $curbb_tpr = 0 ;
                  foreach my $j ( 0 .. ($truepos->{$fpr} - 1)) {
                     $curbb_tpr += $weights->[$j] ; }
                  push @{$bb_tprs->{$fpr}}, $curbb_tpr ;
               }
            }
   
# 4. calc mean/95\% CI for TPR.
            my $bb_tpr_ci ;
            foreach my $fpr ( @{$homolobind_specs->{fpr_cutoffs}}) {
               $bb_tpr_ci->{$fpr} = calc_confidence_interval({
                  data => $bb_tprs->{$fpr}});
               push @outvals, sprintf("%.3f", $bb_tpr_ci->{$fpr}->{mean});
               push @outvals, sprintf("%.3f", $bb_tpr_ci->{$fpr}->{lower});
               push @outvals, sprintf("%.3f", $bb_tpr_ci->{$fpr}->{upper});
            }
   
# in the end: seqid for 1/2/5 FPR for each ligbs; fam-wide TPR/CI at these seqid
            }
            print {$out_fh} join("\t", @outvals)."\n" ;
         }
      }
   }

}


=head2 shuffle_aln_string_keepgapstrx()

   Title:       shuffle_aln_string_keepgapstrx()
   Function:    Shuffle aln string sequence, while maintaining gap structure

   Args:        ->{orig} = original alignmen string
                ->{num_shuffle} = number of shuffles to produce

   Returns:     ->[i] = ith shuffled string

=cut

sub shuffle_aln_string_keepgapstrx {

   my $in = shift ;
   my $orig_string = $in->{orig} ;

# Remove gap (-) or fragment (X) character
   my $seqonly = $in->{orig} ; $seqonly =~ s/[\-X]+//g ;
   my @seqchar = split('', $seqonly) ;

# Keep track of which positions had sequence characters
   my @seqpos ;
   foreach my $j ( 0 .. length($in->{orig})) {
      if (substr($in->{orig}, $j, 1) =~ /[\-X]/) {next;}
      push @seqpos, $j ;
   }

   my $shuffles = [];
   foreach my $j ( 0 .. $in->{num_shuffle}) {
      fy_shuffle(\@seqchar) ;
      $shuffles->[$j] = $orig_string ;
      map {substr($shuffles->[$j], $seqpos[$_], 1) = $seqchar[$_]; }
         (0 .. $#seqchar) ;
#      print STDERR "ORIG: ".$in->{orig}."\n" ;
#      print STDERR " NEW: ".$shuffles->[$j]."\n";
   }

   return $shuffles ;

}


=head2 build_fpcurve()

   Title:       build_fpcurve()
   Function:    Given an array of scores for "negatives", build an FP-vs-score
                  curve - optionally, using weights if provided
   Args:        ->{scores} = [] array of scores for negatives
                ->{goal_fpr} = [] array of goal FPRs (eg, 0.01)
                ->{sort_order} = 'asc' or 'desc'; desc = highest score is best.

   Returns:     ->{fpcpoints}->[i] = {fpr => false positive rate,
                                      last_score => last score}
                ->{cutoffs}->{fprNN} = {fpr => exact false positive rate,
                                        score => score cutoff to yield fpr}
                ->{optimal}->{dist|tpr}

   Displays:    ROC plot properties (AUC, optimal threshold, TP, FP)

=cut

sub build_fpcurve {

   my $in = shift;
   my $scores = $in->{scores} ;

   my @sort_order ;
   if ($in->{sort_order} eq 'desc') {
      @sort_order = sort {$scores->[$b] <=> $scores->[$a]} (0 .. $#{$scores});
   } else {
      @sort_order = sort {$scores->[$a] <=> $scores->[$b]} (0 .. $#{$scores});
   }
   my $n = $#{$scores} + 1 ;

   my $num_neg = $#{$scores} ;
   my $last_score ; my $fp = 0; my $last_fp = 0; my $cumn = 0 ;
   my $fpc ;
   foreach my $ind (@sort_order) {
      if (!defined $last_score || $scores->[$ind] != $last_score) {
         if (!defined $last_score) {$last_score = 'undef';}

         my $curfpr = $fp / $num_neg ;

         push @{$fpc->{fpcpoints}}, {
            fpr => sprintf("%.3f", $curfpr),
            last_score => $last_score
         } ;

         $last_score = $scores->[$ind] ;
         $last_fp = $fp ;
      }

      my $weight = 1 ;
      if (exists $in->{weights}) {
         $weight = $n * $in->{weights}->[$ind] ; }
      $fp += $weight ;
   }
   my $curfpr = $fp / $num_neg ;

   push @{$fpc->{fpcpoints}}, {
      fpr => sprintf("%.3f", $curfpr),
      last_score => $last_score
   } ;

   if (exists $in->{goal_fpr}) {
      $fpc->{cutoffs} = {} ;
      my $findthese ;
      map {$findthese->{$_} = 0;} (@{$in->{goal_fpr}}) ;
      foreach my $fpc_point (@{$fpc->{fpcpoints}}) {
         if ($fpc_point->{last_score} eq 'undef') {next;}
         foreach my $fpr (keys %{$findthese}) {
            if ($findthese->{$fpr} == 1) {next;}
            if ($fpc_point->{fpr} <= $fpr) {
                  $fpc->{cutoffs}->{"fpr$fpr"} = {
                     fpr => $fpc_point->{fpr},
                     score => $fpc_point->{last_score},
                  } ;
            } else {
                  $findthese->{$fpr}++ ;
            }
         }
      }
   }

# if fpr_goal specified, pass through fpc to determine seqid c/o and
#   exact FPR closest to the goal FPR.

   return $fpc ;
}


=head2 build_roc_bayesboot()

   Title:       build_roc_bayesboot()
   Function:    Builds ROC curves with confidence bands using bayesian bootstrap
   Args:        ->{scores}->{p} = array of positive example scores
                ->{scores}->{n} = array of negative examplescores
                ->{sort_order} = 'asc' or 'desc'; desc = highest score is best.
                ->{bb_replicates} = number of bayesian bootstrap replicates

   Returns:     ->{rocpoints}
                ->{auc}->{orig|bb_mean|bb_95_lower|bb_95_upper}
                ->{optimal}->{dist|tpr}

   Displays:    ROC plot properties (AUC, optimal threshold, TP, FP)

=cut

sub build_roc_bayesboot {

   require Math::Random::MT ;

   my $in = shift ;
   my $mtgen ;
   if (!exists $in->{mtgen}) {
      my $mtgen_seed = time ;
      $mtgen = Math::Random::MT->new($mtgen_seed) ;
   } else {
      $mtgen = $in->{mtgen} ;
   }

# Given original scores list, as in for regular build_roc(),
# 1. calculate original ROC score with necessary statistics,
# 2. generate BB replicates of the data points and calc new ROC and get stats
# 3. calc 95\% CI for each particular stat: AUC, TPR/FPR/score thresh @1/2/5/TPR

   my $num_pos = $#{$in->{scores}->{p}} + 1 ;
   my $num_neg = $#{$in->{scores}->{n}} + 1 ;

   my $results ;

# bayesian bootstrap $in->{bb_replicates} times and build ROC curves
   my $replicate_stats; #numbers to generate confidence intervals for
   foreach my $j ( 0 .. $in->{bb_replicates}) {
      my $roc ;
      if ($j == 0) {
         $roc = build_roc({ scores     => $in->{scores},
                            sort_order => $in->{sort_order} }) ;
      } else {
# Sample probability vector from the dirichlet distribution
         my $weights ;
         if ($num_pos == 1) {
            $weights->{p} = [1] ;
         } else {
            $weights->{p} = generaterand_dirichlet({
               mtgen => $mtgen, n => $num_pos }) ;
         }

#         print STDERR "asked for (p) $num_pos-dimension dirichlet draw, but got back: ".($#{$weights->{p}} + 1)."\n" ;

         if ($num_neg == 1) {
            $weights->{n} = [1] ;
         } else {
            $weights->{n} = generaterand_dirichlet({
               mtgen => $mtgen, n => $num_neg }) ;
         }

#         print STDERR "asked for (n) $num_neg-dimension dirichlet draw, but got back: ".($#{$weights->{n}} + 1)."\n" ;

# Use the dirichlet-sampled probabilities as weights in the ROC plot
         $roc = build_roc({ scores     => $in->{scores},
                            weights    => $weights,
                            sort_order => $in->{sort_order} }) ;
      }

# Store statistics for 95\% CI bounds: AUC, TPR/seqid/FPR at 1/2/5TPR
      my $curroc_stats ;
      if (exists $in->{goal_fpr}) {
         my $fpr2tpr = {} ;
         my $findthese ;
         map {$findthese->{$_} = 0;} (@{$in->{goal_fpr}}) ;
         my $pn = 0 ;
         foreach my $roc_point (@{$roc->{rocpoints}}) {
            if ($roc_point->{last_score} eq 'undef') {next;}
            foreach my $fpr (keys %{$findthese}) {
               if ($findthese->{$fpr} == 1) {next;}
               if ($roc_point->{x} <= $fpr) {
                  $fpr2tpr->{$fpr} = {
                     tpr => $roc_point->{y},
                     fpr => $roc_point->{x},
                     score => $roc_point->{last_score},
                  } ;
               } else {
                  $findthese->{$fpr}++ ;
               }
            }
            $pn++ ;
         }
         foreach my $fpr (keys %{$fpr2tpr}) {
            $curroc_stats->{"tpr_at_fpr$fpr"} = $fpr2tpr->{$fpr}->{tpr} ;
            $curroc_stats->{"score_at_fpr$fpr"} = $fpr2tpr->{$fpr}->{score} ;
         }
      }
      $curroc_stats->{auc} = $roc->{auc} ;

      if ($j == 0) {
         $results->{orig} = $curroc_stats ;
      } else {
         map {push @{$replicate_stats->{$_}}, $curroc_stats->{$_};}
            (keys %{$curroc_stats}) ;
      }
   }

   foreach my $statistic (keys %{$replicate_stats}) {
      my $ci = calc_confidence_interval({data => $replicate_stats->{$statistic},
                                         mass => 0.95 }); 
      $results->{bootstrap}->{$statistic} = $ci ;
   }

   return $results ;

}


=head2 fy_shuffle()

   Title:       fy_shuffle()
   Function:    Shuffles an array in place
   Args:        $_ = arrayref
   Returns:     Nothing

=cut

sub fy_shuffle {
   my $deck = shift;  # $deck is a reference to an array
   my $i = @$deck;
   while ($i--) {
      my $j = int rand ($i+1);
      @$deck[$i,$j] = @$deck[$j,$i];
   }
}



=head2 generaterand_dirichlet()

   Title:       generaterand_dirichlet()
   Function:    returns a vector sampled from a flat dirichlet distribution.
                per DB Rubin, 1981. The Bayesian Bootstrap.

   Args:        ->{n} = number of objects
                ->{mtgen} = Math::Random::MT mersene twist generator
   Returns:     $_ = arrayref of sampled objects

=cut

sub generaterand_dirichlet {

   require Math::Random::MT ;

   my $in = shift ;
   my $n = $in->{n} ;
   my $mtgen ;
   if (!exists $in->{mtgen}) {
      my $mtgen_seed = time ;
      $mtgen = Math::Random::MT->new($mtgen_seed) ;
   } else {
      $mtgen = $in->{mtgen} ;
   }

# pick n-1 numbers randomly, sort, calc gaps, return array of gaps
   my @rands ;
   foreach my $j ( 1 .. $n - 1) {
      push @rands, $mtgen->rand() ; }
   my @sort_rands = sort {$a <=> $b} @rands ;

   my @p ;
   push @p, $sort_rands[0] ;
   foreach my $j ( 1 .. $#sort_rands) {
      push @p, ($sort_rands[$j] - $sort_rands[($j -1)]) ; }
   push @p, (1 - $sort_rands[$#sort_rands]) ;

   return \@p ;

}


=head2 generaterand_discrete_pdf()

   Title:       generaterand_discrete_pdf()
   Function:    Randomly draw from a discrete sample
                from fpdgenlib::cisreg.pm; ideas from Perl cookbook 2.10

   Args:        ->{p} = either hash of probabilities {p}->{OBJECT} = PROBABILITY
                           OR array of probabilities {p}->[i] = PROBABILITY
                ->{mtgen} = Math::Random::MT mersene twist generator
   Returns:     $_ = object index (if p is array) or key (if p is hash)

=cut

sub generaterand_discrete_pdf {
# assumes real probabilities, ie total == 1
   require Math::Random::MT ;

   my $in = shift ;

   my $mtgen ;
   if (!exists $in->{mtgen}) {
      my $mtgen_seed = time ;
      $mtgen = Math::Random::MT->new($mtgen_seed) ;
   } else {
      $mtgen = $in->{mtgen} ;
   }

   my $rand = $mtgen->rand() ;

   if (ref($in->{p}) eq 'HASH') {
      foreach my $key ( keys %{$in->{p}}) {
         return $key if (($rand -= $in->{p}->{$key}) < 0) ; }
   } elsif (ref($in->{p}) eq 'ARRAY') {
      foreach my $j ( 0 .. $#{$in->{p}}) {
         return $j if (($rand -= $in->{p}->[$j]) < 0) ; }
   }

}


=head2 calc_confidence_interval()

   Title:       calc_confidence_interval()
   Function:    Given array, return 95\% (or specified) confidence band

   Args:        ->{mass} = probability mass to bound. defaults to 0.95
                ->{mtgen} = Math::Random::MT mersene twist generator

   Returns:     ->{mean} = mean
                ->{lower} = lower confidence interval bound
                ->{upper} = upper confidence interval bound

=cut

sub calc_confidence_interval {

   my $in = shift ;
   my @sorted = sort {$a <=> $b} @{$in->{data}} ;

   my $mass = 0.95 ; #defaults to 95% confidence bands, unless specified
   if (exists $in->{mass}) {
      $mass = $in->{mass} ; }

   my $mean = calc_mean(\@sorted) ;

   my $lower = POSIX::floor(($#sorted + 1) * (1 - $mass) / 2) - 1;
   my $upper = POSIX::ceil(($#sorted + 1) * (1 - ((1 - $mass) / 2))) - 1;
   return {
      mean => $mean,
      lower => $sorted[$lower],
      upper => $sorted[$upper],
   } ;

}


=head2 calc_mean()

   Title:       calc_mean()
   Function:    Calculate mean of values in an arrayref

   Args:        $_ = [] arrayref
   Returns:     $_ = mean

=cut

sub calc_mean {

   my $in = shift ;
   my $total = 0 ;
   my $num = 0 ;

   foreach  my $j ( 0 .. $#{$in}) {
      $total += $in->[$j] ;
      $num++ ;
   }

   my $mean = $total / $num ;
   return $mean ;

}

1 ;
