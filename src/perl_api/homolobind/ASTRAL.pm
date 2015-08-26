=head1 NAME

homolobind::ASTRAL - module to access ASTRAL data

=head1 DESCRIPTION

Contains routines for accessing ASTRAL data for HOMOLOBIND use.
Code is a trimmed down version of pibase::ASTRAL()

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

package homolobind::ASTRAL ;
use strict;
use warnings;
use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw// ;

use homolobind;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;
use POSIX qw/ceil/ ;


=head2 load_asteroids_aln()

   Title:       load_asteroids_aln()
   Function:    Loads an ASTRAL ASTEROIDS alignment file
   Args:        $_->{aln_fn}	name of ASTERIODS alignment file
                $_->{seq_fn}	name of corresponding ASTEROIDS sequence file
                $_->{allchains}	
                $_->{gdseqh}	data structure holding contents of gdseqh file
                $_->{seqclcont100}
                $_->{seqcl100}
                $_->{doms}	optional hash list of domains to load.
   Returns:     parse_aln_raf() alignment structure

=cut

sub load_asteroids_aln {

   my $in = shift ;

   my $seqclcont100 = $in->{seqclcont100} ;
   my $seqcl100 = $in->{seqcl100} ;

   my $doms_enriched ; #supplement with the domain ids of seqcl100 representatives if any of the domains are not themselves representatives
   if (exists $in->{doms}) {
      foreach my $dom (sort keys %{$in->{doms}}) {
         $doms_enriched->{$dom}++ ;
         if (!exists $seqcl100->{$dom}) {next;}
         my $seqcl_no = $seqcl100->{$dom} ;
         my $cluster_rep= $seqclcont100->{$seqcl_no}->[0] ;
         $doms_enriched->{$cluster_rep}++ ; #reprsentative
      }
   }

   my $read_asteroids_aln_options = {
      aln_fn => $in->{aln_fn},
      seq_fn => $in->{seq_fn},
      allchains => $in->{allchains}
   } ;

   if (exists $in->{doms}) {
      $read_asteroids_aln_options->{doms} = $doms_enriched ; }

   my $data = read_asteroids_aln($read_asteroids_aln_options) ;

   my $added ;

   my @origdoms ;
   if (exists $in->{doms}) {
      @origdoms = keys %{$doms_enriched} ;
   } else { #used to always expect $in->{doms} - changed for pilig.pm code
      @origdoms = keys %{$data->{aln}} ;
   }

   foreach my $rep_scopid (@origdoms) {
      if (!exists $seqcl100->{$rep_scopid}) {next;}
      my $clno = $seqcl100->{$rep_scopid} ;
      if ($#{$seqclcont100->{$clno}} > 0) {
         foreach my $j ( 1 .. $#{$seqclcont100->{$clno}}) {
            my $scopid = $seqclcont100->{$clno}->[$j] ;
            if (exists $added->{$scopid}) {next;}

            $data->{seq}->{$scopid} = $data->{seq}->{$rep_scopid} ;
            $data->{aln}->{$scopid} = $data->{aln}->{$rep_scopid} ;
         }
      }
   }

   my $aln ;
   $aln->{alnlength} = $data->{alnlength} ;
   $aln->{aln} = $data->{aln} ;

   return $aln ;
}


=head2 read_asteroids_aln()

   Title:       read_asteroids_aln()
   Function:    Reads an ASTEROIDS alignment file
   Args:        $_->{aln_fn} ASTEROIDS alignment file name
                $_->{seq_fn} corresponding ASTEROIDS sequence file name
                $_->{allchains}->{domain} = pdb chain

   Returns:     ->{seq}->{domain} = 'DOMAINSEQVENCE';
                ->{defstring}->{domain} = definition line from alignment
                ->{class}->{domain} = SCOP class
                ->{aln}->{domain} = domain sequence froma alignment
                ->{pdb}->{domain} = PDB code for the domain
                ->{frags}->{domain} = [{b => startresidue, e => endresidue},...]
                ->{alnlength} = alignment length

=cut

sub read_asteroids_aln {

   my $in = shift;
   my $fn;
   $fn->{aln} = $in->{aln_fn} ;
   $fn->{seq} = $in->{seq_fn} ;
   my $allchains = $in->{allchains} ;


   my $standard_20aa = {
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

   my $data;
   my $cur_dom = '' ;
   my $seq ;
   open(SEQF, $fn->{seq}) ;
   while (my $line = <SEQF>) {
      chomp $line;
      if ($line =~ /^>/) {
         $line =~ s/^>// ;
         my ($t_dom, undef, undef, undef) = split(/\s/, $line) ;
         if (!exists $in->{doms} ||
             (exists $in->{doms} && exists $in->{doms}->{$t_dom})) {
            $data->{seq}->{$t_dom} = ''; }
         $cur_dom = $t_dom ;
      } elsif (!exists $in->{doms} ||
               (exists $in->{doms} && exists $in->{doms}->{$cur_dom})) {
            my $newseq = $line ;
# DONT CHANGE CASE OF x; x = unknown res, X = fragment break.
            $newseq =~ tr/[a-w]/[A-W]/ ;
            $newseq =~ tr/[y-z]/[Y-Z]/ ;
            $data->{seq}->{$cur_dom} .= $newseq ;
#uc to handle 1-seq aln's
      }
   }
   close(SEQF) ;

   $cur_dom = '' ;
   open(ALNF, $fn->{aln}) ;
   while (my $line = <ALNF>) {
      chomp $line;
      if ($line =~ /^>/) {
         $line =~ s/^>// ;
         my ($t_dom, $t_class, $t_def, undef) = split(/\s/, $line) ;
         $cur_dom = $t_dom ;
         if (exists $in->{doms} && !exists $in->{doms}->{$cur_dom}) {next;}

         $t_def =~ s/^\(// ; $t_def =~ s/\)$// ;
         $data->{defstring}->{$t_dom} = $t_def ;
         $data->{class}->{$t_dom} = $t_class;
         $data->{aln}->{$t_dom} = '';
         $data->{pdb}->{$t_dom} = substr($t_dom, 1, 4) ;
         {
            my @t = split(/\,/, $t_def) ;
            my (@ch, @b, @e) ;
            foreach my $j ( 0 .. $#t) {
               my $t_frag ;
               $t_frag->{chain} = '_' ;
               if ($t[$j] =~ /:/) {
                  $t_frag->{chain} = substr($t[$j],0,1) ; }
               $t[$j] =~ s/^.\:// ;
               my ($b, $e) = (' ', ' ');
               if ($t[$j] =~ /.+\-.+/) {
                  ($t_frag->{b}, $t_frag->{e}) =
                  ($t[$j] =~ /(.+)\-(.+)/) ; }
               push @{$data->{frags}->{$t_dom}}, $t_frag ;

               if ($t_frag->{chain} eq '-' || $t_frag->{chain} eq '_') {
                  $t_frag->{chain} = ' ' ;
               } elsif (!exists $allchains->{$data->{pdb}->{$t_dom}}->{$t_frag->{chain}}) {
                  if ($t_frag->{chain} eq uc($t_frag->{chain})) {
                     my $lc = lc($t_frag->{chain}) ;
                     if (exists  $allchains->{$data->{pdb}->{$t_dom}}->{$lc}) {
                        print STDERR "WARNING defstring change: dom $t_dom old chain $t_frag->{chain} changed to $lc\n" ;
                        $t_frag->{chain} = $lc ; }
                  } elsif ($t_frag->{chain} eq lc($t_frag->{chain})) {
                     my $uc = uc($t_frag->{chain}) ;
                     if (exists $allchains->{$data->{pdb}->{$t_dom}}->{$uc}) {
                        print STDERR "WARNING defstring change: dom $t_dom old chain $t_frag->{chain} changed to $uc\n" ;
                        $t_frag->{chain} = $uc ; }
                  }
               }
            }
         }
      } elsif ((!exists $in->{doms}) ||
                (exists $in->{doms} && exists $in->{doms}->{$cur_dom})) {
         my $newseq = $line ;
         $newseq =~ tr/[a-w]/[A-W]/ ; #for 1-seq alns eg a.205.1.1 d1us7b_
         $newseq =~ tr/[y-z]/[Y-Z]/ ;
         $data->{aln}->{$cur_dom} .= $newseq ;
      }
   }
   $data->{alnlength} = length($data->{aln}->{$cur_dom}) ;
   close(ALNF) ;

   return $data ;

}


=head2 load_astral_clusters()

   Title:       load_astral_clusters()
   Function:    Loads ASTRAL sequence cluster definitions
   Args:        $_->{out} - pointer to hash to hold output
                $_->{pibase_specs} -  pibase_specs structure

   Returns:     Nothing - populates the specified $_->{out}
                {out}->{seqcl}->{seq identity}->{scop identifier} = cluster num
                {out}->{seqcl2cont}->{seq identity}->{cluster num}= [scop id,..]

=cut

sub load_astral_clusters {

   my $in = shift ;
   my $out = $in->{out} ;
   my $specs = $in->{specs} ;

   my $seqcl ;
   my $seqcl2cont ;
   {
      foreach my $seqid (keys %{$specs->{astral}->{seqcl}}) {
         open (ASTRALCL, $specs->{astral}->{seqcl}->{$seqid}) ;
         my $clusnum = 0 ;
         while (my $line = <ASTRALCL>) {
            chomp $line;
            if ($line =~ /^Rep/) {
               $line =~ s/^Rep: //g;
               $clusnum++ ; }

            if ($line =~ /score/) {
               $line =~ s/^\s+// ;
               my $scopid = substr($line, 0, 7) ;
               $seqcl->{$seqid}->{$scopid} = $clusnum ;
               push @{$seqcl2cont->{$seqid}->{$clusnum}}, $scopid ;
            }
         }
         close(ASTRALCL) ;
      }
   }

   $out->{seqcl} = $seqcl ;
   $out->{seqcl2cont} = $seqcl2cont ;

}


=head2 get_astral_classlist()

   Title:       get_astral_classlist()
   Function:    Get list of SCOP classes in the ASTRAL compendium
   Args:        $_->{pibase_specs} -  pibase_specs structure
   Returns:     Nothing - populates the specified $_->{out}
                ->{fam}->{scop_family} = number of domains in the family
                ->{sf}->{scop_superfamily} = number of domains in the superfamily

=cut

sub get_astral_classlist {

   my $in = shift;
   my $specs = $in->{specs} ;

   my $classlist ;
   foreach my $type (qw/fam sf/) {
      my @files=glob($specs->{asteroids}->{$type.'_aln'}."/*.fasta_aln");
      for (@files) {
         my $a = $_ ;
         $a =~ s/\.fasta_aln$// ;
         $a =~ s/.*\/// ;
         $classlist->{$type}->{$a}++ ;
       }
   }

   return $classlist ;
  
}


1 ;
