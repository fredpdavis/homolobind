=head1 NAME

homolobind::pilig - perl module for pibase-ligbase overlap calculations

=head1 DESCRIPTION

Perl module with routines to cross-query pibase and ligbase
to get small molecule - protein interaction site overlap
statistics

=head1 AUTHOR

Fred P. Davis, HHMI-JFRC (davisf@janelia.hhmi.org)

=head1 LICENCE AND COPYRIGHT

Copyright 2005,2010 Fred P. Davis (davisf@janelia.hhmi.org).
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

=head2 LOGIC

Note: Code is a streamlined version of pibase::pilig()

=cut

package homolobind::pilig ;
use strict;
use warnings;


use Exporter;
our @ISA = qw/Exporter/ ;
our @EXPORT_OK = qw// ;

use homolobind::pibase ;
use homolobind::ASTRAL ;
use homolobind::SGE ;
use Cwd qw/getcwd/ ;
use File::Temp qw/tempfile tempdir/ ;
use File::Path qw/mkpath/ ;

use File::Basename qw/basename/ ;
use Sys::Hostname qw/hostname/ ;
use File::Temp qw/tempfile tempdir/ ;


sub set_pilig_specs {

   my $params ;
   $params->{PARAM_MIN_LIGMW} = 250;
   $params->{PARAM_MAX_LIGMW} = 1000;

   $params->{PARAM_INTERFACE_DIST_THRESH} = 5.0 ;

   $params->{PARAM_MIN_PEPTIDE_LENGTH} = 5 ;
   $params->{PARAM_MIN_NUMRES_INTERACTING_WITH_PEPTIDE} = 5 ;
   $params->{PARAM_MIN_NUMCONTACTS} = 500 ; #change to minimum SASA cutoff
   $params->{PARAM_MIN_BS_SIMILARITY} = 0.9 ; #threshold to merge binding sites.

   return $params ;

}


sub readin_ligassignments {

   my $in = shift ;

   my $fn = $in->{fn} ;

   my $class2alnlength = $in->{class2alnlength};
   my $standardres= $in->{standardres};
   my $liginfo = $in->{liginfo};

# NOTE: different frmo the in->{cluster_fl} flag that actually _does_ the
#       clustering and prints out memberships to the outfiles specified file

# if this is the clustered assignment list, only read in cluster representative
# (arbitrarily chosen - first one in the list)
   my $clustrep_fl = 0;
   my $clusters_seen ;  # keep track of what clusters have been seen already
   if (exists $in->{clustrep_fl} && $in->{clustrep_fl} == 1) {
      $clustrep_fl = 1 ; }


   my $class2sid ;

   my $ligbits ;
   print STDERR "NOW reading LIG assignment\n" ;
   open(LIGFH, $fn) ;
   while (my $line = <LIGFH>) {
      if ($line =~ /^#/ ||
          $line =~ /^Warning: no access to/ ||
          $line =~ /^Thus no job control/ ) {next;}

      chomp $line;
      my @fields = split(/\t/, $line) ;

      my $cluster_num ;
      if ($clustrep_fl) {
         $cluster_num = shift @fields ; }

      my ($pdb, $sid, $osid, $classtype, $class, $btype, $alnposstring,
          $alnlength, $ligcod, $ligid ) = @fields ;
      $class2sid->{$classtype}->{$class}->{$sid}++ ;

      $ligcod =~ s/ //g ;
      my $ligsig = $pdb."\t".$ligcod."\t".$ligid ;

      if (exists $standardres->{$ligcod}) {next;}
      if ($ligcod eq 'DUM' || $ligcod eq 'UNX' ||
             $ligcod eq 'UNK' || $ligcod eq 'UNL') {next;}

      if (!exists $liginfo->{mw}->{$ligcod}) {
         next; }

      if (!exists $liginfo->{mwinrange}->{$ligcod}) {
         next;}


      if ($clustrep_fl) {
         if (!exists $clusters_seen->{$classtype}->{$cluster_num}) {
            $clusters_seen->{$classtype}->{$cluster_num}++ ;
# don't have to keep track of what was chosen, since only 1 line is needed
# if a line has been read, that cluster is done...
         } else {
            next;
         }
      }

      $class2alnlength->{$classtype}->{$class} = $alnlength ;

      $ligbits->{$classtype}->{$sid}->{$ligsig} = Bit::Vector->new($alnlength) ;

      my @alnpos = split(/\,/, $alnposstring) ;
      foreach my $alnpos (@alnpos) {
         $ligbits->{$classtype}->{$sid}->{$ligsig}->Bit_On($alnpos); }


      if (!exists $ligbits->{$classtype}->{$sid}->{cumulative}) {
         $ligbits->{$classtype}->{$sid}->{cumulative} =
            $ligbits->{$classtype}->{$sid}->{$ligsig}->Clone();
      } else {
         $ligbits->{$classtype}->{$sid}->{cumulative}->Union(
            $ligbits->{$classtype}->{$sid}->{cumulative},
            $ligbits->{$classtype}->{$sid}->{$ligsig}) ;
      }
   }
   close(LIGFH) ;

   return $ligbits ;
}


sub readin_pepnuciassignments {

   my $in = shift ;

   my $fn = $in->{fn} ;

   my $class2alnlength = $in->{class2alnlength};
   my $standardres= $in->{standardres};
   my $pilig_specs = set_pilig_specs() ;

# if this is the clustered assignment list, only read in cluster representative
# (arbitrarily chosen - first one in the list)
   my $clustrep_fl = 0;
   my $clusters_seen ;  # keep track of what clusters have been seen already
   if (exists $in->{clustrep_fl} && $in->{clustrep_fl} == 1) {
      $clustrep_fl = 1 ; }

   my $targetch_info ;
   my $pep_bits = {};
   my $nuc_bits = {};
   print STDERR "NOW reading PEPNUCI assignment\n" ;
   open(PEPNUCIF, $fn) ;
   my $class2sid ;
   my $pepnuci_info ;
   while (my $line = <PEPNUCIF>) {
      if ($line =~ /^#/ ||
          $line =~ /^Warning: no access to/ ||
          $line =~ /^Thus no job control/ ) {next;}

      chomp $line;
      my @fields = split(/\t/, $line) ;

      my $cluster_num ;
      if ($clustrep_fl) {
         $cluster_num = shift @fields ; }

      my ( $pdb, $sid, $osid, $classtype, $class, 
           $targetch_type, $targetch_sid,
           $targetch_len, $alnposstring, $alnlength ) = @fields ;
      $class2sid->{$classtype}->{$class}->{$sid}++ ;

      my $t_alnposstring = $alnposstring ;
      $t_alnposstring =~ s/\,//g ; $t_alnposstring =~ s/undef//g ;
      if ($t_alnposstring eq '') {next;}

      if ($clustrep_fl) {
         if (!exists $clusters_seen->{$classtype}->{$cluster_num}) {
            $clusters_seen->{$classtype}->{$cluster_num}++ ;
         } else {
            next ;
         }
      }


      if ($targetch_type eq 'p' &&
          $targetch_len < $pilig_specs->{PARAM_MIN_PEPTIDE_LENGTH}) {
         next;
      }
      my @t_alnpos = split(/\,/, $alnposstring) ;
      if (($#t_alnpos + 1) <
          $pilig_specs->{PARAM_MIN_NUMRES_INTERACTING_WITH_PEPTIDE}) {
         next; }

      $class2alnlength->{$classtype}->{$class} = $alnlength ;
      $targetch_info->{$targetch_sid}->{chain_type} = $targetch_type ;
      $targetch_info->{$targetch_sid}->{chain_length} = $targetch_len ;

      my $bitref = $pep_bits ;
      if ($targetch_type eq 'n') {
         $bitref = $nuc_bits ; }

      $bitref->{$classtype}->{$sid}->{$targetch_sid} =
         Bit::Vector->new($alnlength) ;

      my @alnpos = split(/\,/, $alnposstring) ;
      foreach my $alnpos (@alnpos) {
         if ($alnpos eq 'undef') {next;}
         $bitref->{$classtype}->{$sid}->{$targetch_sid}->Bit_On($alnpos); }


      if (!exists $bitref->{$classtype}->{$sid}->{cumulative}) {
         $bitref->{$classtype}->{$sid}->{cumulative} =
            $bitref->{$classtype}->{$sid}->{$targetch_sid}->Clone();
      } else {
         $bitref->{$classtype}->{$sid}->{cumulative}->Union(
            $bitref->{$classtype}->{$sid}->{cumulative},
            $bitref->{$classtype}->{$sid}->{$targetch_sid}) ;
      }
   }
   close(PEPNUCIF) ;

   return {
      pep => $pep_bits,
      nuc => $nuc_bits,
      chain_info => $targetch_info
   } ;

}


sub readin_piassignments {

   my $in = shift;
   my $fn = $in->{fn} ;
   my $class2alnlength = $in->{class2alnlength};
   my $pb = $in->{pb} ;
   my $pibits ;
   my $interfaces ;

# if this is the clustered assignment list, only read in cluster representative
# (arbitrarily chosen - first one in the list)
   my $clustrep_fl = 0;
   my $clusters_seen ;  # keep track of what clusters have been seen already
   if (exists $in->{clustrep_fl} && $in->{clustrep_fl} == 1) {
      $clustrep_fl = 1 ; }

   print STDERR "NOW reading PI assignment\n" ;
   my $classes2int ;
   my $class2chains2intside;
   open(PIFH, $fn) ;
   while (my $line = <PIFH>) {
      if ($line =~ /^#/ || $line =~ /^Warning: no access to/ ||
         $line =~ /^Thus no job control/ ) {next;}

      chomp $line;

      my @t = split(/\t/, $line) ;
      my $cluster_num ;
      if ($clustrep_fl) {
         $cluster_num = shift @t ; }

      my ($pdb, $sid, $osid, $classtype, $class, $obtype, $alnposstring,
          $alnlength, $sid1, $sid2, $fam1, $fam2, $chains) = @t ;
      my $sid12 = $sid1."\t".$sid2 ;

#if flag specified, don't read in data for j.* peptide domains
      if (exists $in->{dont_read_jdomains_fl} &&
          ($pb->{sid2class}->{fam}->{$sid1} =~ /^j/ ||
           $pb->{sid2class}->{fam}->{$sid2} =~ /^j/)) {
         next;
      }


      my @alnpos = split(/\,/, $alnposstring) ;
      {
         my @t = ();
         foreach my $p ( @alnpos) {
            if ($p ne 'undef') {push @t, $p;} }
         @alnpos = @t ;
      }
      if ($#alnpos < 0) {next;}

      my $side = 1; if ($sid eq $sid2) {$side = 2;}


      if ($clustrep_fl) {
         if ( !exists $clusters_seen->{$classtype}->{$cluster_num} ||
              $clusters_seen->{$classtype}->{$cluster_num} eq $sid12 ) {
            $clusters_seen->{$classtype}->{$cluster_num} = $sid12 ;
# have to keep track of what was chosen, since need to read in both sides
# of the interface...
         } else {
            next;
         }
      }

      $class2chains2intside->{$classtype}->{$class}->{$chains}->{$sid12."\n".$side}++ ;
      {
         my $temp_class1 = $pb->{sid2class}->{fam}->{$sid1} ;
         my $temp_class2 = $pb->{sid2class}->{fam}->{$sid2} ;
         my $sf1 = $pb->{sid2class}->{sf}->{$sid1} ;
         my $sf2 = $pb->{sid2class}->{sf}->{$sid2} ;
         my $temp_revfl = 0 ;
         if ($temp_class2 lt $temp_class1) {
            $temp_class1 = $pb->{sid2class}->{fam}->{$sid2} ;
            $temp_class2 = $pb->{sid2class}->{fam}->{$sid1} ;
            $sf1 = $pb->{sid2class}->{sf}->{$sid2} ;
            $sf2 = $pb->{sid2class}->{sf}->{$sid1} ;
            $temp_revfl = 1 ;
         }
         $classes2int->{sf}->{$sf1."\t".$sf2}->{$sid12} = $temp_revfl ;
         $classes2int->{fam}->{$temp_class1."\t".$temp_class2}->{$sid12} =
            $temp_revfl ;
      }

      $interfaces->{$sid12}->{$side}->{sid} = $sid ;
      $interfaces->{$sid12}->{chains} = $chains ;
      $interfaces->{$sid12}->{pdb} = $pdb;

      $class2alnlength->{$classtype}->{$class} = $alnlength ;

#      my @btypes = ($obtype) ;
# fpd090104_1840  - changed this so that when P bits are set in
#                   combine_piligpepexpbits() it includes intra,inter,and peptide
      my $btype = 'Pinter';
      if ($chains eq 'same') {
         $btype = "Pintra" ; }

      $interfaces->{$sid12}->{$side}->{pibits}->{$classtype} =
         Bit::Vector->new($alnlength) ;

      foreach my $alnpos (@alnpos)  {
         if ($alnpos eq 'undef') {next;}
         $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Bit_On($alnpos);}


      if (!exists $pibits->{$classtype}->{$sid}->{$btype}) {
         $pibits->{$classtype}->{$sid}->{$btype} =
            $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}->Clone();
      } else {
         $pibits->{$classtype}->{$sid}->{$btype}->Union(
            $pibits->{$classtype}->{$sid}->{$btype},
            $interfaces->{$sid12}->{$side}->{pibits}->{$classtype}) ;
      }
   }
   close(PIFH) ;


   return {
      pibits => $pibits,
      interfaces => $interfaces
   } ;

}


sub readin_expassignments {

   my $in = shift ;
   my $fn = $in->{fn};
   my $class2alnlength = $in->{class2alnlength};

   print STDERR "NOW reading EXP assignment\n" ;
   open(EXPFH, $fn) ;
   my $expbits ;
   while (my $line = <EXPFH>) {
      if ($line =~ /^#/ ||
            $line =~ /^Warning: no access to/ ||
            $line =~ /^Thus no job control/ ) {next;}
      chomp $line;

      my ( $pdb, $sid, $osid, $classtype,
           $class, $btype, $alnposstring, $alnlength) = split(/\t/, $line) ;

      $class2alnlength->{$classtype}->{$class} = $alnlength ;

      if (!exists $expbits->{$classtype}->{$class}) {
         $expbits->{$classtype}->{$class} =
            Bit::Vector->new($alnlength) ;
      }

      my @alnpos = split(/\,/, $alnposstring) ;
      foreach my $alnpos (@alnpos)  {
         if ($alnpos eq "undef") {next;}
         $expbits->{$classtype}->{$class}->Bit_On($alnpos);
      }
   }
   close(EXPFH) ;

   return $expbits ;

}


sub _list_standardres {

   return {
   ALA => 'A' ,
   ARG => 'R' ,
   ASN => 'N' ,
   ASP => 'D' ,
   CYS => 'C' ,
   GLN => 'Q' ,
   GLU => 'E' ,
   GLY => 'G' ,
   HIS => 'H' ,
   HSD => 'H' ,
   HSE => 'H' ,
   ILE => 'I' ,
   LEU => 'L' ,
   LYS => 'K' ,
   MET => 'M' ,
   PHE => 'F' ,
   PRO => 'P' ,
   SER => 'S' ,
   THR => 'T' ,
   TRP => 'W' ,
   TYR => 'Y' ,
   VAL => 'V',
   UNK => 'X',
   UNX => 'X',
   '  C' => 'c',
   '  G' => 'g',
   '  A' => 'a',
   '  T' => 't',
   '  U' => 'u',
   '  I' => 'i',
   'C' => 'c',
   'G' => 'g',
   'A' => 'a',
   'T' => 't',
   'U' => 'u',
   'I' => 'i',
   '+C' => 'c',
   '+G' => 'g',
   '+A' => 'a',
   '+T' => 't',
   '+U' => 'u',
   '+I' => 'i'
   } ;

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


sub _pilig_astral_preload {

   my $in = shift;
   my $specs = $in->{specs} ;

   my $astral ;
   

   $astral->{classes} = homolobind::ASTRAL::get_astral_classlist({
      specs => $specs}) ;

   homolobind::ASTRAL::load_astral_clusters({
      specs => $specs,
      out => $astral
   }) ;

   return $astral ;
}


sub _pilig_ligbase_preload {

   my $in = shift ;
   my $dbh = $in->{dbh} ;

   my $query="SELECT pdb_code, a.ligand_code, a.ligand_idnum, a.ligand_chain, ".
         "a.residue_num, a.ligand_chain FROM active_residues as a";

   my ($pdb, $ligcod, $ligid, $ligchain, $resno, $pdbchain) =
      pibase::mysql_fetchcols($dbh,$query) ;

   my ($pdb2ligid, $pdb2res2ligid) ;
   foreach my $j ( 0 .. $#{$ligcod}) {
      if ($pdbchain->[$j] eq '') {
         $pdbchain->[$j] = ' '; }

      $pdb2ligid->{$pdb->[$j]}->{$ligcod->[$j]."\n".$ligid->[$j]}++ ;
      $pdb2res2ligid->{$pdb->[$j]}->{$resno->[$j]."\n".$pdbchain->[$j]}->{$ligcod->[$j]."\n".$ligid->[$j]}++ ;
   }

   my $lb = {
      pdb2ligid => $pdb2ligid,
      pdb2res2ligid => $pdb2res2ligid
   } ;

   return $lb ;

}


sub _pilig_tod_pibase_preload {

   my $in = shift ;
   my $specs = $in->{specs} ;
   my $astral = $in->{astral} ;
   my $pilig_specs = set_pilig_specs() ;

   my ($bdp_2_raw, $bdp2pdb) ;
   {
      my $tablespecs = homolobind::pibase::table_spec("bdp_files") ;
      my $f2i;
      foreach my $j (0 .. $#{$tablespecs->[0]->{specs}->{field_name}}) {
         $f2i->{$tablespecs->[0]->{specs}->{field_name}->[$j]} = $j ; }

      open(TABLEF, $specs->{pibase}->{tod_dir}."/bdp_files") ;
      while (my $line = <TABLEF>) {
         chomp $line;
         my @t = split(/\t/, $line) ;
         $bdp_2_raw->{$t[$f2i->{'bdp_id'}]} = $t[$f2i->{'raw_pdb'}] ;
         $bdp2pdb->{$t[$f2i->{'bdp_id'}]} = $t[$f2i->{'pdb_id'}] ;
      }
      close(TABLEF) ;
   }

   my $subsetsource_name2id ;
   {
      my $tablespecs = homolobind::pibase::table_spec("subsets_source") ;
      my $f2i;
      foreach my $j (0 .. $#{$tablespecs->[0]->{specs}->{field_name}}) {
         $f2i->{$tablespecs->[0]->{specs}->{field_name}->[$j]} = $j ; }

      open(TABLEF, $specs->{pibase}->{tod_dir}."/subsets_source") ;
      while (my $line = <TABLEF>) {
         chomp $line;
         my @t = split(/\t/, $line) ;
         $subsetsource_name2id->{$t[$f2i->{'subset_source'}]} =
            $t[$f2i->{'subset_source_id'}] ;
      }
      close(TABLEF) ;
   }

   my $sid2pdb ; #new100412_0950 
   my $sid2class ;
   {
      my $tablespecs = homolobind::pibase::table_spec("subsets") ;
      my $f2i;
      foreach my $j (0 .. $#{$tablespecs->[0]->{specs}->{field_name}}) {
         $f2i->{$tablespecs->[0]->{specs}->{field_name}->[$j]} = $j ; }

      open(TABLEF, $specs->{pibase}->{tod_dir}."/subsets") ;
      while (my $line = <TABLEF>) {
         chomp $line;
         my @t = split(/\t/, $line) ;
         if ($t[$f2i->{'subset_source_id'}] != $subsetsource_name2id->{'scop'}){
            next;}

         if (!defined $t[$f2i->{'bdp_id'}] || $t[$f2i->{'bdp_id'}] eq '' ||
             $t[$f2i->{'bdp_id'}] eq 'NULL') {
            next;}

         $sid2class->{fam}->{$t[$f2i->{'subset_id'}]} = $t[$f2i->{'class'}];
         my $fam = $t[$f2i->{'class'}] ;
         my ($sf) = ($fam =~ /([a-z]\.[0-9]+\.[0-9]+)/) ;
         $sid2class->{sf}->{$t[$f2i->{'subset_id'}]} = $sf ;

         $sid2pdb->{$t[$f2i->{'subset_id'}]} = $bdp2pdb->{$t[$f2i->{'bdp_id'}]};
      }
      close(TABLEF) ;
   }

   my $pdbchains ;
   {
      my $tablespecs = homolobind::pibase::table_spec("bdp_chains") ;
      my $f2i;
      foreach my $j (0 .. $#{$tablespecs->[0]->{specs}->{field_name}}) {
         $f2i->{$tablespecs->[0]->{specs}->{field_name}->[$j]} = $j ; }

      open(TABLEF, $specs->{pibase}->{tod_dir}."/bdp_chains") ;
      while (my $line = <TABLEF>) {
         chomp $line;
         my @t = split(/\t/, $line) ;
            if ($bdp_2_raw->{$t[$f2i->{'bdp_id'}]} != 1) {next;}
            if (!defined $t[$f2i->{'pdb_chain_id'}] ||
                $t[$f2i->{'pdb_chain_id'}] eq '') {next;}

            if ($t[$f2i->{'chain_type'}] ne 'p') {next;}

            $pdbchains->{$bdp2pdb->{$t[$f2i->{'bdp_id'}]}}->{$t[$f2i->{'real_chain_id'}]}++ ;
      }
      close(TABLEF) ;
   }


   my $pb = {
      pdbchains => $pdbchains,
      sid2class => $sid2class,
      sid2pdb => $sid2pdb,
   } ;

   if (exists $in->{onlythese}) {
      my $realpb ;
      foreach my $type (keys %{$in->{onlythese}}) {
         $realpb->{$type} = $pb->{$type}; }
      return $realpb ;
   } else {
      return $pb ;
   }

}


sub _pilig_load_liginfo {

   my $in = shift ;
   my $fn = $in->{fn};
   my $pilig_specs = set_pilig_specs() ;

   my $liginfo ;
   open(LIGINFOF, $fn) ;
   while (my $line = <LIGINFOF>) {
      if ($line =~ /^\#/) {next;}
      chomp $line;

      my @t = split(/\t/, $line) ;
      my $lig = shift @t ;
      $liginfo->{ligs}->{$lig}++ ;

      ($liginfo->{name}->{$lig},
       $liginfo->{formula}->{$lig},
       $liginfo->{mw}->{$lig},
       $liginfo->{numatoms}->{$lig},
       $liginfo->{numatoms_nonh}->{$lig} ) = @t ;

      if ($liginfo->{mw}->{$lig} eq '') {
         delete $liginfo->{mw}->{$lig};
      } elsif ($liginfo->{mw}->{$lig} >= $pilig_specs->{PARAM_MIN_LIGMW} &&
               $liginfo->{mw}->{$lig} <= $pilig_specs->{PARAM_MAX_LIGMW}) {
         $liginfo->{mwinrange}->{$lig}++ ;
      }
   }
   close(LIGINFOF) ;

   return $liginfo ;

}


sub combine_piligpepbits {

   my $in = shift ;
   my $pb = $in->{pb} ;
   my $ligbits = $in->{ligbits} ;
   my $liginfo = $in->{liginfo} ;
   my $pibits = $in->{pibits} ;
   my $pepbits = $in->{pepbits} ;

   my $class2bits ;
   my $class2ligs ;
   my $class2allligbits ;

   foreach my $classtype (qw/fam sf/) {
      foreach my $sid (keys %{$ligbits->{$classtype}}) {
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;

         foreach my $ligsig (keys %{$ligbits->{$classtype}->{$sid}}) {
            if ($ligsig eq 'cumulative') {next;}
            $class2ligs->{$classtype}->{$class}->{$ligsig}++;

# fpd090303_2034: changed the elements of this array from (ligsig, ligbits) to
# (ligsig, ligbits, sid) - have to change collate_perfam and collate_perinst
# to respect this

            push @{$class2allligbits->{$classtype}->{$class}},
               [$ligsig,
                $ligbits->{$classtype}->{$sid}->{$ligsig}, $sid] ;

            my ($pdb, $ligcod, $ligid) = split(/\t/, $ligsig) ;

            $class2ligs->{$classtype}->{$class}->{$ligsig}++;
            if (!exists $class2bits->{$classtype}->{$class}->{'L'}) {
               $class2bits->{$classtype}->{$class}->{'L'} =
                  $ligbits->{$classtype}->{$sid}->{$ligsig}->Clone() ;
            } else {
               $class2bits->{$classtype}->{$class}->{'L'}->Union(
                  $class2bits->{$classtype}->{$class}->{'L'},
                  $ligbits->{$classtype}->{$sid}->{$ligsig}) ;
            }
         }
      }

      foreach my $sid (keys %{$pibits->{$classtype}}) {
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;
         foreach my $btype (keys %{$pibits->{$classtype}->{$sid}}) {
            if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
               $class2bits->{$classtype}->{$class}->{$btype} =
                  $pibits->{$classtype}->{$sid}->{$btype}->Clone() ;
            } else {
               $class2bits->{$classtype}->{$class}->{$btype}->Union(
                  $class2bits->{$classtype}->{$class}->{$btype},
                  $pibits->{$classtype}->{$sid}->{$btype}) ;
            }

            if (!exists $class2bits->{$classtype}->{$class}->{'P'}) {
               $class2bits->{$classtype}->{$class}->{'P'} =
                  $pibits->{$classtype}->{$sid}->{$btype}->Clone() ;
            } else {
               $class2bits->{$classtype}->{$class}->{'P'}->Union(
                  $class2bits->{$classtype}->{$class}->{'P'},
                  $pibits->{$classtype}->{$sid}->{$btype}) ;
            }
         }
      }

      foreach my $sid (keys %{$pepbits->{$classtype}}) {
         my $class = $pb->{sid2class}->{$classtype}->{$sid} ;
         foreach my $targetch (keys %{$pepbits->{$classtype}->{$sid}}) {
            foreach my $btype ('p', 'P') {
               if (!exists $class2bits->{$classtype}->{$class}->{$btype}) {
                  $class2bits->{$classtype}->{$class}->{$btype} = 
                     $pepbits->{$classtype}->{$sid}->{$targetch}->Clone() ;
               } else {
                  $class2bits->{$classtype}->{$class}->{$btype}->Union(
                     $class2bits->{$classtype}->{$class}->{$btype},
                     $pepbits->{$classtype}->{$sid}->{$targetch}) ;
               }
            }
         }
      }
   }

   return {
      class2bits => $class2bits,
      class2ligs => $class2ligs,
      class2allligbits => $class2allligbits,
   } ;
   
}

1 ;
