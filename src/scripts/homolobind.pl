#!/usr/local/bin/perl
=head1 NAME

homolobind.pl - predict binding sites by homology transfer

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

=head1 IMPLEMENTATION

=cut

use warnings;
use strict;
use homolobind ;

main() ;

=head2 SUB main()

=cut

sub main {

   homolobind::run_homolobind({ARGV => \@ARGV}) ;

}
