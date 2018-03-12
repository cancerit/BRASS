#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014-2018 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of BRASS.
#
# BRASS is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’."
########## LICENCE ##########


use ExtUtils::MakeMaker;
use File::ShareDir::Install;

install_share module => 'Sanger::CGP::Brass::Implement', 'share';

WriteMakefile(
  NAME          => 'brass',
  LICENSE       => 'agpl_3', # http://search.cpan.org/~dagolden/CPAN-Meta-2.142690/lib/CPAN/Meta/Spec.pm#license
  VERSION_FROM  => 'lib/Bio/Brass.pm',
  EXE_FILES     => [qw(
                        bin/brass-assemble
                        bin/brass_foldback_reads.pl
                        bin/brassI_filter.pl
                        bin/brassI_np_in.pl
                        bin/brassI_pre_filter.pl
                        bin/brassI_prep_bam.pl
                        bin/brass.pl
                        bin/collate_rg_regions.pl
                        bin/combineResults.pl
                        bin/compute_coverage.pl
                        bin/corrected_insertsize.pl
                        bin/coverage_merge.pl
                        bin/filter_with_microbes_and_remapping.pl
                        bin/findbp
                        bin/findbp-all
                        bin/get_abs_bkpts_from_clipped_reads.pl
                        bin/make-repeat-file
                        bin/match_rg_patterns_to_library.pl
                        bin/merge_double_rgs.pl
                        bin/rearr2bedpe
                        )],
  PREREQ_PM     => {
                     'Test::Fatal' => 0.013,
                     'Const::Fast' => 0.014,
                     'Devel::Cover' => 1.09,
                     'Pod::Coverage' => 0.23,
                     'Graph' => 0.96,
                     'Graph::Writer::Dot' => 2.06,
                     'Template' => 2.26,
                     'IPC::Run' => 0.92,
                     'Statistics::Basic' => 1.6611,
                     'Number::Format' => 1.75,
                     'Math::Combinatorics' => 0.09,
                     'Try::Tiny' => 0.22,
                     'Capture::Tiny' => 0.30,
                   }
);

package MY;
use File::ShareDir::Install qw( postamble );
