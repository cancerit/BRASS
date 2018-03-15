/**   LICENCE
* Copyright (c) 2014-2018 Genome Research Ltd.
*
* Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
*
* This file is part of BRASS.
*
* BRASS is free software: you can redistribute it and/or modify it under
* the terms of the GNU Affero General Public License as published by the Free
* Software Foundation; either version 3 of the License, or (at your option) any
* later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
* details.
*
* You should have received a copy of the GNU Affero General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*
* 1. The usage of a range of years within a copyright statement contained within
* this distribution should be interpreted as being equivalent to a list of years
* including the first and last year specified and all consecutive years between
* them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
* 2009, 2011-2012’ should be interpreted as being identical to a statement that
* reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
* statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
* identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
* 2009, 2010, 2011, 2012’."
*/

// Author: John Marshall

#include <iostream>
#include <iomanip>
#include <set>
#include <string>
#include <cstdlib>
#include <cerrno>

#include <unistd.h>

#include "cansam/sam/alignment.h"
#include "cansam/sam/header.h"
#include "cansam/sam/stream.h"
#include "cansam/exception.h"

#include "version.inc"

using std::string;
using namespace sam;

typedef std::set<string> string_set;

void parse_flag_filter(const string& text, int& pos, int& neg) {
  char* eos;
  long value = strtol(text.c_str(), &eos, 0);

  if (*eos != '\0')
    throw std::runtime_error("invalid flag filter '-f " + text + "'");

  if (value > 0)  pos |=  value;
  else            neg |= -value;
}

int main(int argc, char** argv)
try {
  static const char copyright[] =
"Copyright (c) 2014-2018 Genome Research Ltd.\n"
"This is free software: you are free to change and redistribute it.\n"
"There is NO WARRANTY, to the extent permitted by law.\n"
"";

  static const char usage[] =
"Usage: filterout-bam [OPTION]... FILE FILTERFILE...\n"
"Options:\n"
"  -c        Write output BAM compressed [uncompressed]\n"
"  -f FLAGS  Select and discard alignment records matching FLAGS\n"
"  -o FILE   Write output BAM file to FILE rather than standard output\n"
"  -q NUM    Select and discard records with mapping quality less than NUM\n"
"  -Q NUM      ...or with original mapping quality (ZM:i) less than NUM\n"
"  -s NUM    Select and discard records with more than NUM suboptimal hits\n"
"  -S NUM      ...or with more than NUM original suboptimal hits (Z1:i)\n"
"  -v        Display file information and statistics\n"
"\n"
"Copies alignment records from FILE, discarding records whose qname is the\n"
"same as that of a selected record from any of the FILTERFILEs.  Records are\n"
"selected if they match any of the specified criteria.  By default, when no\n"
"criteria are specified, no records are selected or discarded.\n"
"";

  if (argc >= 2) {
    string arg = argv[1];
    if (arg == "--help") { std::cout << usage; return EXIT_SUCCESS; }
    else if (arg == "--version") {
      std::cout << "filterout-bam (Brass) " BRASS_VERSION "\n" << copyright;
      return EXIT_SUCCESS;
    }
  }

  int pos_flags = 0, neg_flags = 0;
  string output_filename = "-";
  int min_quality = 0;
  int min_orig_quality = 0;
  int max_suboptimal = -1;
  int max_orig_suboptimal = -1;
  bool verbose = false;
  bool compress_bam = false;

  int c;
  while ((c = getopt(argc, argv, ":c::f:o:q:Q:s:S:v")) >= 0)
    switch (c) {
    case 'c':  compress_bam = true;  break;
    case 'f':  parse_flag_filter(optarg, pos_flags, neg_flags);  break;
    case 'o':  output_filename = optarg;  break;
    case 'q':  min_quality = atoi(optarg);  break;
    case 'Q':  min_orig_quality = atoi(optarg);  break;
    case 's':  max_suboptimal = atoi(optarg);  break;
    case 'S':  max_orig_suboptimal = atoi(optarg);  break;
    case 'v':  verbose = true;  break;
    default:   std::cerr << usage; return EXIT_FAILURE;
  }

  if (argc - optind < 2)
    { std::cerr  << usage; return EXIT_FAILURE; }

  isamstream in(argv[optind++]);
  if (! in.is_open())
    throw sam::system_error("can't open ", in.filename(), errno);
  in.exceptions(std::ios::failbit | std::ios::badbit);

  string_set discard_qnames;

  while (optind < argc) {
    isamstream filter(argv[optind++]);
    if (! filter.is_open())
      throw sam::system_error("can't open ", filter.filename(), errno);
    filter.exceptions(std::ios::failbit | std::ios::badbit);

    collection headers;
    filter >> headers;

    unsigned long selected = 0;
    alignment aln;
    while (filter >> aln)
      if (((aln.flags() & pos_flags) == pos_flags &&
				(aln.flags() & neg_flags) == 0) ||
	  aln.mapq() < min_quality ||
	  (min_orig_quality && aln.aux("ZM", 255) < min_orig_quality) ||
	  (max_suboptimal >= 0 && aln.aux("X1", 0) > max_suboptimal) ||
	  (max_orig_suboptimal>=0 && aln.aux("Z1", 0) > max_orig_suboptimal) ||
	  false /* (this just regularises the lines above) */) {
	discard_qnames.insert(aln.qname());
	selected++;
      }

    if (verbose)
      std::clog << std::setw(14) << selected
		<< '\t' << filter.filename() << '\n';
  }

  if (verbose)
    std::clog << std::setw(14) << discard_qnames.size()
	      << "\tTotal read names selected to be discarded\n";

  std::ios_base::openmode bam_type = bam_format & ~compressed;
  if(compress_bam) bam_type = bam_format;

  osamstream out(output_filename, bam_type);
  if (! out.is_open())
    throw sam::system_error("can't write to ", output_filename, errno);
  out.exceptions(std::ios::failbit | std::ios::badbit);

  collection headers;
  in >> headers;
  out << headers;

  unsigned long discarded = 0;
  string qname_buffer;
  alignment aln;
  while (in >> aln)
    if (discard_qnames.find(aln.qname(qname_buffer)) == discard_qnames.end())
      out << aln;
    else
      discarded++;

  if (verbose)
    std::clog << std::setw(14) << discarded
	      << "\tAlignment records discarded\n";

  return EXIT_SUCCESS;
}
catch (const std::exception& e) {
  std::cout << std::flush;
  std::cerr << "filterout-bam: " << e.what() << std::endl;
  return EXIT_FAILURE;
}
