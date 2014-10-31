/**   LICENCE
* Copyright (c) 2014 Genome Research Ltd.
*
* Author: Cancer Genome Project <cgpit@sanger.ac.uk>
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
*/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cerrno>

#include <unistd.h>

#ifndef WITHOUT_PSTREAMS
#include <pstreams/pstream.h>
#endif

#include "cansam/sam/alignment.h"
#include "cansam/sam/header.h"
#include "cansam/sam/stream.h"
#include "cansam/exception.h"

#include "version.inc"

using std::string;
using namespace sam;

void copy(osamstream& out, isamstream& in, isamstream& augment, bool pairwise) {
  collection headers, augment_headers;
  in >> headers;
  augment >> augment_headers;

  out << headers;

  const int pairing = FIRST_IN_PAIR | SECOND_IN_PAIR;
  alignment aln1, aug1;

  if (pairwise) {
    alignment aln1, aln2, aug1, aug2;
    while (in >> aln1) {
      if (! (in >> aln2))
	throw sam::exception(in.filename() + " ended with an orphaned read ('" +
			     aln1.qname() + "')");
      if (cmp_by_qname(aln1, aln2) != 0)
	throw sam::exception(in.filename() + " desynchronised ('" +
			     aln1.qname() + "', '" +  aln2.qname() + "')");

      if (! (augment >> aug1 >> aug2))
	throw sam::exception(augment.filename() + " at EOF instead of '" +
			     aln1.qname() + "' pair");

      if (aug2.flags() & FIRST_IN_PAIR)  swap(aug1, aug2);

      if (cmp_by_qname(aln1, aug1) != 0 || cmp_by_qname(aln2, aug2) != 0 ||
	  (aln1.flags() & pairing) != (aug1.flags() & pairing) ||
	  (aln2.flags() & pairing) != (aug2.flags() & pairing)) {
	std::ostringstream message;
	message << "Files out of sync"
		<< " ('" << aln1.qname() << "'/" << aln1.flags()
		<< ", '" << aln2.qname() << "'/" << aln2.flags()
		<< "; '" << aug1.qname() << "'/" << aug1.flags()
		<< ", '" << aug2.qname() << "'/" << aug2.flags() << ")";
	throw sam::exception(message.str());
      }

      aln1.push_back("ZM", aug1.mapq());
      aln2.push_back("ZM", aug2.mapq());

      alignment::iterator it;
      if ((it = aug1.find("X1")) != aug1.end())
	aln1.push_back("Z1", it->value<int>());
      if ((it = aug2.find("X1")) != aug2.end())
	aln2.push_back("Z1", it->value<int>());

      out << aln1 << aln2;
    }
  }
  else {
    alignment aln, aug;
    while (in >> aln) {
      if (! (augment >> aug))
	throw sam::exception(augment.filename() + " at EOF instead of '" +
			     aln.qname() + "'");
      if (cmp_by_qname(aln, aug) != 0 ||
	  (aln.flags() & pairing) != (aug.flags() & pairing)) {
	std::ostringstream message;
	message << "Files out of sync ('" << aln.qname() << "'/" << aln.flags()
		<< "; '" << aug.qname() << "'/" << aug.flags() << ")";
	throw sam::exception(message.str());
      }

      aln.push_back("ZM", aug.mapq());
      alignment::iterator it;
      if ((it = aug.find("X1")) != aug.end())
	aln.push_back("Z1", it->value<int>());

      out << aln;
    }
  }
}

int main(int argc, char** argv)
try {
  static const char copyright[] =
"Copyright (C) 2013 Genome Research Ltd.\n"
"This is free software: you are free to change and redistribute it.\n"
"There is NO WARRANTY, to the extent permitted by law.\n"
"";

  static const char usage[] =
"Usage: augment-bam -a FILE [-p] [-o FILE] [FILE | -c COMMAND...]\n"
"Options:\n"
"  -a FILE  Augment output records with original fields from FILE (required)\n"
"  -c       Invoke a command for its output rather than reading FILE\n"
"  -o FILE  Write output BAM file to FILE rather than standard output\n"
"  -p       Augment records in pairs (input files must be grouped by name)\n"
"\n"
"Copies FILE or the output of COMMAND, converting to BAM if necessary,\n"
"and augmenting alignment records with fields from the corresponding record\n"
"from the augmenting file, which must contain records in the same order:\n"
"  ZM:i  Original MAPQ field\n"
"  Z1:i  Original X1:i field (BWA's suboptimal hit count), if present\n"
"";

  if (argc >= 2) {
    string arg = argv[1];
    if (arg == "--help") { std::cout << usage; return EXIT_SUCCESS; }
    else if (arg == "--version") {
      std::cout << "augment-bam (Brass) " << BRASS_VERSION << "\n" << copyright;
      return EXIT_SUCCESS;
    }
  }

  string augment_filename;
  string output_filename = "-";
  bool invoke_command = false;
  bool pairwise = false;

  opterr = 0;
  int c;
  while ((c = getopt(argc, argv, "+a:co:p")) >= 0)
    switch (c) {
    case 'a':  augment_filename = optarg;  break;
    case 'c':  invoke_command = true;  break;
    case 'o':  output_filename = optarg;  break;
    case 'p':  pairwise = true;  break;
    default:   std::cerr << usage; return EXIT_FAILURE;
  }

  if (augment_filename.empty() || (argc - optind > 1 && ! invoke_command))
    { std::cerr << usage; return EXIT_FAILURE; }

  isamstream augment(augment_filename);
  if (! augment.is_open())
    throw sam::system_error("can't open ", augment_filename, errno);

  int status = EXIT_SUCCESS;

  if (invoke_command) {
#ifdef WITHOUT_PSTREAMS
    throw std::logic_error("-c option not available (rebuild with PSTREAMS)");
#else
    osamstream out(output_filename, bam_format & ~compressed);
    if (! out.is_open())
      throw sam::system_error("can't write to ", output_filename, errno);

    std::vector<string> argvector(&argv[optind], &argv[argc]);
    redi::pstreambuf processbuf(argv[optind], argvector, std::ios::in);
    isamstream in(&processbuf);
    in.set_filename('`' + string(argv[optind]) + '`');

    copy(out, in, augment, pairwise);

    if (! processbuf.close())
      throw sam::system_error("can't close pipe ", in.filename(), errno);

    int child_status = processbuf.status();
    if (WIFSIGNALED(child_status)) {
      std::ostringstream message;
      message << argv[optind] << " killed by signal " << WTERMSIG(child_status);
      throw sam::exception(message.str());
    }
    else if (WIFEXITED(child_status))
      status = WEXITSTATUS(child_status);
#endif
  }
  else {
    isamstream in((optind < argc)? argv[optind] : "-");
    if (! in.is_open())
      throw sam::system_error("can't open ", in.filename(), errno);

    osamstream out(output_filename, bam_format & ~compressed);
    if (! out.is_open())
      throw sam::system_error("can't write to ", output_filename, errno);

    copy(out, in, augment, pairwise);
  }

  return status;
}
catch (const std::exception& e) {
  std::cout << std::flush;
  std::cerr << "augment-bam: " << e.what() << std::endl;
  return EXIT_FAILURE;
}
