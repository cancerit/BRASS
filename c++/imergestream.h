/**   LICENCE
* Copyright (c) 2014-2017 Genome Research Ltd.
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

#ifndef IMERGESTREAM_H
#define IMERGESTREAM_H

#include "cansam/sam/alignment.h"
#include "cansam/sam/header.h"
#include "cansam/sam/stream.h"

namespace sam {

// FIXME This class assumes the headers are "similar enough".  When merging
// them has been implemented, it should be lifted to the Cansam library.

/// Input stream merging two SAM/BAM input streams
template <typename ISamStream2 = isamstream,
	  typename Compare = std::less<alignment> >
class imergestream {
public:
  /// Construct a merged stream from two opened input streams
  imergestream(isamstream& s1, ISamStream2& s2) : stream1(s1), stream2(s2) { }
  ~imergestream() { }

  /// Read the collection of headers
  bool operator>> (collection& headers) {
    stream1 >> headers;
    stream2 >> h2;
    // FIXME Merge the headers; need to rebase the alignments onto merged hdrs
    //headers.merge(h2);
    if (headers.ref_size() != h2.ref_size())
      throw "files have different numbers of references";
    for (sam::collection::ref_iterator it = headers.ref_begin(), it2 = h2.ref_begin(); it != headers.ref_end(); ++it, ++it2)
      if (it->name() != it2->name() || it->length() != it2->length())
	throw "files have a differing reference";

    // FIXME HACK Merge in the RG headers
    // FIXME should be const_iterator but weird pit thing
    for (collection::iterator it = h2.begin(); it != h2.end(); ++it)
      if (it->type_equals("RG")) {
	const readgroup& rg = dynamic_cast<const readgroup&>(*it);
	if (! has_group(headers, rg.id()))  headers.push_back(rg.str());
      }

    if (stream1 >> a1)  state = (stream2 >> a2)? bothpending : pending1;
    else  state = only2;

    return true;
  }

  /// Read an alignment record
  bool operator>> (alignment& aln) {
    switch (state) {
    case bothpending:
      if (less(a1,a2)) { aln.swap(a1); if (!(stream1 >> a1)) state = pending2; }
      else             { aln.swap(a2); if (!(stream2 >> a2)) state = pending1; }
      return true;

    case pending1:  aln.swap(a1); state = only1; return true;
    case pending2:  aln.swap(a2); state = only2; return true;

    case only1:  return stream1 >> aln;
    case only2:  return stream2 >> aln;

    default:  return false;  // Can't happen
    }
  }

private:
  Compare less;
  isamstream& stream1;
  ISamStream2& stream2;
  enum { bothpending, pending1, pending2, only1, only2 } state;
  alignment a1, a2;

  // HACK for the RG header merging
  collection h2;
  static bool has_group(const collection& headers, const std::string& key) {
    for (collection::const_iterator it = headers.begin(); it != headers.end(); ++it)
      if (it->type_equals("RG")) {
	const readgroup& rg = dynamic_cast<const readgroup&>(*it);
	if (rg.id() == key)  return true;
      }

    return false;
  }
};

} // namespace sam

#endif
