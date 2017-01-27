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

// rearrgroup.cpp -- Classes for rearrangement groups and sets thereof.

#include "rearrgroup.h"

#include <cstdlib>

#include "cansam/exception.h"

using std::string;
using namespace sam;


// This is an invalid read group name, used as a "no RG: field" indicator.
const char* const readgroup_set::NO_RG = "\t";

readgroup_set::readgroup_set(const collection& hdrs, scoord_t default_max,
			     string default_sample) {
  std::map<string, int> sample_index;

  for (collection::const_iterator it = hdrs.begin(); it != hdrs.end(); ++it)
    if (it->type_equals("RG")) {
      const readgroup& rg = dynamic_cast<const readgroup&>(*it);

      scoord_t max;
      header::const_iterator mi = rg.find("MI");
      if (mi != rg.end()) {
	const char* text = mi->value<const char*>();
	if (text[0] == 'Z' && text[1] == ':')  text += 2;
	max = atoi(text);
      }
      else {
	if (default_max < 0)
	  throw sam::exception("Read group " + rg.id() +
			       " has no MI: field (-m option required)");
	max = default_max;
      }

      string sample = rg.sample();
      readgroups[rg.id_c_str()] = readgroup_info(sample, max);
      sample_index.insert(make_pair(sample, 0));
    }

  if (! default_sample.empty()) {
    if (default_max < 0)
      throw sam::exception("No maximum insert size given for read pairs"
			   " without a read group (-m option required");

    readgroups[NO_RG] = readgroup_info(default_sample, default_max);
    sample_index.insert(make_pair(default_sample, 0));
  }
  else if (readgroups.empty())
    throw sam::exception("No read groups listed (-s/-m/etc options required)");

  for (std::map<string,int>::iterator it = sample_index.begin();
       it != sample_index.end(); ++it) {
    it->second = samples_.size();
    samples_.push_back(it->first);
  }

  for (readgroup_map::iterator it = readgroups.begin();
       it != readgroups.end(); ++it)
    it->second.sample_index = sample_index[it->second.sample];
}

const readgroup_info& readgroup_set::find(const alignment& aln) const {
  const char* rg = aln.aux("RG", NO_RG);
  readgroup_map::const_iterator it = readgroups.find(rg);
  if (it == readgroups.end()) {
    string rgstr = rg;
    if (rgstr == NO_RG)
      throw sam::exception("Read " + aln.qname() +
			   " has no RG: field (consider -s option)");
    else
      throw sam::exception("Read " + aln.qname() +
			   " has an unknown read group ('" + rgstr + "')");
  }

  return it->second;
}


interval::interval(const alignment& aln, coord_t pos, int strand,
		   coord_t ref_length, const readgroup_info& info) {
    // FIXME This is for short-insert-solexa only

    scoord_t pos0 = pos;

    if (strand == +1) {
      pos5 = pos0, pos3 = pos0 + info.max_insert;
      if (pos3 > ref_length)  pos3 = ref_length;
    }
    else {
      pos0 += aln.length();
      pos5 = pos0 - info.max_insert, pos3 = pos0;
      if (pos5 < 1)  pos5 = 1;
      if (pos3 > ref_length)  pos3 = ref_length;
    }
}


rearr_group::rearr_group(alignment& aln, const interval& alnL,
			 const interval& alnH, const readgroup_info& info,
			 const readgroup_set& readgroups)
  : overlapL(alnL), overlapH(alnH), max_insert(info.max_insert) {

  samples.assign(readgroups.samples().size(), per_sample());
  per_sample& sample = samples[info.sample_index];
  total_count = sample.count = 1;
  sample.readnames.assign(aln.qname_c_str(), aln.qname_length());

  canonical.swap(aln);
}

void rearr_group::insert(const alignment& aln, const interval& alnL,
			 const interval& alnH, const readgroup_info& info) {
#ifndef NOT_PARANOID
  if (! (aln.rindex() == canonical.rindex() &&
	 aln.strand() == canonical.strand() &&
	 aln.mate_rindex() == canonical.mate_rindex() &&
	 aln.mate_strand() == canonical.mate_strand()))
    throw std::logic_error("inserted alignment does not match group");
#endif

  overlapL *= alnL;
  overlapH *= alnH;

  if (max_insert < info.max_insert)  max_insert = info.max_insert;

  per_sample& sample = samples[info.sample_index];
  total_count++;
  sample.count++;
  if (! sample.readnames.empty())  sample.readnames += ';';
  sample.readnames.append(aln.qname_c_str(), aln.qname_length());
}

std::ostream& operator<< (std::ostream& out, const rearr_group& group) {
  const alignment& aln = group.canonical;

  out << aln.rname_c_str() << '\t' << aln.strand_char() << '\t'
      << group.overlapL.pos5 << '\t' << group.overlapL.pos3 << '\t'
      << aln.mate_rname_c_str() << '\t' << aln.mate_strand_char() << '\t'
      << group.overlapH.pos5 << '\t' << group.overlapH.pos3;

  rearr_group::const_sample_iterator it;
  for (it = group.samples.begin(); it != group.samples.end(); ++it)
    out << '\t' << it->count;

  if (! group.notes.empty())  out << '\t' << group.notes;
  else  out << "\t.";

  for (it = group.samples.begin(); it != group.samples.end(); ++it)
    if (it->count > 0)  out << '\t' << it->readnames;
    else  out << "\t.";

  return out;
}


rearr_group_set::rearr_group_set(const collection& refseqs)
  : lists(refseqs.ref_size()), rindex_(-1) {
}

rearr_group_set::iterator_pair rearr_group_set::complete_range() {
  for (size_t i = 1; i < lists.size(); i++)
    lists[0].splice(lists[0].end(), lists[i]);

  return make_pair(lists[0].begin(), lists[0].end());
}

void rearr_group_set::clear() {
  lists[0].clear();
  rindex_ = -1;
}
