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

// rearrgroup.h -- Classes for rearrangement groups and sets thereof.

#ifndef REARRGROUP_H
#define REARRGROUP_H

#include <algorithm>
#include <list>
#include <string>
#include <utility>
#include <vector>

#include "cansam/sam/alignment.h"
#include "cansam/sam/header.h"

using sam::coord_t;
using sam::scoord_t;
using sam::alignment;
using sam::collection;

// Information extracted from @RG headers

struct readgroup_info {
  std::string sample;
  int sample_index;
  scoord_t max_insert;

  readgroup_info() { }
  readgroup_info(const std::string& s, scoord_t mi)
    : sample(s), sample_index(0), max_insert(mi) { }
};

class readgroup_set {
public:
  readgroup_set(const collection& headers, scoord_t default_max,
		std::string default_sample);
  ~readgroup_set() { }

  const readgroup_info& find(const alignment& aln) const;

  const std::vector<std::string>& samples() const { return samples_; }

private:
  struct c_str_less {
    bool operator() (const char* a, const char* b) const
      { return strcmp(a, b) < 0; }
  };

  typedef std::map<const char*, readgroup_info, c_str_less> readgroup_map;
  readgroup_map readgroups;
  std::vector<std::string> samples_;

  static const char* const NO_RG;
};


struct interval {
  scoord_t pos5, pos3;

  interval() { }
  interval(scoord_t p5, scoord_t p3) : pos5(p5), pos3(p3) { }
  interval(const interval& rhs) : pos5(rhs.pos5), pos3(rhs.pos3) { }
  ~interval() { }

  interval(const alignment& aln, coord_t pos, int strand,
	   coord_t ref_length, const readgroup_info& info);

  // Assigns the result of set intersection with RHS.  If this and RHS
  // do not intersect, the result is an empty interval with pos5 > pos3.
  interval& operator*= (const interval& rhs) {
    if (rhs.pos5 > pos5)  pos5 = rhs.pos5;
    if (rhs.pos3 < pos3)  pos3 = rhs.pos3;
    return *this;
  }
};

// Returns whether the set intersection of LHS and RHS is non-empty.
inline bool intersect(const interval& lhs, const interval& rhs) {
  return ! (lhs.pos3 < rhs.pos5 || lhs.pos5 > rhs.pos3);
}


class rearr_group {
public:
  rearr_group() { }
  rearr_group(alignment& aln, const interval& alnL, const interval& alnH,
	      const readgroup_info& info, const readgroup_set& readgroups);
  ~rearr_group() { }

  friend std::ostream& operator<< (std::ostream&, const rearr_group&);

  void insert(const alignment& aln, const interval& alnL, const interval& alnH,
	      const readgroup_info& info);

  bool
  matches(const alignment& aln, const interval& alnL, const interval& alnH) {
    return aln.strand() == canonical.strand() &&
	   aln.mate_strand() == canonical.mate_strand() &&
	   intersect(overlapL, alnL) && intersect(overlapH, alnH) &&
	   1 /* complicated thing involving $dir */;
  }

  int rindex() const { return canonical.rindex(); }
  int mate_rindex() const { return canonical.mate_rindex(); }

//private: // FIXME  Decide on encapsulation of these fields
  alignment canonical;
  interval overlapL, overlapH;
  std::string notes;
  scoord_t max_insert;

  struct per_sample {
    int count;
    std::string readnames;
    per_sample() : count(0) { }
  };
  std::vector<per_sample> samples;
  int total_count;

  typedef std::vector<per_sample>::iterator sample_iterator;
  typedef std::vector<per_sample>::const_iterator const_sample_iterator;
};

std::ostream& operator<< (std::ostream& stream, const rearr_group& group);


class rearr_group_set {
public:
  typedef std::list<rearr_group>::iterator iterator;
  typedef std::pair<iterator, iterator> iterator_pair;

  rearr_group_set(const collection& refseqs);
  ~rearr_group_set() { }

  int rindex() const { return rindex_; }

  iterator_pair mate_range(const alignment& aln)
    { std::list<rearr_group>& list = lists[aln.mate_rindex()];
      return make_pair(list.begin(), list.end()); }

  iterator_pair complete_range();
  void clear();  // To be used only after complete_range()

  void insert(const rearr_group& group)
    { lists[group.mate_rindex()].push_back(group); rindex_ = group.rindex(); }

  iterator erase(iterator it) { return lists[it->mate_rindex()].erase(it); }

private:
  std::vector<std::list<rearr_group> > lists;
  int rindex_;
};

#endif
