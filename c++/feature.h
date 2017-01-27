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

#ifndef FEATURE_H
#define FEATURE_H

#include <iosfwd>
#include <string>
#include <vector>

#include "cansam/intervalmap.h"

class string_cache {
public:
  string_cache() : share(NULL) { }
  ~string_cache();

  std::string operator[] (unsigned index) const { return strings[index]; }
  unsigned find(const std::string& text);
  void clear_lookup(); // FIXME perhap rename

private:
  std::vector<std::string> strings;
  struct sharing_map;
  sharing_map* share;
};


class feature {
public:
  feature(const std::string& name) : name_index(cache.find(name)) { }
  feature(const feature& f) : name_index(f.name_index) { }
  ~feature() { }

  feature& operator= (const feature& f)
    { name_index = f.name_index; return *this; }

  std::string name() const { return cache[name_index]; }

private:
  unsigned name_index;

  static string_cache cache;
};

enum feature_action { none, ignore_reads, filter_reads, transposon };

// This is a bit of a hack.  Features are read from STREAM or FILENAME, and go
// into one of the feature maps depending on action metadata from "track" lines:
//    filter:reads          filter_features
//    annotate:transposon   transposon_features
//    ignore                ignore_features
int insert(sam::interval_multimap<feature>& filter_features,
	   sam::interval_multimap<feature>& transposon_features,
	   sam::interval_multimap<feature>& ignore_features,
	   std::istream& stream, feature_action fixed_action = none);

int insert(sam::interval_multimap<feature>& filter_features,
	   sam::interval_multimap<feature>& transposon_features,
	   sam::interval_multimap<feature>& ignore_features,
	   const std::string& filename, feature_action fixed_action = none);

#endif
