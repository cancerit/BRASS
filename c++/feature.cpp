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

#include "feature.h"

#include <fstream>
#include <map>
#include <sstream>
#include <cerrno>

#ifndef WITHOUT_BOOST
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#endif

#include "cansam/exception.h"

using std::string;
using namespace sam;

struct string_cache::sharing_map {
  std::map<string, unsigned> index;
};

string_cache::~string_cache() {
  delete share;
}

void string_cache::clear_lookup() {
  delete share;
  share = NULL;
}

unsigned string_cache::find(const string& text) {
  if (share == NULL)  share = new sharing_map;

  std::map<string, unsigned>::iterator it = share->index.find(text);
  if (it != share->index.end())  return it->second;

  unsigned ndx = strings.size();
  strings.push_back(text);
  share->index[text] = ndx;
  return ndx;
}


string_cache feature::cache;


typedef std::map<std::string, std::string> dictionary;

// Parse a BED metadata line into DICT.  Syntax is essentially
//    METAKEY KEY=VALUE KEY="QUOTED \"VALUE\"" [...]
void insert_metadata(dictionary& dict, const string& s) {
  static const char whitespace[] = " \t";

  // Insert "track"/"browser"/etc keyword keyed by a sentinel
  size_t pos = s.find_first_of(whitespace);
  dict[""].assign(s, 0, pos);

  string key;
  while ((pos = s.find_first_not_of(whitespace, pos)) != string::npos) {
    size_t eqpos = s.find('=', pos);
    if (eqpos == string::npos)  throw sam::bad_format("No equals");

    key.assign(s, pos, eqpos - pos);
    pos = eqpos+1;

    if (s[pos] == '\'' || s[pos] == '\"') {
      const char quote_backslash[] = { s[pos++], '\\', '\0' };
      string& value = dict[key];
      value.clear();

      while (true) {
	size_t epos = s.find_first_of(quote_backslash, pos);
	if (epos == string::npos)  throw sam::bad_format("Unterminated quotes");

	value.append(s, pos, epos - pos);
	pos = epos+1;
	if (s[epos] == '\\') {
	  if (pos >= s.length())  throw sam::bad_format("Invalid escape");
	  value += s[pos++];
	}
	else
	  break;
      }
    }
    else {
      size_t wspos = s.find_first_of(whitespace, pos);
      dict[key].assign(s, pos, wspos - pos);
      pos = wspos;
    }
  }
}

int insert(sam::interval_multimap<feature>& filters,
	   sam::interval_multimap<feature>& transposons,
	   sam::interval_multimap<feature>& ignores,
	   std::istream& in, feature_action fixed_action) {
  feature_action action = fixed_action;
  int count = 0;

  string line, rname, name;
  seqinterval si;
  for (int lineno = 1; getline(in, line); lineno++)
    if (line.empty() || line[0] == '#' || line.compare(0, 7, "browser") == 0) {
      // Ignore blank lines, comments, and unused BED meta-data
    }
    else if (line.compare(0, 5, "track") == 0) {
      dictionary dict;
      insert_metadata(dict, line);

      dictionary::iterator it = dict.find("brass_action");
      if (it == dict.end())  it = dict.find("action");

      if (it == dict.end())  action = none;
      else if (it->second == "ignore")  action = ignore_reads;
      else if (it->second == "filter:reads")  action = filter_reads;
      else if (it->second == "annotate:transposon")  action = transposon;
      else
	throw sam::bad_format("Invalid action ('" + it->second + "')", lineno);

      if (fixed_action != none && action != fixed_action)
	throw sam::bad_format("Overridden action", lineno);
    }
    else if (action != none) {
      if (line.find('\t') != string::npos) {
	std::istringstream bedline(line);
	coord_t zstart, end;
	if (! (bedline >> rname >> zstart >> end))
	  throw sam::bad_format("too few BED fields", lineno);
	if (! (bedline >> name))  name.clear();

	si.assign(rname, zstart, end);
      }
      else {
	si.assign(line);
	name.clear();
      }

      count++;
      switch (action) {
      case ignore_reads:  ignores.insert(std::make_pair(si, string()));  break;
      case filter_reads:  filters.insert(std::make_pair(si, name));      break;
      case transposon:    transposons.insert(std::make_pair(si, name));  break;
      case none:          break;
      }
    }

  return count;
}

int insert(sam::interval_multimap<feature>& filters,
	   sam::interval_multimap<feature>& transposons,
	   sam::interval_multimap<feature>& ignores,
	   const std::string& filename, feature_action fixed_action)
try {
  std::ios::openmode mode = std::ios::in;
  if (filename.length() >= 3 &&
      filename.compare(filename.length() - 3, 3, ".gz") == 0)
    mode |= std::ios::binary;

  std::ifstream file(filename.c_str(), mode);
  if (! file.is_open())
    throw sam::system_error("can't open ", filename, errno);

  if (mode & std::ios::binary) {
#ifndef WITHOUT_BOOST
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);
    return insert(filters, transposons, ignores, in, fixed_action);
#else
    throw std::logic_error("decompression not available (rebuild with BOOST)");
#endif
  }
  else
    return insert(filters, transposons, ignores, file, fixed_action);
}
catch (sam::exception& e) {
  e.set_filename(filename);
  throw;
}
