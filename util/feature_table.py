#!/usr/bin/env python

import re
import os

import util.feature_table_types

class SeqPosition(object):
    def __init__(self, position, location_operator=None, allow_fuzzy=True):
        self.position = int(position)
        self.location_operator = location_operator
        self.allow_fuzzy = allow_fuzzy

    def __str__(self):
        if self.allow_fuzzy:
            return "{op}{pos}".format(op=self.location_operator if self.location_operator else "", pos=str(self.position))
        else:
            return str(self.position)

    def __int__(self):
        return self.position

    def is_fuzzy(self):
        return self.locaton_operator

    def __eq__(self, other):
        return self.position == other.position

    def __ne__(self, other):
        return self.position != other.position

    def __lt__(self, other):
        return self.position < other.position

    def __le__(self, other):
        return self.position <= other.position

    def __gt__(self, other):
        return self.position > other.position

    def __ge__(self, other):
        return self.position >= other.position

class SeqQualifier(object):
    def __init__(self, qualifier_key, qualifier_value):
        self.qualifier_key   = qualifier_key
        self.qualifier_value = qualifier_value

    def __str__(self):
        if self.qualifier_value:
            return "\t\t\t{k}\t{v}".format(k=self.qualifier_key,v=self.qualifier_value)
        else:
            return "\t\t\t{k}".format(k=self.qualifier_key)

class SeqLocation(object):
    def __init__(self, start_pos, end_pos, feature_type=None):
        self.start = start_pos
        self.end = end_pos
        self.feature_type = feature_type

    def __str__(self):
        return "{start}\t{end}\t{type}".format(start=self.start,end=self.end,type=self.feature_type if self.feature_type else "")

    def __eq__(self, other):
        return ((self.start, self.end) == (other.start, other.end))

    def __ne__(self, other):
        return ((self.start, self.end) != (other.start, other.end))

    def __lt__(self, other):
        return self.start < other.start or (self.start==other.start and self.end < other.end)

    def __le__(self, other):
        return self.start <= other.start or (self.start==other.start and self.end <= other.end)

    def __gt__(self, other):
        return self.start > other.start or (self.start==other.start and self.end > other.end)

    def __ge__(self, other):
        return self.start >= other.start or (self.start==other.start and self.end >= other.end)

class SeqFeature(object):
    def __init__(self, locations=None, feature_type=None):
        self.locations = locations or []
        self.type = feature_type
        self.qualifiers = []

    def add_location(self, location):
        self.locations.append(location)

    def add_location(self, start, location_operator_start, end, location_operator_end, allow_fuzzy=True, *args, **kwargs):
        self.locations.append(SeqLocation(SeqPosition(start, location_operator_start, allow_fuzzy=allow_fuzzy), SeqPosition(end, location_operator_end, allow_fuzzy=allow_fuzzy)))
        # sort with each insertion
        self.sort()

    def sort(self):
        #pass
        self.locations = sorted(self.locations)

    def add_qualifier(self, qualifier_key, qualifier_value, *args, **kwargs):
        self.qualifiers.append(SeqQualifier(qualifier_key, qualifier_value))

    def add_note(self, note_text):
        self.add_qualifier("note", note_text)

    @property
    def lines(self):
        if len(self.locations):
            first_loc = self.locations.pop(0)
            first_loc.feature_type = self.type
            yield first_loc
            for l in self.locations:
                yield l
            # yield notes first
            for q in self.qualifiers:
                if q.qualifier_key == "note":
                    yield q
            for q in self.qualifiers:
                if q.qualifier_key != "note":
                    yield q

class AttrDict(dict):
    """
        Expose (string) dictionary keys as attributes.
        Danger: can cause collisions between dict item keys
                and dict built-ins (the former will override).
                Only to be used with trusted input for dict keys
        Note: memory leak Python < 2.7.4 / Python3 < 3.2.3
              see: https://bugs.python.org/issue1469629
    """
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

class FeatureTable(object):
    # get feature table keys from extensive map
    _default_feature_types = sorted(util.feature_table_types.ft_types.keys())

    def __init__(self, filepath=None, valid_feature_types=None):
        self.refID = None
        self._features = []
        self.valid_feature_types = valid_feature_types or self._default_feature_types

        self.feature_line_regex_map = {        
            "feature_table_header"             : re.compile(r"^>Feature (gb\||ref\|)(?P<refID>.*)\|.*$"),
            "feature_first_location_line"      : re.compile(r"^(?P<location_operator_start>[<>])?(?P<start>\d+)\t(?P<location_operator_end>[<>])?(?P<end>\d+)\t(?P<feature_type>" + "|".join(self.valid_feature_types) + ")$"),
            "feature_subsequent_location_line" : re.compile(r"^(?P<location_operator_start>[<>])?(?P<start>\d+)\t(?P<location_operator_end>[<>])?(?P<end>\d+)\t*$"),
            "offset_line"                      : re.compile(r"^(?:\[offset=(?P<offset>-?\d+)\])$"),
            "feature_qualifier_line"           : re.compile(r"^\t{3}(?P<qualifier_key>[^\t]*)(?:\t(?P<qualifier_value>[^\t]*))?$")
        }

        if filepath:
            self.read_feature_table(filepath)

    def _parse_line(self, line):
        for k,r in self.feature_line_regex_map.items():
            m = r.match(line)
            if m:
                return_dict = AttrDict(m.groupdict())
                return_dict["line_type"] = k
                #return_dict["raw_line"] = line
                return return_dict
        raise LookupError(r"Error parsing feature table line: '%s'" % line)

    @property
    def features(self):
        return self._features

    @property
    def default_feature_types(self):
        return self._default_feature_values.keys()
    

    def add_feature(self, feature):
        self._features.append(feature)

    def read_feature_table(self, filepath, map_function=None, allow_fuzzy=True):
        """
            Read a Genbank-style feature table in tsv/Sequin format
            This applies any offsets present in the input to the positions seen (no offsets in output).
            The map_function accepts two ints: a position interval, and returns two ints
            representing remapped feature coordinates. 
        """
        with open(os.path.realpath(os.path.normpath(filepath)), 'rt') as inf:
            feature_in_progress = None
            offset = 0
            for line in inf:
                line = line.rstrip('\r\n')
                if not line:
                    continue

                l = self._parse_line(line)
                if not self.refID and l.line_type == "feature_table_header":
                    self.refID = l.refID
                elif l.line_type == "feature_first_location_line":
                    if feature_in_progress:
                        self.features.append(feature_in_progress)
                        feature_in_progress = None
                    
                    feature_in_progress = SeqFeature(feature_type=l.feature_type)
                    l.start = int(l.start) + int(offset)
                    l.end = int(l.end) + int(offset)

                    if map_function:
                        l.start, l.end = map_function(l.start, l.end)

                    feature_in_progress.add_location(allow_fuzzy=allow_fuzzy,**l)
                elif l.line_type == "feature_subsequent_location_line":
                    l.start = int(l.start) + int(offset)
                    l.end = int(l.end) + int(offset)

                    if map_function:
                        l.start, l.end = map_function(l.start, l.end)

                    feature_in_progress.add_location(allow_fuzzy=allow_fuzzy,**l)
                elif l.line_type == "feature_qualifier_line":
                    feature_in_progress.add_qualifier(**l)
                elif l.line_type == "offset_line":
                    offset = int(l.offset)

            if feature_in_progress:
                self.features.append(feature_in_progress)
                feature_in_progress = None

        if not self.refID:
            raise Exception("The feature table file did not have a header line in the format '%s'" % self.feature_line_regex_map["feature_table_header"].pattern)

    def remap_locations(self, map_function=None):
        """
            map_function should accept two SeqLocation objects and the containing feature, and return two (mapped) SeqLocation objects
                         if (None,None) is returned, the location will not be included in the remapped location list
        """
        if map_function:
            for feature in self.features:
                remapped_locations = []
                for location in feature.locations:
                    
                    location.start, location.end = map_function(location.start,location.end, feature)
                    # only include locations that are not totally null
                    if location.start is not None and location.end is not None:
                        remapped_locations.append(location)
                    else:
                        feature.add_note("sequencing did not capture all intervals comprising CDS")
                feature.locations = remapped_locations


    def lines(self, exclude_patterns=None):
        yield ">Feature {refID}".format(refID=self.refID)
        for feature in self.features:
            for l in feature.lines:
                if not exclude_patterns:
                    yield l
                else:                    
                    if not any(re.search(pattern, str(l)) for pattern in exclude_patterns):
                       yield l
