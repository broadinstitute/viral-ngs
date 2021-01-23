#!/usr/bin/env python

"""
    Related only to sequence data within this file: 
        "Oligonucleotide sequences Copyright 2016 Illumina, Inc. All rights reserved."
        See: https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pdf
    Python code below is under the license of this repository:
        https://raw.githubusercontent.com/broadinstitute/viral-ngs/master/LICENSE
"""

import re, functools
import csv
import copy
import math
import logging
from collections import OrderedDict, defaultdict

import util.file

log = logging.getLogger(__name__)

def memoize(obj):
    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = "".join([str(args),str(kwargs)])
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer

class IlluminaIndexReference(object):

    def __init__(self, kit=None, instrument=None):
        self.kit = kit
        self.instrument = instrument

    _kits = {
        "TruSight Amplicon Panels":{
            "_short_name":"trusight_amp",
            "i7":{
                "A701":[{"seq":"ATCACGAC","instruments":[]}],
                "A702":[{"seq":"ACAGTGGT","instruments":[]}],
                "A703":[{"seq":"CAGATCCA","instruments":[]}],
                "A704":[{"seq":"ACAAACGG","instruments":[]}],
                "A705":[{"seq":"ACCCAGCA","instruments":[]}],
                "A706":[{"seq":"AACCCCTC","instruments":[]}],
                "A707":[{"seq":"CCCAACCT","instruments":[]}],
                "A708":[{"seq":"CACCACAC","instruments":[]}],
                "A709":[{"seq":"GAAACCCA","instruments":[]}],
                "A710":[{"seq":"TGTGACCA","instruments":[]}],
                "A711":[{"seq":"AGGGTCAA","instruments":[]}],
                "A712":[{"seq":"AGGAGTGG","instruments":[]}]
            },
            "i5":{
                "A501": [{"seq":"TGAACCTT", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"AAGGTTCA", "instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A502": [{"seq":"TGCTAAGT", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"ACTTAGCA", "instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A503": [{"seq":"TGTTCTCT", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"AGAGAACA", "instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A504": [{"seq":"TAAGACAC", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"GTGTCTTA", "instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A505": [{"seq":"CTAATCGA", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TCGATTAG", "instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A506": [{"seq":"CTAGAACA", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TGTTCTAG", "instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A507": [{"seq":"TAAGTTCC", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"GGAACTTA", "instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A508": [{"seq":"TAGACCTA", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TAGGTCTA", "instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}]
            }
        },
        "TruSight Cardio":{
            "_short_name":"trusight_cardio",
            "i7":{
                "N701": [{"seq":"TAAGGCGA","instruments":[]}],
                "N702": [{"seq":"CGTACTAG","instruments":[]}],
                "N703": [{"seq":"AGGCAGAA","instruments":[]}],
                "N704": [{"seq":"TCCTGAGC","instruments":[]}],
                "N705": [{"seq":"GGACTCCT","instruments":[]}],
                "N706": [{"seq":"TAGGCATG","instruments":[]}],
                "N707": [{"seq":"CTCTCTAC","instruments":[]}],
                "N708": [{"seq":"CAGAGAGG","instruments":[]}],
                "N709": [{"seq":"GCTACGCT","instruments":[]}],
                "N710": [{"seq":"CGAGGCTG","instruments":[]}],
                "N711": [{"seq":"AAGAGGCA","instruments":[]}],
                "N712": [{"seq":"GTAGAGGA","instruments":[]}]
            }    ,
            "i5":{
                "E505": [{"seq":"GTAAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"CTCCTTAC","instruments":["NextSeq","HiSeq 3000","HiSeq 4000"]}]
            }
        },
        "TruSight One":{
            "_short_name":"trusight_one",
            "i7":{
                "N701":[{"seq":"TAAGGCGA","instruments":[]}],
                "N702":[{"seq":"CGTACTAG","instruments":[]}],
                "N703":[{"seq":"AGGCAGAA","instruments":[]}],
                "N704":[{"seq":"TCCTGAGC","instruments":[]}],
                "N705":[{"seq":"GGACTCCT","instruments":[]}],
                "N706":[{"seq":"TAGGCATG","instruments":[]}],
                "N707":[{"seq":"CTCTCTAC","instruments":[]}],
                "N708":[{"seq":"CAGAGAGG","instruments":[]}],
                "N709":[{"seq":"GCTACGCT","instruments":[]}],
                "N710":[{"seq":"CGAGGCTG","instruments":[]}],
                "N711":[{"seq":"AAGAGGCA","instruments":[]}],
                "N712":[{"seq":"GTAGAGGA","instruments":[]}]
            },
            "i5":{
                "E502": [{"seq":"CTCTCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"ATAGAGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "E503": [{"seq":"TATCCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"AGAGGATA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "E504": [{"seq":"AGAGTAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TCTACTCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "E505": [{"seq":"GTAAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"CTCCTTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}]
            }
        },
        "TruSight Rapid Capture":{
            "_short_name":"trusight_rapid_capture",
            "i7":{
                "N701":[{"seq":"TAAGGCGA","instruments":[]}],
                "N702":[{"seq":"CGTACTAG","instruments":[]}],
                "N703":[{"seq":"AGGCAGAA","instruments":[]}],
                "N704":[{"seq":"TCCTGAGC","instruments":[]}],
                "N705":[{"seq":"GGACTCCT","instruments":[]}],
                "N706":[{"seq":"TAGGCATG","instruments":[]}],
                "N707":[{"seq":"CTCTCTAC","instruments":[]}],
                "N708":[{"seq":"CAGAGAGG","instruments":[]}],
                "N709":[{"seq":"GCTACGCT","instruments":[]}],
                "N710":[{"seq":"CGAGGCTG","instruments":[]}],
                "N711":[{"seq":"AAGAGGCA","instruments":[]}],
                "N712":[{"seq":"GTAGAGGA","instruments":[]}]
            },
            "i5":{
                "E502": [{"seq":"CTCTCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"ATAGAGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "E505": [{"seq":"GTAAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"CTCCTTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "E506": [{"seq":"ACTGCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TATGCAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "E517": [{"seq":"GCGTAAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TCTTACGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}]
            }
        },
        "TruSight Tumor 15":{
            "_short_name":"trusight_tumor_15",
            "i7":{
                "R701": [{"seq":"ATCACG","instruments":[]}],
                "R702": [{"seq":"CGATGT","instruments":[]}],
                "R703": [{"seq":"TTAGGC","instruments":[]}],
                "R704": [{"seq":"TGACCA","instruments":[]}],
                "R705": [{"seq":"ACAGTG","instruments":[]}],
                "R706": [{"seq":"GCCAAT","instruments":[]}],
                "R707": [{"seq":"CAGATC","instruments":[]}],
                "R708": [{"seq":"ACTTGA","instruments":[]}],
                "R709": [{"seq":"GATCAG","instruments":[]}],
                "R711": [{"seq":"GGCTAC","instruments":[]}],
                "R712": [{"seq":"CTTGTA","instruments":[]}],
                "R725": [{"seq":"ACTGAT","instruments":[]}],
                "R726": [{"seq":"ATGAGC","instruments":[]}],
                "R727": [{"seq":"ATTCCT","instruments":[]}],
                "R728": [{"seq":"CAAAAG","instruments":[]}],
                "R729": [{"seq":"CAACTA","instruments":[]}],
                "R730": [{"seq":"CACCGG","instruments":[]}],
                "R731": [{"seq":"CACGAT","instruments":[]}],
                "R732": [{"seq":"CACTCA","instruments":[]}],
                "R733": [{"seq":"CAGGCG","instruments":[]}],
                "R734": [{"seq":"CATGGC","instruments":[]}],
                "R735": [{"seq":"CATTTT","instruments":[]}],
                "R736": [{"seq":"CCAACA","instruments":[]}],
                "R749": [{"seq":"GATGCT","instruments":[]}]
            },
            "i5":{
                "A501": [{"seq":"TGAACCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"AAGGTTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A502": [{"seq":"TGCTAAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"ACTTAGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}]
            }
        },
        "TruSight RNA Pan-Cancer Panel":{
            "_short_name":"trusight_rna_pan_cancer",
            "indices":{
                "AR001": [{"seq":"ATCACG", "instruments":[], "LT_set":"B"}],
                "AR002": [{"seq":"CGATGT", "instruments":[], "LT_set":"A"}],
                "AR003": [{"seq":"TTAGGC", "instruments":[], "LT_set":"B"}],
                "AR004": [{"seq":"TGACCA", "instruments":[], "LT_set":"A"}],
                "AR005": [{"seq":"ACAGTG", "instruments":[], "LT_set":"A"}],
                "AR006": [{"seq":"GCCAAT", "instruments":[], "LT_set":"A"}],
                "AR007": [{"seq":"CAGATC", "instruments":[], "LT_set":"A"}],
                "AR008": [{"seq":"ACTTGA", "instruments":[], "LT_set":"B"}],
                "AR009": [{"seq":"GATCAG", "instruments":[], "LT_set":"B"}],
                "AR010": [{"seq":"TAGCTT", "instruments":[], "LT_set":"B"}],
                "AR011": [{"seq":"GGCTAC", "instruments":[], "LT_set":"B"}],
                "AR012": [{"seq":"CTTGTA", "instruments":[], "LT_set":"A"}],
                "AR013": [{"seq":"AGTCAA", "instruments":[], "LT_set":"A"}],
                "AR014": [{"seq":"AGTTCC", "instruments":[], "LT_set":"A"}],
                "AR015": [{"seq":"ATGTCA", "instruments":[], "LT_set":"A"}],
                "AR016": [{"seq":"CCGTCC", "instruments":[], "LT_set":"A"}],
                "AR018": [{"seq":"GTCCGC", "instruments":[], "LT_set":"A"}],
                "AR019": [{"seq":"GTGAAA", "instruments":[], "LT_set":"A"}],
                "AR020": [{"seq":"GTGGCC", "instruments":[], "LT_set":"B"}],
                "AR021": [{"seq":"GTTTCG", "instruments":[], "LT_set":"B"}],
                "AR022": [{"seq":"CGTACG", "instruments":[], "LT_set":"B"}],
                "AR023": [{"seq":"GAGTGG", "instruments":[], "LT_set":"B"}],
                "AR025": [{"seq":"ACTGAT", "instruments":[], "LT_set":"B"}],
                "AR027": [{"seq":"ATTCCT", "instruments":[], "LT_set":"B"}]
            },
        },
        "Nextera Index Kit":{
            "_short_name":"nextera",
            "i7":{
                "N701": [{"seq":"TAAGGCGA","instruments":[],"seq_in_adapter":"TCGCCTTA"}],
                "N702": [{"seq":"CGTACTAG","instruments":[],"seq_in_adapter":"CTAGTACG"}],
                "N703": [{"seq":"AGGCAGAA","instruments":[],"seq_in_adapter":"TTCTGCCT"}],
                "N704": [{"seq":"TCCTGAGC","instruments":[],"seq_in_adapter":"GCTCAGGA"}],
                "N705": [{"seq":"GGACTCCT","instruments":[],"seq_in_adapter":"AGGAGTCC"}],
                "N706": [{"seq":"TAGGCATG","instruments":[],"seq_in_adapter":"CATGCCTA"}],
                "N707": [{"seq":"CTCTCTAC","instruments":[],"seq_in_adapter":"GTAGAGAG"}],
                "N708": [{"seq":"CAGAGAGG","instruments":[],"seq_in_adapter":"CCTCTCTG"}],
                "N709": [{"seq":"GCTACGCT","instruments":[],"seq_in_adapter":"AGCGTAGC"}],
                "N710": [{"seq":"CGAGGCTG","instruments":[],"seq_in_adapter":"AAGAGGCA"}],
                "N711": [{"seq":"AAGAGGCA","instruments":[],"seq_in_adapter":"TGCCTCTT"}],
                "N712": [{"seq":"GTAGAGGA","instruments":[],"seq_in_adapter":"TCCTCTAC"}]
            },
            "i5":{
                "[N|S|E]501":[{"seq":"TAGATCGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGATCGC"}, {"seq":"GCGATCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGATCGC"}],
                "[N|S|E]502":[{"seq":"CTCTCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCTCTAT"}, {"seq":"ATAGAGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCTCTAT"}],
                "[N|S|E]503":[{"seq":"TATCCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TATCCTCT"}, {"seq":"AGAGGATA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TATCCTCT"}],
                "[N|S|E]504":[{"seq":"AGAGTAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGAGTAGA"}, {"seq":"TCTACTCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGAGTAGA"}],
                "[N|S|E]505":[{"seq":"GTAAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTAAGGAG"}, {"seq":"CTCCTTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTAAGGAG"}],
                "[N|S|E]506":[{"seq":"ACTGCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTGCATA"}, {"seq":"TATGCAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTGCATA"}],
                "[N|S|E]507":[{"seq":"AAGGAGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGGAGTA"}, {"seq":"TACTCCTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGAGTA"}],
                "[N|S|E]508":[{"seq":"CTAAGCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTAAGCCT"}, {"seq":"AGGCTTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTAAGCCT"}],
                "[N|S|E]517":[{"seq":"GCGTAAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCGTAAGA"}, {"seq":"TCTTACGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCGTAAGA"}],
            }
        },
        "Nextera XT Index Kit v2":{
            "_short_name":"nextera_v2",
            "i7":{
                "N701": [{"seq":"TAAGGCGA","instruments":[],"seq_in_adapter":"TCGCCTTA"}],
                "N702": [{"seq":"CGTACTAG","instruments":[],"seq_in_adapter":"CTAGTACG"}],
                "N703": [{"seq":"AGGCAGAA","instruments":[],"seq_in_adapter":"TTCTGCCT"}],
                "N704": [{"seq":"TCCTGAGC","instruments":[],"seq_in_adapter":"GCTCAGGA"}],
                "N705": [{"seq":"GGACTCCT","instruments":[],"seq_in_adapter":"AGGAGTCC"}],
                "N706": [{"seq":"TAGGCATG","instruments":[],"seq_in_adapter":"CATGCCTA"}],
                "N707": [{"seq":"CTCTCTAC","instruments":[],"seq_in_adapter":"GTAGAGAG"}],
                "N710": [{"seq":"CGAGGCTG","instruments":[],"seq_in_adapter":"CAGCCTCG"}],
                "N711": [{"seq":"AAGAGGCA","instruments":[],"seq_in_adapter":"TGCCTCTT"}],
                "N712": [{"seq":"GTAGAGGA","instruments":[],"seq_in_adapter":"TCCTCTAC"}],
                "N714": [{"seq":"GCTCATGA","instruments":[],"seq_in_adapter":"TCATGAGC"}],
                "N715": [{"seq":"ATCTCAGG","instruments":[],"seq_in_adapter":"CCTGAGAT"}],
                "N716": [{"seq":"ACTCGCTA","instruments":[],"seq_in_adapter":"TAGCGAGT"}],
                "N718": [{"seq":"GGAGCTAC","instruments":[],"seq_in_adapter":"GTAGCTCC"}],
                "N719": [{"seq":"GCGTAGTA","instruments":[],"seq_in_adapter":"TACTACGC"}],
                "N720": [{"seq":"CGGAGCCT","instruments":[],"seq_in_adapter":"AGGCTCCG"}],
                "N721": [{"seq":"TACGCTGC","instruments":[],"seq_in_adapter":"GCAGCGTA"}],
                "N722": [{"seq":"ATGCGCAG","instruments":[],"seq_in_adapter":"CTGCGCAT"}],
                "N723": [{"seq":"TAGCGCTC","instruments":[],"seq_in_adapter":"GAGCGCTA"}],
                "N724": [{"seq":"ACTGAGCG","instruments":[],"seq_in_adapter":"CGCTCAGT"}],
                "N726": [{"seq":"CCTAAGAC","instruments":[],"seq_in_adapter":"GTCTTAGG"}],
                "N727": [{"seq":"CGATCAGT","instruments":[],"seq_in_adapter":"ACTGATCG"}],
                "N728": [{"seq":"TGCAGCTA","instruments":[],"seq_in_adapter":"TAGCTGCA"}],
                "N729": [{"seq":"TCGACGTC","instruments":[],"seq_in_adapter":"GACGTCGA"}]
            },
            "i5":{
                "S502":[{"seq":"CTCTCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCTCTAT"}, {"seq":"ATAGAGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCTCTAT"}],
                "S503":[{"seq":"TATCCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TATCCTCT"}, {"seq":"AGAGGATA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TATCCTCT"}],
                "S505":[{"seq":"GTAAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTAAGGAG"}, {"seq":"CTCCTTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTAAGGAG"}],
                "S506":[{"seq":"ACTGCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTGCATA"}, {"seq":"TATGCAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTGCATA"}],
                "S507":[{"seq":"AAGGAGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGGAGTA"}, {"seq":"TACTCCTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGAGTA"}],
                "S508":[{"seq":"CTAAGCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTAAGCCT"}, {"seq":"AGGCTTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTAAGCCT"}],
                "S510":[{"seq":"CGTCTAAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTCTAAT"}, {"seq":"ATTAGACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTCTAAT"}],
                "S511":[{"seq":"TCTCTCCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTCTCCG"}, {"seq":"CGGAGAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTCTCCG"}],
                "S513":[{"seq":"TCGACTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGACTAG"}, {"seq":"CTAGTCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGACTAG"}],
                "S515":[{"seq":"TTCTAGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCTAGCT"}, {"seq":"AGCTAGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCTAGCT"}],
                "S516":[{"seq":"CCTAGAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTAGAGT"}, {"seq":"ACTCTAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTAGAGT"}],
                "S517":[{"seq":"GCGTAAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCGTAAGA"}, {"seq":"TCTTACGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCGTAAGA"}],
                "S518":[{"seq":"CTATTAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTATTAAG"}, {"seq":"CTTAATAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTATTAAG"}],
                "S520":[{"seq":"AAGGCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGGCTAT"}, {"seq":"ATAGCCTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGCTAT"}],
                "S521":[{"seq":"GAGCCTTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGCCTTA"}, {"seq":"TAAGGCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGCCTTA"}],
                "S522":[{"seq":"TTATGCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTATGCGA"}, {"seq":"TCGCATAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTATGCGA"}],
            }
        },
        "TruSeq Amplicon Kits":{
            "_short_name":"truseq_amplicon",
            "i7":{
                "A701":[{"seq":"ATCACGAC","instruments":[]}],
                "A702":[{"seq":"ACAGTGGT","instruments":[]}],
                "A703":[{"seq":"CAGATCCA","instruments":[]}],
                "A704":[{"seq":"ACAAACGG","instruments":[]}],
                "A705":[{"seq":"ACCCAGCA","instruments":[]}],
                "A706":[{"seq":"AACCCCTC","instruments":[]}],
                "A707":[{"seq":"CCCAACCT","instruments":[]}],
                "A708":[{"seq":"CACCACAC","instruments":[]}],
                "A709":[{"seq":"GAAACCCA","instruments":[]}],
                "A710":[{"seq":"TGTGACCA","instruments":[]}],
                "A711":[{"seq":"AGGGTCAA","instruments":[]}],
                "A712":[{"seq":"AGGAGTGG","instruments":[]}],
            },
            "i5":{
                "A501":[{"seq":"TGAACCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"AAGGTTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A502":[{"seq":"TGCTAAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"ACTTAGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A503":[{"seq":"TGTTCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"AGAGAACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A504":[{"seq":"TAAGACAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"GTGTCTTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A505":[{"seq":"CTAATCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TCGATTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A506":[{"seq":"CTAGAACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TGTTCTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A507":[{"seq":"TAAGTTCC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"GGAACTTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A508":[{"seq":"TAGACCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TAGGTCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
            }
        },
        "TruSeq DNA Methylation":{
            "_short_name":"truseq_methylation",
            "indices":{
                "Index 1": [{"seq":"ATCACG","instruments":[]}],
                "Index 2": [{"seq":"CGATGT","instruments":[]}],
                "Index 3": [{"seq":"TTAGGC","instruments":[]}],
                "Index 4": [{"seq":"TGACCA","instruments":[]}],
                "Index 5": [{"seq":"ACAGTG","instruments":[]}],
                "Index 6": [{"seq":"GCCAAT","instruments":[]}],
                "Index 7": [{"seq":"CAGATC","instruments":[]}],
                "Index 8": [{"seq":"ACTTGA","instruments":[]}],
                "Index 9": [{"seq":"GATCAG","instruments":[]}],
                "Index 10": [{"seq":"TAGCTT","instruments":[]}],
                "Index 11": [{"seq":"GGCTAC","instruments":[]}],
                "Index 12": [{"seq":"CTTGTA","instruments":[]}],
            },
        },
        "TruSeq HT Kits":{
            "_short_name":"truseq_ht",
            "i7":{
                "D701":[{"seq":"ATTACTCG","instruments":[]}],
                "D702":[{"seq":"TCCGGAGA","instruments":[]}],
                "D703":[{"seq":"CGCTCATT","instruments":[]}],
                "D704":[{"seq":"GAGATTCC","instruments":[]}],
                "D705":[{"seq":"ATTCAGAA","instruments":[]}],
                "D706":[{"seq":"GAATTCGT","instruments":[]}],
                "D707":[{"seq":"CTGAAGCT","instruments":[]}],
                "D708":[{"seq":"TAATGCGC","instruments":[]}],
                "D709":[{"seq":"CGGCTATG","instruments":[]}],
                "D710":[{"seq":"TCCGCGAA","instruments":[]}],
                "D711":[{"seq":"TCTCGCGC","instruments":[]}],
                "D712":[{"seq":"AGCGATAG","instruments":[]}],
            },
            "i5":{
                "D501":[{"seq":"TATAGCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"AGGCTATA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D502":[{"seq":"ATAGAGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"GCCTCTAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D503":[{"seq":"CCTATCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"AGGATAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D504":[{"seq":"GGCTCTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TCAGAGCC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D505":[{"seq":"AGGCGAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"CTTCGCCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D506":[{"seq":"TAATCTTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TAAGATTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D507":[{"seq":"CAGGACGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"ACGTCCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D508":[{"seq":"GTACTGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"GTCAGTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
            }
        },
        # "TruSeq LT Kits":{
        #     "_short_name":"",
        # },
        # "TruSeq v1/v2 Kits":{
        #     "_short_name":"",
        # },
        "TruSeq Ribo Profile":{
            "_short_name":"truseq_ribo",
            "indices":{
                "A001": [{"seq":"ATCACG","instruments":[]}],
                "A002": [{"seq":"CGATGT","instruments":[]}],
                "A003": [{"seq":"TTAGGC","instruments":[]}],
                "A004": [{"seq":"TGACCA","instruments":[]}],
                "A005": [{"seq":"ACAGTG","instruments":[]}],
                "A006": [{"seq":"GCCAAT","instruments":[]}],
                "A007": [{"seq":"CAGATC","instruments":[]}],
                "A008": [{"seq":"ACTTGA","instruments":[]}],
                "A009": [{"seq":"GATCAG","instruments":[]}],
                "A010": [{"seq":"TAGCTT","instruments":[]}],
                "A011": [{"seq":"GGCTAC","instruments":[]}],
                "A012": [{"seq":"CTTGTA","instruments":[]}],
            },
        },
        # "TruSeq Synthetic Long-Read DNA":{
        #     "_short_name":"",
        # },
        # "TruSeq Small RNA":{
        #     "_short_name":"",
        # },
        "TruSeq Targeted RNA Expression":{
            "_short_name":"truseq_rna_expression",
            "i7":{
                "R701": [{"seq":"ATCACG","instruments":[]}],
                "R702": [{"seq":"CGATGT","instruments":[]}],
                "R703": [{"seq":"TTAGGC","instruments":[]}],
                "R704": [{"seq":"TGACCA","instruments":[]}],
                "R705": [{"seq":"ACAGTG","instruments":[]}],
                "R706": [{"seq":"GCCAAT","instruments":[]}],
                "R707": [{"seq":"CAGATC","instruments":[]}],
                "R708": [{"seq":"ACTTGA","instruments":[]}],
                "R709": [{"seq":"GATCAG","instruments":[]}],
                "R710": [{"seq":"TAGCTT","instruments":[]}],
                "R711": [{"seq":"GGCTAC","instruments":[]}],
                "R712": [{"seq":"CTTGTA","instruments":[]}],
                "R713": [{"seq":"AGTCAA","instruments":[]}],
                "R714": [{"seq":"AGTTCC","instruments":[]}],
                "R715": [{"seq":"ATGTCA","instruments":[]}],
                "R716": [{"seq":"CCGTCC","instruments":[]}],
                "R717": [{"seq":"GTAGAG","instruments":[]}],
                "R718": [{"seq":"GTCCGC","instruments":[]}],
                "R719": [{"seq":"GTGAAA","instruments":[]}],
                "R720": [{"seq":"GTGGCC","instruments":[]}],
                "R721": [{"seq":"GTTTCG","instruments":[]}],
                "R722": [{"seq":"CGTACG","instruments":[]}],
                "R723": [{"seq":"GAGTGG","instruments":[]}],
                "R724": [{"seq":"GGTAGC","instruments":[]}],
                "R725": [{"seq":"ACTGAT","instruments":[]}],
                "R726": [{"seq":"ATGAGC","instruments":[]}],
                "R727": [{"seq":"ATTCCT","instruments":[]}],
                "R728": [{"seq":"CAAAAG","instruments":[]}],
                "R729": [{"seq":"CAACTA","instruments":[]}],
                "R730": [{"seq":"CACCGG","instruments":[]}],
                "R731": [{"seq":"CACGAT","instruments":[]}],
                "R732": [{"seq":"CACTCA","instruments":[]}],
                "R733": [{"seq":"CAGGCG","instruments":[]}],
                "R734": [{"seq":"CATGGC","instruments":[]}],
                "R735": [{"seq":"CATTTT","instruments":[]}],
                "R736": [{"seq":"CCAACA","instruments":[]}],
                "R737": [{"seq":"CGGAAT","instruments":[]}],
                "R738": [{"seq":"CTAGCT","instruments":[]}],
                "R739": [{"seq":"CTATAC","instruments":[]}],
                "R740": [{"seq":"CTCAGA","instruments":[]}],
                "R741": [{"seq":"GACGAC","instruments":[]}],
                "R742": [{"seq":"TAATCG","instruments":[]}],
                "R743": [{"seq":"TACAGC","instruments":[]}],
                "R744": [{"seq":"TATAAT","instruments":[]}],
                "R745": [{"seq":"TCATTC","instruments":[]}],
                "R746": [{"seq":"TCCCGA","instruments":[]}],
                "R747": [{"seq":"TCGAAG","instruments":[]}],
                "R748": [{"seq":"TCGGCA","instruments":[]}],
            },
            "i5":{
                "A501":[{"seq":"TGAACCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"AAGGTTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A502":[{"seq":"TGCTAAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"ACTTAGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A503":[{"seq":"TGTTCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"AGAGAACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A504":[{"seq":"TAAGACAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"GTGTCTTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A505":[{"seq":"CTAATCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TCGATTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A506":[{"seq":"CTAGAACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TGTTCTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A507":[{"seq":"TAAGTTCC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"GGAACTTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A508":[{"seq":"TAGACCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"]}, {"seq":"TAGGTCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
            }
        },
        # "Process Controls for TruSeq Kits":{
        #     "_short_name":"",
        # },
        # "Nextera DNA Sample Prep Kit":{
        #     "_short_name":"",
        # },
        # "Oligonucleotide Sequences for Genomic DNA":{
        #     "_short_name":"",
        # },
        # "Oligonucleotide Sequences for Paired End DNA":{
        #     "_short_name":"",
        # },
        # "Oligonucleotide Sequences for the Multiplexing Sample Prep Oligo Only Kit":{
        #     "_short_name":"",
        # },
        # "Oligonucleotide Sequences for the v1 and v1.5 Small RNA Kits":{
        #     "_short_name":"",
        # }
        "Broad Institute NEXTERAINDEXPL":{
            "_short_name":"BI_nextera",
            "i7":{
                "BI7A01":[{"seq":"AAGTAGAG","instruments":[],"seq_in_adapter":"CTCTACTT"}],
                "BI7A02":[{"seq":"ACACGATC","instruments":[],"seq_in_adapter":"GATCGTGT"}],
                "BI7A03":[{"seq":"TGTTCCGA","instruments":[],"seq_in_adapter":"TCGGAACA"}],
                "BI7A04":[{"seq":"CATGATCG","instruments":[],"seq_in_adapter":"CGATCATG"}],
                "BI7A05":[{"seq":"CGTTACCA","instruments":[],"seq_in_adapter":"TGGTAACG"}],
                "BI7A06":[{"seq":"CAGAGAGG","instruments":[],"seq_in_adapter":"CCTCTCTG"}],
                "BI7A07":[{"seq":"AACGCATT","instruments":[],"seq_in_adapter":"AATGCGTT"}],
                "BI7A08":[{"seq":"ACAGGTAT","instruments":[],"seq_in_adapter":"ATACCTGT"}],
                "BI7A09":[{"seq":"AGGTAAGG","instruments":[],"seq_in_adapter":"CCTTACCT"}],
                "BI7A10":[{"seq":"AACAATGG","instruments":[],"seq_in_adapter":"CCATTGTT"}],
                "BI7A11":[{"seq":"ACTGTATC","instruments":[],"seq_in_adapter":"GATACAGT"}],
                "BI7A12":[{"seq":"AGGTCGCA","instruments":[],"seq_in_adapter":"TGCGACCT"}],
                "BI7B01":[{"seq":"GGTCCAGA","instruments":[],"seq_in_adapter":"TCTGGACC"}],
                "BI7B02":[{"seq":"CATGCTTA","instruments":[],"seq_in_adapter":"TAAGCATG"}],
                "BI7B03":[{"seq":"AGGATCTA","instruments":[],"seq_in_adapter":"TAGATCCT"}],
                "BI7B04":[{"seq":"TCTGGCGA","instruments":[],"seq_in_adapter":"TCGCCAGA"}],
                "BI7B05":[{"seq":"AGGTTATC","instruments":[],"seq_in_adapter":"GATAACCT"}],
                "BI7B06":[{"seq":"GTCTGATG","instruments":[],"seq_in_adapter":"CATCAGAC"}],
                "BI7B07":[{"seq":"CCAACATT","instruments":[],"seq_in_adapter":"AATGTTGG"}],
                "BI7B08":[{"seq":"CAACTCTC","instruments":[],"seq_in_adapter":"GAGAGTTG"}],
                "BI7B09":[{"seq":"ATTCCTCT","instruments":[],"seq_in_adapter":"AGAGGAAT"}],
                "BI7B10":[{"seq":"CTAACTCG","instruments":[],"seq_in_adapter":"CGAGTTAG"}],
                "BI7B11":[{"seq":"CTGCGGAT","instruments":[],"seq_in_adapter":"ATCCGCAG"}],
                "BI7B12":[{"seq":"CTACCAGG","instruments":[],"seq_in_adapter":"CCTGGTAG"}],
                "BI7C01":[{"seq":"GCACATCT","instruments":[],"seq_in_adapter":"AGATGTGC"}],
                "BI7C02":[{"seq":"GTATAACA","instruments":[],"seq_in_adapter":"TGTTATAC"}],
                "BI7C03":[{"seq":"CATAGCGA","instruments":[],"seq_in_adapter":"TCGCTATG"}],
                "BI7C04":[{"seq":"GACAGTAA","instruments":[],"seq_in_adapter":"TTACTGTC"}],
                "BI7C05":[{"seq":"TTACGCAC","instruments":[],"seq_in_adapter":"GTGCGTAA"}],
                "BI7C06":[{"seq":"GTCATCTA","instruments":[],"seq_in_adapter":"TAGATGAC"}],
                "BI7C07":[{"seq":"CTGTAATC","instruments":[],"seq_in_adapter":"GATTACAG"}],
                "BI7C08":[{"seq":"GCCGTCGA","instruments":[],"seq_in_adapter":"TCGACGGC"}],
                "BI7C09":[{"seq":"GTAACATC","instruments":[],"seq_in_adapter":"GATGTTAC"}],
                "BI7C10":[{"seq":"GAAGGAAG","instruments":[],"seq_in_adapter":"CTTCCTTC"}],
                "BI7C11":[{"seq":"GACCTAAC","instruments":[],"seq_in_adapter":"GTTAGGTC"}],
                "BI7C12":[{"seq":"ACCAACTG","instruments":[],"seq_in_adapter":"CAGTTGGT"}],
                "BI7D01":[{"seq":"TTCGCTGA","instruments":[],"seq_in_adapter":"TCAGCGAA"}],
                "BI7D02":[{"seq":"TGCTCGAC","instruments":[],"seq_in_adapter":"GTCGAGCA"}],
                "BI7D03":[{"seq":"GGACTCCT","instruments":[],"seq_in_adapter":"AGGAGTCC"}],
                "BI7D04":[{"seq":"CAGGAGCC","instruments":[],"seq_in_adapter":"GGCTCCTG"}],
                "BI7D05":[{"seq":"CCTTCGCA","instruments":[],"seq_in_adapter":"TGCGAAGG"}],
                "BI7D06":[{"seq":"TTGAATAG","instruments":[],"seq_in_adapter":"CTATTCAA"}],
                "BI7D07":[{"seq":"TATCTGCC","instruments":[],"seq_in_adapter":"GGCAGATA"}],
                "BI7D08":[{"seq":"TAAGCACA","instruments":[],"seq_in_adapter":"TGTGCTTA"}],
                "BI7D09":[{"seq":"TCGCTAGA","instruments":[],"seq_in_adapter":"TCTAGCGA"}],
                "BI7D10":[{"seq":"TGTAATCA","instruments":[],"seq_in_adapter":"TGATTACA"}],
                "BI7D11":[{"seq":"TTAATCAG","instruments":[],"seq_in_adapter":"CTGATTAA"}],
                "BI7D12":[{"seq":"TGCAAGTA","instruments":[],"seq_in_adapter":"TACTTGCA"}],
                "BI7E01":[{"seq":"AGCAATTC","instruments":[],"seq_in_adapter":"GAATTGCT"}],
                "BI7E02":[{"seq":"AACTTGAC","instruments":[],"seq_in_adapter":"GTCAAGTT"}],
                "BI7E03":[{"seq":"TGTCGGAT","instruments":[],"seq_in_adapter":"ATCCGACA"}],
                "BI7E04":[{"seq":"TCGCCTTG","instruments":[],"seq_in_adapter":"CAAGGCGA"}],
                "BI7E05":[{"seq":"AAGACACT","instruments":[],"seq_in_adapter":"AGTGTCTT"}],
                "BI7E06":[{"seq":"TCTCGGTC","instruments":[],"seq_in_adapter":"GACCGAGA"}],
                "BI7E07":[{"seq":"AATGTTCT","instruments":[],"seq_in_adapter":"AGAACATT"}],
                "BI7E08":[{"seq":"ACTAAGAC","instruments":[],"seq_in_adapter":"GTCTTAGT"}],
                "BI7E09":[{"seq":"ATTATCAA","instruments":[],"seq_in_adapter":"TTGATAAT"}],
                "BI7E10":[{"seq":"ACAGTTGA","instruments":[],"seq_in_adapter":"TCAACTGT"}],
                "BI7E11":[{"seq":"AGCATGGA","instruments":[],"seq_in_adapter":"TCCATGCT"}],
                "BI7E12":[{"seq":"AGGTGCGA","instruments":[],"seq_in_adapter":"TCGCACCT"}],
                "BI7F01":[{"seq":"CACATCCT","instruments":[],"seq_in_adapter":"AGGATGTG"}],
                "BI7F02":[{"seq":"AGTTGCTT","instruments":[],"seq_in_adapter":"AAGCAACT"}],
                "BI7F03":[{"seq":"ATAGCGTC","instruments":[],"seq_in_adapter":"GACGCTAT"}],
                "BI7F04":[{"seq":"ATTATGTT","instruments":[],"seq_in_adapter":"AACATAAT"}],
                "BI7F05":[{"seq":"ATTGTCTG","instruments":[],"seq_in_adapter":"CAGACAAT"}],
                "BI7F06":[{"seq":"CAGCAAGG","instruments":[],"seq_in_adapter":"CCTTGCTG"}],
                "BI7F07":[{"seq":"CGCCTTCC","instruments":[],"seq_in_adapter":"GGAAGGCG"}],
                "BI7F08":[{"seq":"CAGCGGTA","instruments":[],"seq_in_adapter":"TACCGCTG"}],
                "BI7F09":[{"seq":"CAATAGTC","instruments":[],"seq_in_adapter":"GACTATTG"}],
                "BI7F10":[{"seq":"CTATGCGT","instruments":[],"seq_in_adapter":"ACGCATAG"}],
                "BI7F11":[{"seq":"CTGTGGCG","instruments":[],"seq_in_adapter":"CGCCACAG"}],
                "BI7F12":[{"seq":"CGCTATGT","instruments":[],"seq_in_adapter":"ACATAGCG"}],
                "BI7G01":[{"seq":"CCAGTTAG","instruments":[],"seq_in_adapter":"CTAACTGG"}],
                "BI7G02":[{"seq":"CTCTCTAC","instruments":[],"seq_in_adapter":"GTAGAGAG"}],
                "BI7G03":[{"seq":"CCTACCAT","instruments":[],"seq_in_adapter":"ATGGTAGG"}],
                "BI7G04":[{"seq":"GAAGAAGT","instruments":[],"seq_in_adapter":"ACTTCTTC"}],
                "BI7G05":[{"seq":"TCCAGCAA","instruments":[],"seq_in_adapter":"TTGCTGGA"}],
                "BI7G06":[{"seq":"AAGAGGCA","instruments":[],"seq_in_adapter":"TGCCTCTT"}],
                "BI7G07":[{"seq":"GACCAGGA","instruments":[],"seq_in_adapter":"TCCTGGTC"}],
                "BI7G08":[{"seq":"GCCTAGCC","instruments":[],"seq_in_adapter":"GGCTAGGC"}],
                "BI7G09":[{"seq":"GTCCACAG","instruments":[],"seq_in_adapter":"CTGTGGAC"}],
                "BI7G10":[{"seq":"GACCGTTG","instruments":[],"seq_in_adapter":"CAACGGTC"}],
                "BI7G11":[{"seq":"GATATCCA","instruments":[],"seq_in_adapter":"TGGATATC"}],
                "BI7G12":[{"seq":"GCCGCAAC","instruments":[],"seq_in_adapter":"GTTGCGGC"}],
                "BI7H01":[{"seq":"AAGGATGT","instruments":[],"seq_in_adapter":"ACATCCTT"}],
                "BI7H02":[{"seq":"AGGCAGAA","instruments":[],"seq_in_adapter":"TTCTGCCT"}],
                "BI7H03":[{"seq":"ATTCTAGG","instruments":[],"seq_in_adapter":"CCTAGAAT"}],
                "BI7H04":[{"seq":"TCCTGAGC","instruments":[],"seq_in_adapter":"GCTCAGGA"}],
                "BI7H05":[{"seq":"TAATGAAC","instruments":[],"seq_in_adapter":"GTTCATTA"}],
                "BI7H06":[{"seq":"CCAGAGCT","instruments":[],"seq_in_adapter":"AGCTCTGG"}],
                "BI7H07":[{"seq":"CGTACTAG","instruments":[],"seq_in_adapter":"CTAGTACG"}],
                "BI7H08":[{"seq":"TATCCAGG","instruments":[],"seq_in_adapter":"CCTGGATA"}],
                "BI7H09":[{"seq":"TCTGCAAG","instruments":[],"seq_in_adapter":"CTTGCAGA"}],
                "BI7H10":[{"seq":"TTGTCTAT","instruments":[],"seq_in_adapter":"ATAGACAA"}],
                "BI7H11":[{"seq":"TTATATCT","instruments":[],"seq_in_adapter":"AGATATAA"}],
                "BI7H12":[{"seq":"TAGGCATG","instruments":[],"seq_in_adapter":"CATGCCTA"}],
            },
            "i5":{
                "BI5A01":[{"seq":"ATCGACTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCGACTG"}, {"seq":"CAGTCGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCGACTG"}],
                "BI5A02":[{"seq":"CGGTTCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGGTTCTT"}, {"seq":"AAGAACCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGGTTCTT"}],
                "BI5A03":[{"seq":"AACCTCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACCTCTT"}, {"seq":"AAGAGGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACCTCTT"}],
                "BI5A04":[{"seq":"CGCATATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGCATATT"}, {"seq":"AATATGCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGCATATT"}],
                "BI5A05":[{"seq":"CTGCTCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGCTCCT"}, {"seq":"AGGAGCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGCTCCT"}],
                "BI5A06":[{"seq":"TAGATCGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGATCGC"}, {"seq":"GCGATCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGATCGC"}],
                "BI5A07":[{"seq":"CCTGTCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTGTCAT"}, {"seq":"ATGACAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTGTCAT"}],
                "BI5A08":[{"seq":"CACTTCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACTTCAT"}, {"seq":"ATGAAGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACTTCAT"}],
                "BI5A09":[{"seq":"AATACCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AATACCAT"}, {"seq":"ATGGTATT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AATACCAT"}],
                "BI5A10":[{"seq":"ATGAATTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATGAATTA"}, {"seq":"TAATTCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATGAATTA"}],
                "BI5A11":[{"seq":"TGCTTCAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCTTCAC"}, {"seq":"GTGAAGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCTTCAC"}],
                "BI5A12":[{"seq":"ATCCTTAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCCTTAA"}, {"seq":"TTAAGGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCCTTAA"}],
                "BI5B01":[{"seq":"GCTAGCAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTAGCAG"}, {"seq":"CTGCTAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTAGCAG"}],
                "BI5B02":[{"seq":"TAGCATTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGCATTG"}, {"seq":"CAATGCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGCATTG"}],
                "BI5B03":[{"seq":"CTACATTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTACATTG"}, {"seq":"CAATGTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTACATTG"}],
                "BI5B04":[{"seq":"TCATTCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCATTCGA"}, {"seq":"TCGAATGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCATTCGA"}],
                "BI5B05":[{"seq":"TTAGCCAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTAGCCAG"}, {"seq":"CTGGCTAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTAGCCAG"}],
                "BI5B06":[{"seq":"GAACTTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAACTTCG"}, {"seq":"CGAAGTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAACTTCG"}],
                "BI5B07":[{"seq":"GACGGTTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACGGTTA"}, {"seq":"TAACCGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACGGTTA"}],
                "BI5B08":[{"seq":"CAAGCTTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAAGCTTA"}, {"seq":"TAAGCTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAAGCTTA"}],
                "BI5B09":[{"seq":"TCAGGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCAGGCTT"}, {"seq":"AAGCCTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCAGGCTT"}],
                "BI5B10":[{"seq":"TACTCCAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACTCCAG"}, {"seq":"CTGGAGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACTCCAG"}],
                "BI5B11":[{"seq":"GCTTCCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTTCCTA"}, {"seq":"TAGGAAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTTCCTA"}],
                "BI5B12":[{"seq":"TTCTTGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCTTGGC"}, {"seq":"GCCAAGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCTTGGC"}],
                "BI5C01":[{"seq":"TACTCTCC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACTCTCC"}, {"seq":"GGAGAGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACTCTCC"}],
                "BI5C02":[{"seq":"AATTCAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AATTCAAC"}, {"seq":"GTTGAATT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AATTCAAC"}],
                "BI5C03":[{"seq":"GCGATTAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCGATTAC"}, {"seq":"GTAATCGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCGATTAC"}],
                "BI5C04":[{"seq":"GTCCAATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCCAATC"}, {"seq":"GATTGGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCCAATC"}],
                "BI5C05":[{"seq":"GCTGATTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTGATTC"}, {"seq":"GAATCAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTGATTC"}],
                "BI5C06":[{"seq":"CTGTATTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGTATTC"}, {"seq":"GAATACAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGTATTC"}],
                "BI5C07":[{"seq":"CTATTAGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTATTAGC"}, {"seq":"GCTAATAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTATTAGC"}],
                "BI5C08":[{"seq":"AGGTACCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGTACCA"}, {"seq":"TGGTACCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGTACCA"}],
                "BI5C09":[{"seq":"GAACGCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAACGCTA"}, {"seq":"TAGCGTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAACGCTA"}],
                "BI5C10":[{"seq":"ATCATACC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCATACC"}, {"seq":"GGTATGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCATACC"}],
                "BI5C11":[{"seq":"GACCATCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACCATCT"}, {"seq":"AGATGGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACCATCT"}],
                "BI5C12":[{"seq":"CATCACTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATCACTT"}, {"seq":"AAGTGATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATCACTT"}],
                "BI5D01":[{"seq":"TGACAGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGACAGCA"}, {"seq":"TGCTGTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGACAGCA"}],
                "BI5D02":[{"seq":"TTCACAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCACAGA"}, {"seq":"TCTGTGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCACAGA"}],
                "BI5D03":[{"seq":"CTCTCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCTCTAT"}, {"seq":"ATAGAGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCTCTAT"}],
                "BI5D04":[{"seq":"CTTGGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTGGCTT"}, {"seq":"AAGCCAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTGGCTT"}],
                "BI5D05":[{"seq":"GAATCGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAATCGAC"}, {"seq":"GTCGATTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAATCGAC"}],
                "BI5D06":[{"seq":"ATATCCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATATCCGA"}, {"seq":"TCGGATAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATATCCGA"}],
                "BI5D07":[{"seq":"TCCAACCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCAACCA"}, {"seq":"TGGTTGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCAACCA"}],
                "BI5D08":[{"seq":"TCCATAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCATAAC"}, {"seq":"GTTATGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCATAAC"}],
                "BI5D09":[{"seq":"CTGACATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGACATC"}, {"seq":"GATGTCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGACATC"}],
                "BI5D10":[{"seq":"CCTCTAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTCTAAC"}, {"seq":"GTTAGAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTCTAAC"}],
                "BI5D11":[{"seq":"CTGGTATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGGTATT"}, {"seq":"AATACCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGGTATT"}],
                "BI5D12":[{"seq":"CGAACTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGAACTTC"}, {"seq":"GAAGTTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGAACTTC"}],
                "BI5E01":[{"seq":"GCAGGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCAGGTTG"}, {"seq":"CAACCTGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCAGGTTG"}],
                "BI5E02":[{"seq":"GCTCTCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTCTCTT"}, {"seq":"AAGAGAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTCTCTT"}],
                "BI5E03":[{"seq":"AATTGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AATTGCTT"}, {"seq":"AAGCAATT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AATTGCTT"}],
                "BI5E04":[{"seq":"CCAACGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAACGCT"}, {"seq":"AGCGTTGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAACGCT"}],
                "BI5E05":[{"seq":"AGTCACCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTCACCT"}, {"seq":"AGGTGACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTCACCT"}],
                "BI5E06":[{"seq":"GCGGACTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCGGACTT"}, {"seq":"AAGTCCGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCGGACTT"}],
                "BI5E07":[{"seq":"CTGGCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGGCTAT"}, {"seq":"ATAGCCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGGCTAT"}],
                "BI5E08":[{"seq":"GTCCTCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCCTCAT"}, {"seq":"ATGAGGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCCTCAT"}],
                "BI5E09":[{"seq":"GCCACCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCCACCAT"}, {"seq":"ATGGTGGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCCACCAT"}],
                "BI5E10":[{"seq":"ATCTTCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCTTCTC"}, {"seq":"GAGAAGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCTTCTC"}],
                "BI5E11":[{"seq":"TTAATCAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTAATCAC"}, {"seq":"GTGATTAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTAATCAC"}],
                "BI5E12":[{"seq":"GACATTAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACATTAA"}, {"seq":"TTAATGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACATTAA"}],
                "BI5F01":[{"seq":"TTCCAGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCCAGCT"}, {"seq":"AGCTGGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCCAGCT"}],
                "BI5F02":[{"seq":"TGACTTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGACTTGG"}, {"seq":"CCAAGTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGACTTGG"}],
                "BI5F03":[{"seq":"TTGGTCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGGTCTG"}, {"seq":"CAGACCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGGTCTG"}],
                "BI5F04":[{"seq":"TCCACTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCACTTC"}, {"seq":"GAAGTGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCACTTC"}],
                "BI5F05":[{"seq":"CACGATTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACGATTC"}, {"seq":"GAATCGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACGATTC"}],
                "BI5F06":[{"seq":"GCGATATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCGATATT"}, {"seq":"AATATCGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCGATATT"}],
                "BI5F07":[{"seq":"CAGCGATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGCGATT"}, {"seq":"AATCGCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGCGATT"}],
                "BI5F08":[{"seq":"AGTACTGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTACTGC"}, {"seq":"GCAGTACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTACTGC"}],
                "BI5F09":[{"seq":"CGACTCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGACTCTC"}, {"seq":"GAGAGTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGACTCTC"}],
                "BI5F10":[{"seq":"CAGCTCAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGCTCAC"}, {"seq":"GTGAGCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGCTCAC"}],
                "BI5F11":[{"seq":"CGCGAATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGCGAATA"}, {"seq":"TATTCGCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGCGAATA"}],
                "BI5F12":[{"seq":"TTCACCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCACCTT"}, {"seq":"AAGGTGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCACCTT"}],
                "BI5G01":[{"seq":"TAGTTAGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGTTAGC"}, {"seq":"GCTAACTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGTTAGC"}],
                "BI5G02":[{"seq":"TATCCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TATCCTCT"}, {"seq":"AGAGGATA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TATCCTCT"}],
                "BI5G03":[{"seq":"CATCCTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATCCTGG"}, {"seq":"CCAGGATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATCCTGG"}],
                "BI5G04":[{"seq":"AATCTCCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AATCTCCA"}, {"seq":"TGGAGATT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AATCTCCA"}],
                "BI5G05":[{"seq":"GCTCCGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTCCGAT"}, {"seq":"ATCGGAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTCCGAT"}],
                "BI5G06":[{"seq":"AGAGTAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGAGTAGA"}, {"seq":"TCTACTCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGAGTAGA"}],
                "BI5G07":[{"seq":"CCATCACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCATCACA"}, {"seq":"TGTGATGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCATCACA"}],
                "BI5G08":[{"seq":"CTTGAATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTGAATC"}, {"seq":"GATTCAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTGAATC"}],
                "BI5G09":[{"seq":"TGCTATTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCTATTA"}, {"seq":"TAATAGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCTATTA"}],
                "BI5G10":[{"seq":"GGTTATCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GGTTATCT"}, {"seq":"AGATAACC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GGTTATCT"}],
                "BI5G11":[{"seq":"GCTCACCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTCACCA"}, {"seq":"TGGTGAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTCACCA"}],
                "BI5G12":[{"seq":"CCAATCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAATCTG"}, {"seq":"CAGATTGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAATCTG"}],
                "BI5H01":[{"seq":"AGCGCTAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCGCTAA"}, {"seq":"TTAGCGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCGCTAA"}],
                "BI5H02":[{"seq":"GTAAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTAAGGAG"}, {"seq":"CTCCTTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTAAGGAG"}],
                "BI5H03":[{"seq":"GGATTAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GGATTAAC"}, {"seq":"GTTAATCC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GGATTAAC"}],
                "BI5H04":[{"seq":"ACTGCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTGCATA"}, {"seq":"TATGCAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTGCATA"}],
                "BI5H05":[{"seq":"GCACAATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCACAATT"}, {"seq":"AATTGTGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCACAATT"}],
                "BI5H06":[{"seq":"CAACTGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAACTGAT"}, {"seq":"ATCAGTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAACTGAT"}],
                "BI5H07":[{"seq":"AAGGAGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGGAGTA"}, {"seq":"TACTCCTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGAGTA"}],
                "BI5H08":[{"seq":"CCAACTAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAACTAA"}, {"seq":"TTAGTTGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAACTAA"}],
                "BI5H09":[{"seq":"CTTCTGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTCTGGC"}, {"seq":"GCCAGAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTCTGGC"}],
                "BI5H10":[{"seq":"TCCGCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCGCATA"}, {"seq":"TATGCGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCGCATA"}],
                "BI5H11":[{"seq":"TCATGTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCATGTCT"}, {"seq":"AGACATGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCATGTCT"}],
                "BI5H12":[{"seq":"CTAAGCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTAAGCCT"}, {"seq":"AGGCTTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTAAGCCT"}],
            }
        },
        "Additional Pardis Sabeti Lab NEXTERAINDEXPL":{
            "_short_name":"PS_nextera",
            "i7":{
                "PS7A01":[{"seq":"CTCTACTT","instruments":[],"seq_in_adapter":"AAGTAGAG"}],
                "PS7A02":[{"seq":"GATCGTGT","instruments":[],"seq_in_adapter":"ACACGATC"}],
                "PS7B01":[{"seq":"TCTGGACC","instruments":[],"seq_in_adapter":"GGTCCAGA"}],
                "PS7B02":[{"seq":"TAAGCATG","instruments":[],"seq_in_adapter":"CATGCTTA"}],
                "PS7C01":[{"seq":"AGATGTGC","instruments":[],"seq_in_adapter":"GCACATCT"}],
                "PS7C02":[{"seq":"TGTTATAC","instruments":[],"seq_in_adapter":"GTATAACA"}],
                "PS7D01":[{"seq":"TCAGCGAA","instruments":[],"seq_in_adapter":"TTCGCTGA"}],
                "PS7D02":[{"seq":"GTCGAGCA","instruments":[],"seq_in_adapter":"TGCTCGAC"}],
                "PS7E01":[{"seq":"GAATTGCT","instruments":[],"seq_in_adapter":"AGCAATTC"}],
                "PS7F01":[{"seq":"AGGATGTG","instruments":[],"seq_in_adapter":"CACATCCT"}],
                "PS7G01":[{"seq":"CTAACTGG","instruments":[],"seq_in_adapter":"CCAGTTAG"}],
                "PS7H01":[{"seq":"ACATCCTT","instruments":[],"seq_in_adapter":"AAGGATGT"}],
            }
        },
        "IDT xGen UDI-UMI":{
            "_short_name":"xGen_UDI-UMI",
            "i7":{
                "X7001":[{"seq":"CTGATCGT","instruments":[],"seq_in_adapter":"ACGATCAG"}],
                "X7002":[{"seq":"ACTCTCGA","instruments":[],"seq_in_adapter":"TCGAGAGT"}],
                "X7003":[{"seq":"TGAGCTAG","instruments":[],"seq_in_adapter":"CTAGCTCA"}],
                "X7004":[{"seq":"GAGACGAT","instruments":[],"seq_in_adapter":"ATCGTCTC"}],
                "X7005":[{"seq":"CTTGTCGA","instruments":[],"seq_in_adapter":"TCGACAAG"}],
                "X7006":[{"seq":"TTCCAAGG","instruments":[],"seq_in_adapter":"CCTTGGAA"}],
                "X7007":[{"seq":"CGCATGAT","instruments":[],"seq_in_adapter":"ATCATGCG"}],
                "X7008":[{"seq":"ACGGAACA","instruments":[],"seq_in_adapter":"TGTTCCGT"}],
                "X7009":[{"seq":"CGGCTAAT","instruments":[],"seq_in_adapter":"ATTAGCCG"}],
                "X7010":[{"seq":"ATCGATCG","instruments":[],"seq_in_adapter":"CGATCGAT"}],
                "X7011":[{"seq":"GCAAGATC","instruments":[],"seq_in_adapter":"GATCTTGC"}],
                "X7012":[{"seq":"GCTATCCT","instruments":[],"seq_in_adapter":"AGGATAGC"}],
                "X7013":[{"seq":"TACGCTAC","instruments":[],"seq_in_adapter":"GTAGCGTA"}],
                "X7014":[{"seq":"TGGACTCT","instruments":[],"seq_in_adapter":"AGAGTCCA"}],
                "X7015":[{"seq":"AGAGTAGC","instruments":[],"seq_in_adapter":"GCTACTCT"}],
                "X7016":[{"seq":"ATCCAGAG","instruments":[],"seq_in_adapter":"CTCTGGAT"}],
                "X7017":[{"seq":"GACGATCT","instruments":[],"seq_in_adapter":"AGATCGTC"}],
                "X7018":[{"seq":"AACTGAGC","instruments":[],"seq_in_adapter":"GCTCAGTT"}],
                "X7019":[{"seq":"CTTAGGAC","instruments":[],"seq_in_adapter":"GTCCTAAG"}],
                "X7020":[{"seq":"GTGCCATA","instruments":[],"seq_in_adapter":"TATGGCAC"}],
                "X7021":[{"seq":"GAATCCGA","instruments":[],"seq_in_adapter":"TCGGATTC"}],
                "X7022":[{"seq":"TCGCTGTT","instruments":[],"seq_in_adapter":"AACAGCGA"}],
                "X7023":[{"seq":"TTCGTTGG","instruments":[],"seq_in_adapter":"CCAACGAA"}],
                "X7024":[{"seq":"AAGCACTG","instruments":[],"seq_in_adapter":"CAGTGCTT"}],
                "X7025":[{"seq":"CCTTGATC","instruments":[],"seq_in_adapter":"GATCAAGG"}],
                "X7026":[{"seq":"GTCGAAGA","instruments":[],"seq_in_adapter":"TCTTCGAC"}],
                "X7027":[{"seq":"ACCACGAT","instruments":[],"seq_in_adapter":"ATCGTGGT"}],
                "X7028":[{"seq":"GATTACCG","instruments":[],"seq_in_adapter":"CGGTAATC"}],
                "X7029":[{"seq":"GCACAACT","instruments":[],"seq_in_adapter":"AGTTGTGC"}],
                "X7030":[{"seq":"GCGTCATT","instruments":[],"seq_in_adapter":"AATGACGC"}],
                "X7031":[{"seq":"ATCCGGTA","instruments":[],"seq_in_adapter":"TACCGGAT"}],
                "X7032":[{"seq":"CGTTGCAA","instruments":[],"seq_in_adapter":"TTGCAACG"}],
                "X7033":[{"seq":"GTGAAGTG","instruments":[],"seq_in_adapter":"CACTTCAC"}],
                "X7034":[{"seq":"CATGGCTA","instruments":[],"seq_in_adapter":"TAGCCATG"}],
                "X7035":[{"seq":"ATGCCTGT","instruments":[],"seq_in_adapter":"ACAGGCAT"}],
                "X7036":[{"seq":"CAACACCT","instruments":[],"seq_in_adapter":"AGGTGTTG"}],
                "X7037":[{"seq":"TGTGACTG","instruments":[],"seq_in_adapter":"CAGTCACA"}],
                "X7038":[{"seq":"GTCATCGA","instruments":[],"seq_in_adapter":"TCGATGAC"}],
                "X7039":[{"seq":"AGCACTTC","instruments":[],"seq_in_adapter":"GAAGTGCT"}],
                "X7040":[{"seq":"GAAGGAAG","instruments":[],"seq_in_adapter":"CTTCCTTC"}],
                "X7041":[{"seq":"GTTGTTCG","instruments":[],"seq_in_adapter":"CGAACAAC"}],
                "X7042":[{"seq":"CGGTTGTT","instruments":[],"seq_in_adapter":"AACAACCG"}],
                "X7043":[{"seq":"ACTGAGGT","instruments":[],"seq_in_adapter":"ACCTCAGT"}],
                "X7044":[{"seq":"TGAAGACG","instruments":[],"seq_in_adapter":"CGTCTTCA"}],
                "X7045":[{"seq":"GTTACGCA","instruments":[],"seq_in_adapter":"TGCGTAAC"}],
                "X7046":[{"seq":"AGCGTGTT","instruments":[],"seq_in_adapter":"AACACGCT"}],
                "X7047":[{"seq":"GATCGAGT","instruments":[],"seq_in_adapter":"ACTCGATC"}],
                "X7048":[{"seq":"ACAGCTCA","instruments":[],"seq_in_adapter":"TGAGCTGT"}],
                "X7049":[{"seq":"GAGCAGTA","instruments":[],"seq_in_adapter":"TACTGCTC"}],
                "X7050":[{"seq":"AGTTCGTC","instruments":[],"seq_in_adapter":"GACGAACT"}],
                "X7051":[{"seq":"TTGCGAAG","instruments":[],"seq_in_adapter":"CTTCGCAA"}],
                "X7052":[{"seq":"ATCGCCAT","instruments":[],"seq_in_adapter":"ATGGCGAT"}],
                "X7053":[{"seq":"TGGCATGT","instruments":[],"seq_in_adapter":"ACATGCCA"}],
                "X7054":[{"seq":"CTGTTGAC","instruments":[],"seq_in_adapter":"GTCAACAG"}],
                "X7055":[{"seq":"CATACCAC","instruments":[],"seq_in_adapter":"GTGGTATG"}],
                "X7056":[{"seq":"GAAGTTGG","instruments":[],"seq_in_adapter":"CCAACTTC"}],
                "X7057":[{"seq":"ATGACGTC","instruments":[],"seq_in_adapter":"GACGTCAT"}],
                "X7058":[{"seq":"TTGGACGT","instruments":[],"seq_in_adapter":"ACGTCCAA"}],
                "X7059":[{"seq":"AGTGGATC","instruments":[],"seq_in_adapter":"GATCCACT"}],
                "X7060":[{"seq":"GATAGGCT","instruments":[],"seq_in_adapter":"AGCCTATC"}],
                "X7061":[{"seq":"TGGTAGCT","instruments":[],"seq_in_adapter":"AGCTACCA"}],
                "X7062":[{"seq":"CGCAATCT","instruments":[],"seq_in_adapter":"AGATTGCG"}],
                "X7063":[{"seq":"GATGTGTG","instruments":[],"seq_in_adapter":"CACACATC"}],
                "X7064":[{"seq":"GATTGCTC","instruments":[],"seq_in_adapter":"GAGCAATC"}],
                "X7065":[{"seq":"CGCTCTAT","instruments":[],"seq_in_adapter":"ATAGAGCG"}],
                "X7066":[{"seq":"TATCGGTC","instruments":[],"seq_in_adapter":"GACCGATA"}],
                "X7067":[{"seq":"AACGTCTG","instruments":[],"seq_in_adapter":"CAGACGTT"}],
                "X7068":[{"seq":"ACGTTCAG","instruments":[],"seq_in_adapter":"CTGAACGT"}],
                "X7069":[{"seq":"CAGTCCAA","instruments":[],"seq_in_adapter":"TTGGACTG"}],
                "X7070":[{"seq":"TTGCAGAC","instruments":[],"seq_in_adapter":"GTCTGCAA"}],
                "X7071":[{"seq":"CAATGTGG","instruments":[],"seq_in_adapter":"CCACATTG"}],
                "X7072":[{"seq":"ACTCCATC","instruments":[],"seq_in_adapter":"GATGGAGT"}],
                "X7073":[{"seq":"GTTGACCT","instruments":[],"seq_in_adapter":"AGGTCAAC"}],
                "X7074":[{"seq":"CGTGTGTA","instruments":[],"seq_in_adapter":"TACACACG"}],
                "X7075":[{"seq":"ACGACTTG","instruments":[],"seq_in_adapter":"CAAGTCGT"}],
                "X7076":[{"seq":"CACTAGCT","instruments":[],"seq_in_adapter":"AGCTAGTG"}],
                "X7077":[{"seq":"ACTAGGAG","instruments":[],"seq_in_adapter":"CTCCTAGT"}],
                "X7078":[{"seq":"GTAGGAGT","instruments":[],"seq_in_adapter":"ACTCCTAC"}],
                "X7079":[{"seq":"CCTGATTG","instruments":[],"seq_in_adapter":"CAATCAGG"}],
                "X7080":[{"seq":"ATGCACGA","instruments":[],"seq_in_adapter":"TCGTGCAT"}],
                "X7081":[{"seq":"CGACGTTA","instruments":[],"seq_in_adapter":"TAACGTCG"}],
                "X7082":[{"seq":"TACGCCTT","instruments":[],"seq_in_adapter":"AAGGCGTA"}],
                "X7083":[{"seq":"CCGTAAGA","instruments":[],"seq_in_adapter":"TCTTACGG"}],
                "X7084":[{"seq":"ATCACACG","instruments":[],"seq_in_adapter":"CGTGTGAT"}],
                "X7085":[{"seq":"CACCTGTT","instruments":[],"seq_in_adapter":"AACAGGTG"}],
                "X7086":[{"seq":"CTTCGACT","instruments":[],"seq_in_adapter":"AGTCGAAG"}],
                "X7087":[{"seq":"TGCTTCCA","instruments":[],"seq_in_adapter":"TGGAAGCA"}],
                "X7088":[{"seq":"AGAACGAG","instruments":[],"seq_in_adapter":"CTCGTTCT"}],
                "X7089":[{"seq":"GTTCTCGT","instruments":[],"seq_in_adapter":"ACGAGAAC"}],
                "X7090":[{"seq":"TCAGGCTT","instruments":[],"seq_in_adapter":"AAGCCTGA"}],
                "X7091":[{"seq":"CCTTGTAG","instruments":[],"seq_in_adapter":"CTACAAGG"}],
                "X7092":[{"seq":"GAACATCG","instruments":[],"seq_in_adapter":"CGATGTTC"}],
                "X7093":[{"seq":"TAACCGGT","instruments":[],"seq_in_adapter":"ACCGGTTA"}],
                "X7094":[{"seq":"AACCGTTC","instruments":[],"seq_in_adapter":"GAACGGTT"}],
                "X7095":[{"seq":"TGGTACAG","instruments":[],"seq_in_adapter":"CTGTACCA"}],
                "X7096":[{"seq":"ATATGCGC","instruments":[],"seq_in_adapter":"GCGCATAT"}],
                "X7097":[{"seq":"GCCTATCA","instruments":[],"seq_in_adapter":"TGATAGGC"}],
                "X7098":[{"seq":"CTTGGATG","instruments":[],"seq_in_adapter":"CATCCAAG"}],
                "X7099":[{"seq":"AGTCTCAC","instruments":[],"seq_in_adapter":"GTGAGACT"}],
                "X7100":[{"seq":"CTCATCAG","instruments":[],"seq_in_adapter":"CTGATGAG"}],
                "X7101":[{"seq":"TGTACCGT","instruments":[],"seq_in_adapter":"ACGGTACA"}],
                "X7102":[{"seq":"AAGTCGAG","instruments":[],"seq_in_adapter":"CTCGACTT"}],
                "X7103":[{"seq":"CACGTTGT","instruments":[],"seq_in_adapter":"ACAACGTG"}],
                "X7104":[{"seq":"TCACAGCA","instruments":[],"seq_in_adapter":"TGCTGTGA"}],
                "X7105":[{"seq":"CTACTTGG","instruments":[],"seq_in_adapter":"CCAAGTAG"}],
                "X7106":[{"seq":"CCTCAGTT","instruments":[],"seq_in_adapter":"AACTGAGG"}],
                "X7107":[{"seq":"TCCTACCT","instruments":[],"seq_in_adapter":"AGGTAGGA"}],
                "X7108":[{"seq":"ATGGCGAA","instruments":[],"seq_in_adapter":"TTCGCCAT"}],
                "X7109":[{"seq":"CTTACCTG","instruments":[],"seq_in_adapter":"CAGGTAAG"}],
                "X7110":[{"seq":"CTCGATAC","instruments":[],"seq_in_adapter":"GTATCGAG"}],
                "X7111":[{"seq":"TCCGTGAA","instruments":[],"seq_in_adapter":"TTCACGGA"}],
                "X7112":[{"seq":"TAGAGCTC","instruments":[],"seq_in_adapter":"GAGCTCTA"}],
                "X7113":[{"seq":"TGACTGAC","instruments":[],"seq_in_adapter":"GTCAGTCA"}],
                "X7114":[{"seq":"TAGACGTG","instruments":[],"seq_in_adapter":"CACGTCTA"}],
                "X7115":[{"seq":"CCGGAATT","instruments":[],"seq_in_adapter":"AATTCCGG"}],
                "X7116":[{"seq":"CTCCTAGA","instruments":[],"seq_in_adapter":"TCTAGGAG"}],
                "X7117":[{"seq":"CAACGGAT","instruments":[],"seq_in_adapter":"ATCCGTTG"}],
                "X7118":[{"seq":"TGGCTATC","instruments":[],"seq_in_adapter":"GATAGCCA"}],
                "X7119":[{"seq":"CGGTCATA","instruments":[],"seq_in_adapter":"TATGACCG"}],
                "X7120":[{"seq":"TCCAATCG","instruments":[],"seq_in_adapter":"CGATTGGA"}],
                "X7121":[{"seq":"GAGCTTGT","instruments":[],"seq_in_adapter":"ACAAGCTC"}],
                "X7122":[{"seq":"GAAGGTTC","instruments":[],"seq_in_adapter":"GAACCTTC"}],
                "X7123":[{"seq":"ATCTCGCT","instruments":[],"seq_in_adapter":"AGCGAGAT"}],
                "X7124":[{"seq":"AGTTACGG","instruments":[],"seq_in_adapter":"CCGTAACT"}],
                "X7125":[{"seq":"GTGTCTGA","instruments":[],"seq_in_adapter":"TCAGACAC"}],
                "X7126":[{"seq":"TGACTTCG","instruments":[],"seq_in_adapter":"CGAAGTCA"}],
                "X7127":[{"seq":"TGGATCAC","instruments":[],"seq_in_adapter":"GTGATCCA"}],
                "X7128":[{"seq":"ACACCAGT","instruments":[],"seq_in_adapter":"ACTGGTGT"}],
                "X7129":[{"seq":"CAGGTTAG","instruments":[],"seq_in_adapter":"CTAACCTG"}],
                "X7130":[{"seq":"AGTTGGCT","instruments":[],"seq_in_adapter":"AGCCAACT"}],
                "X7131":[{"seq":"TCAACTGG","instruments":[],"seq_in_adapter":"CCAGTTGA"}],
                "X7132":[{"seq":"CTGCACTT","instruments":[],"seq_in_adapter":"AAGTGCAG"}],
                "X7133":[{"seq":"ACACGGTT","instruments":[],"seq_in_adapter":"AACCGTGT"}],
                "X7134":[{"seq":"AATACGCG","instruments":[],"seq_in_adapter":"CGCGTATT"}],
                "X7135":[{"seq":"TGCGAACT","instruments":[],"seq_in_adapter":"AGTTCGCA"}],
                "X7136":[{"seq":"GCTGACTA","instruments":[],"seq_in_adapter":"TAGTCAGC"}],
                "X7137":[{"seq":"GTGGTGTT","instruments":[],"seq_in_adapter":"AACACCAC"}],
                "X7138":[{"seq":"GTGCTTAC","instruments":[],"seq_in_adapter":"GTAAGCAC"}],
                "X7139":[{"seq":"TCAAGGAC","instruments":[],"seq_in_adapter":"GTCCTTGA"}],
                "X7140":[{"seq":"TGAACCTG","instruments":[],"seq_in_adapter":"CAGGTTCA"}],
                "X7141":[{"seq":"AGTGTTGG","instruments":[],"seq_in_adapter":"CCAACACT"}],
                "X7142":[{"seq":"GTACTCTC","instruments":[],"seq_in_adapter":"GAGAGTAC"}],
                "X7143":[{"seq":"CCGTATCT","instruments":[],"seq_in_adapter":"AGATACGG"}],
                "X7144":[{"seq":"CGAAGAAC","instruments":[],"seq_in_adapter":"GTTCTTCG"}],
                "X7145":[{"seq":"AGCGGAAT","instruments":[],"seq_in_adapter":"ATTCCGCT"}],
                "X7146":[{"seq":"GTGAGCTT","instruments":[],"seq_in_adapter":"AAGCTCAC"}],
                "X7147":[{"seq":"CGTGATCA","instruments":[],"seq_in_adapter":"TGATCACG"}],
                "X7148":[{"seq":"TCGCATTG","instruments":[],"seq_in_adapter":"CAATGCGA"}],
                "X7149":[{"seq":"TGACGCAT","instruments":[],"seq_in_adapter":"ATGCGTCA"}],
                "X7150":[{"seq":"CCGATGTA","instruments":[],"seq_in_adapter":"TACATCGG"}],
                "X7151":[{"seq":"TTCGCAGT","instruments":[],"seq_in_adapter":"ACTGCGAA"}],
                "X7152":[{"seq":"ACGACAGA","instruments":[],"seq_in_adapter":"TCTGTCGT"}],
                "X7153":[{"seq":"AGCTTGAG","instruments":[],"seq_in_adapter":"CTCAAGCT"}],
                "X7154":[{"seq":"GAGTGGTT","instruments":[],"seq_in_adapter":"AACCACTC"}],
                "X7155":[{"seq":"GCTGTAAG","instruments":[],"seq_in_adapter":"CTTACAGC"}],
                "X7156":[{"seq":"CCAAGACT","instruments":[],"seq_in_adapter":"AGTCTTGG"}],
                "X7157":[{"seq":"ATTGCGTG","instruments":[],"seq_in_adapter":"CACGCAAT"}],
                "X7158":[{"seq":"CTGAAGCT","instruments":[],"seq_in_adapter":"AGCTTCAG"}],
                "X7159":[{"seq":"TAACGAGG","instruments":[],"seq_in_adapter":"CCTCGTTA"}],
                "X7160":[{"seq":"TCGTCTCA","instruments":[],"seq_in_adapter":"TGAGACGA"}],
                "X7161":[{"seq":"TTCCTGTG","instruments":[],"seq_in_adapter":"CACAGGAA"}],
                "X7162":[{"seq":"CGTTGAGT","instruments":[],"seq_in_adapter":"ACTCAACG"}],
                "X7163":[{"seq":"AGTCGCTT","instruments":[],"seq_in_adapter":"AAGCGACT"}],
                "X7164":[{"seq":"TAGGTAGG","instruments":[],"seq_in_adapter":"CCTACCTA"}],
                "X7165":[{"seq":"CAGGAGAT","instruments":[],"seq_in_adapter":"ATCTCCTG"}],
                "X7166":[{"seq":"CATCGTGA","instruments":[],"seq_in_adapter":"TCACGATG"}],
                "X7167":[{"seq":"TGTTGTGG","instruments":[],"seq_in_adapter":"CCACAACA"}],
                "X7168":[{"seq":"ACAGACCT","instruments":[],"seq_in_adapter":"AGGTCTGT"}],
                "X7169":[{"seq":"GTCCTTCT","instruments":[],"seq_in_adapter":"AGAAGGAC"}],
                "X7170":[{"seq":"TGATACGC","instruments":[],"seq_in_adapter":"GCGTATCA"}],
                "X7171":[{"seq":"CTGTGTTG","instruments":[],"seq_in_adapter":"CAACACAG"}],
                "X7172":[{"seq":"AACGTGGA","instruments":[],"seq_in_adapter":"TCCACGTT"}],
                "X7173":[{"seq":"GTTGCGAT","instruments":[],"seq_in_adapter":"ATCGCAAC"}],
                "X7174":[{"seq":"AACGACGT","instruments":[],"seq_in_adapter":"ACGTCGTT"}],
                "X7175":[{"seq":"CGTATTCG","instruments":[],"seq_in_adapter":"CGAATACG"}],
                "X7176":[{"seq":"AGCAAGCA","instruments":[],"seq_in_adapter":"TGCTTGCT"}],
                "X7177":[{"seq":"TGTTCGAG","instruments":[],"seq_in_adapter":"CTCGAACA"}],
                "X7178":[{"seq":"CTCCATGT","instruments":[],"seq_in_adapter":"ACATGGAG"}],
                "X7179":[{"seq":"CGTCTTGT","instruments":[],"seq_in_adapter":"ACAAGACG"}],
                "X7180":[{"seq":"ATAAGGCG","instruments":[],"seq_in_adapter":"CGCCTTAT"}],
                "X7181":[{"seq":"TGTCTGCT","instruments":[],"seq_in_adapter":"AGCAGACA"}],
                "X7182":[{"seq":"CGCTTAAC","instruments":[],"seq_in_adapter":"GTTAAGCG"}],
                "X7183":[{"seq":"GATCCATG","instruments":[],"seq_in_adapter":"CATGGATC"}],
                "X7184":[{"seq":"ACCTCTGT","instruments":[],"seq_in_adapter":"ACAGAGGT"}],
                "X7185":[{"seq":"GCCACTTA","instruments":[],"seq_in_adapter":"TAAGTGGC"}],
                "X7186":[{"seq":"ACCTGACT","instruments":[],"seq_in_adapter":"AGTCAGGT"}],
                "X7187":[{"seq":"GTTAAGGC","instruments":[],"seq_in_adapter":"GCCTTAAC"}],
                "X7188":[{"seq":"ATGCCAAC","instruments":[],"seq_in_adapter":"GTTGGCAT"}],
                "X7189":[{"seq":"AGAGGTTG","instruments":[],"seq_in_adapter":"CAACCTCT"}],
                "X7190":[{"seq":"ACCATCCA","instruments":[],"seq_in_adapter":"TGGATGGT"}],
                "X7191":[{"seq":"GTGGATAG","instruments":[],"seq_in_adapter":"CTATCCAC"}],
                "X7192":[{"seq":"CTGAGATC","instruments":[],"seq_in_adapter":"GATCTCAG"}],
                "X7193":[{"seq":"CTTCGTTC","instruments":[],"seq_in_adapter":"GAACGAAG"}],
                "X7194":[{"seq":"GTCTAGGT","instruments":[],"seq_in_adapter":"ACCTAGAC"}],
                "X7195":[{"seq":"ACGTCGTA","instruments":[],"seq_in_adapter":"TACGACGT"}],
                "X7196":[{"seq":"GAGCTCAA","instruments":[],"seq_in_adapter":"TTGAGCTC"}],
                "X7197":[{"seq":"CGTGTACT","instruments":[],"seq_in_adapter":"AGTACACG"}],
                "X7198":[{"seq":"CACTGACA","instruments":[],"seq_in_adapter":"TGTCAGTG"}],
                "X7199":[{"seq":"TCGTAGTC","instruments":[],"seq_in_adapter":"GACTACGA"}],
                "X7200":[{"seq":"GCACGTAA","instruments":[],"seq_in_adapter":"TTACGTGC"}],
                "X7201":[{"seq":"CAAGCAGT","instruments":[],"seq_in_adapter":"ACTGCTTG"}],
                "X7202":[{"seq":"ACATAGGC","instruments":[],"seq_in_adapter":"GCCTATGT"}],
                "X7203":[{"seq":"TGTGGTAC","instruments":[],"seq_in_adapter":"GTACCACA"}],
                "X7204":[{"seq":"CACCACTA","instruments":[],"seq_in_adapter":"TAGTGGTG"}],
                "X7205":[{"seq":"CTGCGTAT","instruments":[],"seq_in_adapter":"ATACGCAG"}],
                "X7206":[{"seq":"ACGGTCTT","instruments":[],"seq_in_adapter":"AAGACCGT"}],
                "X7207":[{"seq":"GATTGGAG","instruments":[],"seq_in_adapter":"CTCCAATC"}],
                "X7208":[{"seq":"TGTCCAGA","instruments":[],"seq_in_adapter":"TCTGGACA"}],
                "X7209":[{"seq":"CCAGTGTT","instruments":[],"seq_in_adapter":"AACACTGG"}],
                "X7210":[{"seq":"TGCACCAA","instruments":[],"seq_in_adapter":"TTGGTGCA"}],
                "X7211":[{"seq":"TTGACAGG","instruments":[],"seq_in_adapter":"CCTGTCAA"}],
                "X7212":[{"seq":"AGGCATAG","instruments":[],"seq_in_adapter":"CTATGCCT"}],
                "X7213":[{"seq":"TAGCCGAA","instruments":[],"seq_in_adapter":"TTCGGCTA"}],
                "X7214":[{"seq":"TTGTCGGT","instruments":[],"seq_in_adapter":"ACCGACAA"}],
                "X7215":[{"seq":"CATCTACG","instruments":[],"seq_in_adapter":"CGTAGATG"}],
                "X7216":[{"seq":"GCATACAG","instruments":[],"seq_in_adapter":"CTGTATGC"}],
                "X7217":[{"seq":"ACAGCAAC","instruments":[],"seq_in_adapter":"GTTGCTGT"}],
                "X7218":[{"seq":"CTGGTTCT","instruments":[],"seq_in_adapter":"AGAACCAG"}],
                "X7219":[{"seq":"TCGACATC","instruments":[],"seq_in_adapter":"GATGTCGA"}],
                "X7220":[{"seq":"AACCTCCT","instruments":[],"seq_in_adapter":"AGGAGGTT"}],
                "X7221":[{"seq":"CAGCGATT","instruments":[],"seq_in_adapter":"AATCGCTG"}],
                "X7222":[{"seq":"AGGTCACT","instruments":[],"seq_in_adapter":"AGTGACCT"}],
                "X7223":[{"seq":"GCAATTCG","instruments":[],"seq_in_adapter":"CGAATTGC"}],
                "X7224":[{"seq":"GCTTCTTG","instruments":[],"seq_in_adapter":"CAAGAAGC"}],
                "X7225":[{"seq":"AACTGGTG","instruments":[],"seq_in_adapter":"CACCAGTT"}],
                "X7226":[{"seq":"CGGAATAC","instruments":[],"seq_in_adapter":"GTATTCCG"}],
                "X7227":[{"seq":"GCTTCGAA","instruments":[],"seq_in_adapter":"TTCGAAGC"}],
                "X7228":[{"seq":"CAAGGTCT","instruments":[],"seq_in_adapter":"AGACCTTG"}],
                "X7229":[{"seq":"AACCTTGG","instruments":[],"seq_in_adapter":"CCAAGGTT"}],
                "X7230":[{"seq":"CCATACGT","instruments":[],"seq_in_adapter":"ACGTATGG"}],
                "X7231":[{"seq":"TGGTCCTT","instruments":[],"seq_in_adapter":"AAGGACCA"}],
                "X7232":[{"seq":"ACCGCATA","instruments":[],"seq_in_adapter":"TATGCGGT"}],
                "X7233":[{"seq":"CCTTCCTT","instruments":[],"seq_in_adapter":"AAGGAAGG"}],
                "X7234":[{"seq":"TACACGCT","instruments":[],"seq_in_adapter":"AGCGTGTA"}],
                "X7235":[{"seq":"TGCGTAGA","instruments":[],"seq_in_adapter":"TCTACGCA"}],
                "X7236":[{"seq":"AAGAGCCA","instruments":[],"seq_in_adapter":"TGGCTCTT"}],
                "X7237":[{"seq":"ATGGAAGG","instruments":[],"seq_in_adapter":"CCTTCCAT"}],
                "X7238":[{"seq":"GCCAGTAT","instruments":[],"seq_in_adapter":"ATACTGGC"}],
                "X7239":[{"seq":"CGTAGGTT","instruments":[],"seq_in_adapter":"AACCTACG"}],
                "X7240":[{"seq":"CGAGTATG","instruments":[],"seq_in_adapter":"CATACTCG"}],
                "X7241":[{"seq":"CAAGTGCA","instruments":[],"seq_in_adapter":"TGCACTTG"}],
                "X7242":[{"seq":"TCGAGTGA","instruments":[],"seq_in_adapter":"TCACTCGA"}],
                "X7243":[{"seq":"CTACAGTG","instruments":[],"seq_in_adapter":"CACTGTAG"}],
                "X7244":[{"seq":"GATCGTAC","instruments":[],"seq_in_adapter":"GTACGATC"}],
                "X7245":[{"seq":"CTTCACCA","instruments":[],"seq_in_adapter":"TGGTGAAG"}],
                "X7246":[{"seq":"CTCAGCTA","instruments":[],"seq_in_adapter":"TAGCTGAG"}],
                "X7247":[{"seq":"TCTGCTCT","instruments":[],"seq_in_adapter":"AGAGCAGA"}],
                "X7248":[{"seq":"AACCGAAG","instruments":[],"seq_in_adapter":"CTTCGGTT"}],
                "X7249":[{"seq":"GCTGTTGT","instruments":[],"seq_in_adapter":"ACAACAGC"}],
                "X7250":[{"seq":"TTACGGCT","instruments":[],"seq_in_adapter":"AGCCGTAA"}],
                "X7251":[{"seq":"GACAAGAG","instruments":[],"seq_in_adapter":"CTCTTGTC"}],
                "X7252":[{"seq":"AGGATCTG","instruments":[],"seq_in_adapter":"CAGATCCT"}],
                "X7253":[{"seq":"GTAGCATC","instruments":[],"seq_in_adapter":"GATGCTAC"}],
                "X7254":[{"seq":"GTGTTCCT","instruments":[],"seq_in_adapter":"AGGAACAC"}],
                "X7255":[{"seq":"AGGATGGT","instruments":[],"seq_in_adapter":"ACCATCCT"}],
                "X7256":[{"seq":"TCACGTTC","instruments":[],"seq_in_adapter":"GAACGTGA"}],
                "X7257":[{"seq":"GCGTTCTA","instruments":[],"seq_in_adapter":"TAGAACGC"}],
                "X7258":[{"seq":"CTCTGGTT","instruments":[],"seq_in_adapter":"AACCAGAG"}],
                "X7259":[{"seq":"TTAGGTCG","instruments":[],"seq_in_adapter":"CGACCTAA"}],
                "X7260":[{"seq":"TCTGAGAG","instruments":[],"seq_in_adapter":"CTCTCAGA"}],
                "X7261":[{"seq":"TTCAGCCT","instruments":[],"seq_in_adapter":"AGGCTGAA"}],
                "X7262":[{"seq":"TCTCCGAT","instruments":[],"seq_in_adapter":"ATCGGAGA"}],
                "X7263":[{"seq":"CAGGTATC","instruments":[],"seq_in_adapter":"GATACCTG"}],
                "X7264":[{"seq":"AGTCAGGA","instruments":[],"seq_in_adapter":"TCCTGACT"}],
                "X7265":[{"seq":"AAGGCTGA","instruments":[],"seq_in_adapter":"TCAGCCTT"}],
                "X7266":[{"seq":"CGATGCTT","instruments":[],"seq_in_adapter":"AAGCATCG"}],
                "X7267":[{"seq":"GTATTGGC","instruments":[],"seq_in_adapter":"GCCAATAC"}],
                "X7268":[{"seq":"ACTGTGTC","instruments":[],"seq_in_adapter":"GACACAGT"}],
                "X7269":[{"seq":"TGCCTCTT","instruments":[],"seq_in_adapter":"AAGAGGCA"}],
                "X7270":[{"seq":"CAGTCTTC","instruments":[],"seq_in_adapter":"GAAGACTG"}],
                "X7271":[{"seq":"CATAACGG","instruments":[],"seq_in_adapter":"CCGTTATG"}],
                "X7272":[{"seq":"ACTGCTAG","instruments":[],"seq_in_adapter":"CTAGCAGT"}],
                "X7273":[{"seq":"ATTCTGGC","instruments":[],"seq_in_adapter":"GCCAGAAT"}],
                "X7274":[{"seq":"TTCTCTCG","instruments":[],"seq_in_adapter":"CGAGAGAA"}],
                "X7275":[{"seq":"TCCGAGTT","instruments":[],"seq_in_adapter":"AACTCGGA"}],
                "X7276":[{"seq":"CGAACTGT","instruments":[],"seq_in_adapter":"ACAGTTCG"}],
                "X7277":[{"seq":"AACGGTCA","instruments":[],"seq_in_adapter":"TGACCGTT"}],
                "X7278":[{"seq":"AGCAGATG","instruments":[],"seq_in_adapter":"CATCTGCT"}],
                "X7279":[{"seq":"TATCAGCG","instruments":[],"seq_in_adapter":"CGCTGATA"}],
                "X7280":[{"seq":"TCAGACGA","instruments":[],"seq_in_adapter":"TCGTCTGA"}],
                "X7281":[{"seq":"ACCATGTG","instruments":[],"seq_in_adapter":"CACATGGT"}],
                "X7282":[{"seq":"CTAACTCG","instruments":[],"seq_in_adapter":"CGAGTTAG"}],
                "X7283":[{"seq":"GCTTAGCT","instruments":[],"seq_in_adapter":"AGCTAAGC"}],
                "X7284":[{"seq":"CATGGAAC","instruments":[],"seq_in_adapter":"GTTCCATG"}],
                "X7285":[{"seq":"TAGGATGC","instruments":[],"seq_in_adapter":"GCATCCTA"}],
                "X7286":[{"seq":"GTTCATGG","instruments":[],"seq_in_adapter":"CCATGAAC"}],
                "X7287":[{"seq":"TCGTGGAT","instruments":[],"seq_in_adapter":"ATCCACGA"}],
                "X7288":[{"seq":"ACCTTCTC","instruments":[],"seq_in_adapter":"GAGAAGGT"}],
                "X7289":[{"seq":"CATTGCCT","instruments":[],"seq_in_adapter":"AGGCAATG"}],
                "X7290":[{"seq":"CTAGGTGA","instruments":[],"seq_in_adapter":"TCACCTAG"}],
                "X7291":[{"seq":"TCCGTATG","instruments":[],"seq_in_adapter":"CATACGGA"}],
                "X7292":[{"seq":"ACGATGAC","instruments":[],"seq_in_adapter":"GTCATCGT"}],
                "X7293":[{"seq":"GTCGGTAA","instruments":[],"seq_in_adapter":"TTACCGAC"}],
                "X7294":[{"seq":"TCGAAGGT","instruments":[],"seq_in_adapter":"ACCTTCGA"}],
                "X7295":[{"seq":"AGAAGCGT","instruments":[],"seq_in_adapter":"ACGCTTCT"}],
                "X7296":[{"seq":"CTCTACTC","instruments":[],"seq_in_adapter":"GAGTAGAG"}],
                "X7297":[{"seq":"CTAGGCAT","instruments":[],"seq_in_adapter":"ATGCCTAG"}],
                "X7298":[{"seq":"TGGAGTTG","instruments":[],"seq_in_adapter":"CAACTCCA"}],
                "X7299":[{"seq":"GAGGACTT","instruments":[],"seq_in_adapter":"AAGTCCTC"}],
                "X7300":[{"seq":"CAATCGAC","instruments":[],"seq_in_adapter":"GTCGATTG"}],
                "X7301":[{"seq":"TCTAACGC","instruments":[],"seq_in_adapter":"GCGTTAGA"}],
                "X7302":[{"seq":"TCTCGCAA","instruments":[],"seq_in_adapter":"TTGCGAGA"}],
                "X7303":[{"seq":"ATCGGTGT","instruments":[],"seq_in_adapter":"ACACCGAT"}],
                "X7304":[{"seq":"GAGATACG","instruments":[],"seq_in_adapter":"CGTATCTC"}],
                "X7305":[{"seq":"GTCTCCTT","instruments":[],"seq_in_adapter":"AAGGAGAC"}],
                "X7306":[{"seq":"AGTCGACA","instruments":[],"seq_in_adapter":"TGTCGACT"}],
                "X7307":[{"seq":"CGGATTGA","instruments":[],"seq_in_adapter":"TCAATCCG"}],
                "X7308":[{"seq":"CACAAGTC","instruments":[],"seq_in_adapter":"GACTTGTG"}],
                "X7309":[{"seq":"TACATCGG","instruments":[],"seq_in_adapter":"CCGATGTA"}],
                "X7310":[{"seq":"AGCTCCTA","instruments":[],"seq_in_adapter":"TAGGAGCT"}],
                "X7311":[{"seq":"ACTCGTTG","instruments":[],"seq_in_adapter":"CAACGAGT"}],
                "X7312":[{"seq":"CTGACACA","instruments":[],"seq_in_adapter":"TGTGTCAG"}],
                "X7313":[{"seq":"CAACCTAG","instruments":[],"seq_in_adapter":"CTAGGTTG"}],
                "X7314":[{"seq":"AAGGACAC","instruments":[],"seq_in_adapter":"GTGTCCTT"}],
                "X7315":[{"seq":"TGCAGGTA","instruments":[],"seq_in_adapter":"TACCTGCA"}],
                "X7316":[{"seq":"ACCTAAGG","instruments":[],"seq_in_adapter":"CCTTAGGT"}],
                "X7317":[{"seq":"AGTCTGTG","instruments":[],"seq_in_adapter":"CACAGACT"}],
                "X7318":[{"seq":"AGGTTCGA","instruments":[],"seq_in_adapter":"TCGAACCT"}],
                "X7319":[{"seq":"GACTATGC","instruments":[],"seq_in_adapter":"GCATAGTC"}],
                "X7320":[{"seq":"TTCAGGAG","instruments":[],"seq_in_adapter":"CTCCTGAA"}],
                "X7321":[{"seq":"TGTGCGTT","instruments":[],"seq_in_adapter":"AACGCACA"}],
                "X7322":[{"seq":"CGAGACTA","instruments":[],"seq_in_adapter":"TAGTCTCG"}],
                "X7323":[{"seq":"CTCAGAGT","instruments":[],"seq_in_adapter":"ACTCTGAG"}],
                "X7324":[{"seq":"GCCATAAC","instruments":[],"seq_in_adapter":"GTTATGGC"}],
                "X7325":[{"seq":"TTACCGAG","instruments":[],"seq_in_adapter":"CTCGGTAA"}],
                "X7326":[{"seq":"GCTCTGTA","instruments":[],"seq_in_adapter":"TACAGAGC"}],
                "X7327":[{"seq":"CGTTATGC","instruments":[],"seq_in_adapter":"GCATAACG"}],
                "X7328":[{"seq":"GTCTGATC","instruments":[],"seq_in_adapter":"GATCAGAC"}],
                "X7329":[{"seq":"TAGTTGCG","instruments":[],"seq_in_adapter":"CGCAACTA"}],
                "X7330":[{"seq":"TGATCGGA","instruments":[],"seq_in_adapter":"TCCGATCA"}],
                "X7331":[{"seq":"CCAAGTTG","instruments":[],"seq_in_adapter":"CAACTTGG"}],
                "X7332":[{"seq":"CCTACTGA","instruments":[],"seq_in_adapter":"TCAGTAGG"}],
                "X7333":[{"seq":"CTTGCTGT","instruments":[],"seq_in_adapter":"ACAGCAAG"}],
                "X7334":[{"seq":"TGCCATTC","instruments":[],"seq_in_adapter":"GAATGGCA"}],
                "X7335":[{"seq":"TTGATCCG","instruments":[],"seq_in_adapter":"CGGATCAA"}],
                "X7336":[{"seq":"AGTGCAGT","instruments":[],"seq_in_adapter":"ACTGCACT"}],
                "X7337":[{"seq":"GACTTAGG","instruments":[],"seq_in_adapter":"CCTAAGTC"}],
                "X7338":[{"seq":"CGTACGAA","instruments":[],"seq_in_adapter":"TTCGTACG"}],
                "X7339":[{"seq":"TACCAGGA","instruments":[],"seq_in_adapter":"TCCTGGTA"}],
                "X7340":[{"seq":"CGTCAATG","instruments":[],"seq_in_adapter":"CATTGACG"}],
                "X7341":[{"seq":"GAAGAGGT","instruments":[],"seq_in_adapter":"ACCTCTTC"}],
                "X7342":[{"seq":"GACGAATG","instruments":[],"seq_in_adapter":"CATTCGTC"}],
                "X7343":[{"seq":"AGGAGGAA","instruments":[],"seq_in_adapter":"TTCCTCCT"}],
                "X7344":[{"seq":"CTTACAGC","instruments":[],"seq_in_adapter":"GCTGTAAG"}],
                "X7345":[{"seq":"GAGATGTC","instruments":[],"seq_in_adapter":"GACATCTC"}],
                "X7346":[{"seq":"TACGGTTG","instruments":[],"seq_in_adapter":"CAACCGTA"}],
                "X7347":[{"seq":"CTATCGCA","instruments":[],"seq_in_adapter":"TGCGATAG"}],
                "X7348":[{"seq":"TCGAACCA","instruments":[],"seq_in_adapter":"TGGTTCGA"}],
                "X7349":[{"seq":"GAACGCTT","instruments":[],"seq_in_adapter":"AAGCGTTC"}],
                "X7350":[{"seq":"CAGAATCG","instruments":[],"seq_in_adapter":"CGATTCTG"}],
                "X7351":[{"seq":"ATGGTTGC","instruments":[],"seq_in_adapter":"GCAACCAT"}],
                "X7352":[{"seq":"GCTGGATT","instruments":[],"seq_in_adapter":"AATCCAGC"}],
                "X7353":[{"seq":"GATGCACT","instruments":[],"seq_in_adapter":"AGTGCATC"}],
                "X7354":[{"seq":"ACCAATGC","instruments":[],"seq_in_adapter":"GCATTGGT"}],
                "X7355":[{"seq":"GTCCTAAG","instruments":[],"seq_in_adapter":"CTTAGGAC"}],
                "X7356":[{"seq":"CCGACTAT","instruments":[],"seq_in_adapter":"ATAGTCGG"}],
                "X7357":[{"seq":"TTGGTCTC","instruments":[],"seq_in_adapter":"GAGACCAA"}],
                "X7358":[{"seq":"GCCTTGTT","instruments":[],"seq_in_adapter":"AACAAGGC"}],
                "X7359":[{"seq":"GATACTGG","instruments":[],"seq_in_adapter":"CCAGTATC"}],
                "X7360":[{"seq":"ATTCGAGG","instruments":[],"seq_in_adapter":"CCTCGAAT"}],
                "X7361":[{"seq":"GTCAGTTG","instruments":[],"seq_in_adapter":"CAACTGAC"}],
                "X7362":[{"seq":"GTAGAGCA","instruments":[],"seq_in_adapter":"TGCTCTAC"}],
                "X7363":[{"seq":"ACGTGATG","instruments":[],"seq_in_adapter":"CATCACGT"}],
                "X7364":[{"seq":"TAAGTGGC","instruments":[],"seq_in_adapter":"GCCACTTA"}],
                "X7365":[{"seq":"TGTGAAGC","instruments":[],"seq_in_adapter":"GCTTCACA"}],
                "X7366":[{"seq":"CATTCGGT","instruments":[],"seq_in_adapter":"ACCGAATG"}],
                "X7367":[{"seq":"TTGGTGAG","instruments":[],"seq_in_adapter":"CTCACCAA"}],
                "X7368":[{"seq":"CAGTTCTG","instruments":[],"seq_in_adapter":"CAGAACTG"}],
                "X7369":[{"seq":"AGGCTTCT","instruments":[],"seq_in_adapter":"AGAAGCCT"}],
                "X7370":[{"seq":"GAATCGTG","instruments":[],"seq_in_adapter":"CACGATTC"}],
                "X7371":[{"seq":"ACCAGCTT","instruments":[],"seq_in_adapter":"AAGCTGGT"}],
                "X7372":[{"seq":"CTCATTGC","instruments":[],"seq_in_adapter":"GCAATGAG"}],
                "X7373":[{"seq":"CGATAGAG","instruments":[],"seq_in_adapter":"CTCTATCG"}],
                "X7374":[{"seq":"TGGAGAGT","instruments":[],"seq_in_adapter":"ACTCTCCA"}],
                "X7375":[{"seq":"GTATGCTG","instruments":[],"seq_in_adapter":"CAGCATAC"}],
                "X7376":[{"seq":"CTGGAGTA","instruments":[],"seq_in_adapter":"TACTCCAG"}],
                "X7377":[{"seq":"AATGCCTC","instruments":[],"seq_in_adapter":"GAGGCATT"}],
                "X7378":[{"seq":"TGAGGTGT","instruments":[],"seq_in_adapter":"ACACCTCA"}],
                "X7379":[{"seq":"ACATTGCG","instruments":[],"seq_in_adapter":"CGCAATGT"}],
                "X7380":[{"seq":"TCTCTAGG","instruments":[],"seq_in_adapter":"CCTAGAGA"}],
                "X7381":[{"seq":"CGCTAGTA","instruments":[],"seq_in_adapter":"TACTAGCG"}],
                "X7382":[{"seq":"AATGGACG","instruments":[],"seq_in_adapter":"CGTCCATT"}],
                "X7383":[{"seq":"GATAGCGA","instruments":[],"seq_in_adapter":"TCGCTATC"}],
                "X7384":[{"seq":"CGACCATT","instruments":[],"seq_in_adapter":"AATGGTCG"}],
            },
            "i5":{
                "X5001":[{"seq":"ATATGCGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCGCATAT"}, {"seq":"GCGCATAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCGCATAT"}],
                "X5002":[{"seq":"TGGTACAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGTACCA"}, {"seq":"CTGTACCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGTACCA"}],
                "X5003":[{"seq":"AACCGTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAACGGTT"}, {"seq":"GAACGGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAACGGTT"}],
                "X5004":[{"seq":"TAACCGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCGGTTA"}, {"seq":"ACCGGTTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCGGTTA"}],
                "X5005":[{"seq":"GAACATCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGATGTTC"}, {"seq":"CGATGTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGATGTTC"}],
                "X5006":[{"seq":"CCTTGTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTACAAGG"}, {"seq":"CTACAAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTACAAGG"}],
                "X5007":[{"seq":"TCAGGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGCCTGA"}, {"seq":"AAGCCTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGCCTGA"}],
                "X5008":[{"seq":"GTTCTCGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGAGAAC"}, {"seq":"ACGAGAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGAGAAC"}],
                "X5009":[{"seq":"AGAACGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCGTTCT"}, {"seq":"CTCGTTCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCGTTCT"}],
                "X5010":[{"seq":"TGCTTCCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGAAGCA"}, {"seq":"TGGAAGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGAAGCA"}],
                "X5011":[{"seq":"CTTCGACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTCGAAG"}, {"seq":"AGTCGAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTCGAAG"}],
                "X5012":[{"seq":"CACCTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACAGGTG"}, {"seq":"AACAGGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACAGGTG"}],
                "X5013":[{"seq":"ATCACACG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTGTGAT"}, {"seq":"CGTGTGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTGTGAT"}],
                "X5014":[{"seq":"CCGTAAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTTACGG"}, {"seq":"TCTTACGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTTACGG"}],
                "X5015":[{"seq":"TACGCCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGGCGTA"}, {"seq":"AAGGCGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGCGTA"}],
                "X5016":[{"seq":"CGACGTTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAACGTCG"}, {"seq":"TAACGTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAACGTCG"}],
                "X5017":[{"seq":"ATGCACGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGTGCAT"}, {"seq":"TCGTGCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGTGCAT"}],
                "X5018":[{"seq":"CCTGATTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAATCAGG"}, {"seq":"CAATCAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAATCAGG"}],
                "X5019":[{"seq":"GTAGGAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTCCTAC"}, {"seq":"ACTCCTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTCCTAC"}],
                "X5020":[{"seq":"ACTAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCCTAGT"}, {"seq":"CTCCTAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCCTAGT"}],
                "X5021":[{"seq":"CACTAGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCTAGTG"}, {"seq":"AGCTAGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCTAGTG"}],
                "X5022":[{"seq":"ACGACTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAAGTCGT"}, {"seq":"CAAGTCGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAAGTCGT"}],
                "X5023":[{"seq":"CGTGTGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACACG"}, {"seq":"TACACACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACACG"}],
                "X5024":[{"seq":"GTTGACCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGTCAAC"}, {"seq":"AGGTCAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGTCAAC"}],
                "X5025":[{"seq":"ACTCCATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATGGAGT"}, {"seq":"GATGGAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATGGAGT"}],
                "X5026":[{"seq":"CAATGTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCACATTG"}, {"seq":"CCACATTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCACATTG"}],
                "X5027":[{"seq":"TTGCAGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCTGCAA"}, {"seq":"GTCTGCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCTGCAA"}],
                "X5028":[{"seq":"CAGTCCAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGGACTG"}, {"seq":"TTGGACTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGGACTG"}],
                "X5029":[{"seq":"ACGTTCAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGAACGT"}, {"seq":"CTGAACGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGAACGT"}],
                "X5030":[{"seq":"AACGTCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGACGTT"}, {"seq":"CAGACGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGACGTT"}],
                "X5031":[{"seq":"TATCGGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACCGATA"}, {"seq":"GACCGATA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACCGATA"}],
                "X5032":[{"seq":"CGCTCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATAGAGCG"}, {"seq":"ATAGAGCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATAGAGCG"}],
                "X5033":[{"seq":"GATTGCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGCAATC"}, {"seq":"GAGCAATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGCAATC"}],
                "X5034":[{"seq":"GATGTGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACACATC"}, {"seq":"CACACATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACACATC"}],
                "X5035":[{"seq":"CGCAATCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGATTGCG"}, {"seq":"AGATTGCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGATTGCG"}],
                "X5036":[{"seq":"TGGTAGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCTACCA"}, {"seq":"AGCTACCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCTACCA"}],
                "X5037":[{"seq":"GATAGGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCCTATC"}, {"seq":"AGCCTATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCCTATC"}],
                "X5038":[{"seq":"AGTGGATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATCCACT"}, {"seq":"GATCCACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATCCACT"}],
                "X5039":[{"seq":"TTGGACGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGTCCAA"}, {"seq":"ACGTCCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGTCCAA"}],
                "X5040":[{"seq":"ATGACGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACGTCAT"}, {"seq":"GACGTCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACGTCAT"}],
                "X5041":[{"seq":"GAAGTTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAACTTC"}, {"seq":"CCAACTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAACTTC"}],
                "X5042":[{"seq":"CATACCAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTGGTATG"}, {"seq":"GTGGTATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTGGTATG"}],
                "X5043":[{"seq":"CTGTTGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCAACAG"}, {"seq":"GTCAACAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCAACAG"}],
                "X5044":[{"seq":"TGGCATGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACATGCCA"}, {"seq":"ACATGCCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACATGCCA"}],
                "X5045":[{"seq":"ATCGCCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATGGCGAT"}, {"seq":"ATGGCGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATGGCGAT"}],
                "X5046":[{"seq":"TTGCGAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTCGCAA"}, {"seq":"CTTCGCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTCGCAA"}],
                "X5047":[{"seq":"AGTTCGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACGAACT"}, {"seq":"GACGAACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACGAACT"}],
                "X5048":[{"seq":"GAGCAGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACTGCTC"}, {"seq":"TACTGCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACTGCTC"}],
                "X5049":[{"seq":"ACAGCTCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGAGCTGT"}, {"seq":"TGAGCTGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGAGCTGT"}],
                "X5050":[{"seq":"GATCGAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTCGATC"}, {"seq":"ACTCGATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTCGATC"}],
                "X5051":[{"seq":"AGCGTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACACGCT"}, {"seq":"AACACGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACACGCT"}],
                "X5052":[{"seq":"GTTACGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCGTAAC"}, {"seq":"TGCGTAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCGTAAC"}],
                "X5053":[{"seq":"TGAAGACG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTCTTCA"}, {"seq":"CGTCTTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTCTTCA"}],
                "X5054":[{"seq":"ACTGAGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCTCAGT"}, {"seq":"ACCTCAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCTCAGT"}],
                "X5055":[{"seq":"CGGTTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACAACCG"}, {"seq":"AACAACCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACAACCG"}],
                "X5056":[{"seq":"GTTGTTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGAACAAC"}, {"seq":"CGAACAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGAACAAC"}],
                "X5057":[{"seq":"GAAGGAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTCCTTC"}, {"seq":"CTTCCTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTCCTTC"}],
                "X5058":[{"seq":"AGCACTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAAGTGCT"}, {"seq":"GAAGTGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAAGTGCT"}],
                "X5059":[{"seq":"GTCATCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGATGAC"}, {"seq":"TCGATGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGATGAC"}],
                "X5060":[{"seq":"TGTGACTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGTCACA"}, {"seq":"CAGTCACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGTCACA"}],
                "X5061":[{"seq":"CAACACCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGTGTTG"}, {"seq":"AGGTGTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGTGTTG"}],
                "X5062":[{"seq":"ATGCCTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACAGGCAT"}, {"seq":"ACAGGCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACAGGCAT"}],
                "X5063":[{"seq":"CATGGCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGCCATG"}, {"seq":"TAGCCATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGCCATG"}],
                "X5064":[{"seq":"GTGAAGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACTTCAC"}, {"seq":"CACTTCAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACTTCAC"}],
                "X5065":[{"seq":"CGTTGCAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGCAACG"}, {"seq":"TTGCAACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGCAACG"}],
                "X5066":[{"seq":"ATCCGGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACCGGAT"}, {"seq":"TACCGGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACCGGAT"}],
                "X5067":[{"seq":"GCGTCATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AATGACGC"}, {"seq":"AATGACGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AATGACGC"}],
                "X5068":[{"seq":"GCACAACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTTGTGC"}, {"seq":"AGTTGTGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTTGTGC"}],
                "X5069":[{"seq":"GATTACCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGGTAATC"}, {"seq":"CGGTAATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGGTAATC"}],
                "X5070":[{"seq":"ACCACGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCGTGGT"}, {"seq":"ATCGTGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCGTGGT"}],
                "X5071":[{"seq":"GTCGAAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTTCGAC"}, {"seq":"TCTTCGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTTCGAC"}],
                "X5072":[{"seq":"CCTTGATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATCAAGG"}, {"seq":"GATCAAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATCAAGG"}],
                "X5073":[{"seq":"AAGCACTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGTGCTT"}, {"seq":"CAGTGCTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGTGCTT"}],
                "X5074":[{"seq":"TTCGTTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAACGAA"}, {"seq":"CCAACGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAACGAA"}],
                "X5075":[{"seq":"TCGCTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACAGCGA"}, {"seq":"AACAGCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACAGCGA"}],
                "X5076":[{"seq":"GAATCCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGGATTC"}, {"seq":"TCGGATTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGGATTC"}],
                "X5077":[{"seq":"GTGCCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TATGGCAC"}, {"seq":"TATGGCAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TATGGCAC"}],
                "X5078":[{"seq":"CTTAGGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCCTAAG"}, {"seq":"GTCCTAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCCTAAG"}],
                "X5079":[{"seq":"AACTGAGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTCAGTT"}, {"seq":"GCTCAGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTCAGTT"}],
                "X5080":[{"seq":"GACGATCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGATCGTC"}, {"seq":"AGATCGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGATCGTC"}],
                "X5081":[{"seq":"ATCCAGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCTGGAT"}, {"seq":"CTCTGGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCTGGAT"}],
                "X5082":[{"seq":"AGAGTAGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTACTCT"}, {"seq":"GCTACTCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTACTCT"}],
                "X5083":[{"seq":"TGGACTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGAGTCCA"}, {"seq":"AGAGTCCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGAGTCCA"}],
                "X5084":[{"seq":"TACGCTAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTAGCGTA"}, {"seq":"GTAGCGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTAGCGTA"}],
                "X5085":[{"seq":"GCTATCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGATAGC"}, {"seq":"AGGATAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGATAGC"}],
                "X5086":[{"seq":"GCAAGATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATCTTGC"}, {"seq":"GATCTTGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATCTTGC"}],
                "X5087":[{"seq":"ATCGATCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGATCGAT"}, {"seq":"CGATCGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGATCGAT"}],
                "X5088":[{"seq":"CGGCTAAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATTAGCCG"}, {"seq":"ATTAGCCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATTAGCCG"}],
                "X5089":[{"seq":"ACGGAACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTTCCGT"}, {"seq":"TGTTCCGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTTCCGT"}],
                "X5090":[{"seq":"CGCATGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCATGCG"}, {"seq":"ATCATGCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCATGCG"}],
                "X5091":[{"seq":"TTCCAAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTTGGAA"}, {"seq":"CCTTGGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTTGGAA"}],
                "X5092":[{"seq":"CTTGTCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGACAAG"}, {"seq":"TCGACAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGACAAG"}],
                "X5093":[{"seq":"GAGACGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCGTCTC"}, {"seq":"ATCGTCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCGTCTC"}],
                "X5094":[{"seq":"TGAGCTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTAGCTCA"}, {"seq":"CTAGCTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTAGCTCA"}],
                "X5095":[{"seq":"ACTCTCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGAGAGT"}, {"seq":"TCGAGAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGAGAGT"}],
                "X5096":[{"seq":"CTGATCGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGATCAG"}, {"seq":"ACGATCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGATCAG"}],
                "X5097":[{"seq":"CGACCATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AATGGTCG"}, {"seq":"AATGGTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AATGGTCG"}],
                "X5098":[{"seq":"GATAGCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGCTATC"}, {"seq":"TCGCTATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGCTATC"}],
                "X5099":[{"seq":"AATGGACG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTCCATT"}, {"seq":"CGTCCATT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTCCATT"}],
                "X5100":[{"seq":"CGCTAGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACTAGCG"}, {"seq":"TACTAGCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACTAGCG"}],
                "X5101":[{"seq":"TCTCTAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTAGAGA"}, {"seq":"CCTAGAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTAGAGA"}],
                "X5102":[{"seq":"ACATTGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGCAATGT"}, {"seq":"CGCAATGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGCAATGT"}],
                "X5103":[{"seq":"TGAGGTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACACCTCA"}, {"seq":"ACACCTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACACCTCA"}],
                "X5104":[{"seq":"AATGCCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGGCATT"}, {"seq":"GAGGCATT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGGCATT"}],
                "X5105":[{"seq":"CTGGAGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACTCCAG"}, {"seq":"TACTCCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACTCCAG"}],
                "X5106":[{"seq":"GTATGCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGCATAC"}, {"seq":"CAGCATAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGCATAC"}],
                "X5107":[{"seq":"TGGAGAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTCTCCA"}, {"seq":"ACTCTCCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTCTCCA"}],
                "X5108":[{"seq":"CGATAGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCTATCG"}, {"seq":"CTCTATCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCTATCG"}],
                "X5109":[{"seq":"CTCATTGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCAATGAG"}, {"seq":"GCAATGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCAATGAG"}],
                "X5110":[{"seq":"ACCAGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGCTGGT"}, {"seq":"AAGCTGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGCTGGT"}],
                "X5111":[{"seq":"GAATCGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACGATTC"}, {"seq":"CACGATTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACGATTC"}],
                "X5112":[{"seq":"AGGCTTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGAAGCCT"}, {"seq":"AGAAGCCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGAAGCCT"}],
                "X5113":[{"seq":"CAGTTCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGAACTG"}, {"seq":"CAGAACTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGAACTG"}],
                "X5114":[{"seq":"TTGGTGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCACCAA"}, {"seq":"CTCACCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCACCAA"}],
                "X5115":[{"seq":"CATTCGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCGAATG"}, {"seq":"ACCGAATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCGAATG"}],
                "X5116":[{"seq":"TGTGAAGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTTCACA"}, {"seq":"GCTTCACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTTCACA"}],
                "X5117":[{"seq":"TAAGTGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCCACTTA"}, {"seq":"GCCACTTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCCACTTA"}],
                "X5118":[{"seq":"ACGTGATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATCACGT"}, {"seq":"CATCACGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATCACGT"}],
                "X5119":[{"seq":"GTAGAGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCTCTAC"}, {"seq":"TGCTCTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCTCTAC"}],
                "X5120":[{"seq":"GTCAGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAACTGAC"}, {"seq":"CAACTGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAACTGAC"}],
                "X5121":[{"seq":"ATTCGAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTCGAAT"}, {"seq":"CCTCGAAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTCGAAT"}],
                "X5122":[{"seq":"GATACTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAGTATC"}, {"seq":"CCAGTATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAGTATC"}],
                "X5123":[{"seq":"GCCTTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACAAGGC"}, {"seq":"AACAAGGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACAAGGC"}],
                "X5124":[{"seq":"TTGGTCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGACCAA"}, {"seq":"GAGACCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGACCAA"}],
                "X5125":[{"seq":"CCGACTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATAGTCGG"}, {"seq":"ATAGTCGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATAGTCGG"}],
                "X5126":[{"seq":"GTCCTAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTAGGAC"}, {"seq":"CTTAGGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTAGGAC"}],
                "X5127":[{"seq":"ACCAATGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCATTGGT"}, {"seq":"GCATTGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCATTGGT"}],
                "X5128":[{"seq":"GATGCACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTGCATC"}, {"seq":"AGTGCATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTGCATC"}],
                "X5129":[{"seq":"GCTGGATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AATCCAGC"}, {"seq":"AATCCAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AATCCAGC"}],
                "X5130":[{"seq":"ATGGTTGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCAACCAT"}, {"seq":"GCAACCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCAACCAT"}],
                "X5131":[{"seq":"CAGAATCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGATTCTG"}, {"seq":"CGATTCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGATTCTG"}],
                "X5132":[{"seq":"GAACGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGCGTTC"}, {"seq":"AAGCGTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGCGTTC"}],
                "X5133":[{"seq":"TCGAACCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGTTCGA"}, {"seq":"TGGTTCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGTTCGA"}],
                "X5134":[{"seq":"CTATCGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCGATAG"}, {"seq":"TGCGATAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCGATAG"}],
                "X5135":[{"seq":"TACGGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAACCGTA"}, {"seq":"CAACCGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAACCGTA"}],
                "X5136":[{"seq":"GAGATGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACATCTC"}, {"seq":"GACATCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACATCTC"}],
                "X5137":[{"seq":"CTTACAGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTGTAAG"}, {"seq":"GCTGTAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTGTAAG"}],
                "X5138":[{"seq":"AGGAGGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCCTCCT"}, {"seq":"TTCCTCCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCCTCCT"}],
                "X5139":[{"seq":"GACGAATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATTCGTC"}, {"seq":"CATTCGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATTCGTC"}],
                "X5140":[{"seq":"GAAGAGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCTCTTC"}, {"seq":"ACCTCTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCTCTTC"}],
                "X5141":[{"seq":"CGTCAATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATTGACG"}, {"seq":"CATTGACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATTGACG"}],
                "X5142":[{"seq":"TACCAGGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCTGGTA"}, {"seq":"TCCTGGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCTGGTA"}],
                "X5143":[{"seq":"CGTACGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCGTACG"}, {"seq":"TTCGTACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCGTACG"}],
                "X5144":[{"seq":"GACTTAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTAAGTC"}, {"seq":"CCTAAGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTAAGTC"}],
                "X5145":[{"seq":"AGTGCAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTGCACT"}, {"seq":"ACTGCACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTGCACT"}],
                "X5146":[{"seq":"TTGATCCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGGATCAA"}, {"seq":"CGGATCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGGATCAA"}],
                "X5147":[{"seq":"TGCCATTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAATGGCA"}, {"seq":"GAATGGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAATGGCA"}],
                "X5148":[{"seq":"CTTGCTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACAGCAAG"}, {"seq":"ACAGCAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACAGCAAG"}],
                "X5149":[{"seq":"CCTACTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCAGTAGG"}, {"seq":"TCAGTAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCAGTAGG"}],
                "X5150":[{"seq":"CCAAGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAACTTGG"}, {"seq":"CAACTTGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAACTTGG"}],
                "X5151":[{"seq":"TGATCGGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCGATCA"}, {"seq":"TCCGATCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCGATCA"}],
                "X5152":[{"seq":"TAGTTGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGCAACTA"}, {"seq":"CGCAACTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGCAACTA"}],
                "X5153":[{"seq":"GTCTGATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATCAGAC"}, {"seq":"GATCAGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATCAGAC"}],
                "X5154":[{"seq":"CGTTATGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCATAACG"}, {"seq":"GCATAACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCATAACG"}],
                "X5155":[{"seq":"GCTCTGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACAGAGC"}, {"seq":"TACAGAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACAGAGC"}],
                "X5156":[{"seq":"TTACCGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCGGTAA"}, {"seq":"CTCGGTAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCGGTAA"}],
                "X5157":[{"seq":"GCCATAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTATGGC"}, {"seq":"GTTATGGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTATGGC"}],
                "X5158":[{"seq":"CTCAGAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTCTGAG"}, {"seq":"ACTCTGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTCTGAG"}],
                "X5159":[{"seq":"CGAGACTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGTCTCG"}, {"seq":"TAGTCTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGTCTCG"}],
                "X5160":[{"seq":"TGTGCGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACGCACA"}, {"seq":"AACGCACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACGCACA"}],
                "X5161":[{"seq":"TTCAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCCTGAA"}, {"seq":"CTCCTGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCCTGAA"}],
                "X5162":[{"seq":"GACTATGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCATAGTC"}, {"seq":"GCATAGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCATAGTC"}],
                "X5163":[{"seq":"AGGTTCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGAACCT"}, {"seq":"TCGAACCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGAACCT"}],
                "X5164":[{"seq":"AGTCTGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACAGACT"}, {"seq":"CACAGACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACAGACT"}],
                "X5165":[{"seq":"ACCTAAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTTAGGT"}, {"seq":"CCTTAGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTTAGGT"}],
                "X5166":[{"seq":"TGCAGGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACCTGCA"}, {"seq":"TACCTGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACCTGCA"}],
                "X5167":[{"seq":"AAGGACAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTGTCCTT"}, {"seq":"GTGTCCTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTGTCCTT"}],
                "X5168":[{"seq":"CAACCTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTAGGTTG"}, {"seq":"CTAGGTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTAGGTTG"}],
                "X5169":[{"seq":"CTGACACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTGTCAG"}, {"seq":"TGTGTCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTGTCAG"}],
                "X5170":[{"seq":"ACTCGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAACGAGT"}, {"seq":"CAACGAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAACGAGT"}],
                "X5171":[{"seq":"AGCTCCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGGAGCT"}, {"seq":"TAGGAGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGGAGCT"}],
                "X5172":[{"seq":"TACATCGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCGATGTA"}, {"seq":"CCGATGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCGATGTA"}],
                "X5173":[{"seq":"CACAAGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACTTGTG"}, {"seq":"GACTTGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACTTGTG"}],
                "X5174":[{"seq":"CGGATTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCAATCCG"}, {"seq":"TCAATCCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCAATCCG"}],
                "X5175":[{"seq":"AGTCGACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTCGACT"}, {"seq":"TGTCGACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTCGACT"}],
                "X5176":[{"seq":"GTCTCCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGGAGAC"}, {"seq":"AAGGAGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGAGAC"}],
                "X5177":[{"seq":"GAGATACG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTATCTC"}, {"seq":"CGTATCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTATCTC"}],
                "X5178":[{"seq":"ATCGGTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACACCGAT"}, {"seq":"ACACCGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACACCGAT"}],
                "X5179":[{"seq":"TCTCGCAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGCGAGA"}, {"seq":"TTGCGAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGCGAGA"}],
                "X5180":[{"seq":"TCTAACGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCGTTAGA"}, {"seq":"GCGTTAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCGTTAGA"}],
                "X5181":[{"seq":"CAATCGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCGATTG"}, {"seq":"GTCGATTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCGATTG"}],
                "X5182":[{"seq":"GAGGACTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGTCCTC"}, {"seq":"AAGTCCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGTCCTC"}],
                "X5183":[{"seq":"TGGAGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAACTCCA"}, {"seq":"CAACTCCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAACTCCA"}],
                "X5184":[{"seq":"CTAGGCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATGCCTAG"}, {"seq":"ATGCCTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATGCCTAG"}],
                "X5185":[{"seq":"CTCTACTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGTAGAG"}, {"seq":"GAGTAGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGTAGAG"}],
                "X5186":[{"seq":"AGAAGCGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGCTTCT"}, {"seq":"ACGCTTCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGCTTCT"}],
                "X5187":[{"seq":"TCGAAGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCTTCGA"}, {"seq":"ACCTTCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCTTCGA"}],
                "X5188":[{"seq":"GTCGGTAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTACCGAC"}, {"seq":"TTACCGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTACCGAC"}],
                "X5189":[{"seq":"ACGATGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCATCGT"}, {"seq":"GTCATCGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCATCGT"}],
                "X5190":[{"seq":"TCCGTATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATACGGA"}, {"seq":"CATACGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATACGGA"}],
                "X5191":[{"seq":"CTAGGTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCACCTAG"}, {"seq":"TCACCTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCACCTAG"}],
                "X5192":[{"seq":"CATTGCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGCAATG"}, {"seq":"AGGCAATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGCAATG"}],
                "X5193":[{"seq":"ACCTTCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGAAGGT"}, {"seq":"GAGAAGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGAAGGT"}],
                "X5194":[{"seq":"TCGTGGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCCACGA"}, {"seq":"ATCCACGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCCACGA"}],
                "X5195":[{"seq":"GTTCATGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCATGAAC"}, {"seq":"CCATGAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCATGAAC"}],
                "X5196":[{"seq":"TAGGATGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCATCCTA"}, {"seq":"GCATCCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCATCCTA"}],
                "X5197":[{"seq":"CATGGAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTCCATG"}, {"seq":"GTTCCATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTCCATG"}],
                "X5198":[{"seq":"GCTTAGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCTAAGC"}, {"seq":"AGCTAAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCTAAGC"}],
                "X5199":[{"seq":"CTAACTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGAGTTAG"}, {"seq":"CGAGTTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGAGTTAG"}],
                "X5200":[{"seq":"ACCATGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACATGGT"}, {"seq":"CACATGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACATGGT"}],
                "X5201":[{"seq":"TCAGACGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGTCTGA"}, {"seq":"TCGTCTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGTCTGA"}],
                "X5202":[{"seq":"TATCAGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGCTGATA"}, {"seq":"CGCTGATA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGCTGATA"}],
                "X5203":[{"seq":"AGCAGATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATCTGCT"}, {"seq":"CATCTGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATCTGCT"}],
                "X5204":[{"seq":"AACGGTCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGACCGTT"}, {"seq":"TGACCGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGACCGTT"}],
                "X5205":[{"seq":"CGAACTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACAGTTCG"}, {"seq":"ACAGTTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACAGTTCG"}],
                "X5206":[{"seq":"TCCGAGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACTCGGA"}, {"seq":"AACTCGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACTCGGA"}],
                "X5207":[{"seq":"TTCTCTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGAGAGAA"}, {"seq":"CGAGAGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGAGAGAA"}],
                "X5208":[{"seq":"ATTCTGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCCAGAAT"}, {"seq":"GCCAGAAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCCAGAAT"}],
                "X5209":[{"seq":"ACTGCTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTAGCAGT"}, {"seq":"CTAGCAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTAGCAGT"}],
                "X5210":[{"seq":"CATAACGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCGTTATG"}, {"seq":"CCGTTATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCGTTATG"}],
                "X5211":[{"seq":"CAGTCTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAAGACTG"}, {"seq":"GAAGACTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAAGACTG"}],
                "X5212":[{"seq":"TGCCTCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGAGGCA"}, {"seq":"AAGAGGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGAGGCA"}],
                "X5213":[{"seq":"ACTGTGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACACAGT"}, {"seq":"GACACAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACACAGT"}],
                "X5214":[{"seq":"GTATTGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCCAATAC"}, {"seq":"GCCAATAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCCAATAC"}],
                "X5215":[{"seq":"CGATGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGCATCG"}, {"seq":"AAGCATCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGCATCG"}],
                "X5216":[{"seq":"AAGGCTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCAGCCTT"}, {"seq":"TCAGCCTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCAGCCTT"}],
                "X5217":[{"seq":"AGTCAGGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCTGACT"}, {"seq":"TCCTGACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCTGACT"}],
                "X5218":[{"seq":"CAGGTATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATACCTG"}, {"seq":"GATACCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATACCTG"}],
                "X5219":[{"seq":"TCTCCGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCGGAGA"}, {"seq":"ATCGGAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCGGAGA"}],
                "X5220":[{"seq":"TTCAGCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGCTGAA"}, {"seq":"AGGCTGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGCTGAA"}],
                "X5221":[{"seq":"TCTGAGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCTCAGA"}, {"seq":"CTCTCAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCTCAGA"}],
                "X5222":[{"seq":"TTAGGTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGACCTAA"}, {"seq":"CGACCTAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGACCTAA"}],
                "X5223":[{"seq":"CTCTGGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACCAGAG"}, {"seq":"AACCAGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACCAGAG"}],
                "X5224":[{"seq":"GCGTTCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGAACGC"}, {"seq":"TAGAACGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGAACGC"}],
                "X5225":[{"seq":"TCACGTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAACGTGA"}, {"seq":"GAACGTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAACGTGA"}],
                "X5226":[{"seq":"AGGATGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCATCCT"}, {"seq":"ACCATCCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCATCCT"}],
                "X5227":[{"seq":"GTGTTCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGAACAC"}, {"seq":"AGGAACAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGAACAC"}],
                "X5228":[{"seq":"GTAGCATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATGCTAC"}, {"seq":"GATGCTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATGCTAC"}],
                "X5229":[{"seq":"AGGATCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGATCCT"}, {"seq":"CAGATCCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGATCCT"}],
                "X5230":[{"seq":"GACAAGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCTTGTC"}, {"seq":"CTCTTGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCTTGTC"}],
                "X5231":[{"seq":"TTACGGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCCGTAA"}, {"seq":"AGCCGTAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCCGTAA"}],
                "X5232":[{"seq":"GCTGTTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACAACAGC"}, {"seq":"ACAACAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACAACAGC"}],
                "X5233":[{"seq":"AACCGAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTCGGTT"}, {"seq":"CTTCGGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTCGGTT"}],
                "X5234":[{"seq":"TCTGCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGAGCAGA"}, {"seq":"AGAGCAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGAGCAGA"}],
                "X5235":[{"seq":"CTCAGCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGCTGAG"}, {"seq":"TAGCTGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGCTGAG"}],
                "X5236":[{"seq":"CTTCACCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGTGAAG"}, {"seq":"TGGTGAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGTGAAG"}],
                "X5237":[{"seq":"GATCGTAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTACGATC"}, {"seq":"GTACGATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTACGATC"}],
                "X5238":[{"seq":"CTACAGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACTGTAG"}, {"seq":"CACTGTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACTGTAG"}],
                "X5239":[{"seq":"TCGAGTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCACTCGA"}, {"seq":"TCACTCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCACTCGA"}],
                "X5240":[{"seq":"CAAGTGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCACTTG"}, {"seq":"TGCACTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCACTTG"}],
                "X5241":[{"seq":"CGAGTATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATACTCG"}, {"seq":"CATACTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATACTCG"}],
                "X5242":[{"seq":"CGTAGGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACCTACG"}, {"seq":"AACCTACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACCTACG"}],
                "X5243":[{"seq":"GCCAGTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATACTGGC"}, {"seq":"ATACTGGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATACTGGC"}],
                "X5244":[{"seq":"ATGGAAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTTCCAT"}, {"seq":"CCTTCCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTTCCAT"}],
                "X5245":[{"seq":"AAGAGCCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGCTCTT"}, {"seq":"TGGCTCTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGCTCTT"}],
                "X5246":[{"seq":"TGCGTAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTACGCA"}, {"seq":"TCTACGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTACGCA"}],
                "X5247":[{"seq":"TACACGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCGTGTA"}, {"seq":"AGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCGTGTA"}],
                "X5248":[{"seq":"CCTTCCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGGAAGG"}, {"seq":"AAGGAAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGAAGG"}],
                "X5249":[{"seq":"ACCGCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TATGCGGT"}, {"seq":"TATGCGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TATGCGGT"}],
                "X5250":[{"seq":"TGGTCCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGGACCA"}, {"seq":"AAGGACCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGACCA"}],
                "X5251":[{"seq":"CCATACGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGTATGG"}, {"seq":"ACGTATGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGTATGG"}],
                "X5252":[{"seq":"AACCTTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAAGGTT"}, {"seq":"CCAAGGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAAGGTT"}],
                "X5253":[{"seq":"CAAGGTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGACCTTG"}, {"seq":"AGACCTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGACCTTG"}],
                "X5254":[{"seq":"GCTTCGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCGAAGC"}, {"seq":"TTCGAAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCGAAGC"}],
                "X5255":[{"seq":"CGGAATAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTATTCCG"}, {"seq":"GTATTCCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTATTCCG"}],
                "X5256":[{"seq":"AACTGGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACCAGTT"}, {"seq":"CACCAGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACCAGTT"}],
                "X5257":[{"seq":"GCTTCTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAAGAAGC"}, {"seq":"CAAGAAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAAGAAGC"}],
                "X5258":[{"seq":"GCAATTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGAATTGC"}, {"seq":"CGAATTGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGAATTGC"}],
                "X5259":[{"seq":"AGGTCACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTGACCT"}, {"seq":"AGTGACCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTGACCT"}],
                "X5260":[{"seq":"CAGCGATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AATCGCTG"}, {"seq":"AATCGCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AATCGCTG"}],
                "X5261":[{"seq":"AACCTCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGAGGTT"}, {"seq":"AGGAGGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGAGGTT"}],
                "X5262":[{"seq":"TCGACATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATGTCGA"}, {"seq":"GATGTCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATGTCGA"}],
                "X5263":[{"seq":"CTGGTTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGAACCAG"}, {"seq":"AGAACCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGAACCAG"}],
                "X5264":[{"seq":"ACAGCAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTGCTGT"}, {"seq":"GTTGCTGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTGCTGT"}],
                "X5265":[{"seq":"GCATACAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGTATGC"}, {"seq":"CTGTATGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGTATGC"}],
                "X5266":[{"seq":"CATCTACG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTAGATG"}, {"seq":"CGTAGATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTAGATG"}],
                "X5267":[{"seq":"TTGTCGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCGACAA"}, {"seq":"ACCGACAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCGACAA"}],
                "X5268":[{"seq":"TAGCCGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCGGCTA"}, {"seq":"TTCGGCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCGGCTA"}],
                "X5269":[{"seq":"AGGCATAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTATGCCT"}, {"seq":"CTATGCCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTATGCCT"}],
                "X5270":[{"seq":"TTGACAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTGTCAA"}, {"seq":"CCTGTCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTGTCAA"}],
                "X5271":[{"seq":"TGCACCAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGGTGCA"}, {"seq":"TTGGTGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGGTGCA"}],
                "X5272":[{"seq":"CCAGTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACACTGG"}, {"seq":"AACACTGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACACTGG"}],
                "X5273":[{"seq":"TGTCCAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTGGACA"}, {"seq":"TCTGGACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTGGACA"}],
                "X5274":[{"seq":"GATTGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCCAATC"}, {"seq":"CTCCAATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCCAATC"}],
                "X5275":[{"seq":"ACGGTCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGACCGT"}, {"seq":"AAGACCGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGACCGT"}],
                "X5276":[{"seq":"CTGCGTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATACGCAG"}, {"seq":"ATACGCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATACGCAG"}],
                "X5277":[{"seq":"CACCACTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGTGGTG"}, {"seq":"TAGTGGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGTGGTG"}],
                "X5278":[{"seq":"TGTGGTAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTACCACA"}, {"seq":"GTACCACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTACCACA"}],
                "X5279":[{"seq":"ACATAGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCCTATGT"}, {"seq":"GCCTATGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCCTATGT"}],
                "X5280":[{"seq":"CAAGCAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTGCTTG"}, {"seq":"ACTGCTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTGCTTG"}],
                "X5281":[{"seq":"GCACGTAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTACGTGC"}, {"seq":"TTACGTGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTACGTGC"}],
                "X5282":[{"seq":"TCGTAGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACTACGA"}, {"seq":"GACTACGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACTACGA"}],
                "X5283":[{"seq":"CACTGACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTCAGTG"}, {"seq":"TGTCAGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTCAGTG"}],
                "X5284":[{"seq":"CGTGTACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTACACG"}, {"seq":"AGTACACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTACACG"}],
                "X5285":[{"seq":"GAGCTCAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGAGCTC"}, {"seq":"TTGAGCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGAGCTC"}],
                "X5286":[{"seq":"ACGTCGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACGACGT"}, {"seq":"TACGACGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACGACGT"}],
                "X5287":[{"seq":"GTCTAGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCTAGAC"}, {"seq":"ACCTAGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCTAGAC"}],
                "X5288":[{"seq":"CTTCGTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAACGAAG"}, {"seq":"GAACGAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAACGAAG"}],
                "X5289":[{"seq":"CTGAGATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATCTCAG"}, {"seq":"GATCTCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATCTCAG"}],
                "X5290":[{"seq":"GTGGATAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTATCCAC"}, {"seq":"CTATCCAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTATCCAC"}],
                "X5291":[{"seq":"ACCATCCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGATGGT"}, {"seq":"TGGATGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGATGGT"}],
                "X5292":[{"seq":"AGAGGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAACCTCT"}, {"seq":"CAACCTCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAACCTCT"}],
                "X5293":[{"seq":"ATGCCAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTGGCAT"}, {"seq":"GTTGGCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTGGCAT"}],
                "X5294":[{"seq":"GTTAAGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCCTTAAC"}, {"seq":"GCCTTAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCCTTAAC"}],
                "X5295":[{"seq":"ACCTGACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTCAGGT"}, {"seq":"AGTCAGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTCAGGT"}],
                "X5296":[{"seq":"GCCACTTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAAGTGGC"}, {"seq":"TAAGTGGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAAGTGGC"}],
                "X5297":[{"seq":"ACCTCTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACAGAGGT"}, {"seq":"ACAGAGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACAGAGGT"}],
                "X5298":[{"seq":"GATCCATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATGGATC"}, {"seq":"CATGGATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATGGATC"}],
                "X5299":[{"seq":"CGCTTAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTAAGCG"}, {"seq":"GTTAAGCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTAAGCG"}],
                "X5300":[{"seq":"TGTCTGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCAGACA"}, {"seq":"AGCAGACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCAGACA"}],
                "X5301":[{"seq":"ATAAGGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGCCTTAT"}, {"seq":"CGCCTTAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGCCTTAT"}],
                "X5302":[{"seq":"CGTCTTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACAAGACG"}, {"seq":"ACAAGACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACAAGACG"}],
                "X5303":[{"seq":"CTCCATGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACATGGAG"}, {"seq":"ACATGGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACATGGAG"}],
                "X5304":[{"seq":"TGTTCGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCGAACA"}, {"seq":"CTCGAACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCGAACA"}],
                "X5305":[{"seq":"AGCAAGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCTTGCT"}, {"seq":"TGCTTGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCTTGCT"}],
                "X5306":[{"seq":"CGTATTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGAATACG"}, {"seq":"CGAATACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGAATACG"}],
                "X5307":[{"seq":"AACGACGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGTCGTT"}, {"seq":"ACGTCGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGTCGTT"}],
                "X5308":[{"seq":"GTTGCGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCGCAAC"}, {"seq":"ATCGCAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCGCAAC"}],
                "X5309":[{"seq":"AACGTGGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCACGTT"}, {"seq":"TCCACGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCACGTT"}],
                "X5310":[{"seq":"CTGTGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAACACAG"}, {"seq":"CAACACAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAACACAG"}],
                "X5311":[{"seq":"TGATACGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCGTATCA"}, {"seq":"GCGTATCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCGTATCA"}],
                "X5312":[{"seq":"GTCCTTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGAAGGAC"}, {"seq":"AGAAGGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGAAGGAC"}],
                "X5313":[{"seq":"ACAGACCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGTCTGT"}, {"seq":"AGGTCTGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGTCTGT"}],
                "X5314":[{"seq":"TGTTGTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCACAACA"}, {"seq":"CCACAACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCACAACA"}],
                "X5315":[{"seq":"CATCGTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCACGATG"}, {"seq":"TCACGATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCACGATG"}],
                "X5316":[{"seq":"CAGGAGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCTCCTG"}, {"seq":"ATCTCCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCTCCTG"}],
                "X5317":[{"seq":"TAGGTAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTACCTA"}, {"seq":"CCTACCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTACCTA"}],
                "X5318":[{"seq":"AGTCGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGCGACT"}, {"seq":"AAGCGACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGCGACT"}],
                "X5319":[{"seq":"CGTTGAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTCAACG"}, {"seq":"ACTCAACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTCAACG"}],
                "X5320":[{"seq":"TTCCTGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACAGGAA"}, {"seq":"CACAGGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACAGGAA"}],
                "X5321":[{"seq":"TCGTCTCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGAGACGA"}, {"seq":"TGAGACGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGAGACGA"}],
                "X5322":[{"seq":"TAACGAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTCGTTA"}, {"seq":"CCTCGTTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTCGTTA"}],
                "X5323":[{"seq":"CTGAAGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCTTCAG"}, {"seq":"AGCTTCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCTTCAG"}],
                "X5324":[{"seq":"ATTGCGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACGCAAT"}, {"seq":"CACGCAAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACGCAAT"}],
                "X5325":[{"seq":"CCAAGACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTCTTGG"}, {"seq":"AGTCTTGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTCTTGG"}],
                "X5326":[{"seq":"GCTGTAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTACAGC"}, {"seq":"CTTACAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTACAGC"}],
                "X5327":[{"seq":"GAGTGGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACCACTC"}, {"seq":"AACCACTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACCACTC"}],
                "X5328":[{"seq":"AGCTTGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCAAGCT"}, {"seq":"CTCAAGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCAAGCT"}],
                "X5329":[{"seq":"ACGACAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTGTCGT"}, {"seq":"TCTGTCGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTGTCGT"}],
                "X5330":[{"seq":"TTCGCAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTGCGAA"}, {"seq":"ACTGCGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTGCGAA"}],
                "X5331":[{"seq":"CCGATGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACATCGG"}, {"seq":"TACATCGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACATCGG"}],
                "X5332":[{"seq":"TGACGCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATGCGTCA"}, {"seq":"ATGCGTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATGCGTCA"}],
                "X5333":[{"seq":"TCGCATTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAATGCGA"}, {"seq":"CAATGCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAATGCGA"}],
                "X5334":[{"seq":"CGTGATCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGATCACG"}, {"seq":"TGATCACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGATCACG"}],
                "X5335":[{"seq":"GTGAGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGCTCAC"}, {"seq":"AAGCTCAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGCTCAC"}],
                "X5336":[{"seq":"AGCGGAAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATTCCGCT"}, {"seq":"ATTCCGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATTCCGCT"}],
                "X5337":[{"seq":"CGAAGAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTCTTCG"}, {"seq":"GTTCTTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTCTTCG"}],
                "X5338":[{"seq":"CCGTATCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGATACGG"}, {"seq":"AGATACGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGATACGG"}],
                "X5339":[{"seq":"GTACTCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGAGTAC"}, {"seq":"GAGAGTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGAGTAC"}],
                "X5340":[{"seq":"AGTGTTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAACACT"}, {"seq":"CCAACACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAACACT"}],
                "X5341":[{"seq":"TGAACCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGGTTCA"}, {"seq":"CAGGTTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGGTTCA"}],
                "X5342":[{"seq":"TCAAGGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCCTTGA"}, {"seq":"GTCCTTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCCTTGA"}],
                "X5343":[{"seq":"GTGCTTAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTAAGCAC"}, {"seq":"GTAAGCAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTAAGCAC"}],
                "X5344":[{"seq":"GTGGTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACACCAC"}, {"seq":"AACACCAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACACCAC"}],
                "X5345":[{"seq":"GCTGACTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGTCAGC"}, {"seq":"TAGTCAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGTCAGC"}],
                "X5346":[{"seq":"TGCGAACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTTCGCA"}, {"seq":"AGTTCGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTTCGCA"}],
                "X5347":[{"seq":"AATACGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGCGTATT"}, {"seq":"CGCGTATT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGCGTATT"}],
                "X5348":[{"seq":"ACACGGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACCGTGT"}, {"seq":"AACCGTGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACCGTGT"}],
                "X5349":[{"seq":"CTGCACTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGTGCAG"}, {"seq":"AAGTGCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGTGCAG"}],
                "X5350":[{"seq":"TCAACTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAGTTGA"}, {"seq":"CCAGTTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAGTTGA"}],
                "X5351":[{"seq":"AGTTGGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCCAACT"}, {"seq":"AGCCAACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCCAACT"}],
                "X5352":[{"seq":"CAGGTTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTAACCTG"}, {"seq":"CTAACCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTAACCTG"}],
                "X5353":[{"seq":"ACACCAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTGGTGT"}, {"seq":"ACTGGTGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTGGTGT"}],
                "X5354":[{"seq":"TGGATCAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTGATCCA"}, {"seq":"GTGATCCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTGATCCA"}],
                "X5355":[{"seq":"TGACTTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGAAGTCA"}, {"seq":"CGAAGTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGAAGTCA"}],
                "X5356":[{"seq":"GTGTCTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCAGACAC"}, {"seq":"TCAGACAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCAGACAC"}],
                "X5357":[{"seq":"AGTTACGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCGTAACT"}, {"seq":"CCGTAACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCGTAACT"}],
                "X5358":[{"seq":"ATCTCGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCGAGAT"}, {"seq":"AGCGAGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCGAGAT"}],
                "X5359":[{"seq":"GAAGGTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAACCTTC"}, {"seq":"GAACCTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAACCTTC"}],
                "X5360":[{"seq":"GAGCTTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACAAGCTC"}, {"seq":"ACAAGCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACAAGCTC"}],
                "X5361":[{"seq":"TCCAATCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGATTGGA"}, {"seq":"CGATTGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGATTGGA"}],
                "X5362":[{"seq":"CGGTCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TATGACCG"}, {"seq":"TATGACCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TATGACCG"}],
                "X5363":[{"seq":"TGGCTATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATAGCCA"}, {"seq":"GATAGCCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATAGCCA"}],
                "X5364":[{"seq":"CAACGGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCCGTTG"}, {"seq":"ATCCGTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCCGTTG"}],
                "X5365":[{"seq":"CTCCTAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTAGGAG"}, {"seq":"TCTAGGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTAGGAG"}],
                "X5366":[{"seq":"CCGGAATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AATTCCGG"}, {"seq":"AATTCCGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AATTCCGG"}],
                "X5367":[{"seq":"TAGACGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACGTCTA"}, {"seq":"CACGTCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACGTCTA"}],
                "X5368":[{"seq":"TGACTGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCAGTCA"}, {"seq":"GTCAGTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCAGTCA"}],
                "X5369":[{"seq":"TAGAGCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGCTCTA"}, {"seq":"GAGCTCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGCTCTA"}],
                "X5370":[{"seq":"TCCGTGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCACGGA"}, {"seq":"TTCACGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCACGGA"}],
                "X5371":[{"seq":"CTCGATAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTATCGAG"}, {"seq":"GTATCGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTATCGAG"}],
                "X5372":[{"seq":"CTTACCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGGTAAG"}, {"seq":"CAGGTAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGGTAAG"}],
                "X5373":[{"seq":"ATGGCGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCGCCAT"}, {"seq":"TTCGCCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCGCCAT"}],
                "X5374":[{"seq":"TCCTACCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGTAGGA"}, {"seq":"AGGTAGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGTAGGA"}],
                "X5375":[{"seq":"CCTCAGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACTGAGG"}, {"seq":"AACTGAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACTGAGG"}],
                "X5376":[{"seq":"CTACTTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAAGTAG"}, {"seq":"CCAAGTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAAGTAG"}],
                "X5377":[{"seq":"TCACAGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCTGTGA"}, {"seq":"TGCTGTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCTGTGA"}],
                "X5378":[{"seq":"CACGTTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACAACGTG"}, {"seq":"ACAACGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACAACGTG"}],
                "X5379":[{"seq":"AAGTCGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCGACTT"}, {"seq":"CTCGACTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCGACTT"}],
                "X5380":[{"seq":"TGTACCGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGGTACA"}, {"seq":"ACGGTACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGGTACA"}],
                "X5381":[{"seq":"CTCATCAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGATGAG"}, {"seq":"CTGATGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGATGAG"}],
                "X5382":[{"seq":"AGTCTCAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTGAGACT"}, {"seq":"GTGAGACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTGAGACT"}],
                "X5383":[{"seq":"CTTGGATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATCCAAG"}, {"seq":"CATCCAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATCCAAG"}],
                "X5384":[{"seq":"GCCTATCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGATAGGC"}, {"seq":"TGATAGGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGATAGGC"}],
            }
        },
        "Sabeti Lab Custom RNA-seq v5":{
            "_short_name":"nextera",
            "i7":{
                "RPI1":[{"seq":"ATCACGAT","instruments":[],"seq_in_adapter":"ATCGTGAT"}],
                "RPI2":[{"seq":"CGATGTAT","instruments":[],"seq_in_adapter":"ATACATCG"}],
                "RPI3":[{"seq":"TTAGGCAT","instruments":[],"seq_in_adapter":"ATGCCTAA"}],
                "RPI4":[{"seq":"TGACCAAT","instruments":[],"seq_in_adapter":"ATTGGTCA"}],
                "RPI5":[{"seq":"ACAGTGAT","instruments":[],"seq_in_adapter":"ATCACTGT"}],
                "RPI6":[{"seq":"GCCAATAT","instruments":[],"seq_in_adapter":"ATATTGGC"}],
                "RPI7":[{"seq":"CAGATCAT","instruments":[],"seq_in_adapter":"ATGATCTG"}],
                "RPI8":[{"seq":"ACTTGAAT","instruments":[],"seq_in_adapter":"ATTCAAGT"}],
                "RPI9":[{"seq":"GATCAGAT","instruments":[],"seq_in_adapter":"ATCTGATC"}],
                "RPI10":[{"seq":"TAGCTTAT","instruments":[],"seq_in_adapter":"ATAAGCTA"}],
                "RPI11":[{"seq":"GGCTACAT","instruments":[],"seq_in_adapter":"ATGTAGCC"}],
                "RPI12":[{"seq":"CTTGTAAT","instruments":[],"seq_in_adapter":"ATTACAAG"}],
                "RPI13":[{"seq":"AGTCAAAT","instruments":[],"seq_in_adapter":"ATTTGACT"}],
                "RPI14":[{"seq":"AGTTCCAT","instruments":[],"seq_in_adapter":"ATGGAACT"}],
                "RPI15":[{"seq":"ATGTCAAT","instruments":[],"seq_in_adapter":"ATTGACAT"}],
                "RPI16":[{"seq":"CCGTCCAT","instruments":[],"seq_in_adapter":"ATGGACGG"}],
                "RPI17":[{"seq":"GTAGAGAT","instruments":[],"seq_in_adapter":"ATCTCTAC"}],
                "RPI18":[{"seq":"GTCCGCAT","instruments":[],"seq_in_adapter":"ATGCGGAC"}],
                "RPI19":[{"seq":"GTGAAAAT","instruments":[],"seq_in_adapter":"ATTTTCAC"}],
                "RPI20":[{"seq":"GTGGCCAT","instruments":[],"seq_in_adapter":"ATGGCCAC"}],
                "RPI21":[{"seq":"GTTTCGAT","instruments":[],"seq_in_adapter":"ATCGAAAC"}],
                "RPI22":[{"seq":"CGTACGAT","instruments":[],"seq_in_adapter":"ATCGTACG"}],
                "RPI23":[{"seq":"GAGTGGAT","instruments":[],"seq_in_adapter":"ATCCACTC"}],
                "RPI24":[{"seq":"GGTAGCAT","instruments":[],"seq_in_adapter":"ATGCTACC"}],
                "RPI25":[{"seq":"ACTGATAT","instruments":[],"seq_in_adapter":"ATATCAGT"}],
                "RPI26":[{"seq":"ATGAGCAT","instruments":[],"seq_in_adapter":"ATGCTCAT"}],
                "RPI27":[{"seq":"ATTCCTAT","instruments":[],"seq_in_adapter":"ATAGGAAT"}],
                "RPI28":[{"seq":"CAAAAGAT","instruments":[],"seq_in_adapter":"ATCTTTTG"}],
                "RPI29":[{"seq":"CAACTAAT","instruments":[],"seq_in_adapter":"ATTAGTTG"}],
                "RPI30":[{"seq":"CACCGGAT","instruments":[],"seq_in_adapter":"ATCCGGTG"}],
                "RPI31":[{"seq":"CACGATAT","instruments":[],"seq_in_adapter":"ATATCGTG"}],
                "RPI32":[{"seq":"CACTCAAT","instruments":[],"seq_in_adapter":"ATTGAGTG"}],
                "RPI33":[{"seq":"CAGGCGAT","instruments":[],"seq_in_adapter":"ATCGCCTG"}],
                "RPI34":[{"seq":"CATGGCAT","instruments":[],"seq_in_adapter":"ATGCCATG"}],
                "RPI35":[{"seq":"CATTTTAT","instruments":[],"seq_in_adapter":"ATAAAATG"}],
                "RPI36":[{"seq":"CCAACAAT","instruments":[],"seq_in_adapter":"ATTGTTGG"}],
                "RPI37":[{"seq":"CGGAATAT","instruments":[],"seq_in_adapter":"ATATTCCG"}],
                "RPI38":[{"seq":"CTAGCTAT","instruments":[],"seq_in_adapter":"ATAGCTAG"}],
                "RPI39":[{"seq":"CTATACAT","instruments":[],"seq_in_adapter":"ATGTATAG"}],
                "RPI40":[{"seq":"CTCAGAAT","instruments":[],"seq_in_adapter":"ATTCTGAG"}],
                "RPI41":[{"seq":"GACGACAT","instruments":[],"seq_in_adapter":"ATGTCGTC"}],
                "RPI42":[{"seq":"TAATCGAT","instruments":[],"seq_in_adapter":"ATCGATTA"}],
                "RPI43":[{"seq":"TACAGCAT","instruments":[],"seq_in_adapter":"ATGCTGTA"}],
                "RPI44":[{"seq":"TATAATAT","instruments":[],"seq_in_adapter":"ATATTATA"}],
                "RPI45":[{"seq":"TCATTCAT","instruments":[],"seq_in_adapter":"ATGAATGA"}],
                "RPI46":[{"seq":"TCCCGAAT","instruments":[],"seq_in_adapter":"ATTCGGGA"}],
                "RPI47":[{"seq":"TCGAAGAT","instruments":[],"seq_in_adapter":"ATCTTCGA"}],
                "RPI48":[{"seq":"TCGGCAAT","instruments":[],"seq_in_adapter":"ATTGCCGA"}],
            },
            "i5":{
                "RP1":[{"seq":"GTTCAGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTCAGAG"}, {"seq":"CTCTGAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTCAGAG"}],
                "RP1-v4":[{"seq":"AAGTGCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGTGCTA"}, {"seq":"CTCTGAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGTGCTA"}],
                "RP2-v4":[{"seq":"CTCGAAGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCGAAGC"}, {"seq":"CTCTGAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCGAAGC"}],
                "RP3-v4":[{"seq":"TCACTGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCACTGCT"}, {"seq":"CTCTGAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCACTGCT"}],
                "RP4-v4":[{"seq":"GGTACTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GGTACTAG"}, {"seq":"CTCTGAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GGTACTAG"}],
            }
        },
    }

    @classmethod
    @memoize
    def neighbors(cls, seq, distance=1):
        ''' Returns a list Hamming-distance of 1 from the input seq'''
        neighbor_list = []

        seq = seq.rstrip()
        for i in range(0, len(seq.rstrip())):
            for j in ["A", "C", "G", "T", "N"]:
                if j == seq[i].upper():
                    continue    #  Do not include the seq itself
                q = seq[0:i] + j + seq[i+1:]
                neighbor_list.append(q)

        if distance==1:
            return sorted(neighbor_list)
        else:
            return sorted(list(set([item for sublist in (cls.neighbors(s, distance=distance-1) for s in neighbor_list) for item in sublist if item!=seq])))

    @staticmethod
    def reverse_complement(seq):
        """
            Returns the reverse complement using string.maketrans
        """
        table = bytearray.maketrans(b"ACTGN",b"TGACN")
        return bytearray(seq.encode("UTF8")).translate(table)[::-1].decode("UTF8")

    @property
    @memoize
    def _barcodes_meta_all(self):
        barcodes = []
        for kit_name,kit_value in self._kits.items():
            for item_key,item_value in kit_value.items():

                if type(item_value) == dict:
                    for index_name,index_value in item_value.items():
                        for barcode_meta in index_value:
                            barcodes.append(barcode_meta)
        return barcodes

    @classmethod
    @memoize
    def kits(cls):
        return sorted(list(cls._kits.keys()))

    @property
    def instruments(self):
        instruments = set()
        for barcode_meta in self._barcodes_meta_all:
            ins_for_barcode = barcode_meta.get("instruments","")
            if len(ins_for_barcode):
            
                instruments |= set(ins_for_barcode)
        return sorted(list(instruments))

    @memoize
    def index_for_seq(self, seq, kit=None, instrument=None):
        # use kit/instrument passed in if present
        # value forinstrument, if present, if neither return list of options
        possible_indices = set()
        for kit_name,kit_value in self._kits.items():
            if not kit or (kit_name == kit) or (kit_value["_short_name"]==kit):
                for item_key,item_value in kit_value.items():
                    if type(item_value) == dict:
                        for index_name,index_value in item_value.items():
                            for barcode_meta in index_value:
                                if barcode_meta["seq"] == seq:
                                    if not instrument or instrument in barcode_meta["instruments"]:
                                        possible_indices |= set([index_name])
        return sorted(list(possible_indices))

    @memoize
    def seq_for_index(self, index, kit=None, instrument=None):
        # use kit/instrument passed in if present
        # value forinstrument, if present, if neither return list of options
        possible_seqs = set()
        for kit_name,kit_value in self._kits.items():
            if not kit or (kit_name == kit) or (kit_value["_short_name"]==kit):
                for item_key,item_value in kit_value.items():
                    if type(item_value) == dict:
                        for index_name,index_value in item_value.items():
                            matches = re.match("(?P<index>"+index_name+")", index)
                            if matches:
                                for barcode_meta in index_value:
                                    if not instrument or instrument in barcode_meta["instruments"]:
                                        possible_seqs |= set([barcode_meta["seq"]])
        return sorted(list(possible_seqs))

    @memoize
    def guess_index(self, seq, distance=1, kit=None, instrument=None):
        # using neighbors()...
        # use kit/instrument if present
        possible_indices = set()
        # first include exact matches
        exact_matches = self.index_for_seq(seq, kit=kit, instrument=instrument)
        if len(exact_matches):
            possible_indices |= set(exact_matches)
        else:
            # then include barcodes _distance_ away
            for neighbor_seq in self.neighbors(seq, distance=distance):
                possible_indices |= set(self.index_for_seq(neighbor_seq, kit=kit, instrument=instrument))
        return sorted(list(possible_indices))

class UncertainSamplesheetError(Exception):
    pass

class IlluminaBarcodeHelper(object):
    index_reference = IlluminaIndexReference()
    def __init__(self, barcode_counts, picard_metrics, sample_name, rows_limit=1000):
        # (barcode1,barcode2): count
        self.barcodes_seen = OrderedDict()
        # barcode: illumina index name (N507, etc.)
        self.barcode_name_map = OrderedDict()
        # sample_name: (barcode1,barcode2)
        self.sample_to_barcodes = OrderedDict()
        # list of all sample names (barcode_name from Picard metrics file)
        self.samples = []
        # sample_name: count
        self.sample_to_read_counts = OrderedDict()
        # unresolved barcode pair count
        self.unassigned_read_count = 0

        # read Picard demux metrics file
        for row in util.file.read_tabfile_dict(picard_metrics, skip_prefix="#"):
            barcodes = tuple(row["BARCODE"].split("-"))
            if "BARCODE_NAME" in row:
                self.sample_to_barcodes[row["BARCODE_NAME"]] = barcodes
                self.samples.append(row["BARCODE_NAME"])
                self.sample_to_read_counts[row["BARCODE_NAME"]] = int(row["READS"])
            elif all(re.match(r'^N+$',barcode) for barcode in barcodes):
                self.unassigned_read_count = int(row["READS"])
                continue

        # read barcodes seen in the file, in the format:
        #Barcode1   Likely_Index_Names1 Barcode2    Likely_Index_Names2 Count
        #CTCTCTAC   N707    AAGGAGTA    S507,[N|S|E]507 40324834
        for row in util.file.read_tabfile_dict(barcode_counts, rowcount_limit=rows_limit):
            if (row["Barcode1"],row.get("Barcode2",None)) not in self.barcodes_seen:
                self.barcodes_seen[(row["Barcode1"],row.get("Barcode2",None))] = int(row["Count"])
            self.barcode_name_map[row["Barcode1"]] = row["Likely_Index_Names1"]
            if "Barcode2" in row and row["Barcode2"]:
                self.barcode_name_map[row["Barcode2"]] = row["Likely_Index_Names2"]

    def outlier_barcodes(self, outlier_threshold=0.775, expected_assigned_fraction=0.7, number_of_negative_controls=1):
        """
            This identifies samples listed in the Picard metrics (derived from 
            the sample sheet) that are believed to have an anomalously low 
            number of reads after demultiplexing, based on the following assumptions:
             - The pool of samples will have (at least) one negative control that should yield zero reads
             - The total number of reads should be evenly distributed among the rest of the samples
               since the nucleic acid concentration should be roughly balanced in the pool
             - The variation among what should be evenly distributed is gaussian
            

            A warning is raised if <70% of the total reads are assigned (likely indicating
            multiple problems with the samplesheet).

            Loosely based on the method for outlier detection specified in
                Boris Iglewicz and David Hoaglin (1993), 
                "Volume 16: How to Detect and Handle Outliers", 
                The ASQC Basic References in Quality Control: Statistical Techniques

            Params:
            outlier_threshold (float): 0.675, corresponds to 75th percentile
            expected_assigned_fraction (float): fraction of reads assigned to samples
            number_of_negative_controls (int): the number of samples in the pool expected to have zero reads
        """
        ##log.debug(f"outlier_threshold {outlier_threshold}")
        assigned_read_count = sum(self.sample_to_read_counts.values())
        total_read_count = assigned_read_count+self.unassigned_read_count
        fraction_assigned = float(assigned_read_count)/float(total_read_count)
        log.info("fraction_assigned %s", fraction_assigned)
        if fraction_assigned < expected_assigned_fraction:
            raise UncertainSamplesheetError("Only {:.0%} of reads were assigned to barcode pairs listed in the samplesheet. Check the sample sheet for errors.".format(fraction_assigned))

        num_samples = len(self.sample_to_read_counts)
        assert number_of_negative_controls<num_samples, "number_of_negative_controls must be < num_samples"

        log_obs_fractions_of_pool = [ -math.log(float(x)/float(total_read_count),10) if x>0 
                                        else 0
                                        for x in 
                                        list(self.sample_to_read_counts.values())
                                    ]
        ##log.debug(f"log_obs_fractions_of_pool {log_obs_fractions_of_pool}")
        log_exp_fractions_of_pool = [-math.log(1.0/float(num_samples-number_of_negative_controls),10)]*num_samples
        ##log.warning(f"log_exp_fractions_of_pool {log_exp_fractions_of_pool}")
        residuals = [obs-exp for (obs,exp) in zip(log_obs_fractions_of_pool,log_exp_fractions_of_pool)]
        resid_stdev = self.stddevp(residuals) #essentially RMSE
        resid_mean = self.mean(residuals) # mean error
        resid_median = self.median(residuals) # median error

        ##log.debug(f"residuals {residuals}")
        ##log.debug(f"resid_stdev {resid_stdev}")
        ##log.debug(f"resid_mean {resid_mean}")
        ##log.debug(f"resid_median {resid_median}")

        # modifed zscore using median to reduce influence of outliers
        zscores_residual_relative_to_median = [float(1.0 * (x-resid_median))/resid_stdev for x in residuals]
        ##log.debug(f"zscores_residual_relative_to_median {zscores_residual_relative_to_median}")
        # only consider LOW variance
        indices_of_outliers = [i for i,v in enumerate(zscores_residual_relative_to_median) if v > outlier_threshold]
        indices_of_outliers += [i for i,v in enumerate(list(self.sample_to_read_counts.values())) if v==0]
        ##log.debug("self.sample_to_read_counts.keys() {}".format(self.sample_to_read_counts.keys()))
        names_of_outlier_samples = [(list(self.sample_to_read_counts.keys()))[i] for i in indices_of_outliers]
        log.warning("outlier samples")
        for s in names_of_outlier_samples:
            log.warning("\t%s (%s reads)", s,self.sample_to_read_counts[s])

        return names_of_outlier_samples

    @classmethod
    def mean(cls, nums):
        return sum(nums) / len(nums)

    @classmethod
    def stddevp(cls, nums):
        """population standard deviation of nums"""
        mn = cls.mean(nums)
        variance = sum([(e-mn)**2 for e in nums]) / len(nums)
        return math.sqrt(variance)

    @classmethod
    def median(cls, nums):
        length = len(nums)
        if length < 1:
            return None
        if length % 2 == 1:
            # if this is a list with an odd number of elements, 
            # return the middle element from a sorted version of the list
            return sorted(nums)[length//2]
        else:
            # otherwise if the list has an even number of elements
            # return the interpolated mean of the middle two numbers
            return sum(sorted(nums)[length//2-1:length//2+1])/2.0


    def guess_barcodes_for_sample(self, sample_name):
        """
            This guesses the barcode value for a sample name,
            based on the following:
             - a list is made of novel barcode pairs seen in the data, but not in the picard metrics
             - for the sample in question, get the most abundant novel barcode pair where one of the 
               barcodes seen in the data matches one of the barcodes in the picard metrics (partial match)
             - if there are no partial matches, get the most abundant novel barcode pair 

            Limitations:
             - If multiple samples share a barcode with multiple novel barcodes, disentangling them
               is difficult or impossible
        """
        # output_header = ["sample_name",
        #                     "expected_barcode_1","expected_barcode_1_name",
        #                     "expected_barcode_2","expected_barcode_2_name"
        #                     "expected_barcodes_read_count",
        #                     "guessed_barcode_1","guessed_barcode_1_name",
        #                     "guessed_barcode_2","guessed_barcode_2_name"
        #                     "guessed_barcodes_read_count",
        #                     "match_type"
        #                 ]

        out_dict = OrderedDict()
        out_dict["sample_name"] = sample_name
        out_dict["match_type"] = "not_found"
        
        barcodes_seen_novel = copy.deepcopy(self.barcodes_seen)

        # From barcodes seen in data, collect barcode pairs not expected based on sample sheet
        novel_barcode_pairs = []

        for barcode_pair in self.barcodes_seen.keys():
            barcode = tuple((b for b in barcode_pair if b is not None))
            if barcode not in self.sample_to_barcodes.values():
                novel_barcode_pairs.append(barcode_pair)
            else:
                del barcodes_seen_novel[barcode_pair]

        is_dual_index = len(self.sample_to_barcodes[sample_name]) > 1

        out_dict["expected_barcode_1"]           = self.sample_to_barcodes[sample_name][0]
        out_dict["expected_barcode_1_name"]      = ",".join(self.index_reference.guess_index(self.sample_to_barcodes[sample_name][0]))
        if is_dual_index:
            out_dict["expected_barcode_2"]           = self.sample_to_barcodes[sample_name][1]
            out_dict["expected_barcode_2_name"]      = ",".join(self.index_reference.guess_index(self.sample_to_barcodes[sample_name][1]))
        out_dict["expected_barcodes_read_count"] = self.sample_to_read_counts[sample_name]

        claimed_barcodes = self.sample_to_barcodes[sample_name]

        found_partial_match = False
        putative_match = None

        if is_dual_index:
            # barcodes_seen_novel is sorted by read count, desc
            for (barcode_pair,count) in barcodes_seen_novel.items():
                if barcode_pair[0]==claimed_barcodes[0] or barcode_pair[1]==claimed_barcodes[1]:
                    found_partial_match=True
                    putative_match = barcode_pair
                    out_dict["match_type"] = "one_barcode_match"
                    break
        else:
            # barcodes_seen_novel is sorted by read count, desc
            for (barcode_pair,count) in barcodes_seen_novel.items():
                # get barcode with greatest count that isn't the one we're looking for
                if barcode_pair[0]!=claimed_barcodes[0]:
                    found_partial_match=True
                    putative_match = barcode_pair
                    out_dict["match_type"] = "one_barcode_match"
                    break

        # find index of match to help determine if it is a reasonable guess
        idx_of_match = -1
        for (idx,(barcode_pair,count)) in enumerate(self.barcodes_seen.items()):
            if barcode_pair==putative_match:
                idx_of_match=idx
                break

        # if the one-barcode match is too far down the list of barcode pairs seen
        # (farther down than 1.5x the number of samples)
        if not found_partial_match or idx_of_match>(len(self.sample_to_read_counts)*1.5):
            for (barcode_pair,count) in barcodes_seen_novel.items():
                # return the most pair with greatest count
                putative_match = barcode_pair
                out_dict["match_type"] = "high_count_novel_barcode"
                break
        
        out_dict["guessed_barcode_1"]           = "unknown"
        out_dict["guessed_barcode_1_name"]      = "unknown"
        #out_dict["guessed_barcodes_read_count"] = 0
        if is_dual_index and (putative_match is not None and putative_match[1] != None):
            out_dict["guessed_barcode_2"]           = "unknown"
            out_dict["guessed_barcode_2_name"]      = "unknown"
        

        if putative_match is not None and putative_match[0]!=None:
            out_dict["guessed_barcode_1"]           = putative_match[0]
            out_dict["guessed_barcode_1_name"]      = self.barcode_name_map[putative_match[0]]
            if is_dual_index and putative_match[1] != None:
                out_dict["guessed_barcode_2"]           = putative_match[1]
                out_dict["guessed_barcode_2_name"]      = self.barcode_name_map[putative_match[1]]
                out_dict["guessed_barcodes_read_count"] = self.barcodes_seen[(putative_match[0],putative_match[1])]

        return out_dict


    def find_uncertain_barcodes(self, sample_names=None, outlier_threshold=0.675, expected_assigned_fraction=0.7, number_of_negative_controls=1, readcount_threshold=None):
        """
            If sample_names is specified, barcodes for the named samples 
            will be re-examined. 

            If readcount_threshold is specified, this will be used as the
            cutoff above which sample barcodes will not be re-examined.
            If this is not specified, outliers will be detected automatically.
        """
        samples_to_examine = sample_names or []

        guessed_barcodes = []

        if not samples_to_examine:
            if readcount_threshold:
                for sample_name,count in self.sample_to_read_counts.values():
                    if count < readcount_threshold:
                        samples_to_examine.append(sample_name)
            else:
                samples_to_examine = self.outlier_barcodes(outlier_threshold=outlier_threshold, expected_assigned_fraction=expected_assigned_fraction, number_of_negative_controls=number_of_negative_controls)

        for s in samples_to_examine:
            guessed_barcodes.append(self.guess_barcodes_for_sample(s))

        consolidated_guesses = defaultdict(list)

        for row in guessed_barcodes:
            consolidated_guesses[(row["guessed_barcode_1"],row.get("guessed_barcode_2",None))].append(row)

        final_guesses = []

        def clear_guessed_fields(sample,match_reason,fields_to_clear=None):
            fields_to_clear = fields_to_clear or ["guessed_barcode_1", "guessed_barcode_2",
                                                    "guessed_barcode_1_name", "guessed_barcode_2_name",
                                                    "guessed_barcodes_read_count"]
            for field in fields_to_clear:
                if field in sample:
                    sample[field] = None
            sample["match_type"] = match_reason
            return sample

        for barcode_pair,samples in consolidated_guesses.items():
            if len(samples) > 1:
                log.warning("Ambiguous! Multiple samples corresponding to guessed barcodes %s:", barcode_pair)
                for sample in samples:
                    log.warning("\t%s expected (%s,%s); -> Guessed (%s,%s); match type: %s", sample["sample_name"], sample["expected_barcode_1"],sample.get("expected_barcode_2",""),sample["guessed_barcode_1"],sample.get("guessed_barcode_2",""),sample["match_type"])

                    final_guesses.append(clear_guessed_fields(sample, "alternative_indices_uncertain"))
            else:
                for sample in samples:
                    if sample["guessed_barcodes_read_count"] < sample["expected_barcodes_read_count"]:
                        final_guesses.append(clear_guessed_fields(sample, "alternatives_have_lower_read_counts"))
                    else:
                        final_guesses.append(sample)
        return final_guesses

    def write_guessed_barcodes(self, out_tsv, guessed_barcodes):
        possible_header_cols = ["sample_name",
                            "expected_barcode_1","expected_barcode_2",
                            "expected_barcode_1_name","expected_barcode_2_name",
                            "expected_barcodes_read_count",
                            "guessed_barcode_1", "guessed_barcode_2",
                            "guessed_barcode_1_name", "guessed_barcode_2_name",
                            "guessed_barcodes_read_count",
                            "match_type"
                        ]
        output_header_columns = []
        for header_key in possible_header_cols:
            for row in guessed_barcodes:
                if header_key in row.keys() and header_key not in output_header_columns:
                    output_header_columns.append(header_key)

        with open(out_tsv, 'w') as tsvfile:
            csv.register_dialect('dict_tsv', quoting=csv.QUOTE_MINIMAL, delimiter="\t")
            writer = csv.DictWriter(tsvfile, fieldnames=output_header_columns, dialect="dict_tsv")

            writer.writeheader()
            writer.writerows(guessed_barcodes)
