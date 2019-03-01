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
                "A501": [{"seq":"TGAACCTT", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"AAGGTTCA", "instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A502": [{"seq":"TGCTAAGT", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"ACTTAGCA", "instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A503": [{"seq":"TGTTCTCT", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"AGAGAACA", "instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A504": [{"seq":"TAAGACAC", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"GTGTCTTA", "instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A505": [{"seq":"CTAATCGA", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TCGATTAG", "instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A506": [{"seq":"CTAGAACA", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TGTTCTAG", "instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A507": [{"seq":"TAAGTTCC", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"GGAACTTA", "instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A508": [{"seq":"TAGACCTA", "instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TAGGTCTA", "instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}]
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
                "E505": [{"seq":"GTAAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"CTCCTTAC","instruments":["NextSeq","HiSeq 3000","HiSeq 4000"]}]
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
                "E502": [{"seq":"CTCTCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"ATAGAGAG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "E503": [{"seq":"TATCCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"AGAGGATA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "E504": [{"seq":"AGAGTAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TCTACTCT","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "E505": [{"seq":"GTAAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"CTCCTTAC","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}]
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
                "E502": [{"seq":"CTCTCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"ATAGAGAG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "E505": [{"seq":"GTAAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"CTCCTTAC","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "E506": [{"seq":"ACTGCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TATGCAGT","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "E517": [{"seq":"GCGTAAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TCTTACGC","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}]
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
                "A501": [{"seq":"TGAACCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"AAGGTTCA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A502": [{"seq":"TGCTAAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"ACTTAGCA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}]
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
                "[N|S|E]501":[{"seq":"TAGATCGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"TAGATCGC"}, {"seq":"GCGATCTA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGATCGC"}],
                "[N|S|E]502":[{"seq":"CTCTCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"CTCTCTAT"}, {"seq":"ATAGAGAG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCTCTAT"}],
                "[N|S|E]503":[{"seq":"TATCCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"TATCCTCT"}, {"seq":"AGAGGATA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TATCCTCT"}],
                "[N|S|E]504":[{"seq":"AGAGTAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"AGAGTAGA"}, {"seq":"TCTACTCT","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGAGTAGA"}],
                "[N|S|E]505":[{"seq":"GTAAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"GTAAGGAG"}, {"seq":"CTCCTTAC","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTAAGGAG"}],
                "[N|S|E]506":[{"seq":"ACTGCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"ACTGCATA"}, {"seq":"TATGCAGT","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTGCATA"}],
                "[N|S|E]507":[{"seq":"AAGGAGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"AAGGAGTA"}, {"seq":"TACTCCTT","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGAGTA"}],
                "[N|S|E]508":[{"seq":"CTAAGCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"CTAAGCCT"}, {"seq":"AGGCTTAG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTAAGCCT"}],
                "[N|S|E]517":[{"seq":"GCGTAAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"GCGTAAGA"}, {"seq":"TCTTACGC","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCGTAAGA"}],
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
                "S502":[{"seq":"CTCTCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"CTCTCTAT"}, {"seq":"ATAGAGAG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCTCTAT"}],
                "S503":[{"seq":"TATCCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"TATCCTCT"}, {"seq":"AGAGGATA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TATCCTCT"}],
                "S505":[{"seq":"GTAAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"GTAAGGAG"}, {"seq":"CTCCTTAC","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTAAGGAG"}],
                "S506":[{"seq":"ACTGCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"ACTGCATA"}, {"seq":"TATGCAGT","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTGCATA"}],
                "S507":[{"seq":"AAGGAGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"AAGGAGTA"}, {"seq":"TACTCCTT","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGAGTA"}],
                "S508":[{"seq":"CTAAGCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"CTAAGCCT"}, {"seq":"AGGCTTAG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTAAGCCT"}],
                "S510":[{"seq":"CGTCTAAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"CGTCTAAT"}, {"seq":"ATTAGACG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTCTAAT"}],
                "S511":[{"seq":"TCTCTCCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"TCTCTCCG"}, {"seq":"CGGAGAGA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTCTCCG"}],
                "S513":[{"seq":"TCGACTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"TCGACTAG"}, {"seq":"CTAGTCGA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGACTAG"}],
                "S515":[{"seq":"TTCTAGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"TTCTAGCT"}, {"seq":"AGCTAGAA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCTAGCT"}],
                "S516":[{"seq":"CCTAGAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"CCTAGAGT"}, {"seq":"ACTCTAGG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTAGAGT"}],
                "S517":[{"seq":"GCGTAAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"GCGTAAGA"}, {"seq":"TCTTACGC","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCGTAAGA"}],
                "S518":[{"seq":"CTATTAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"CTATTAAG"}, {"seq":"CTTAATAG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTATTAAG"}],
                "S520":[{"seq":"AAGGCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"AAGGCTAT"}, {"seq":"ATAGCCTT","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGCTAT"}],
                "S521":[{"seq":"GAGCCTTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"GAGCCTTA"}, {"seq":"TAAGGCTC","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGCCTTA"}],
                "S522":[{"seq":"TTATGCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"],"seq_in_adapter":"TTATGCGA"}, {"seq":"TCGCATAA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTATGCGA"}],
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
                "A501":[{"seq":"TGAACCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"AAGGTTCA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A502":[{"seq":"TGCTAAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"ACTTAGCA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A503":[{"seq":"TGTTCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"AGAGAACA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A504":[{"seq":"TAAGACAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"GTGTCTTA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A505":[{"seq":"CTAATCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TCGATTAG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A506":[{"seq":"CTAGAACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TGTTCTAG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A507":[{"seq":"TAAGTTCC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"GGAACTTA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A508":[{"seq":"TAGACCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TAGGTCTA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
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
                "D501":[{"seq":"TATAGCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"AGGCTATA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D502":[{"seq":"ATAGAGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"GCCTCTAT","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D503":[{"seq":"CCTATCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"AGGATAGG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D504":[{"seq":"GGCTCTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TCAGAGCC","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D505":[{"seq":"AGGCGAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"CTTCGCCT","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D506":[{"seq":"TAATCTTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TAAGATTA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D507":[{"seq":"CAGGACGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"ACGTCCTG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "D508":[{"seq":"GTACTGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"GTCAGTAC","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
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
                "A501":[{"seq":"TGAACCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"AAGGTTCA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A502":[{"seq":"TGCTAAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"ACTTAGCA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A503":[{"seq":"TGTTCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"AGAGAACA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A504":[{"seq":"TAAGACAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"GTGTCTTA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A505":[{"seq":"CTAATCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TCGATTAG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A506":[{"seq":"CTAGAACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TGTTCTAG","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A507":[{"seq":"TAAGTTCC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"GGAACTTA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
                "A508":[{"seq":"TAGACCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500"]}, {"seq":"TAGGTCTA","instruments":["MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"]}],
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

    def outlier_barcodes(self, outlier_threshold=0.675, expected_assigned_fraction=0.7, number_of_negative_controls=1):
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
            outlier_threshold (float): 0.675 corresponds to 75th percentile
            expected_assigned_fraction (float): fraction of reads assigned to samples
            number_of_negative_controls (int): the number of samples in the pool expected to have zero reads
        """
        assigned_read_count = sum(self.sample_to_read_counts.values())
        total_read_count = assigned_read_count+self.unassigned_read_count
        fraction_assigned = float(assigned_read_count)/float(total_read_count)
        log.info("fraction_assigned %s", fraction_assigned)
        if fraction_assigned < expected_assigned_fraction:
            raise UncertainSamplesheetError("Only {:.0%} of reads were assigned to barcode pairs listed in the samplesheet. Check the sample sheet for errors.".format(fraction_assigned))

        num_samples = len(self.sample_to_read_counts)
        log_obs_fractions_of_pool = [ -math.log(float(x)/float(total_read_count),10) if x>0 
                                        else 0 
                                        for x in 
                                        list(self.sample_to_read_counts.values())+[self.unassigned_read_count]
                                    ]

        log_exp_fractions_of_pool = [-math.log(1.0/float(num_samples-number_of_negative_controls),10)]*num_samples + [0]

        residuals = [obs-exp for (obs,exp) in zip(log_obs_fractions_of_pool,log_exp_fractions_of_pool)]
        resid_stdev = self.stddevp(residuals) #essentially RMSE
        resid_mean = self.mean(residuals) # mean error
        resid_median = self.median(residuals) # median error

        # modifed zscore using median to reduce influence of outliers
        zscores_residual_relative_to_median = [float(1.0 * (x-resid_median))/resid_stdev for x in residuals]
        # only consider LOW variance
        indices_of_outliers = [i for i,v in enumerate(zscores_residual_relative_to_median[:-1]) if v > outlier_threshold]
        indices_of_outliers += [i for i,v in enumerate(list(self.sample_to_read_counts.values())[:-1]) if v==0]
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

        barcodes_seen_novel = copy.deepcopy(self.barcodes_seen)

        # From barcodes seen in data, collect barcode pairs not expected based on sample sheet
        novel_barcode_pairs = []

        for barcode_pair in self.barcodes_seen.keys():
            if barcode_pair not in self.sample_to_barcodes.values():
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
                out_dict["match_type"] = "high_count_novel_pair"
                break
        
        out_dict["guessed_barcode_1"]           = putative_match[0]
        out_dict["guessed_barcode_1_name"]      = self.barcode_name_map[putative_match[0]]
        if is_dual_index:
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
