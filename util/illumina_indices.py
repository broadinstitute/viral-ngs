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
                "BI5A04":[{"seq":"TACACCGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCGC"}, {"seq":"GCGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCGC"}],
                "BI5A05":[{"seq":"TACACCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCTG"}, {"seq":"CAGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCTG"}],
                "BI5A06":[{"seq":"TACACTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTAG"}, {"seq":"CTAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTAG"}],
                "BI5A07":[{"seq":"TACACCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCCT"}, {"seq":"AGGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCCT"}],
                "BI5A08":[{"seq":"TACACCAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCAC"}, {"seq":"GTGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCAC"}],
                "BI5A09":[{"seq":"TACACAAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACAAT"}, {"seq":"ATTGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACAAT"}],
                "BI5A10":[{"seq":"TACACATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACATG"}, {"seq":"CATGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACATG"}],
                "BI5A11":[{"seq":"TACACTGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTGC"}, {"seq":"GCAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTGC"}],
                "BI5A12":[{"seq":"TACACATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACATC"}, {"seq":"GATGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACATC"}],
                "BI5B01":[{"seq":"TACACGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCT"}, {"seq":"AGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCT"}],
                "BI5B02":[{"seq":"TACACTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTAG"}, {"seq":"CTAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTAG"}],
                "BI5B03":[{"seq":"TACACCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCTA"}, {"seq":"TAGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCTA"}],
                "BI5B04":[{"seq":"TACACTCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTCA"}, {"seq":"TGAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTCA"}],
                "BI5B05":[{"seq":"TACACTTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTTA"}, {"seq":"TAAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTTA"}],
                "BI5B06":[{"seq":"TACACGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGAA"}, {"seq":"TTCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGAA"}],
                "BI5B07":[{"seq":"TACACGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGAC"}, {"seq":"GTCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGAC"}],
                "BI5B08":[{"seq":"TACACCAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCAA"}, {"seq":"TTGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCAA"}],
                "BI5B09":[{"seq":"TACACTCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTCA"}, {"seq":"TGAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTCA"}],
                "BI5B10":[{"seq":"TACACTAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTAC"}, {"seq":"GTAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTAC"}],
                "BI5B11":[{"seq":"TACACGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCT"}, {"seq":"AGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCT"}],
                "BI5B12":[{"seq":"TACACTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTTC"}, {"seq":"GAAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTTC"}],
                "BI5C01":[{"seq":"TACACTAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTAC"}, {"seq":"GTAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTAC"}],
                "BI5C02":[{"seq":"TACACAAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACAAT"}, {"seq":"ATTGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACAAT"}],
                "BI5C03":[{"seq":"TACACGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCG"}, {"seq":"CGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCG"}],
                "BI5C04":[{"seq":"TACACGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGTC"}, {"seq":"GACGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGTC"}],
                "BI5C05":[{"seq":"TACACGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCT"}, {"seq":"AGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCT"}],
                "BI5C06":[{"seq":"TACACCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCTG"}, {"seq":"CAGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCTG"}],
                "BI5C07":[{"seq":"TACACCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCTA"}, {"seq":"TAGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCTA"}],
                "BI5C08":[{"seq":"TACACAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACAGG"}, {"seq":"CCTGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACAGG"}],
                "BI5C09":[{"seq":"TACACGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGAA"}, {"seq":"TTCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGAA"}],
                "BI5C10":[{"seq":"TACACATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACATC"}, {"seq":"GATGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACATC"}],
                "BI5C11":[{"seq":"TACACGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGAC"}, {"seq":"GTCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGAC"}],
                "BI5C12":[{"seq":"TACACCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCAT"}, {"seq":"ATGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCAT"}],
                "BI5D01":[{"seq":"TACACTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTGA"}, {"seq":"TCAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTGA"}],
                "BI5D02":[{"seq":"TACACTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTTC"}, {"seq":"GAAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTTC"}],
                "BI5D03":[{"seq":"TACACCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCTC"}, {"seq":"GAGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCTC"}],
                "BI5D04":[{"seq":"TACACCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCTT"}, {"seq":"AAGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCTT"}],
                "BI5D05":[{"seq":"TACACGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGAA"}, {"seq":"TTCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGAA"}],
                "BI5D06":[{"seq":"TACACATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACATA"}, {"seq":"TATGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACATA"}],
                "BI5D07":[{"seq":"TACACTCC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTCC"}, {"seq":"GGAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTCC"}],
                "BI5D08":[{"seq":"TACACTCC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTCC"}, {"seq":"GGAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTCC"}],
                "BI5D09":[{"seq":"TACACCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCTG"}, {"seq":"CAGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCTG"}],
                "BI5D10":[{"seq":"TACACCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCCT"}, {"seq":"AGGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCCT"}],
                "BI5D11":[{"seq":"TACACCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCTG"}, {"seq":"CAGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCTG"}],
                "BI5D12":[{"seq":"TACACCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCGA"}, {"seq":"TCGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCGA"}],
                "BI5E01":[{"seq":"TACACGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCA"}, {"seq":"TGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCA"}],
                "BI5E02":[{"seq":"TACACGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCT"}, {"seq":"AGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCT"}],
                "BI5E03":[{"seq":"TACACAAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACAAT"}, {"seq":"ATTGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACAAT"}],
                "BI5E04":[{"seq":"TACACCCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCCA"}, {"seq":"TGGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCCA"}],
                "BI5E05":[{"seq":"TACACAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACAGT"}, {"seq":"ACTGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACAGT"}],
                "BI5E06":[{"seq":"TACACGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCG"}, {"seq":"CGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCG"}],
                "BI5E07":[{"seq":"TACACCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCTG"}, {"seq":"CAGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCTG"}],
                "BI5E08":[{"seq":"TACACGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGTC"}, {"seq":"GACGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGTC"}],
                "BI5E09":[{"seq":"TACACGCC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCC"}, {"seq":"GGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCC"}],
                "BI5E10":[{"seq":"TACACATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACATC"}, {"seq":"GATGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACATC"}],
                "BI5E11":[{"seq":"TACACTTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTTA"}, {"seq":"TAAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTTA"}],
                "BI5E12":[{"seq":"TACACGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGAC"}, {"seq":"GTCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGAC"}],
                "BI5F01":[{"seq":"TACACTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTTC"}, {"seq":"GAAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTTC"}],
                "BI5F02":[{"seq":"TACACTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTGA"}, {"seq":"TCAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTGA"}],
                "BI5F03":[{"seq":"TACACTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTTG"}, {"seq":"CAAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTTG"}],
                "BI5F04":[{"seq":"TACACTCC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTCC"}, {"seq":"GGAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTCC"}],
                "BI5F05":[{"seq":"TACACCAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCAC"}, {"seq":"GTGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCAC"}],
                "BI5F06":[{"seq":"TACACGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCG"}, {"seq":"CGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCG"}],
                "BI5F07":[{"seq":"TACACCAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCAG"}, {"seq":"CTGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCAG"}],
                "BI5F08":[{"seq":"TACACAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACAGT"}, {"seq":"ACTGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACAGT"}],
                "BI5F09":[{"seq":"TACACCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCGA"}, {"seq":"TCGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCGA"}],
                "BI5F10":[{"seq":"TACACCAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCAG"}, {"seq":"CTGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCAG"}],
                "BI5F11":[{"seq":"TACACCGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCGC"}, {"seq":"GCGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCGC"}],
                "BI5F12":[{"seq":"TACACTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTTC"}, {"seq":"GAAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTTC"}],
                "BI5G01":[{"seq":"TACACTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTAG"}, {"seq":"CTAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTAG"}],
                "BI5G02":[{"seq":"TACACTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTAT"}, {"seq":"ATAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTAT"}],
                "BI5G03":[{"seq":"TACACCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCAT"}, {"seq":"ATGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCAT"}],
                "BI5G04":[{"seq":"TACACAAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACAAT"}, {"seq":"ATTGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACAAT"}],
                "BI5G05":[{"seq":"TACACGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCT"}, {"seq":"AGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCT"}],
                "BI5G06":[{"seq":"TACACAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACAGA"}, {"seq":"TCTGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACAGA"}],
                "BI5G07":[{"seq":"TACACCCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCCA"}, {"seq":"TGGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCCA"}],
                "BI5G08":[{"seq":"TACACCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCTT"}, {"seq":"AAGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCTT"}],
                "BI5G09":[{"seq":"TACACTGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTGC"}, {"seq":"GCAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTGC"}],
                "BI5G10":[{"seq":"TACACGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGGT"}, {"seq":"ACCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGGT"}],
                "BI5G11":[{"seq":"TACACGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCT"}, {"seq":"AGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCT"}],
                "BI5G12":[{"seq":"TACACCCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCCA"}, {"seq":"TGGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCCA"}],
                "BI5H01":[{"seq":"TACACAGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACAGC"}, {"seq":"GCTGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACAGC"}],
                "BI5H02":[{"seq":"TACACGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGTA"}, {"seq":"TACGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGTA"}],
                "BI5H03":[{"seq":"TACACGGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGGA"}, {"seq":"TCCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGGA"}],
                "BI5H04":[{"seq":"TACACACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACACT"}, {"seq":"AGTGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACACT"}],
                "BI5H05":[{"seq":"TACACGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCA"}, {"seq":"TGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCA"}],
                "BI5H06":[{"seq":"TACACCAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCAA"}, {"seq":"TTGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCAA"}],
                "BI5H07":[{"seq":"TACACAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACAAG"}, {"seq":"CTTGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACAAG"}],
                "BI5H08":[{"seq":"TACACCCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCCA"}, {"seq":"TGGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCCA"}],
                "BI5H09":[{"seq":"TACACCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCTT"}, {"seq":"AAGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCTT"}],
                "BI5H10":[{"seq":"TACACTCC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTCC"}, {"seq":"GGAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTCC"}],
                "BI5H11":[{"seq":"TACACTCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACTCA"}, {"seq":"TGAGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACTCA"}],
                "BI5H12":[{"seq":"TACACCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACCTA"}, {"seq":"TAGGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACCTA"}],
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
                "X7001":[{"seq":"CTGATCGT","instruments":[],"seq_in_adapter":"CTGATCGT"}],
                "X7002":[{"seq":"ACTCTCGA","instruments":[],"seq_in_adapter":"ACTCTCGA"}],
                "X7003":[{"seq":"TGAGCTAG","instruments":[],"seq_in_adapter":"TGAGCTAG"}],
                "X7004":[{"seq":"GAGACGAT","instruments":[],"seq_in_adapter":"GAGACGAT"}],
                "X7005":[{"seq":"CTTGTCGA","instruments":[],"seq_in_adapter":"CTTGTCGA"}],
                "X7006":[{"seq":"TTCCAAGG","instruments":[],"seq_in_adapter":"TTCCAAGG"}],
                "X7007":[{"seq":"CGCATGAT","instruments":[],"seq_in_adapter":"CGCATGAT"}],
                "X7008":[{"seq":"ACGGAACA","instruments":[],"seq_in_adapter":"ACGGAACA"}],
                "X7009":[{"seq":"CGGCTAAT","instruments":[],"seq_in_adapter":"CGGCTAAT"}],
                "X7010":[{"seq":"ATCGATCG","instruments":[],"seq_in_adapter":"ATCGATCG"}],
                "X7011":[{"seq":"GCAAGATC","instruments":[],"seq_in_adapter":"GCAAGATC"}],
                "X7012":[{"seq":"GCTATCCT","instruments":[],"seq_in_adapter":"GCTATCCT"}],
                "X7013":[{"seq":"TACGCTAC","instruments":[],"seq_in_adapter":"TACGCTAC"}],
                "X7014":[{"seq":"TGGACTCT","instruments":[],"seq_in_adapter":"TGGACTCT"}],
                "X7015":[{"seq":"AGAGTAGC","instruments":[],"seq_in_adapter":"AGAGTAGC"}],
                "X7016":[{"seq":"ATCCAGAG","instruments":[],"seq_in_adapter":"ATCCAGAG"}],
                "X7017":[{"seq":"GACGATCT","instruments":[],"seq_in_adapter":"GACGATCT"}],
                "X7018":[{"seq":"AACTGAGC","instruments":[],"seq_in_adapter":"AACTGAGC"}],
                "X7019":[{"seq":"CTTAGGAC","instruments":[],"seq_in_adapter":"CTTAGGAC"}],
                "X7020":[{"seq":"GTGCCATA","instruments":[],"seq_in_adapter":"GTGCCATA"}],
                "X7021":[{"seq":"GAATCCGA","instruments":[],"seq_in_adapter":"GAATCCGA"}],
                "X7022":[{"seq":"TCGCTGTT","instruments":[],"seq_in_adapter":"TCGCTGTT"}],
                "X7023":[{"seq":"TTCGTTGG","instruments":[],"seq_in_adapter":"TTCGTTGG"}],
                "X7024":[{"seq":"AAGCACTG","instruments":[],"seq_in_adapter":"AAGCACTG"}],
                "X7025":[{"seq":"CCTTGATC","instruments":[],"seq_in_adapter":"CCTTGATC"}],
                "X7026":[{"seq":"GTCGAAGA","instruments":[],"seq_in_adapter":"GTCGAAGA"}],
                "X7027":[{"seq":"ACCACGAT","instruments":[],"seq_in_adapter":"ACCACGAT"}],
                "X7028":[{"seq":"GATTACCG","instruments":[],"seq_in_adapter":"GATTACCG"}],
                "X7029":[{"seq":"GCACAACT","instruments":[],"seq_in_adapter":"GCACAACT"}],
                "X7030":[{"seq":"GCGTCATT","instruments":[],"seq_in_adapter":"GCGTCATT"}],
                "X7031":[{"seq":"ATCCGGTA","instruments":[],"seq_in_adapter":"ATCCGGTA"}],
                "X7032":[{"seq":"CGTTGCAA","instruments":[],"seq_in_adapter":"CGTTGCAA"}],
                "X7033":[{"seq":"GTGAAGTG","instruments":[],"seq_in_adapter":"GTGAAGTG"}],
                "X7034":[{"seq":"CATGGCTA","instruments":[],"seq_in_adapter":"CATGGCTA"}],
                "X7035":[{"seq":"ATGCCTGT","instruments":[],"seq_in_adapter":"ATGCCTGT"}],
                "X7036":[{"seq":"CAACACCT","instruments":[],"seq_in_adapter":"CAACACCT"}],
                "X7037":[{"seq":"TGTGACTG","instruments":[],"seq_in_adapter":"TGTGACTG"}],
                "X7038":[{"seq":"GTCATCGA","instruments":[],"seq_in_adapter":"GTCATCGA"}],
                "X7039":[{"seq":"AGCACTTC","instruments":[],"seq_in_adapter":"AGCACTTC"}],
                "X7040":[{"seq":"GAAGGAAG","instruments":[],"seq_in_adapter":"GAAGGAAG"}],
                "X7041":[{"seq":"GTTGTTCG","instruments":[],"seq_in_adapter":"GTTGTTCG"}],
                "X7042":[{"seq":"CGGTTGTT","instruments":[],"seq_in_adapter":"CGGTTGTT"}],
                "X7043":[{"seq":"ACTGAGGT","instruments":[],"seq_in_adapter":"ACTGAGGT"}],
                "X7044":[{"seq":"TGAAGACG","instruments":[],"seq_in_adapter":"TGAAGACG"}],
                "X7045":[{"seq":"GTTACGCA","instruments":[],"seq_in_adapter":"GTTACGCA"}],
                "X7046":[{"seq":"AGCGTGTT","instruments":[],"seq_in_adapter":"AGCGTGTT"}],
                "X7047":[{"seq":"GATCGAGT","instruments":[],"seq_in_adapter":"GATCGAGT"}],
                "X7048":[{"seq":"ACAGCTCA","instruments":[],"seq_in_adapter":"ACAGCTCA"}],
                "X7049":[{"seq":"GAGCAGTA","instruments":[],"seq_in_adapter":"GAGCAGTA"}],
                "X7050":[{"seq":"AGTTCGTC","instruments":[],"seq_in_adapter":"AGTTCGTC"}],
                "X7051":[{"seq":"TTGCGAAG","instruments":[],"seq_in_adapter":"TTGCGAAG"}],
                "X7052":[{"seq":"ATCGCCAT","instruments":[],"seq_in_adapter":"ATCGCCAT"}],
                "X7053":[{"seq":"TGGCATGT","instruments":[],"seq_in_adapter":"TGGCATGT"}],
                "X7054":[{"seq":"CTGTTGAC","instruments":[],"seq_in_adapter":"CTGTTGAC"}],
                "X7055":[{"seq":"CATACCAC","instruments":[],"seq_in_adapter":"CATACCAC"}],
                "X7056":[{"seq":"GAAGTTGG","instruments":[],"seq_in_adapter":"GAAGTTGG"}],
                "X7057":[{"seq":"ATGACGTC","instruments":[],"seq_in_adapter":"ATGACGTC"}],
                "X7058":[{"seq":"TTGGACGT","instruments":[],"seq_in_adapter":"TTGGACGT"}],
                "X7059":[{"seq":"AGTGGATC","instruments":[],"seq_in_adapter":"AGTGGATC"}],
                "X7060":[{"seq":"GATAGGCT","instruments":[],"seq_in_adapter":"GATAGGCT"}],
                "X7061":[{"seq":"TGGTAGCT","instruments":[],"seq_in_adapter":"TGGTAGCT"}],
                "X7062":[{"seq":"CGCAATCT","instruments":[],"seq_in_adapter":"CGCAATCT"}],
                "X7063":[{"seq":"GATGTGTG","instruments":[],"seq_in_adapter":"GATGTGTG"}],
                "X7064":[{"seq":"GATTGCTC","instruments":[],"seq_in_adapter":"GATTGCTC"}],
                "X7065":[{"seq":"CGCTCTAT","instruments":[],"seq_in_adapter":"CGCTCTAT"}],
                "X7066":[{"seq":"TATCGGTC","instruments":[],"seq_in_adapter":"TATCGGTC"}],
                "X7067":[{"seq":"AACGTCTG","instruments":[],"seq_in_adapter":"AACGTCTG"}],
                "X7068":[{"seq":"ACGTTCAG","instruments":[],"seq_in_adapter":"ACGTTCAG"}],
                "X7069":[{"seq":"CAGTCCAA","instruments":[],"seq_in_adapter":"CAGTCCAA"}],
                "X7070":[{"seq":"TTGCAGAC","instruments":[],"seq_in_adapter":"TTGCAGAC"}],
                "X7071":[{"seq":"CAATGTGG","instruments":[],"seq_in_adapter":"CAATGTGG"}],
                "X7072":[{"seq":"ACTCCATC","instruments":[],"seq_in_adapter":"ACTCCATC"}],
                "X7073":[{"seq":"GTTGACCT","instruments":[],"seq_in_adapter":"GTTGACCT"}],
                "X7074":[{"seq":"CGTGTGTA","instruments":[],"seq_in_adapter":"CGTGTGTA"}],
                "X7075":[{"seq":"ACGACTTG","instruments":[],"seq_in_adapter":"ACGACTTG"}],
                "X7076":[{"seq":"CACTAGCT","instruments":[],"seq_in_adapter":"CACTAGCT"}],
                "X7077":[{"seq":"ACTAGGAG","instruments":[],"seq_in_adapter":"ACTAGGAG"}],
                "X7078":[{"seq":"GTAGGAGT","instruments":[],"seq_in_adapter":"GTAGGAGT"}],
                "X7079":[{"seq":"CCTGATTG","instruments":[],"seq_in_adapter":"CCTGATTG"}],
                "X7080":[{"seq":"ATGCACGA","instruments":[],"seq_in_adapter":"ATGCACGA"}],
                "X7081":[{"seq":"CGACGTTA","instruments":[],"seq_in_adapter":"CGACGTTA"}],
                "X7082":[{"seq":"TACGCCTT","instruments":[],"seq_in_adapter":"TACGCCTT"}],
                "X7083":[{"seq":"CCGTAAGA","instruments":[],"seq_in_adapter":"CCGTAAGA"}],
                "X7084":[{"seq":"ATCACACG","instruments":[],"seq_in_adapter":"ATCACACG"}],
                "X7085":[{"seq":"CACCTGTT","instruments":[],"seq_in_adapter":"CACCTGTT"}],
                "X7086":[{"seq":"CTTCGACT","instruments":[],"seq_in_adapter":"CTTCGACT"}],
                "X7087":[{"seq":"TGCTTCCA","instruments":[],"seq_in_adapter":"TGCTTCCA"}],
                "X7088":[{"seq":"AGAACGAG","instruments":[],"seq_in_adapter":"AGAACGAG"}],
                "X7089":[{"seq":"GTTCTCGT","instruments":[],"seq_in_adapter":"GTTCTCGT"}],
                "X7090":[{"seq":"TCAGGCTT","instruments":[],"seq_in_adapter":"TCAGGCTT"}],
                "X7091":[{"seq":"CCTTGTAG","instruments":[],"seq_in_adapter":"CCTTGTAG"}],
                "X7092":[{"seq":"GAACATCG","instruments":[],"seq_in_adapter":"GAACATCG"}],
                "X7093":[{"seq":"TAACCGGT","instruments":[],"seq_in_adapter":"TAACCGGT"}],
                "X7094":[{"seq":"AACCGTTC","instruments":[],"seq_in_adapter":"AACCGTTC"}],
                "X7095":[{"seq":"TGGTACAG","instruments":[],"seq_in_adapter":"TGGTACAG"}],
                "X7096":[{"seq":"ATATGCGC","instruments":[],"seq_in_adapter":"ATATGCGC"}],
                "X7097":[{"seq":"GCCTATCA","instruments":[],"seq_in_adapter":"GCCTATCA"}],
                "X7098":[{"seq":"CTTGGATG","instruments":[],"seq_in_adapter":"CTTGGATG"}],
                "X7099":[{"seq":"AGTCTCAC","instruments":[],"seq_in_adapter":"AGTCTCAC"}],
                "X7100":[{"seq":"CTCATCAG","instruments":[],"seq_in_adapter":"CTCATCAG"}],
                "X7101":[{"seq":"TGTACCGT","instruments":[],"seq_in_adapter":"TGTACCGT"}],
                "X7102":[{"seq":"AAGTCGAG","instruments":[],"seq_in_adapter":"AAGTCGAG"}],
                "X7103":[{"seq":"CACGTTGT","instruments":[],"seq_in_adapter":"CACGTTGT"}],
                "X7104":[{"seq":"TCACAGCA","instruments":[],"seq_in_adapter":"TCACAGCA"}],
                "X7105":[{"seq":"CTACTTGG","instruments":[],"seq_in_adapter":"CTACTTGG"}],
                "X7106":[{"seq":"CCTCAGTT","instruments":[],"seq_in_adapter":"CCTCAGTT"}],
                "X7107":[{"seq":"TCCTACCT","instruments":[],"seq_in_adapter":"TCCTACCT"}],
                "X7108":[{"seq":"ATGGCGAA","instruments":[],"seq_in_adapter":"ATGGCGAA"}],
                "X7109":[{"seq":"CTTACCTG","instruments":[],"seq_in_adapter":"CTTACCTG"}],
                "X7110":[{"seq":"CTCGATAC","instruments":[],"seq_in_adapter":"CTCGATAC"}],
                "X7111":[{"seq":"TCCGTGAA","instruments":[],"seq_in_adapter":"TCCGTGAA"}],
                "X7112":[{"seq":"TAGAGCTC","instruments":[],"seq_in_adapter":"TAGAGCTC"}],
                "X7113":[{"seq":"TGACTGAC","instruments":[],"seq_in_adapter":"TGACTGAC"}],
                "X7114":[{"seq":"TAGACGTG","instruments":[],"seq_in_adapter":"TAGACGTG"}],
                "X7115":[{"seq":"CCGGAATT","instruments":[],"seq_in_adapter":"CCGGAATT"}],
                "X7116":[{"seq":"CTCCTAGA","instruments":[],"seq_in_adapter":"CTCCTAGA"}],
                "X7117":[{"seq":"CAACGGAT","instruments":[],"seq_in_adapter":"CAACGGAT"}],
                "X7118":[{"seq":"TGGCTATC","instruments":[],"seq_in_adapter":"TGGCTATC"}],
                "X7119":[{"seq":"CGGTCATA","instruments":[],"seq_in_adapter":"CGGTCATA"}],
                "X7120":[{"seq":"TCCAATCG","instruments":[],"seq_in_adapter":"TCCAATCG"}],
                "X7121":[{"seq":"GAGCTTGT","instruments":[],"seq_in_adapter":"GAGCTTGT"}],
                "X7122":[{"seq":"GAAGGTTC","instruments":[],"seq_in_adapter":"GAAGGTTC"}],
                "X7123":[{"seq":"ATCTCGCT","instruments":[],"seq_in_adapter":"ATCTCGCT"}],
                "X7124":[{"seq":"AGTTACGG","instruments":[],"seq_in_adapter":"AGTTACGG"}],
                "X7125":[{"seq":"GTGTCTGA","instruments":[],"seq_in_adapter":"GTGTCTGA"}],
                "X7126":[{"seq":"TGACTTCG","instruments":[],"seq_in_adapter":"TGACTTCG"}],
                "X7127":[{"seq":"TGGATCAC","instruments":[],"seq_in_adapter":"TGGATCAC"}],
                "X7128":[{"seq":"ACACCAGT","instruments":[],"seq_in_adapter":"ACACCAGT"}],
                "X7129":[{"seq":"CAGGTTAG","instruments":[],"seq_in_adapter":"CAGGTTAG"}],
                "X7130":[{"seq":"AGTTGGCT","instruments":[],"seq_in_adapter":"AGTTGGCT"}],
                "X7131":[{"seq":"TCAACTGG","instruments":[],"seq_in_adapter":"TCAACTGG"}],
                "X7132":[{"seq":"CTGCACTT","instruments":[],"seq_in_adapter":"CTGCACTT"}],
                "X7133":[{"seq":"ACACGGTT","instruments":[],"seq_in_adapter":"ACACGGTT"}],
                "X7134":[{"seq":"AATACGCG","instruments":[],"seq_in_adapter":"AATACGCG"}],
                "X7135":[{"seq":"TGCGAACT","instruments":[],"seq_in_adapter":"TGCGAACT"}],
                "X7136":[{"seq":"GCTGACTA","instruments":[],"seq_in_adapter":"GCTGACTA"}],
                "X7137":[{"seq":"GTGGTGTT","instruments":[],"seq_in_adapter":"GTGGTGTT"}],
                "X7138":[{"seq":"GTGCTTAC","instruments":[],"seq_in_adapter":"GTGCTTAC"}],
                "X7139":[{"seq":"TCAAGGAC","instruments":[],"seq_in_adapter":"TCAAGGAC"}],
                "X7140":[{"seq":"TGAACCTG","instruments":[],"seq_in_adapter":"TGAACCTG"}],
                "X7141":[{"seq":"AGTGTTGG","instruments":[],"seq_in_adapter":"AGTGTTGG"}],
                "X7142":[{"seq":"GTACTCTC","instruments":[],"seq_in_adapter":"GTACTCTC"}],
                "X7143":[{"seq":"CCGTATCT","instruments":[],"seq_in_adapter":"CCGTATCT"}],
                "X7144":[{"seq":"CGAAGAAC","instruments":[],"seq_in_adapter":"CGAAGAAC"}],
                "X7145":[{"seq":"AGCGGAAT","instruments":[],"seq_in_adapter":"AGCGGAAT"}],
                "X7146":[{"seq":"GTGAGCTT","instruments":[],"seq_in_adapter":"GTGAGCTT"}],
                "X7147":[{"seq":"CGTGATCA","instruments":[],"seq_in_adapter":"CGTGATCA"}],
                "X7148":[{"seq":"TCGCATTG","instruments":[],"seq_in_adapter":"TCGCATTG"}],
                "X7149":[{"seq":"TGACGCAT","instruments":[],"seq_in_adapter":"TGACGCAT"}],
                "X7150":[{"seq":"CCGATGTA","instruments":[],"seq_in_adapter":"CCGATGTA"}],
                "X7151":[{"seq":"TTCGCAGT","instruments":[],"seq_in_adapter":"TTCGCAGT"}],
                "X7152":[{"seq":"ACGACAGA","instruments":[],"seq_in_adapter":"ACGACAGA"}],
                "X7153":[{"seq":"AGCTTGAG","instruments":[],"seq_in_adapter":"AGCTTGAG"}],
                "X7154":[{"seq":"GAGTGGTT","instruments":[],"seq_in_adapter":"GAGTGGTT"}],
                "X7155":[{"seq":"GCTGTAAG","instruments":[],"seq_in_adapter":"GCTGTAAG"}],
                "X7156":[{"seq":"CCAAGACT","instruments":[],"seq_in_adapter":"CCAAGACT"}],
                "X7157":[{"seq":"ATTGCGTG","instruments":[],"seq_in_adapter":"ATTGCGTG"}],
                "X7158":[{"seq":"CTGAAGCT","instruments":[],"seq_in_adapter":"CTGAAGCT"}],
                "X7159":[{"seq":"TAACGAGG","instruments":[],"seq_in_adapter":"TAACGAGG"}],
                "X7160":[{"seq":"TCGTCTCA","instruments":[],"seq_in_adapter":"TCGTCTCA"}],
                "X7161":[{"seq":"TTCCTGTG","instruments":[],"seq_in_adapter":"TTCCTGTG"}],
                "X7162":[{"seq":"CGTTGAGT","instruments":[],"seq_in_adapter":"CGTTGAGT"}],
                "X7163":[{"seq":"AGTCGCTT","instruments":[],"seq_in_adapter":"AGTCGCTT"}],
                "X7164":[{"seq":"TAGGTAGG","instruments":[],"seq_in_adapter":"TAGGTAGG"}],
                "X7165":[{"seq":"CAGGAGAT","instruments":[],"seq_in_adapter":"CAGGAGAT"}],
                "X7166":[{"seq":"CATCGTGA","instruments":[],"seq_in_adapter":"CATCGTGA"}],
                "X7167":[{"seq":"TGTTGTGG","instruments":[],"seq_in_adapter":"TGTTGTGG"}],
                "X7168":[{"seq":"ACAGACCT","instruments":[],"seq_in_adapter":"ACAGACCT"}],
                "X7169":[{"seq":"GTCCTTCT","instruments":[],"seq_in_adapter":"GTCCTTCT"}],
                "X7170":[{"seq":"TGATACGC","instruments":[],"seq_in_adapter":"TGATACGC"}],
                "X7171":[{"seq":"CTGTGTTG","instruments":[],"seq_in_adapter":"CTGTGTTG"}],
                "X7172":[{"seq":"AACGTGGA","instruments":[],"seq_in_adapter":"AACGTGGA"}],
                "X7173":[{"seq":"GTTGCGAT","instruments":[],"seq_in_adapter":"GTTGCGAT"}],
                "X7174":[{"seq":"AACGACGT","instruments":[],"seq_in_adapter":"AACGACGT"}],
                "X7175":[{"seq":"CGTATTCG","instruments":[],"seq_in_adapter":"CGTATTCG"}],
                "X7176":[{"seq":"AGCAAGCA","instruments":[],"seq_in_adapter":"AGCAAGCA"}],
                "X7177":[{"seq":"TGTTCGAG","instruments":[],"seq_in_adapter":"TGTTCGAG"}],
                "X7178":[{"seq":"CTCCATGT","instruments":[],"seq_in_adapter":"CTCCATGT"}],
                "X7179":[{"seq":"CGTCTTGT","instruments":[],"seq_in_adapter":"CGTCTTGT"}],
                "X7180":[{"seq":"ATAAGGCG","instruments":[],"seq_in_adapter":"ATAAGGCG"}],
                "X7181":[{"seq":"TGTCTGCT","instruments":[],"seq_in_adapter":"TGTCTGCT"}],
                "X7182":[{"seq":"CGCTTAAC","instruments":[],"seq_in_adapter":"CGCTTAAC"}],
                "X7183":[{"seq":"GATCCATG","instruments":[],"seq_in_adapter":"GATCCATG"}],
                "X7184":[{"seq":"ACCTCTGT","instruments":[],"seq_in_adapter":"ACCTCTGT"}],
                "X7185":[{"seq":"GCCACTTA","instruments":[],"seq_in_adapter":"GCCACTTA"}],
                "X7186":[{"seq":"ACCTGACT","instruments":[],"seq_in_adapter":"ACCTGACT"}],
                "X7187":[{"seq":"GTTAAGGC","instruments":[],"seq_in_adapter":"GTTAAGGC"}],
                "X7188":[{"seq":"ATGCCAAC","instruments":[],"seq_in_adapter":"ATGCCAAC"}],
                "X7189":[{"seq":"AGAGGTTG","instruments":[],"seq_in_adapter":"AGAGGTTG"}],
                "X7190":[{"seq":"ACCATCCA","instruments":[],"seq_in_adapter":"ACCATCCA"}],
                "X7191":[{"seq":"GTGGATAG","instruments":[],"seq_in_adapter":"GTGGATAG"}],
                "X7192":[{"seq":"CTGAGATC","instruments":[],"seq_in_adapter":"CTGAGATC"}],
                "X7193":[{"seq":"CTTCGTTC","instruments":[],"seq_in_adapter":"CTTCGTTC"}],
                "X7194":[{"seq":"GTCTAGGT","instruments":[],"seq_in_adapter":"GTCTAGGT"}],
                "X7195":[{"seq":"ACGTCGTA","instruments":[],"seq_in_adapter":"ACGTCGTA"}],
                "X7196":[{"seq":"GAGCTCAA","instruments":[],"seq_in_adapter":"GAGCTCAA"}],
                "X7197":[{"seq":"CGTGTACT","instruments":[],"seq_in_adapter":"CGTGTACT"}],
                "X7198":[{"seq":"CACTGACA","instruments":[],"seq_in_adapter":"CACTGACA"}],
                "X7199":[{"seq":"TCGTAGTC","instruments":[],"seq_in_adapter":"TCGTAGTC"}],
                "X7200":[{"seq":"GCACGTAA","instruments":[],"seq_in_adapter":"GCACGTAA"}],
                "X7201":[{"seq":"CAAGCAGT","instruments":[],"seq_in_adapter":"CAAGCAGT"}],
                "X7202":[{"seq":"ACATAGGC","instruments":[],"seq_in_adapter":"ACATAGGC"}],
                "X7203":[{"seq":"TGTGGTAC","instruments":[],"seq_in_adapter":"TGTGGTAC"}],
                "X7204":[{"seq":"CACCACTA","instruments":[],"seq_in_adapter":"CACCACTA"}],
                "X7205":[{"seq":"CTGCGTAT","instruments":[],"seq_in_adapter":"CTGCGTAT"}],
                "X7206":[{"seq":"ACGGTCTT","instruments":[],"seq_in_adapter":"ACGGTCTT"}],
                "X7207":[{"seq":"GATTGGAG","instruments":[],"seq_in_adapter":"GATTGGAG"}],
                "X7208":[{"seq":"TGTCCAGA","instruments":[],"seq_in_adapter":"TGTCCAGA"}],
                "X7209":[{"seq":"CCAGTGTT","instruments":[],"seq_in_adapter":"CCAGTGTT"}],
                "X7210":[{"seq":"TGCACCAA","instruments":[],"seq_in_adapter":"TGCACCAA"}],
                "X7211":[{"seq":"TTGACAGG","instruments":[],"seq_in_adapter":"TTGACAGG"}],
                "X7212":[{"seq":"AGGCATAG","instruments":[],"seq_in_adapter":"AGGCATAG"}],
                "X7213":[{"seq":"TAGCCGAA","instruments":[],"seq_in_adapter":"TAGCCGAA"}],
                "X7214":[{"seq":"TTGTCGGT","instruments":[],"seq_in_adapter":"TTGTCGGT"}],
                "X7215":[{"seq":"CATCTACG","instruments":[],"seq_in_adapter":"CATCTACG"}],
                "X7216":[{"seq":"GCATACAG","instruments":[],"seq_in_adapter":"GCATACAG"}],
                "X7217":[{"seq":"ACAGCAAC","instruments":[],"seq_in_adapter":"ACAGCAAC"}],
                "X7218":[{"seq":"CTGGTTCT","instruments":[],"seq_in_adapter":"CTGGTTCT"}],
                "X7219":[{"seq":"TCGACATC","instruments":[],"seq_in_adapter":"TCGACATC"}],
                "X7220":[{"seq":"AACCTCCT","instruments":[],"seq_in_adapter":"AACCTCCT"}],
                "X7221":[{"seq":"CAGCGATT","instruments":[],"seq_in_adapter":"CAGCGATT"}],
                "X7222":[{"seq":"AGGTCACT","instruments":[],"seq_in_adapter":"AGGTCACT"}],
                "X7223":[{"seq":"GCAATTCG","instruments":[],"seq_in_adapter":"GCAATTCG"}],
                "X7224":[{"seq":"GCTTCTTG","instruments":[],"seq_in_adapter":"GCTTCTTG"}],
                "X7225":[{"seq":"AACTGGTG","instruments":[],"seq_in_adapter":"AACTGGTG"}],
                "X7226":[{"seq":"CGGAATAC","instruments":[],"seq_in_adapter":"CGGAATAC"}],
                "X7227":[{"seq":"GCTTCGAA","instruments":[],"seq_in_adapter":"GCTTCGAA"}],
                "X7228":[{"seq":"CAAGGTCT","instruments":[],"seq_in_adapter":"CAAGGTCT"}],
                "X7229":[{"seq":"AACCTTGG","instruments":[],"seq_in_adapter":"AACCTTGG"}],
                "X7230":[{"seq":"CCATACGT","instruments":[],"seq_in_adapter":"CCATACGT"}],
                "X7231":[{"seq":"TGGTCCTT","instruments":[],"seq_in_adapter":"TGGTCCTT"}],
                "X7232":[{"seq":"ACCGCATA","instruments":[],"seq_in_adapter":"ACCGCATA"}],
                "X7233":[{"seq":"CCTTCCTT","instruments":[],"seq_in_adapter":"CCTTCCTT"}],
                "X7234":[{"seq":"TACACGCT","instruments":[],"seq_in_adapter":"TACACGCT"}],
                "X7235":[{"seq":"TGCGTAGA","instruments":[],"seq_in_adapter":"TGCGTAGA"}],
                "X7236":[{"seq":"AAGAGCCA","instruments":[],"seq_in_adapter":"AAGAGCCA"}],
                "X7237":[{"seq":"ATGGAAGG","instruments":[],"seq_in_adapter":"ATGGAAGG"}],
                "X7238":[{"seq":"GCCAGTAT","instruments":[],"seq_in_adapter":"GCCAGTAT"}],
                "X7239":[{"seq":"CGTAGGTT","instruments":[],"seq_in_adapter":"CGTAGGTT"}],
                "X7240":[{"seq":"CGAGTATG","instruments":[],"seq_in_adapter":"CGAGTATG"}],
                "X7241":[{"seq":"CAAGTGCA","instruments":[],"seq_in_adapter":"CAAGTGCA"}],
                "X7242":[{"seq":"TCGAGTGA","instruments":[],"seq_in_adapter":"TCGAGTGA"}],
                "X7243":[{"seq":"CTACAGTG","instruments":[],"seq_in_adapter":"CTACAGTG"}],
                "X7244":[{"seq":"GATCGTAC","instruments":[],"seq_in_adapter":"GATCGTAC"}],
                "X7245":[{"seq":"CTTCACCA","instruments":[],"seq_in_adapter":"CTTCACCA"}],
                "X7246":[{"seq":"CTCAGCTA","instruments":[],"seq_in_adapter":"CTCAGCTA"}],
                "X7247":[{"seq":"TCTGCTCT","instruments":[],"seq_in_adapter":"TCTGCTCT"}],
                "X7248":[{"seq":"AACCGAAG","instruments":[],"seq_in_adapter":"AACCGAAG"}],
                "X7249":[{"seq":"GCTGTTGT","instruments":[],"seq_in_adapter":"GCTGTTGT"}],
                "X7250":[{"seq":"TTACGGCT","instruments":[],"seq_in_adapter":"TTACGGCT"}],
                "X7251":[{"seq":"GACAAGAG","instruments":[],"seq_in_adapter":"GACAAGAG"}],
                "X7252":[{"seq":"AGGATCTG","instruments":[],"seq_in_adapter":"AGGATCTG"}],
                "X7253":[{"seq":"GTAGCATC","instruments":[],"seq_in_adapter":"GTAGCATC"}],
                "X7254":[{"seq":"GTGTTCCT","instruments":[],"seq_in_adapter":"GTGTTCCT"}],
                "X7255":[{"seq":"AGGATGGT","instruments":[],"seq_in_adapter":"AGGATGGT"}],
                "X7256":[{"seq":"TCACGTTC","instruments":[],"seq_in_adapter":"TCACGTTC"}],
                "X7257":[{"seq":"GCGTTCTA","instruments":[],"seq_in_adapter":"GCGTTCTA"}],
                "X7258":[{"seq":"CTCTGGTT","instruments":[],"seq_in_adapter":"CTCTGGTT"}],
                "X7259":[{"seq":"TTAGGTCG","instruments":[],"seq_in_adapter":"TTAGGTCG"}],
                "X7260":[{"seq":"TCTGAGAG","instruments":[],"seq_in_adapter":"TCTGAGAG"}],
                "X7261":[{"seq":"TTCAGCCT","instruments":[],"seq_in_adapter":"TTCAGCCT"}],
                "X7262":[{"seq":"TCTCCGAT","instruments":[],"seq_in_adapter":"TCTCCGAT"}],
                "X7263":[{"seq":"CAGGTATC","instruments":[],"seq_in_adapter":"CAGGTATC"}],
                "X7264":[{"seq":"AGTCAGGA","instruments":[],"seq_in_adapter":"AGTCAGGA"}],
                "X7265":[{"seq":"AAGGCTGA","instruments":[],"seq_in_adapter":"AAGGCTGA"}],
                "X7266":[{"seq":"CGATGCTT","instruments":[],"seq_in_adapter":"CGATGCTT"}],
                "X7267":[{"seq":"GTATTGGC","instruments":[],"seq_in_adapter":"GTATTGGC"}],
                "X7268":[{"seq":"ACTGTGTC","instruments":[],"seq_in_adapter":"ACTGTGTC"}],
                "X7269":[{"seq":"TGCCTCTT","instruments":[],"seq_in_adapter":"TGCCTCTT"}],
                "X7270":[{"seq":"CAGTCTTC","instruments":[],"seq_in_adapter":"CAGTCTTC"}],
                "X7271":[{"seq":"CATAACGG","instruments":[],"seq_in_adapter":"CATAACGG"}],
                "X7272":[{"seq":"ACTGCTAG","instruments":[],"seq_in_adapter":"ACTGCTAG"}],
                "X7273":[{"seq":"ATTCTGGC","instruments":[],"seq_in_adapter":"ATTCTGGC"}],
                "X7274":[{"seq":"TTCTCTCG","instruments":[],"seq_in_adapter":"TTCTCTCG"}],
                "X7275":[{"seq":"TCCGAGTT","instruments":[],"seq_in_adapter":"TCCGAGTT"}],
                "X7276":[{"seq":"CGAACTGT","instruments":[],"seq_in_adapter":"CGAACTGT"}],
                "X7277":[{"seq":"AACGGTCA","instruments":[],"seq_in_adapter":"AACGGTCA"}],
                "X7278":[{"seq":"AGCAGATG","instruments":[],"seq_in_adapter":"AGCAGATG"}],
                "X7279":[{"seq":"TATCAGCG","instruments":[],"seq_in_adapter":"TATCAGCG"}],
                "X7280":[{"seq":"TCAGACGA","instruments":[],"seq_in_adapter":"TCAGACGA"}],
                "X7281":[{"seq":"ACCATGTG","instruments":[],"seq_in_adapter":"ACCATGTG"}],
                "X7282":[{"seq":"CTAACTCG","instruments":[],"seq_in_adapter":"CTAACTCG"}],
                "X7283":[{"seq":"GCTTAGCT","instruments":[],"seq_in_adapter":"GCTTAGCT"}],
                "X7284":[{"seq":"CATGGAAC","instruments":[],"seq_in_adapter":"CATGGAAC"}],
                "X7285":[{"seq":"TAGGATGC","instruments":[],"seq_in_adapter":"TAGGATGC"}],
                "X7286":[{"seq":"GTTCATGG","instruments":[],"seq_in_adapter":"GTTCATGG"}],
                "X7287":[{"seq":"TCGTGGAT","instruments":[],"seq_in_adapter":"TCGTGGAT"}],
                "X7288":[{"seq":"ACCTTCTC","instruments":[],"seq_in_adapter":"ACCTTCTC"}],
                "X7289":[{"seq":"CATTGCCT","instruments":[],"seq_in_adapter":"CATTGCCT"}],
                "X7290":[{"seq":"CTAGGTGA","instruments":[],"seq_in_adapter":"CTAGGTGA"}],
                "X7291":[{"seq":"TCCGTATG","instruments":[],"seq_in_adapter":"TCCGTATG"}],
                "X7292":[{"seq":"ACGATGAC","instruments":[],"seq_in_adapter":"ACGATGAC"}],
                "X7293":[{"seq":"GTCGGTAA","instruments":[],"seq_in_adapter":"GTCGGTAA"}],
                "X7294":[{"seq":"TCGAAGGT","instruments":[],"seq_in_adapter":"TCGAAGGT"}],
                "X7295":[{"seq":"AGAAGCGT","instruments":[],"seq_in_adapter":"AGAAGCGT"}],
                "X7296":[{"seq":"CTCTACTC","instruments":[],"seq_in_adapter":"CTCTACTC"}],
                "X7297":[{"seq":"CTAGGCAT","instruments":[],"seq_in_adapter":"CTAGGCAT"}],
                "X7298":[{"seq":"TGGAGTTG","instruments":[],"seq_in_adapter":"TGGAGTTG"}],
                "X7299":[{"seq":"GAGGACTT","instruments":[],"seq_in_adapter":"GAGGACTT"}],
                "X7300":[{"seq":"CAATCGAC","instruments":[],"seq_in_adapter":"CAATCGAC"}],
                "X7301":[{"seq":"TCTAACGC","instruments":[],"seq_in_adapter":"TCTAACGC"}],
                "X7302":[{"seq":"TCTCGCAA","instruments":[],"seq_in_adapter":"TCTCGCAA"}],
                "X7303":[{"seq":"ATCGGTGT","instruments":[],"seq_in_adapter":"ATCGGTGT"}],
                "X7304":[{"seq":"GAGATACG","instruments":[],"seq_in_adapter":"GAGATACG"}],
                "X7305":[{"seq":"GTCTCCTT","instruments":[],"seq_in_adapter":"GTCTCCTT"}],
                "X7306":[{"seq":"AGTCGACA","instruments":[],"seq_in_adapter":"AGTCGACA"}],
                "X7307":[{"seq":"CGGATTGA","instruments":[],"seq_in_adapter":"CGGATTGA"}],
                "X7308":[{"seq":"CACAAGTC","instruments":[],"seq_in_adapter":"CACAAGTC"}],
                "X7309":[{"seq":"TACATCGG","instruments":[],"seq_in_adapter":"TACATCGG"}],
                "X7310":[{"seq":"AGCTCCTA","instruments":[],"seq_in_adapter":"AGCTCCTA"}],
                "X7311":[{"seq":"ACTCGTTG","instruments":[],"seq_in_adapter":"ACTCGTTG"}],
                "X7312":[{"seq":"CTGACACA","instruments":[],"seq_in_adapter":"CTGACACA"}],
                "X7313":[{"seq":"CAACCTAG","instruments":[],"seq_in_adapter":"CAACCTAG"}],
                "X7314":[{"seq":"AAGGACAC","instruments":[],"seq_in_adapter":"AAGGACAC"}],
                "X7315":[{"seq":"TGCAGGTA","instruments":[],"seq_in_adapter":"TGCAGGTA"}],
                "X7316":[{"seq":"ACCTAAGG","instruments":[],"seq_in_adapter":"ACCTAAGG"}],
                "X7317":[{"seq":"AGTCTGTG","instruments":[],"seq_in_adapter":"AGTCTGTG"}],
                "X7318":[{"seq":"AGGTTCGA","instruments":[],"seq_in_adapter":"AGGTTCGA"}],
                "X7319":[{"seq":"GACTATGC","instruments":[],"seq_in_adapter":"GACTATGC"}],
                "X7320":[{"seq":"TTCAGGAG","instruments":[],"seq_in_adapter":"TTCAGGAG"}],
                "X7321":[{"seq":"TGTGCGTT","instruments":[],"seq_in_adapter":"TGTGCGTT"}],
                "X7322":[{"seq":"CGAGACTA","instruments":[],"seq_in_adapter":"CGAGACTA"}],
                "X7323":[{"seq":"CTCAGAGT","instruments":[],"seq_in_adapter":"CTCAGAGT"}],
                "X7324":[{"seq":"GCCATAAC","instruments":[],"seq_in_adapter":"GCCATAAC"}],
                "X7325":[{"seq":"TTACCGAG","instruments":[],"seq_in_adapter":"TTACCGAG"}],
                "X7326":[{"seq":"GCTCTGTA","instruments":[],"seq_in_adapter":"GCTCTGTA"}],
                "X7327":[{"seq":"CGTTATGC","instruments":[],"seq_in_adapter":"CGTTATGC"}],
                "X7328":[{"seq":"GTCTGATC","instruments":[],"seq_in_adapter":"GTCTGATC"}],
                "X7329":[{"seq":"TAGTTGCG","instruments":[],"seq_in_adapter":"TAGTTGCG"}],
                "X7330":[{"seq":"TGATCGGA","instruments":[],"seq_in_adapter":"TGATCGGA"}],
                "X7331":[{"seq":"CCAAGTTG","instruments":[],"seq_in_adapter":"CCAAGTTG"}],
                "X7332":[{"seq":"CCTACTGA","instruments":[],"seq_in_adapter":"CCTACTGA"}],
                "X7333":[{"seq":"CTTGCTGT","instruments":[],"seq_in_adapter":"CTTGCTGT"}],
                "X7334":[{"seq":"TGCCATTC","instruments":[],"seq_in_adapter":"TGCCATTC"}],
                "X7335":[{"seq":"TTGATCCG","instruments":[],"seq_in_adapter":"TTGATCCG"}],
                "X7336":[{"seq":"AGTGCAGT","instruments":[],"seq_in_adapter":"AGTGCAGT"}],
                "X7337":[{"seq":"GACTTAGG","instruments":[],"seq_in_adapter":"GACTTAGG"}],
                "X7338":[{"seq":"CGTACGAA","instruments":[],"seq_in_adapter":"CGTACGAA"}],
                "X7339":[{"seq":"TACCAGGA","instruments":[],"seq_in_adapter":"TACCAGGA"}],
                "X7340":[{"seq":"CGTCAATG","instruments":[],"seq_in_adapter":"CGTCAATG"}],
                "X7341":[{"seq":"GAAGAGGT","instruments":[],"seq_in_adapter":"GAAGAGGT"}],
                "X7342":[{"seq":"GACGAATG","instruments":[],"seq_in_adapter":"GACGAATG"}],
                "X7343":[{"seq":"AGGAGGAA","instruments":[],"seq_in_adapter":"AGGAGGAA"}],
                "X7344":[{"seq":"CTTACAGC","instruments":[],"seq_in_adapter":"CTTACAGC"}],
                "X7345":[{"seq":"GAGATGTC","instruments":[],"seq_in_adapter":"GAGATGTC"}],
                "X7346":[{"seq":"TACGGTTG","instruments":[],"seq_in_adapter":"TACGGTTG"}],
                "X7347":[{"seq":"CTATCGCA","instruments":[],"seq_in_adapter":"CTATCGCA"}],
                "X7348":[{"seq":"TCGAACCA","instruments":[],"seq_in_adapter":"TCGAACCA"}],
                "X7349":[{"seq":"GAACGCTT","instruments":[],"seq_in_adapter":"GAACGCTT"}],
                "X7350":[{"seq":"CAGAATCG","instruments":[],"seq_in_adapter":"CAGAATCG"}],
                "X7351":[{"seq":"ATGGTTGC","instruments":[],"seq_in_adapter":"ATGGTTGC"}],
                "X7352":[{"seq":"GCTGGATT","instruments":[],"seq_in_adapter":"GCTGGATT"}],
                "X7353":[{"seq":"GATGCACT","instruments":[],"seq_in_adapter":"GATGCACT"}],
                "X7354":[{"seq":"ACCAATGC","instruments":[],"seq_in_adapter":"ACCAATGC"}],
                "X7355":[{"seq":"GTCCTAAG","instruments":[],"seq_in_adapter":"GTCCTAAG"}],
                "X7356":[{"seq":"CCGACTAT","instruments":[],"seq_in_adapter":"CCGACTAT"}],
                "X7357":[{"seq":"TTGGTCTC","instruments":[],"seq_in_adapter":"TTGGTCTC"}],
                "X7358":[{"seq":"GCCTTGTT","instruments":[],"seq_in_adapter":"GCCTTGTT"}],
                "X7359":[{"seq":"GATACTGG","instruments":[],"seq_in_adapter":"GATACTGG"}],
                "X7360":[{"seq":"ATTCGAGG","instruments":[],"seq_in_adapter":"ATTCGAGG"}],
                "X7361":[{"seq":"GTCAGTTG","instruments":[],"seq_in_adapter":"GTCAGTTG"}],
                "X7362":[{"seq":"GTAGAGCA","instruments":[],"seq_in_adapter":"GTAGAGCA"}],
                "X7363":[{"seq":"ACGTGATG","instruments":[],"seq_in_adapter":"ACGTGATG"}],
                "X7364":[{"seq":"TAAGTGGC","instruments":[],"seq_in_adapter":"TAAGTGGC"}],
                "X7365":[{"seq":"TGTGAAGC","instruments":[],"seq_in_adapter":"TGTGAAGC"}],
                "X7366":[{"seq":"CATTCGGT","instruments":[],"seq_in_adapter":"CATTCGGT"}],
                "X7367":[{"seq":"TTGGTGAG","instruments":[],"seq_in_adapter":"TTGGTGAG"}],
                "X7368":[{"seq":"CAGTTCTG","instruments":[],"seq_in_adapter":"CAGTTCTG"}],
                "X7369":[{"seq":"AGGCTTCT","instruments":[],"seq_in_adapter":"AGGCTTCT"}],
                "X7370":[{"seq":"GAATCGTG","instruments":[],"seq_in_adapter":"GAATCGTG"}],
                "X7371":[{"seq":"ACCAGCTT","instruments":[],"seq_in_adapter":"ACCAGCTT"}],
                "X7372":[{"seq":"CTCATTGC","instruments":[],"seq_in_adapter":"CTCATTGC"}],
                "X7373":[{"seq":"CGATAGAG","instruments":[],"seq_in_adapter":"CGATAGAG"}],
                "X7374":[{"seq":"TGGAGAGT","instruments":[],"seq_in_adapter":"TGGAGAGT"}],
                "X7375":[{"seq":"GTATGCTG","instruments":[],"seq_in_adapter":"GTATGCTG"}],
                "X7376":[{"seq":"CTGGAGTA","instruments":[],"seq_in_adapter":"CTGGAGTA"}],
                "X7377":[{"seq":"AATGCCTC","instruments":[],"seq_in_adapter":"AATGCCTC"}],
                "X7378":[{"seq":"TGAGGTGT","instruments":[],"seq_in_adapter":"TGAGGTGT"}],
                "X7379":[{"seq":"ACATTGCG","instruments":[],"seq_in_adapter":"ACATTGCG"}],
                "X7380":[{"seq":"TCTCTAGG","instruments":[],"seq_in_adapter":"TCTCTAGG"}],
                "X7381":[{"seq":"CGCTAGTA","instruments":[],"seq_in_adapter":"CGCTAGTA"}],
                "X7382":[{"seq":"AATGGACG","instruments":[],"seq_in_adapter":"AATGGACG"}],
                "X7383":[{"seq":"GATAGCGA","instruments":[],"seq_in_adapter":"GATAGCGA"}],
                "X7384":[{"seq":"CGACCATT","instruments":[],"seq_in_adapter":"CGACCATT"}],
            },
            "i5":{
                "X5001":[{"seq":"ATATGCGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATATGCGC"}, {"seq":"GCGCATAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATATGCGC"}],
                "X5002":[{"seq":"TGGTACAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGTACAG"}, {"seq":"CTGTACCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGTACAG"}],
                "X5003":[{"seq":"AACCGTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACCGTTC"}, {"seq":"GAACGGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACCGTTC"}],
                "X5004":[{"seq":"TAACCGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAACCGGT"}, {"seq":"ACCGGTTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAACCGGT"}],
                "X5005":[{"seq":"GAACATCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAACATCG"}, {"seq":"CGATGTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAACATCG"}],
                "X5006":[{"seq":"CCTTGTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTTGTAG"}, {"seq":"CTACAAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTTGTAG"}],
                "X5007":[{"seq":"TCAGGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCAGGCTT"}, {"seq":"AAGCCTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCAGGCTT"}],
                "X5008":[{"seq":"GTTCTCGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTCTCGT"}, {"seq":"ACGAGAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTCTCGT"}],
                "X5009":[{"seq":"AGAACGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGAACGAG"}, {"seq":"CTCGTTCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGAACGAG"}],
                "X5010":[{"seq":"TGCTTCCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCTTCCA"}, {"seq":"TGGAAGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCTTCCA"}],
                "X5011":[{"seq":"CTTCGACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTCGACT"}, {"seq":"AGTCGAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTCGACT"}],
                "X5012":[{"seq":"CACCTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACCTGTT"}, {"seq":"AACAGGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACCTGTT"}],
                "X5013":[{"seq":"ATCACACG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCACACG"}, {"seq":"CGTGTGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCACACG"}],
                "X5014":[{"seq":"CCGTAAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCGTAAGA"}, {"seq":"TCTTACGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCGTAAGA"}],
                "X5015":[{"seq":"TACGCCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACGCCTT"}, {"seq":"AAGGCGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACGCCTT"}],
                "X5016":[{"seq":"CGACGTTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGACGTTA"}, {"seq":"TAACGTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGACGTTA"}],
                "X5017":[{"seq":"ATGCACGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATGCACGA"}, {"seq":"TCGTGCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATGCACGA"}],
                "X5018":[{"seq":"CCTGATTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTGATTG"}, {"seq":"CAATCAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTGATTG"}],
                "X5019":[{"seq":"GTAGGAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTAGGAGT"}, {"seq":"ACTCCTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTAGGAGT"}],
                "X5020":[{"seq":"ACTAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTAGGAG"}, {"seq":"CTCCTAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTAGGAG"}],
                "X5021":[{"seq":"CACTAGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACTAGCT"}, {"seq":"AGCTAGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACTAGCT"}],
                "X5022":[{"seq":"ACGACTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGACTTG"}, {"seq":"CAAGTCGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGACTTG"}],
                "X5023":[{"seq":"CGTGTGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTGTGTA"}, {"seq":"TACACACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTGTGTA"}],
                "X5024":[{"seq":"GTTGACCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTGACCT"}, {"seq":"AGGTCAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTGACCT"}],
                "X5025":[{"seq":"ACTCCATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTCCATC"}, {"seq":"GATGGAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTCCATC"}],
                "X5026":[{"seq":"CAATGTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAATGTGG"}, {"seq":"CCACATTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAATGTGG"}],
                "X5027":[{"seq":"TTGCAGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGCAGAC"}, {"seq":"GTCTGCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGCAGAC"}],
                "X5028":[{"seq":"CAGTCCAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGTCCAA"}, {"seq":"TTGGACTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGTCCAA"}],
                "X5029":[{"seq":"ACGTTCAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGTTCAG"}, {"seq":"CTGAACGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGTTCAG"}],
                "X5030":[{"seq":"AACGTCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACGTCTG"}, {"seq":"CAGACGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACGTCTG"}],
                "X5031":[{"seq":"TATCGGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TATCGGTC"}, {"seq":"GACCGATA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TATCGGTC"}],
                "X5032":[{"seq":"CGCTCTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGCTCTAT"}, {"seq":"ATAGAGCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGCTCTAT"}],
                "X5033":[{"seq":"GATTGCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATTGCTC"}, {"seq":"GAGCAATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATTGCTC"}],
                "X5034":[{"seq":"GATGTGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATGTGTG"}, {"seq":"CACACATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATGTGTG"}],
                "X5035":[{"seq":"CGCAATCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGCAATCT"}, {"seq":"AGATTGCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGCAATCT"}],
                "X5036":[{"seq":"TGGTAGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGTAGCT"}, {"seq":"AGCTACCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGTAGCT"}],
                "X5037":[{"seq":"GATAGGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATAGGCT"}, {"seq":"AGCCTATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATAGGCT"}],
                "X5038":[{"seq":"AGTGGATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTGGATC"}, {"seq":"GATCCACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTGGATC"}],
                "X5039":[{"seq":"TTGGACGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGGACGT"}, {"seq":"ACGTCCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGGACGT"}],
                "X5040":[{"seq":"ATGACGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATGACGTC"}, {"seq":"GACGTCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATGACGTC"}],
                "X5041":[{"seq":"GAAGTTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAAGTTGG"}, {"seq":"CCAACTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAAGTTGG"}],
                "X5042":[{"seq":"CATACCAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATACCAC"}, {"seq":"GTGGTATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATACCAC"}],
                "X5043":[{"seq":"CTGTTGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGTTGAC"}, {"seq":"GTCAACAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGTTGAC"}],
                "X5044":[{"seq":"TGGCATGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGCATGT"}, {"seq":"ACATGCCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGCATGT"}],
                "X5045":[{"seq":"ATCGCCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCGCCAT"}, {"seq":"ATGGCGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCGCCAT"}],
                "X5046":[{"seq":"TTGCGAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGCGAAG"}, {"seq":"CTTCGCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGCGAAG"}],
                "X5047":[{"seq":"AGTTCGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTTCGTC"}, {"seq":"GACGAACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTTCGTC"}],
                "X5048":[{"seq":"GAGCAGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGCAGTA"}, {"seq":"TACTGCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGCAGTA"}],
                "X5049":[{"seq":"ACAGCTCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACAGCTCA"}, {"seq":"TGAGCTGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACAGCTCA"}],
                "X5050":[{"seq":"GATCGAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATCGAGT"}, {"seq":"ACTCGATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATCGAGT"}],
                "X5051":[{"seq":"AGCGTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCGTGTT"}, {"seq":"AACACGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCGTGTT"}],
                "X5052":[{"seq":"GTTACGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTACGCA"}, {"seq":"TGCGTAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTACGCA"}],
                "X5053":[{"seq":"TGAAGACG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGAAGACG"}, {"seq":"CGTCTTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGAAGACG"}],
                "X5054":[{"seq":"ACTGAGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTGAGGT"}, {"seq":"ACCTCAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTGAGGT"}],
                "X5055":[{"seq":"CGGTTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGGTTGTT"}, {"seq":"AACAACCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGGTTGTT"}],
                "X5056":[{"seq":"GTTGTTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTGTTCG"}, {"seq":"CGAACAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTGTTCG"}],
                "X5057":[{"seq":"GAAGGAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAAGGAAG"}, {"seq":"CTTCCTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAAGGAAG"}],
                "X5058":[{"seq":"AGCACTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCACTTC"}, {"seq":"GAAGTGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCACTTC"}],
                "X5059":[{"seq":"GTCATCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCATCGA"}, {"seq":"TCGATGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCATCGA"}],
                "X5060":[{"seq":"TGTGACTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTGACTG"}, {"seq":"CAGTCACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTGACTG"}],
                "X5061":[{"seq":"CAACACCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAACACCT"}, {"seq":"AGGTGTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAACACCT"}],
                "X5062":[{"seq":"ATGCCTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATGCCTGT"}, {"seq":"ACAGGCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATGCCTGT"}],
                "X5063":[{"seq":"CATGGCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATGGCTA"}, {"seq":"TAGCCATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATGGCTA"}],
                "X5064":[{"seq":"GTGAAGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTGAAGTG"}, {"seq":"CACTTCAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTGAAGTG"}],
                "X5065":[{"seq":"CGTTGCAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTTGCAA"}, {"seq":"TTGCAACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTTGCAA"}],
                "X5066":[{"seq":"ATCCGGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCCGGTA"}, {"seq":"TACCGGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCCGGTA"}],
                "X5067":[{"seq":"GCGTCATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCGTCATT"}, {"seq":"AATGACGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCGTCATT"}],
                "X5068":[{"seq":"GCACAACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCACAACT"}, {"seq":"AGTTGTGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCACAACT"}],
                "X5069":[{"seq":"GATTACCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATTACCG"}, {"seq":"CGGTAATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATTACCG"}],
                "X5070":[{"seq":"ACCACGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCACGAT"}, {"seq":"ATCGTGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCACGAT"}],
                "X5071":[{"seq":"GTCGAAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCGAAGA"}, {"seq":"TCTTCGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCGAAGA"}],
                "X5072":[{"seq":"CCTTGATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTTGATC"}, {"seq":"GATCAAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTTGATC"}],
                "X5073":[{"seq":"AAGCACTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGCACTG"}, {"seq":"CAGTGCTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGCACTG"}],
                "X5074":[{"seq":"TTCGTTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCGTTGG"}, {"seq":"CCAACGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCGTTGG"}],
                "X5075":[{"seq":"TCGCTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGCTGTT"}, {"seq":"AACAGCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGCTGTT"}],
                "X5076":[{"seq":"GAATCCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAATCCGA"}, {"seq":"TCGGATTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAATCCGA"}],
                "X5077":[{"seq":"GTGCCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTGCCATA"}, {"seq":"TATGGCAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTGCCATA"}],
                "X5078":[{"seq":"CTTAGGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTAGGAC"}, {"seq":"GTCCTAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTAGGAC"}],
                "X5079":[{"seq":"AACTGAGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACTGAGC"}, {"seq":"GCTCAGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACTGAGC"}],
                "X5080":[{"seq":"GACGATCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACGATCT"}, {"seq":"AGATCGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACGATCT"}],
                "X5081":[{"seq":"ATCCAGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCCAGAG"}, {"seq":"CTCTGGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCCAGAG"}],
                "X5082":[{"seq":"AGAGTAGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGAGTAGC"}, {"seq":"GCTACTCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGAGTAGC"}],
                "X5083":[{"seq":"TGGACTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGACTCT"}, {"seq":"AGAGTCCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGACTCT"}],
                "X5084":[{"seq":"TACGCTAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACGCTAC"}, {"seq":"GTAGCGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACGCTAC"}],
                "X5085":[{"seq":"GCTATCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTATCCT"}, {"seq":"AGGATAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTATCCT"}],
                "X5086":[{"seq":"GCAAGATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCAAGATC"}, {"seq":"GATCTTGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCAAGATC"}],
                "X5087":[{"seq":"ATCGATCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCGATCG"}, {"seq":"CGATCGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCGATCG"}],
                "X5088":[{"seq":"CGGCTAAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGGCTAAT"}, {"seq":"ATTAGCCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGGCTAAT"}],
                "X5089":[{"seq":"ACGGAACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGGAACA"}, {"seq":"TGTTCCGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGGAACA"}],
                "X5090":[{"seq":"CGCATGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGCATGAT"}, {"seq":"ATCATGCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGCATGAT"}],
                "X5091":[{"seq":"TTCCAAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCCAAGG"}, {"seq":"CCTTGGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCCAAGG"}],
                "X5092":[{"seq":"CTTGTCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTGTCGA"}, {"seq":"TCGACAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTGTCGA"}],
                "X5093":[{"seq":"GAGACGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGACGAT"}, {"seq":"ATCGTCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGACGAT"}],
                "X5094":[{"seq":"TGAGCTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGAGCTAG"}, {"seq":"CTAGCTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGAGCTAG"}],
                "X5095":[{"seq":"ACTCTCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTCTCGA"}, {"seq":"TCGAGAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTCTCGA"}],
                "X5096":[{"seq":"CTGATCGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGATCGT"}, {"seq":"ACGATCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGATCGT"}],
                "X5097":[{"seq":"CGACCATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGACCATT"}, {"seq":"AATGGTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGACCATT"}],
                "X5098":[{"seq":"GATAGCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATAGCGA"}, {"seq":"TCGCTATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATAGCGA"}],
                "X5099":[{"seq":"AATGGACG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AATGGACG"}, {"seq":"CGTCCATT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AATGGACG"}],
                "X5100":[{"seq":"CGCTAGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGCTAGTA"}, {"seq":"TACTAGCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGCTAGTA"}],
                "X5101":[{"seq":"TCTCTAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTCTAGG"}, {"seq":"CCTAGAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTCTAGG"}],
                "X5102":[{"seq":"ACATTGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACATTGCG"}, {"seq":"CGCAATGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACATTGCG"}],
                "X5103":[{"seq":"TGAGGTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGAGGTGT"}, {"seq":"ACACCTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGAGGTGT"}],
                "X5104":[{"seq":"AATGCCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AATGCCTC"}, {"seq":"GAGGCATT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AATGCCTC"}],
                "X5105":[{"seq":"CTGGAGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGGAGTA"}, {"seq":"TACTCCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGGAGTA"}],
                "X5106":[{"seq":"GTATGCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTATGCTG"}, {"seq":"CAGCATAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTATGCTG"}],
                "X5107":[{"seq":"TGGAGAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGAGAGT"}, {"seq":"ACTCTCCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGAGAGT"}],
                "X5108":[{"seq":"CGATAGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGATAGAG"}, {"seq":"CTCTATCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGATAGAG"}],
                "X5109":[{"seq":"CTCATTGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCATTGC"}, {"seq":"GCAATGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCATTGC"}],
                "X5110":[{"seq":"ACCAGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCAGCTT"}, {"seq":"AAGCTGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCAGCTT"}],
                "X5111":[{"seq":"GAATCGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAATCGTG"}, {"seq":"CACGATTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAATCGTG"}],
                "X5112":[{"seq":"AGGCTTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGCTTCT"}, {"seq":"AGAAGCCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGCTTCT"}],
                "X5113":[{"seq":"CAGTTCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGTTCTG"}, {"seq":"CAGAACTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGTTCTG"}],
                "X5114":[{"seq":"TTGGTGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGGTGAG"}, {"seq":"CTCACCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGGTGAG"}],
                "X5115":[{"seq":"CATTCGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATTCGGT"}, {"seq":"ACCGAATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATTCGGT"}],
                "X5116":[{"seq":"TGTGAAGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTGAAGC"}, {"seq":"GCTTCACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTGAAGC"}],
                "X5117":[{"seq":"TAAGTGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAAGTGGC"}, {"seq":"GCCACTTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAAGTGGC"}],
                "X5118":[{"seq":"ACGTGATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGTGATG"}, {"seq":"CATCACGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGTGATG"}],
                "X5119":[{"seq":"GTAGAGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTAGAGCA"}, {"seq":"TGCTCTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTAGAGCA"}],
                "X5120":[{"seq":"GTCAGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCAGTTG"}, {"seq":"CAACTGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCAGTTG"}],
                "X5121":[{"seq":"ATTCGAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATTCGAGG"}, {"seq":"CCTCGAAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATTCGAGG"}],
                "X5122":[{"seq":"GATACTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATACTGG"}, {"seq":"CCAGTATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATACTGG"}],
                "X5123":[{"seq":"GCCTTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCCTTGTT"}, {"seq":"AACAAGGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCCTTGTT"}],
                "X5124":[{"seq":"TTGGTCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGGTCTC"}, {"seq":"GAGACCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGGTCTC"}],
                "X5125":[{"seq":"CCGACTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCGACTAT"}, {"seq":"ATAGTCGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCGACTAT"}],
                "X5126":[{"seq":"GTCCTAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCCTAAG"}, {"seq":"CTTAGGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCCTAAG"}],
                "X5127":[{"seq":"ACCAATGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCAATGC"}, {"seq":"GCATTGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCAATGC"}],
                "X5128":[{"seq":"GATGCACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATGCACT"}, {"seq":"AGTGCATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATGCACT"}],
                "X5129":[{"seq":"GCTGGATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTGGATT"}, {"seq":"AATCCAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTGGATT"}],
                "X5130":[{"seq":"ATGGTTGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATGGTTGC"}, {"seq":"GCAACCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATGGTTGC"}],
                "X5131":[{"seq":"CAGAATCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGAATCG"}, {"seq":"CGATTCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGAATCG"}],
                "X5132":[{"seq":"GAACGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAACGCTT"}, {"seq":"AAGCGTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAACGCTT"}],
                "X5133":[{"seq":"TCGAACCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGAACCA"}, {"seq":"TGGTTCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGAACCA"}],
                "X5134":[{"seq":"CTATCGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTATCGCA"}, {"seq":"TGCGATAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTATCGCA"}],
                "X5135":[{"seq":"TACGGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACGGTTG"}, {"seq":"CAACCGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACGGTTG"}],
                "X5136":[{"seq":"GAGATGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGATGTC"}, {"seq":"GACATCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGATGTC"}],
                "X5137":[{"seq":"CTTACAGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTACAGC"}, {"seq":"GCTGTAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTACAGC"}],
                "X5138":[{"seq":"AGGAGGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGAGGAA"}, {"seq":"TTCCTCCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGAGGAA"}],
                "X5139":[{"seq":"GACGAATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACGAATG"}, {"seq":"CATTCGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACGAATG"}],
                "X5140":[{"seq":"GAAGAGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAAGAGGT"}, {"seq":"ACCTCTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAAGAGGT"}],
                "X5141":[{"seq":"CGTCAATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTCAATG"}, {"seq":"CATTGACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTCAATG"}],
                "X5142":[{"seq":"TACCAGGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACCAGGA"}, {"seq":"TCCTGGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACCAGGA"}],
                "X5143":[{"seq":"CGTACGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTACGAA"}, {"seq":"TTCGTACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTACGAA"}],
                "X5144":[{"seq":"GACTTAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACTTAGG"}, {"seq":"CCTAAGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACTTAGG"}],
                "X5145":[{"seq":"AGTGCAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTGCAGT"}, {"seq":"ACTGCACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTGCAGT"}],
                "X5146":[{"seq":"TTGATCCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGATCCG"}, {"seq":"CGGATCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGATCCG"}],
                "X5147":[{"seq":"TGCCATTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCCATTC"}, {"seq":"GAATGGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCCATTC"}],
                "X5148":[{"seq":"CTTGCTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTGCTGT"}, {"seq":"ACAGCAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTGCTGT"}],
                "X5149":[{"seq":"CCTACTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTACTGA"}, {"seq":"TCAGTAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTACTGA"}],
                "X5150":[{"seq":"CCAAGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAAGTTG"}, {"seq":"CAACTTGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAAGTTG"}],
                "X5151":[{"seq":"TGATCGGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGATCGGA"}, {"seq":"TCCGATCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGATCGGA"}],
                "X5152":[{"seq":"TAGTTGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGTTGCG"}, {"seq":"CGCAACTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGTTGCG"}],
                "X5153":[{"seq":"GTCTGATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCTGATC"}, {"seq":"GATCAGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCTGATC"}],
                "X5154":[{"seq":"CGTTATGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTTATGC"}, {"seq":"GCATAACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTTATGC"}],
                "X5155":[{"seq":"GCTCTGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTCTGTA"}, {"seq":"TACAGAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTCTGTA"}],
                "X5156":[{"seq":"TTACCGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTACCGAG"}, {"seq":"CTCGGTAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTACCGAG"}],
                "X5157":[{"seq":"GCCATAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCCATAAC"}, {"seq":"GTTATGGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCCATAAC"}],
                "X5158":[{"seq":"CTCAGAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCAGAGT"}, {"seq":"ACTCTGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCAGAGT"}],
                "X5159":[{"seq":"CGAGACTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGAGACTA"}, {"seq":"TAGTCTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGAGACTA"}],
                "X5160":[{"seq":"TGTGCGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTGCGTT"}, {"seq":"AACGCACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTGCGTT"}],
                "X5161":[{"seq":"TTCAGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCAGGAG"}, {"seq":"CTCCTGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCAGGAG"}],
                "X5162":[{"seq":"GACTATGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACTATGC"}, {"seq":"GCATAGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACTATGC"}],
                "X5163":[{"seq":"AGGTTCGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGTTCGA"}, {"seq":"TCGAACCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGTTCGA"}],
                "X5164":[{"seq":"AGTCTGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTCTGTG"}, {"seq":"CACAGACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTCTGTG"}],
                "X5165":[{"seq":"ACCTAAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCTAAGG"}, {"seq":"CCTTAGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCTAAGG"}],
                "X5166":[{"seq":"TGCAGGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCAGGTA"}, {"seq":"TACCTGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCAGGTA"}],
                "X5167":[{"seq":"AAGGACAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGGACAC"}, {"seq":"GTGTCCTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGACAC"}],
                "X5168":[{"seq":"CAACCTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAACCTAG"}, {"seq":"CTAGGTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAACCTAG"}],
                "X5169":[{"seq":"CTGACACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGACACA"}, {"seq":"TGTGTCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGACACA"}],
                "X5170":[{"seq":"ACTCGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTCGTTG"}, {"seq":"CAACGAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTCGTTG"}],
                "X5171":[{"seq":"AGCTCCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCTCCTA"}, {"seq":"TAGGAGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCTCCTA"}],
                "X5172":[{"seq":"TACATCGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACATCGG"}, {"seq":"CCGATGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACATCGG"}],
                "X5173":[{"seq":"CACAAGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACAAGTC"}, {"seq":"GACTTGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACAAGTC"}],
                "X5174":[{"seq":"CGGATTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGGATTGA"}, {"seq":"TCAATCCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGGATTGA"}],
                "X5175":[{"seq":"AGTCGACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTCGACA"}, {"seq":"TGTCGACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTCGACA"}],
                "X5176":[{"seq":"GTCTCCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCTCCTT"}, {"seq":"AAGGAGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCTCCTT"}],
                "X5177":[{"seq":"GAGATACG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGATACG"}, {"seq":"CGTATCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGATACG"}],
                "X5178":[{"seq":"ATCGGTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCGGTGT"}, {"seq":"ACACCGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCGGTGT"}],
                "X5179":[{"seq":"TCTCGCAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTCGCAA"}, {"seq":"TTGCGAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTCGCAA"}],
                "X5180":[{"seq":"TCTAACGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTAACGC"}, {"seq":"GCGTTAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTAACGC"}],
                "X5181":[{"seq":"CAATCGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAATCGAC"}, {"seq":"GTCGATTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAATCGAC"}],
                "X5182":[{"seq":"GAGGACTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGGACTT"}, {"seq":"AAGTCCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGGACTT"}],
                "X5183":[{"seq":"TGGAGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGAGTTG"}, {"seq":"CAACTCCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGAGTTG"}],
                "X5184":[{"seq":"CTAGGCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTAGGCAT"}, {"seq":"ATGCCTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTAGGCAT"}],
                "X5185":[{"seq":"CTCTACTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCTACTC"}, {"seq":"GAGTAGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCTACTC"}],
                "X5186":[{"seq":"AGAAGCGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGAAGCGT"}, {"seq":"ACGCTTCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGAAGCGT"}],
                "X5187":[{"seq":"TCGAAGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGAAGGT"}, {"seq":"ACCTTCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGAAGGT"}],
                "X5188":[{"seq":"GTCGGTAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCGGTAA"}, {"seq":"TTACCGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCGGTAA"}],
                "X5189":[{"seq":"ACGATGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGATGAC"}, {"seq":"GTCATCGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGATGAC"}],
                "X5190":[{"seq":"TCCGTATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCGTATG"}, {"seq":"CATACGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCGTATG"}],
                "X5191":[{"seq":"CTAGGTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTAGGTGA"}, {"seq":"TCACCTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTAGGTGA"}],
                "X5192":[{"seq":"CATTGCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATTGCCT"}, {"seq":"AGGCAATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATTGCCT"}],
                "X5193":[{"seq":"ACCTTCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCTTCTC"}, {"seq":"GAGAAGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCTTCTC"}],
                "X5194":[{"seq":"TCGTGGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGTGGAT"}, {"seq":"ATCCACGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGTGGAT"}],
                "X5195":[{"seq":"GTTCATGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTCATGG"}, {"seq":"CCATGAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTCATGG"}],
                "X5196":[{"seq":"TAGGATGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGGATGC"}, {"seq":"GCATCCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGGATGC"}],
                "X5197":[{"seq":"CATGGAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATGGAAC"}, {"seq":"GTTCCATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATGGAAC"}],
                "X5198":[{"seq":"GCTTAGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTTAGCT"}, {"seq":"AGCTAAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTTAGCT"}],
                "X5199":[{"seq":"CTAACTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTAACTCG"}, {"seq":"CGAGTTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTAACTCG"}],
                "X5200":[{"seq":"ACCATGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCATGTG"}, {"seq":"CACATGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCATGTG"}],
                "X5201":[{"seq":"TCAGACGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCAGACGA"}, {"seq":"TCGTCTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCAGACGA"}],
                "X5202":[{"seq":"TATCAGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TATCAGCG"}, {"seq":"CGCTGATA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TATCAGCG"}],
                "X5203":[{"seq":"AGCAGATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCAGATG"}, {"seq":"CATCTGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCAGATG"}],
                "X5204":[{"seq":"AACGGTCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACGGTCA"}, {"seq":"TGACCGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACGGTCA"}],
                "X5205":[{"seq":"CGAACTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGAACTGT"}, {"seq":"ACAGTTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGAACTGT"}],
                "X5206":[{"seq":"TCCGAGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCGAGTT"}, {"seq":"AACTCGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCGAGTT"}],
                "X5207":[{"seq":"TTCTCTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCTCTCG"}, {"seq":"CGAGAGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCTCTCG"}],
                "X5208":[{"seq":"ATTCTGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATTCTGGC"}, {"seq":"GCCAGAAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATTCTGGC"}],
                "X5209":[{"seq":"ACTGCTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTGCTAG"}, {"seq":"CTAGCAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTGCTAG"}],
                "X5210":[{"seq":"CATAACGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATAACGG"}, {"seq":"CCGTTATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATAACGG"}],
                "X5211":[{"seq":"CAGTCTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGTCTTC"}, {"seq":"GAAGACTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGTCTTC"}],
                "X5212":[{"seq":"TGCCTCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCCTCTT"}, {"seq":"AAGAGGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCCTCTT"}],
                "X5213":[{"seq":"ACTGTGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACTGTGTC"}, {"seq":"GACACAGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACTGTGTC"}],
                "X5214":[{"seq":"GTATTGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTATTGGC"}, {"seq":"GCCAATAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTATTGGC"}],
                "X5215":[{"seq":"CGATGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGATGCTT"}, {"seq":"AAGCATCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGATGCTT"}],
                "X5216":[{"seq":"AAGGCTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGGCTGA"}, {"seq":"TCAGCCTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGGCTGA"}],
                "X5217":[{"seq":"AGTCAGGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTCAGGA"}, {"seq":"TCCTGACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTCAGGA"}],
                "X5218":[{"seq":"CAGGTATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGGTATC"}, {"seq":"GATACCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGGTATC"}],
                "X5219":[{"seq":"TCTCCGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTCCGAT"}, {"seq":"ATCGGAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTCCGAT"}],
                "X5220":[{"seq":"TTCAGCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCAGCCT"}, {"seq":"AGGCTGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCAGCCT"}],
                "X5221":[{"seq":"TCTGAGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTGAGAG"}, {"seq":"CTCTCAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTGAGAG"}],
                "X5222":[{"seq":"TTAGGTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTAGGTCG"}, {"seq":"CGACCTAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTAGGTCG"}],
                "X5223":[{"seq":"CTCTGGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCTGGTT"}, {"seq":"AACCAGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCTGGTT"}],
                "X5224":[{"seq":"GCGTTCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCGTTCTA"}, {"seq":"TAGAACGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCGTTCTA"}],
                "X5225":[{"seq":"TCACGTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCACGTTC"}, {"seq":"GAACGTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCACGTTC"}],
                "X5226":[{"seq":"AGGATGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGATGGT"}, {"seq":"ACCATCCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGATGGT"}],
                "X5227":[{"seq":"GTGTTCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTGTTCCT"}, {"seq":"AGGAACAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTGTTCCT"}],
                "X5228":[{"seq":"GTAGCATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTAGCATC"}, {"seq":"GATGCTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTAGCATC"}],
                "X5229":[{"seq":"AGGATCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGATCTG"}, {"seq":"CAGATCCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGATCTG"}],
                "X5230":[{"seq":"GACAAGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GACAAGAG"}, {"seq":"CTCTTGTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GACAAGAG"}],
                "X5231":[{"seq":"TTACGGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTACGGCT"}, {"seq":"AGCCGTAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTACGGCT"}],
                "X5232":[{"seq":"GCTGTTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTGTTGT"}, {"seq":"ACAACAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTGTTGT"}],
                "X5233":[{"seq":"AACCGAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACCGAAG"}, {"seq":"CTTCGGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACCGAAG"}],
                "X5234":[{"seq":"TCTGCTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCTGCTCT"}, {"seq":"AGAGCAGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCTGCTCT"}],
                "X5235":[{"seq":"CTCAGCTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCAGCTA"}, {"seq":"TAGCTGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCAGCTA"}],
                "X5236":[{"seq":"CTTCACCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTCACCA"}, {"seq":"TGGTGAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTCACCA"}],
                "X5237":[{"seq":"GATCGTAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATCGTAC"}, {"seq":"GTACGATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATCGTAC"}],
                "X5238":[{"seq":"CTACAGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTACAGTG"}, {"seq":"CACTGTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTACAGTG"}],
                "X5239":[{"seq":"TCGAGTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGAGTGA"}, {"seq":"TCACTCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGAGTGA"}],
                "X5240":[{"seq":"CAAGTGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAAGTGCA"}, {"seq":"TGCACTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAAGTGCA"}],
                "X5241":[{"seq":"CGAGTATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGAGTATG"}, {"seq":"CATACTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGAGTATG"}],
                "X5242":[{"seq":"CGTAGGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTAGGTT"}, {"seq":"AACCTACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTAGGTT"}],
                "X5243":[{"seq":"GCCAGTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCCAGTAT"}, {"seq":"ATACTGGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCCAGTAT"}],
                "X5244":[{"seq":"ATGGAAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATGGAAGG"}, {"seq":"CCTTCCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATGGAAGG"}],
                "X5245":[{"seq":"AAGAGCCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGAGCCA"}, {"seq":"TGGCTCTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGAGCCA"}],
                "X5246":[{"seq":"TGCGTAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCGTAGA"}, {"seq":"TCTACGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCGTAGA"}],
                "X5247":[{"seq":"TACACGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TACACGCT"}, {"seq":"AGCGTGTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TACACGCT"}],
                "X5248":[{"seq":"CCTTCCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTTCCTT"}, {"seq":"AAGGAAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTTCCTT"}],
                "X5249":[{"seq":"ACCGCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCGCATA"}, {"seq":"TATGCGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCGCATA"}],
                "X5250":[{"seq":"TGGTCCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGTCCTT"}, {"seq":"AAGGACCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGTCCTT"}],
                "X5251":[{"seq":"CCATACGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCATACGT"}, {"seq":"ACGTATGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCATACGT"}],
                "X5252":[{"seq":"AACCTTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACCTTGG"}, {"seq":"CCAAGGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACCTTGG"}],
                "X5253":[{"seq":"CAAGGTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAAGGTCT"}, {"seq":"AGACCTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAAGGTCT"}],
                "X5254":[{"seq":"GCTTCGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTTCGAA"}, {"seq":"TTCGAAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTTCGAA"}],
                "X5255":[{"seq":"CGGAATAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGGAATAC"}, {"seq":"GTATTCCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGGAATAC"}],
                "X5256":[{"seq":"AACTGGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACTGGTG"}, {"seq":"CACCAGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACTGGTG"}],
                "X5257":[{"seq":"GCTTCTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTTCTTG"}, {"seq":"CAAGAAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTTCTTG"}],
                "X5258":[{"seq":"GCAATTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCAATTCG"}, {"seq":"CGAATTGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCAATTCG"}],
                "X5259":[{"seq":"AGGTCACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGTCACT"}, {"seq":"AGTGACCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGTCACT"}],
                "X5260":[{"seq":"CAGCGATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGCGATT"}, {"seq":"AATCGCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGCGATT"}],
                "X5261":[{"seq":"AACCTCCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACCTCCT"}, {"seq":"AGGAGGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACCTCCT"}],
                "X5262":[{"seq":"TCGACATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGACATC"}, {"seq":"GATGTCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGACATC"}],
                "X5263":[{"seq":"CTGGTTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGGTTCT"}, {"seq":"AGAACCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGGTTCT"}],
                "X5264":[{"seq":"ACAGCAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACAGCAAC"}, {"seq":"GTTGCTGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACAGCAAC"}],
                "X5265":[{"seq":"GCATACAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCATACAG"}, {"seq":"CTGTATGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCATACAG"}],
                "X5266":[{"seq":"CATCTACG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATCTACG"}, {"seq":"CGTAGATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATCTACG"}],
                "X5267":[{"seq":"TTGTCGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGTCGGT"}, {"seq":"ACCGACAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGTCGGT"}],
                "X5268":[{"seq":"TAGCCGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGCCGAA"}, {"seq":"TTCGGCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGCCGAA"}],
                "X5269":[{"seq":"AGGCATAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGGCATAG"}, {"seq":"CTATGCCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGGCATAG"}],
                "X5270":[{"seq":"TTGACAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTGACAGG"}, {"seq":"CCTGTCAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTGACAGG"}],
                "X5271":[{"seq":"TGCACCAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCACCAA"}, {"seq":"TTGGTGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCACCAA"}],
                "X5272":[{"seq":"CCAGTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAGTGTT"}, {"seq":"AACACTGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAGTGTT"}],
                "X5273":[{"seq":"TGTCCAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTCCAGA"}, {"seq":"TCTGGACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTCCAGA"}],
                "X5274":[{"seq":"GATTGGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATTGGAG"}, {"seq":"CTCCAATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATTGGAG"}],
                "X5275":[{"seq":"ACGGTCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGGTCTT"}, {"seq":"AAGACCGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGGTCTT"}],
                "X5276":[{"seq":"CTGCGTAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGCGTAT"}, {"seq":"ATACGCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGCGTAT"}],
                "X5277":[{"seq":"CACCACTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACCACTA"}, {"seq":"TAGTGGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACCACTA"}],
                "X5278":[{"seq":"TGTGGTAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTGGTAC"}, {"seq":"GTACCACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTGGTAC"}],
                "X5279":[{"seq":"ACATAGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACATAGGC"}, {"seq":"GCCTATGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACATAGGC"}],
                "X5280":[{"seq":"CAAGCAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAAGCAGT"}, {"seq":"ACTGCTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAAGCAGT"}],
                "X5281":[{"seq":"GCACGTAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCACGTAA"}, {"seq":"TTACGTGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCACGTAA"}],
                "X5282":[{"seq":"TCGTAGTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGTAGTC"}, {"seq":"GACTACGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGTAGTC"}],
                "X5283":[{"seq":"CACTGACA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACTGACA"}, {"seq":"TGTCAGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACTGACA"}],
                "X5284":[{"seq":"CGTGTACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTGTACT"}, {"seq":"AGTACACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTGTACT"}],
                "X5285":[{"seq":"GAGCTCAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGCTCAA"}, {"seq":"TTGAGCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGCTCAA"}],
                "X5286":[{"seq":"ACGTCGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGTCGTA"}, {"seq":"TACGACGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGTCGTA"}],
                "X5287":[{"seq":"GTCTAGGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCTAGGT"}, {"seq":"ACCTAGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCTAGGT"}],
                "X5288":[{"seq":"CTTCGTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTCGTTC"}, {"seq":"GAACGAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTCGTTC"}],
                "X5289":[{"seq":"CTGAGATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGAGATC"}, {"seq":"GATCTCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGAGATC"}],
                "X5290":[{"seq":"GTGGATAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTGGATAG"}, {"seq":"CTATCCAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTGGATAG"}],
                "X5291":[{"seq":"ACCATCCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCATCCA"}, {"seq":"TGGATGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCATCCA"}],
                "X5292":[{"seq":"AGAGGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGAGGTTG"}, {"seq":"CAACCTCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGAGGTTG"}],
                "X5293":[{"seq":"ATGCCAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATGCCAAC"}, {"seq":"GTTGGCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATGCCAAC"}],
                "X5294":[{"seq":"GTTAAGGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTAAGGC"}, {"seq":"GCCTTAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTAAGGC"}],
                "X5295":[{"seq":"ACCTGACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCTGACT"}, {"seq":"AGTCAGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCTGACT"}],
                "X5296":[{"seq":"GCCACTTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCCACTTA"}, {"seq":"TAAGTGGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCCACTTA"}],
                "X5297":[{"seq":"ACCTCTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACCTCTGT"}, {"seq":"ACAGAGGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACCTCTGT"}],
                "X5298":[{"seq":"GATCCATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GATCCATG"}, {"seq":"CATGGATC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GATCCATG"}],
                "X5299":[{"seq":"CGCTTAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGCTTAAC"}, {"seq":"GTTAAGCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGCTTAAC"}],
                "X5300":[{"seq":"TGTCTGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTCTGCT"}, {"seq":"AGCAGACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTCTGCT"}],
                "X5301":[{"seq":"ATAAGGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATAAGGCG"}, {"seq":"CGCCTTAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATAAGGCG"}],
                "X5302":[{"seq":"CGTCTTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTCTTGT"}, {"seq":"ACAAGACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTCTTGT"}],
                "X5303":[{"seq":"CTCCATGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCCATGT"}, {"seq":"ACATGGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCCATGT"}],
                "X5304":[{"seq":"TGTTCGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTTCGAG"}, {"seq":"CTCGAACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTTCGAG"}],
                "X5305":[{"seq":"AGCAAGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCAAGCA"}, {"seq":"TGCTTGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCAAGCA"}],
                "X5306":[{"seq":"CGTATTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTATTCG"}, {"seq":"CGAATACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTATTCG"}],
                "X5307":[{"seq":"AACGACGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACGACGT"}, {"seq":"ACGTCGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACGACGT"}],
                "X5308":[{"seq":"GTTGCGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTTGCGAT"}, {"seq":"ATCGCAAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTTGCGAT"}],
                "X5309":[{"seq":"AACGTGGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AACGTGGA"}, {"seq":"TCCACGTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AACGTGGA"}],
                "X5310":[{"seq":"CTGTGTTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGTGTTG"}, {"seq":"CAACACAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGTGTTG"}],
                "X5311":[{"seq":"TGATACGC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGATACGC"}, {"seq":"GCGTATCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGATACGC"}],
                "X5312":[{"seq":"GTCCTTCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTCCTTCT"}, {"seq":"AGAAGGAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTCCTTCT"}],
                "X5313":[{"seq":"ACAGACCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACAGACCT"}, {"seq":"AGGTCTGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACAGACCT"}],
                "X5314":[{"seq":"TGTTGTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTTGTGG"}, {"seq":"CCACAACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTTGTGG"}],
                "X5315":[{"seq":"CATCGTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CATCGTGA"}, {"seq":"TCACGATG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CATCGTGA"}],
                "X5316":[{"seq":"CAGGAGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGGAGAT"}, {"seq":"ATCTCCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGGAGAT"}],
                "X5317":[{"seq":"TAGGTAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGGTAGG"}, {"seq":"CCTACCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGGTAGG"}],
                "X5318":[{"seq":"AGTCGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTCGCTT"}, {"seq":"AAGCGACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTCGCTT"}],
                "X5319":[{"seq":"CGTTGAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTTGAGT"}, {"seq":"ACTCAACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTTGAGT"}],
                "X5320":[{"seq":"TTCCTGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCCTGTG"}, {"seq":"CACAGGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCCTGTG"}],
                "X5321":[{"seq":"TCGTCTCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGTCTCA"}, {"seq":"TGAGACGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGTCTCA"}],
                "X5322":[{"seq":"TAACGAGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAACGAGG"}, {"seq":"CCTCGTTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAACGAGG"}],
                "X5323":[{"seq":"CTGAAGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGAAGCT"}, {"seq":"AGCTTCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGAAGCT"}],
                "X5324":[{"seq":"ATTGCGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATTGCGTG"}, {"seq":"CACGCAAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATTGCGTG"}],
                "X5325":[{"seq":"CCAAGACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCAAGACT"}, {"seq":"AGTCTTGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCAAGACT"}],
                "X5326":[{"seq":"GCTGTAAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTGTAAG"}, {"seq":"CTTACAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTGTAAG"}],
                "X5327":[{"seq":"GAGTGGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGTGGTT"}, {"seq":"AACCACTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGTGGTT"}],
                "X5328":[{"seq":"AGCTTGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCTTGAG"}, {"seq":"CTCAAGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCTTGAG"}],
                "X5329":[{"seq":"ACGACAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACGACAGA"}, {"seq":"TCTGTCGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACGACAGA"}],
                "X5330":[{"seq":"TTCGCAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TTCGCAGT"}, {"seq":"ACTGCGAA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TTCGCAGT"}],
                "X5331":[{"seq":"CCGATGTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCGATGTA"}, {"seq":"TACATCGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCGATGTA"}],
                "X5332":[{"seq":"TGACGCAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGACGCAT"}, {"seq":"ATGCGTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGACGCAT"}],
                "X5333":[{"seq":"TCGCATTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCGCATTG"}, {"seq":"CAATGCGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCGCATTG"}],
                "X5334":[{"seq":"CGTGATCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGTGATCA"}, {"seq":"TGATCACG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGTGATCA"}],
                "X5335":[{"seq":"GTGAGCTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTGAGCTT"}, {"seq":"AAGCTCAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTGAGCTT"}],
                "X5336":[{"seq":"AGCGGAAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGCGGAAT"}, {"seq":"ATTCCGCT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGCGGAAT"}],
                "X5337":[{"seq":"CGAAGAAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGAAGAAC"}, {"seq":"GTTCTTCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGAAGAAC"}],
                "X5338":[{"seq":"CCGTATCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCGTATCT"}, {"seq":"AGATACGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCGTATCT"}],
                "X5339":[{"seq":"GTACTCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTACTCTC"}, {"seq":"GAGAGTAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTACTCTC"}],
                "X5340":[{"seq":"AGTGTTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTGTTGG"}, {"seq":"CCAACACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTGTTGG"}],
                "X5341":[{"seq":"TGAACCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGAACCTG"}, {"seq":"CAGGTTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGAACCTG"}],
                "X5342":[{"seq":"TCAAGGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCAAGGAC"}, {"seq":"GTCCTTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCAAGGAC"}],
                "X5343":[{"seq":"GTGCTTAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTGCTTAC"}, {"seq":"GTAAGCAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTGCTTAC"}],
                "X5344":[{"seq":"GTGGTGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTGGTGTT"}, {"seq":"AACACCAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTGGTGTT"}],
                "X5345":[{"seq":"GCTGACTA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCTGACTA"}, {"seq":"TAGTCAGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCTGACTA"}],
                "X5346":[{"seq":"TGCGAACT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGCGAACT"}, {"seq":"AGTTCGCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGCGAACT"}],
                "X5347":[{"seq":"AATACGCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AATACGCG"}, {"seq":"CGCGTATT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AATACGCG"}],
                "X5348":[{"seq":"ACACGGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACACGGTT"}, {"seq":"AACCGTGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACACGGTT"}],
                "X5349":[{"seq":"CTGCACTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTGCACTT"}, {"seq":"AAGTGCAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTGCACTT"}],
                "X5350":[{"seq":"TCAACTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCAACTGG"}, {"seq":"CCAGTTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCAACTGG"}],
                "X5351":[{"seq":"AGTTGGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTTGGCT"}, {"seq":"AGCCAACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTTGGCT"}],
                "X5352":[{"seq":"CAGGTTAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAGGTTAG"}, {"seq":"CTAACCTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAGGTTAG"}],
                "X5353":[{"seq":"ACACCAGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ACACCAGT"}, {"seq":"ACTGGTGT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ACACCAGT"}],
                "X5354":[{"seq":"TGGATCAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGATCAC"}, {"seq":"GTGATCCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGATCAC"}],
                "X5355":[{"seq":"TGACTTCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGACTTCG"}, {"seq":"CGAAGTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGACTTCG"}],
                "X5356":[{"seq":"GTGTCTGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GTGTCTGA"}, {"seq":"TCAGACAC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GTGTCTGA"}],
                "X5357":[{"seq":"AGTTACGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTTACGG"}, {"seq":"CCGTAACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTTACGG"}],
                "X5358":[{"seq":"ATCTCGCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATCTCGCT"}, {"seq":"AGCGAGAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATCTCGCT"}],
                "X5359":[{"seq":"GAAGGTTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAAGGTTC"}, {"seq":"GAACCTTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAAGGTTC"}],
                "X5360":[{"seq":"GAGCTTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GAGCTTGT"}, {"seq":"ACAAGCTC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GAGCTTGT"}],
                "X5361":[{"seq":"TCCAATCG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCAATCG"}, {"seq":"CGATTGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCAATCG"}],
                "X5362":[{"seq":"CGGTCATA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CGGTCATA"}, {"seq":"TATGACCG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CGGTCATA"}],
                "X5363":[{"seq":"TGGCTATC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGGCTATC"}, {"seq":"GATAGCCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGGCTATC"}],
                "X5364":[{"seq":"CAACGGAT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CAACGGAT"}, {"seq":"ATCCGTTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CAACGGAT"}],
                "X5365":[{"seq":"CTCCTAGA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCCTAGA"}, {"seq":"TCTAGGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCCTAGA"}],
                "X5366":[{"seq":"CCGGAATT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCGGAATT"}, {"seq":"AATTCCGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCGGAATT"}],
                "X5367":[{"seq":"TAGACGTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGACGTG"}, {"seq":"CACGTCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGACGTG"}],
                "X5368":[{"seq":"TGACTGAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGACTGAC"}, {"seq":"GTCAGTCA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGACTGAC"}],
                "X5369":[{"seq":"TAGAGCTC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TAGAGCTC"}, {"seq":"GAGCTCTA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TAGAGCTC"}],
                "X5370":[{"seq":"TCCGTGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCGTGAA"}, {"seq":"TTCACGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCGTGAA"}],
                "X5371":[{"seq":"CTCGATAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCGATAC"}, {"seq":"GTATCGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCGATAC"}],
                "X5372":[{"seq":"CTTACCTG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTACCTG"}, {"seq":"CAGGTAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTACCTG"}],
                "X5373":[{"seq":"ATGGCGAA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"ATGGCGAA"}, {"seq":"TTCGCCAT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"ATGGCGAA"}],
                "X5374":[{"seq":"TCCTACCT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCCTACCT"}, {"seq":"AGGTAGGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCCTACCT"}],
                "X5375":[{"seq":"CCTCAGTT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CCTCAGTT"}, {"seq":"AACTGAGG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CCTCAGTT"}],
                "X5376":[{"seq":"CTACTTGG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTACTTGG"}, {"seq":"CCAAGTAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTACTTGG"}],
                "X5377":[{"seq":"TCACAGCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TCACAGCA"}, {"seq":"TGCTGTGA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TCACAGCA"}],
                "X5378":[{"seq":"CACGTTGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CACGTTGT"}, {"seq":"ACAACGTG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CACGTTGT"}],
                "X5379":[{"seq":"AAGTCGAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AAGTCGAG"}, {"seq":"CTCGACTT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AAGTCGAG"}],
                "X5380":[{"seq":"TGTACCGT","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"TGTACCGT"}, {"seq":"ACGGTACA","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"TGTACCGT"}],
                "X5381":[{"seq":"CTCATCAG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTCATCAG"}, {"seq":"CTGATGAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTCATCAG"}],
                "X5382":[{"seq":"AGTCTCAC","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"AGTCTCAC"}, {"seq":"GTGAGACT","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"AGTCTCAC"}],
                "X5383":[{"seq":"CTTGGATG","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"CTTGGATG"}, {"seq":"CATCCAAG","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"CTTGGATG"}],
                "X5384":[{"seq":"GCCTATCA","instruments":["MiSeq","HiSeq 2000","HiSeq 2500","NovaSeq"],"seq_in_adapter":"GCCTATCA"}, {"seq":"TGATAGGC","instruments":["iSeq","MiniSeq","NextSeq","HiSeq 3000","HiSeq 4000"],"seq_in_adapter":"GCCTATCA"}],
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
