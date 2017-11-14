#!/usr/bin/env python

"""
    Related only to sequence data within this file: 
        "Oligonucleotide sequences Copyright 2016 Illumina, Inc. All rights reserved."
        See: https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pdf
    Python code below is under the license of this repository:
        https://raw.githubusercontent.com/broadinstitute/viral-ngs/master/LICENSE
"""

import re, functools

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
            "i7":{
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
