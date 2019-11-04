#!/usr/bin/env python

# This is intended to assist with parsing Sequin-format feature tables.
# 
# Note: Molecule scope, Organism scope, and qualifier conditions are not currently included here
#       Definitions and comments adapted from INSDC documentation.
#
# see the feature table key reference:
#   https://www.insdc.org/documents/feature-table#7.2
#   https://www.ncbi.nlm.nih.gov/Sequin/table.html
ft_types = {
    "assembly_gap":{
        # Definition            gap between two components of a genome or transcriptome assembly;
        "mandatory_qualifiers":["estimated_length","gap_type","linkage_evidence"]
        # Comment               the location span of the assembly_gap feature for an unknown gap has 
        #                       to be specified by the submitter; the specified gap length has to be 
        #                       reasonable (less or = 1000) and will be indicated as "n"'s in sequence.
        #                       However, the value for the estimated_length of assembly_gap features 
        #                       within a single (non-CON) transcriptome record must be an integer 
        #                       and can not be "unknown";
    },
    "C_region":{
        # Definition            constant region of immunoglobulin light and heavy 
        #                       chains, and T-cell receptor alpha, beta, and gamma 
        #                       chains; includes one or more exons depending on the 
        #                       particular chain
        "optional_qualifiers":[
            "allele",
            "citation",
            "db_xref",
            "experiment",
            "gene",
            "gene_synonym",
            "inference",
            "locus_tag",
            "map",
            "note",
            "old_locus_tag",
            "product",
            "pseudo",
            "pseudogene",
            "standard_name"
        ]
    },
    "CDS":{
        # Definition            coding sequence; sequence of nucleotides that
        #                       corresponds with the sequence of amino acids in a
        #                       protein (location includes stop codon); 
        #                       feature includes amino acid conceptual translation.
        "optional_qualifiers":[
            "allele",
            "artificial_location",
            "citation",
            "codon_start",
            "db_xref",
            "EC_number",
            "exception",
            "experiment",
            "function",
            "gene",
            "gene_synonym",
            "inference",
            "locus_tag",
            "map",
            "note",
            "number=unquoted text (single token)",
            "old_locus_tag",
            "operon",
            "product",
            "protein_id",
            "pseudo",
            "pseudogene",
            "ribosomal_slippage",
            "standard_name",
            "translation",
            "transl_except=(pos:<base_range>,aa:<amino_acid>)",
            "transl_table ",
            "trans_splicing"
        ]
        # Comment               codon_start has valid value of 1 or 2 or 3, indicating
        #                       the offset at which the first complete codon of a coding
        #                       feature can be found, relative to the first base of
        #                       that feature;
        #                       "transl_table defines the genetic code table used if",
        #                       other than the universal genetic code table;
        #                       genetic code exceptions outside the range of the specified
        #                       tables is reported in /transl_except qualifier;
        #                       "protein_id consists of a stable ID portion (3+5 format",
        #                       with 3 position letters and 5 numbers) plus a version 
        #                       number after the decimal point; when the protein 
        #                       sequence encoded by the CDS changes, only the version 
        #                       number of the /protein_id value is incremented; the
        #                       stable part of the /protein_id remains unchanged and as 
        #                       a result will permanently be associated with a given 
        #                       protein; 
    },
    "centromere":{
        # Definition            region of biological interest identified as a centromere and
        #                       which has been experimentally characterized;
        "optional_qualifiers":[
            "citation",
            "db_xref",
            "experiment",
            "inference",
            "note",
            "standard_name"
        ]
        # Comment               the centromere feature describes the interval of DNA 
        #                       that corresponds to a region where chromatids are held 
        #                       and a kinetochore is formed
    },
    "D-loop":{
        # Definition            displacement loop; a region within mitochondrial DNA in
        #                       which a short stretch of RNA is paired with one strand
        #                       of DNA, displacing the original partner DNA strand in
        #                       this region; also used to describe the displacement of a
        #                       region of one strand of duplex DNA by a single stranded
        #                       invader in the reaction catalyzed by RecA protein
        "optional_qualifiers":[
            "allele",
            "citation",
            "db_xref",
            "experiment",
            "gene",
            "gene_synonym",
            "inference",
            "locus_tag",
            "map",
            "note",
            "old_locus_tag"
        ]
    },
    "D_segment":{
        # Definition            Diversity segment of immunoglobulin heavy chain, and 
        #                       T-cell receptor beta chain;  
        "optional_qualifiers":[
            "allele",
            "citation",
            "db_xref",
            "experiment",
            "gene",
            "gene_synonym",
            "inference",
            "locus_tag",
            "map",
            "note",
            "old_locus_tag",
            "product",
            "pseudo",
            "pseudogene",
            "standard_name"
        ]
    },
    "exon":{
        # Definition            region of genome that codes for portion of spliced mRNA, 
        #                       rRNA and tRNA; may contain 5'UTR, all CDSs and 3' UTR; 
        "optional_qualifiers":[
            "allele",
            "citation",
            "db_xref",
            "EC_number",
            "experiment",
            "function",
            "gene",
            "gene_synonym",
            "inference",
            "locus_tag",
            "map",
            "note",
            "number=unquoted text (single token)",
            "old_locus_tag",
            "product",
            "pseudo",
            "pseudogene",
            "standard_name",
            "trans_splicing"
        ]
    },
    "gap":{
        # Definition            gap in the sequence
        "mandatory_qualifiers":[
            "estimated_length"
        ],
        "optional_qualifiers":[
            "experiment",
            "inference",
            "map",
            "note"
        ]
        # Comment               the location span of the gap feature for an unknown 
        #                       gap is 100 bp, with the 100 bp indicated as 100 "n"'s in 
        #                       the sequence.  Where estimated length is indicated by 
        #                       an integer, this is indicated by the same number of 
        #                       "n"'s in the sequence. 
        #                       No upper or lower limit is set on the size of the gap.
    },
    "gene":{
        # Definition            region of biological interest identified as a gene 
        #                       and for which a name has been assigned;
        "optional_qualifiers":[
            "allele",
            "citation",
            "db_xref",
            "experiment",
            "function",
            "gene",
            "gene_synonym",
            "inference",
            "locus_tag",
            "map",
            "note",
            "old_locus_tag",
            "operon",
            "product",
            "pseudo",
            "pseudogene",
            "phenotype",
            "standard_name",
            "trans_splicing"
        ]
        # Comment               the gene feature describes the interval of DNA that 
        #                       corresponds to a genetic trait or phenotype; the feature is,
        #                       by definition, not strictly bound to it's positions at the 
        #                       ends;  it is meant to represent a region where the gene is 
        #                       located.
    },
    "iDNA":{
        # Definition            intervening DNA; DNA which is eliminated through any of
        #                       several kinds of recombination;
        "optional_qualifiers":[
            "allele",
            "citation",
            "db_xref",
            "experiment",
            "function",
            "gene",
            "gene_synonym",
            "inference",
            "locus_tag",
            "map",
            "note",
            "number=unquoted text (single token)",
            "old_locus_tag",
            "standard_name"
        ]
        # Comment               e.g., in the somatic processing of immunoglobulin genes.
    },
    "intron":{
        # Definition            a segment of DNA that is transcribed, but removed from
        #                       within the transcript by splicing together the sequences
        #                       (exons) on either side of it;
        "optional_qualifiers":[
            "allele",
            "citation",
            "db_xref",
            "experiment",
            "function",
            "gene",
            "gene_synonym",
            "inference",
            "locus_tag",
            "map",
            "note",
            "number=unquoted text (single token)",
            "old_locus_tag",
            "pseudo",
            "pseudogene",
            "standard_name",
            "trans_splicing"
        ]
    },
    "J_segment":{
        # Definition            joining segment of immunoglobulin light and heavy 
        #                       chains, and T-cell receptor alpha, beta, and gamma 
        #                       chains;  
        "optional_qualifiers":[
            "allele",
            "citation",
            "db_xref",
            "experiment",
            "gene",
            "gene_synonym",
            "inference",
            "locus_tag",
            "map",
            "note",
            "old_locus_tag",
            "product",
            "pseudo",
            "pseudogene",
            "standard_name"
        ]
    },
    "mat_peptide":{
        # Definition            mature peptide or protein coding sequence; coding
        #                       sequence for the mature or final peptide or protein
        #                       product following post-translational modification; the
        #                       location does not include the stop codon (unlike the
        #                       corresponding CDS);
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "EC_number",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name"
        ]
    },
    "misc_binding":{
        # Definition            site in nucleic acid which covalently or non-covalently
        #                       binds another moiety that cannot be described by any
        #                       other binding key (primer_bind or protein_bind);
        "mandatory_qualifiers":["bound_moiety"],
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag"
        ]
        # Comment               note that feature key regulatory with /regulatory_class
        #                       should be used for ribosome binding sites.
    },
    "misc_difference":{
        # Definition            feature sequence is different from that presented 
        #                       in the entry and cannot be described by any other 
        #                       difference key (old_sequence, variation, or modified_base);
        "optional_qualifiers":[
            "allele",
             "citation",
             "clone",
             "compare",
             "db_xref",
             "experiment",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map=",
             "note",
             "old_locus_tag",
             "phenotype",
             "replace",
             "standard_name"
        ]
        # Comment               the misc_difference feature key should be used to 
        #                       describe variability that arises as a result of 
        #                       genetic manipulation (e.g. site directed mutagenesis);
        #                       use /replace="" to annotate deletion, e.g. 
        #                       misc_difference 412..433
        #                       "replace=""  ",
    },
    "misc_feature":{
        # Definition            region of biological interest which cannot be described
        #                       by any other feature key; a new or rare feature;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "number=unquoted text (single token)",
             "old_locus_tag",
             "phenotype",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name"
        ]
        # Comment               this key should not be used when the need is merely to 
        #                       mark a region in order to comment on it or to use it in 
        #                       another feature's location
    },
    "misc_recomb":{
        # Definition            site of any generalized, site-specific or replicative
        #                       recombination event where there is a breakage and
        #                       reunion of duplex DNA that cannot be described by other
        #                       recombination keys or qualifiers of source key 
        #                       (/proviral);
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "recombination_class",
             "standard_name"
        ]
    },
    "misc_RNA":{
        # Definition            any transcript or RNA product that cannot be defined by
        #                       other RNA keys (prim_transcript, precursor_RNA, mRNA,
        #                       5'UTR, 3'UTR, exon, CDS, sig_peptide, transit_peptide,
        #                       mat_peptide, intron, polyA_site, ncRNA, rRNA and tRNA);
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "operon",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name",
             "trans_splicing"
        ]
    },
    "misc_structure":{
        # Definition            any secondary or tertiary nucleotide structure or 
        #                       conformation that cannot be described by other Structure
        #                       keys (stem_loop and D-loop);
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "standard_name"
        ]
    },
    "mobile_element":{
        # Definition            region of genome containing mobile elements;
        "mandatory_qualifiers":["mobile_element_type"],
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "rpt_family",
             "rpt_type",
             "standard_name"
        ]
    },
    "modified_base":{
        # Definition            the indicated nucleotide is a modified nucleotide and
        #                       should be substituted for by the indicated molecule
        #                       (given in the mod_base qualifier value)
        "mandatory_qualifiers":["mod_base"],
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "frequency",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag"
        ]
        # Comment               value is limited to the restricted vocabulary for 
        #                       modified base abbreviations;
    },
    "mRNA":{
        # Definition            messenger RNA; includes 5'untranslated region (5'UTR),
        #                       coding sequences (CDS, exon) and 3'untranslated region
        #                       (3'UTR);
        "optional_qualifiers":[
            "allele",
             "artificial_location",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "operon",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name",
             "trans_splicing"
        ]
    },
    "ncRNA":{
        # Definition            a non-protein-coding gene, other than ribosomal RNA and
        #                       transfer RNA, the functional molecule of which is the RNA
        #                       transcript;
        "mandatory_qualifiers":["ncRNA_class"],
                      
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "operon",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name",
             "trans_splicing"
        ]
        # Comment               the ncRNA feature is not used for ribosomal and transfer
        #                       RNA annotation, for which the rRNA and tRNA feature keys
        #                       should be used, respectively;
    },
    "N_region":{
        # Definition            extra nucleotides inserted between rearranged 
        #                       immunoglobulin segments.
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name"
        ]
    },
    "old_sequence":{
        # Definition            the presented sequence revises a previous version of the
        #                       sequence at this location;
        "mandatory_qualifiers":["citation", "compare"],
        "optional_qualifiers":[
            "allele",
             "db_xref",
             "experiment",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "replace"
        ]
        # Comment               /replace="" is used to annotate deletion, e.g. 
        #                       old_sequence 12..15 "replace="" ",
        #                       NOTE: This feature key is not valid in entries/records
        #                       created from 15-Oct-2007.
    },
    "operon":{
        # Definition            region containing polycistronic transcript including a cluster of
        #                       genes that are under the control of the same regulatory sequences/promoter
        #                       and in the same biological pathway
        "mandatory_qualifiers":["operon"],
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "inference",
             "map",
             "note",
             "phenotype",
             "pseudo",
             "pseudogene",
             "standard_name"
        ]
    },
    "oriT":{
        # Definition            origin of transfer; region of a DNA molecule where transfer is
        #                       initiated during the process of conjugation or mobilization
        "optional_qualifiers":[
            "allele",
             "bound_moiety",
             "citation",
             "db_xref",
             "direction=value",
             "experiment",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "rpt_family",
             "rpt_type",
             "rpt_unit_range",
             "rpt_unit_seq",
             "standard_name"
        ]
        # Comment               rep_origin should be used for origins of replication; 
        #  "direction has legal values RIGHT, LEFT and BOTH, however only                ",
        #                       RIGHT and LEFT are valid when used in conjunction with the oriT  
        #                       feature;
        #                       origins of transfer can be present in the chromosome; 
        #                       plasmids can contain multiple origins of transfer
    },
    "polyA_site":{
        # Definition            site on an RNA transcript to which will be added adenine
        #                       residues by post-transcriptional polyadenylation;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag"
        ]
    },
    "precursor_RNA":{
        # Definition            any RNA species that is not yet the mature RNA product;
        #                       may include ncRNA, rRNA, tRNA, 5' untranslated region
        #                       (5'UTR), coding sequences (CDS, exon), intervening
        #                       sequences (intron) and 3' untranslated region (3'UTR);
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "operon",
             "product",
             "standard_name",
             "trans_splicing"
        ]
        # Comment               used for RNA which may be the result of 
        #                       post-transcriptional processing;  if the RNA in question 
        #                       is known not to have been processed, use the 
        #                       prim_transcript key.
    },
    "prim_transcript":{
        # Definition            primary (initial, unprocessed) transcript;  
        #                       may include ncRNA, rRNA, tRNA, 5' untranslated region
        #                       (5'UTR), coding sequences (CDS, exon), intervening
        #                       sequences (intron) and 3' untranslated region (3'UTR);
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "operon",
             "standard_name"
        ]
    },
    "primer_bind":{
        # Definition            non-covalent primer binding site for initiation of
        #                       replication, transcription, or reverse transcription;
        #                       includes site(s) for synthetic e.g., PCR primer elements;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "standard_name",
             "PCR_conditions"
        ]
        # Comment               used to annotate the site on a given sequence to which a primer 
        #                       molecule binds - not intended to represent the sequence of the
        #                       primer molecule itself; PCR components and reaction times may 
        #                       be stored under the "/PCR_conditions" qualifier; 
        #                       since PCR reactions most often involve pairs of primers,
        #                       a single primer_bind key may use the order() operator
        #                       with two locations, or a pair of primer_bind keys may be
        #                       used.
    },
    "propeptide":{
        # Definition            propeptide coding sequence; coding sequence for the domain of a 
        #                       proprotein that is cleaved to form the mature protein product.
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name"
        ]
    },
    "protein_bind":{
        # Definition            non-covalent protein binding site on nucleic acid;
        "mandatory_qualifiers":["bound_moiety"],
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "operon",
             "standard_name"
        ]
        # Comment               note that feature key regulatory with /regulatory_class
        #                       should be used for ribosome binding sites.
    },
    "regulatory":{
        # Definition            any region of sequence that functions in the regulation of
        #                       transcription, translation, replication or chromatin structure;
        "mandatory_qualifiers":["regulatory_class"],
        "optional_qualifiers":[
            "allele",
             "bound_moiety",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "operon",
             "phenotype",
             "pseudo",
             "pseudogene",
             "standard_name"
        ]
        # Comment               This feature has replaced the following Feature Keys on 15-DEC-2014:
        #                       enhancer, promoter, CAAT_signal, TATA_signal, -35_signal, -10_signal,
        #                       RBS, GC_signal, polyA_signal, attenuator, terminator, misc_signal.
    },              
    "repeat_region":{
        # Definition            region of genome containing repeating units;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "rpt_family",
             "rpt_type",
             "rpt_unit_range",
             "rpt_unit_seq",
             "satellite",
             "standard_name"
        ]
    },
    "rep_origin":{
        # Definition            origin of replication; starting site for duplication of
        #                       nucleic acid to give two identical copies; 
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "direction=value",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "standard_name"
        ]
        # Comment               /direction has valid values: RIGHT, LEFT, or BOTH.
    },
    "rRNA":{
        # Definition            mature ribosomal RNA; RNA component of the
        #                       ribonucleoprotein particle (ribosome) which assembles
        #                       amino acids into proteins.
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "operon",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name"
        ]
        # Comment               rRNA sizes should be annotated with the /product
        #                       qualifier.
    },
    "S_region":{
        # Definition            switch region of immunoglobulin heavy chains;  
        #                       involved in the rearrangement of heavy chain DNA leading 
        #                       to the expression of a different immunoglobulin class 
        #                       from the same B-cell;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name"
        ]
    },
    "sig_peptide":{
        # Definition            signal peptide coding sequence; coding sequence for an
        #                       N-terminal domain of a secreted protein; this domain is
        #                       involved in attaching nascent polypeptide to the
        #                       membrane leader sequence;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name"
        ]
    },
    "source":{
        # Definition            identifies the biological source of the specified span of
        #                       the sequence; this key is mandatory; more than one source
        #                       key per sequence is allowed; every entry/record will have, as a
        #                       minimum, either a single source key spanning the entire
        #                       sequence or multiple source keys, which together, span the
        #                       entire sequence.
        "mandatory_qualifiers":["organism","mol_type"],
        "optional_qualifiers":[
            "altitude",
             "bio_material",
             "cell_line",
             "cell_type",
             "chromosome",
             "citation",
             "clone",
             "clone_lib",
             "collected_by",
             "collection_date",
             "country",
             "cultivar",
             "culture_collection",
             "db_xref",
             "dev_stage",
             "ecotype",
             "environmental_sample",
             "focus",
             "germline",
             "haplogroup",
             "haplotype",
             "host",
             "identified_by",
             "isolate",
             "isolation_source",
             "lab_host",
             "lat_lon",
             "macronuclear",
             "map",
             "mating_type",
             "note",
             "organelle",
             "PCR_primers",
             "plasmid",
             "pop_variant",
             "proviral",
             "rearranged",
             "segment",
             "serotype",
             "serovar",
             "sex",
             "specimen_voucher",
             "strain",
             "sub_clone",
             "submitter_seqid",
             "sub_species",
             "sub_strain",
             "tissue_lib",
             "tissue_type",
             "transgenic",
             "type_material",
             "variety"
        ]
        # Comment               transgenic sequences must have at least two source feature
        #                       keys; in a transgenic sequence the source feature key
        #                       describing the organism that is the recipient of the DNA
        #                       must span the entire sequence;
        #                       see Appendix III /organelle for a list of <organelle_value>
    },
    "stem_loop":{
        # Definition            hairpin; a double-helical region formed by base-pairing
        #                       between adjacent (inverted) complementary sequences in a
        #                       single strand of RNA or DNA. 
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "operon",
             "standard_name"
        ]
    },
    "STS":{
        # Definition            sequence tagged site; short, single-copy DNA sequence
        #                       that characterizes a mapping landmark on the genome and
        #                       can be detected by PCR; a region of the genome can be
        #                       mapped by determining the order of a series of STSs;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "standard_name",
        ]
        # Comment               STS location to include primer(s) in primer_bind key or
        #                       primers.
    },
    "telomere":{
        # Definition            region of biological interest identified as a telomere 
        #                       and which has been experimentally characterized;
        "optional_qualifiers":[
            "citation",
             "db_xref",
             "experiment",
             "inference",
             "note",
             "rpt_type",
             "rpt_unit_range",
             "rpt_unit_seq",
             "standard_name"
        ]
        # Comment               the telomere feature describes the interval of DNA 
        #                       that corresponds to a specific structure at the end of   
        #                       the linear eukaryotic chromosome which is required for                
        #               the integrity and maintenance of the end; this region
        #                       is unique compared to the rest of the chromosome and 
        #                       represent the physical end of the chromosome;
    },
    "tmRNA":{
        # Definition            transfer messenger RNA; tmRNA acts as a tRNA first,
        #                       and then as an mRNA that encodes a peptide tag; the
        #                       ribosome translates this mRNA region of tmRNA and attaches
        #                       the encoded peptide tag to the C-terminus of the
        #                       unfinished protein; this attached tag targets the protein for
        #                       destruction or proteolysis;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name",
             "tag_peptide"
        ]
    },
    "transit_peptide":{
        # Definition            transit peptide coding sequence; coding sequence for an
        #                       N-terminal domain of a nuclear-encoded organellar
        #                       protein; this domain is involved in post-translational
        #                       import of the protein into the organelle;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name"
        ]
    },
    "tRNA":{
        # Definition            mature transfer RNA, a small RNA molecule (75-85 bases
        #                       long) that mediates the translation of a nucleic acid
        #                       sequence into an amino acid sequence;
        "optional_qualifiers":[
            "allele",
             "anticodon=(pos:<location>,aa:<amino_acid>,seq:<text>)",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "operon",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name",
             "trans_splicing"
        ]
    },
    "unsure":{
        # Definition            a small region of sequenced bases, generally 10 or fewer in its length, which 
        #                       could not be confidently identified. Such a region might contain called bases 
        #                       (A, T, G, or C), or a mixture of called-bases and uncalled-bases ('N').
        #                       The unsure feature should not be used when annotating gaps in genome assemblies.
        #                       Please refer to assembly_gap feature for gaps within the sequence of an assembled
        #                       genome. For annotation of gaps in other sequences than assembled genomes use the 
        #                       gap feature.
        "optional_qualifiers":[
            "allele",
             "citation",
             "compare",
             "db_xref",
             "experiment",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "replace"
        ]
        # Comment               use /replace="" to annotate deletion, e.g. 
        #                       unsure      11..15 "replace=""  ",
    },
    "V_region":{
        # Definition            variable region of immunoglobulin light and heavy
        #                       chains, and T-cell receptor alpha, beta, and gamma
        #                       chains;  codes for the variable amino terminal portion;
        #                       can be composed of V_segments, D_segments, N_regions,
        #                       and J_segments;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name",
        ]
    },
    "V_segment":{
        # Definition            variable segment of immunoglobulin light and heavy
        #                       chains, and T-cell receptor alpha, beta, and gamma
        #                       chains; codes for most of the variable region (V_region)
        #                       and the last few amino acids of the leader peptide;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "product",
             "pseudo",
             "pseudogene",
             "standard_name"
        ]
    },
    "variation":{
        # Definition            a related strain contains stable mutations from the same
        #                       gene (e.g., RFLPs, polymorphisms, etc.) which differ
        #                       from the presented sequence at this location (and
        #                       possibly others);
        "optional_qualifiers":[
            "allele",
             "citation",
             "compare",
             "db_xref",
             "experiment",
             "frequency",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "phenotype",
             "product",
             "replace",
             "standard_name"
        ]
        # Comment               used to describe alleles, RFLP's,and other naturally occurring 
        #                       mutations and  polymorphisms; variability arising as a result 
        #                       of genetic manipulation (e.g. site directed mutagenesis) should 
        #                       be described with the misc_difference feature;
        #                       use /replace="" to annotate deletion, e.g. 
    },
    "3'UTR":{
        # Definition            1) region at the 3' end of a mature transcript (following 
        #                       the stop codon) that is not translated into a protein;
        #                       2) region at the 3' end of an RNA virus (following the last stop
        #                       codon) that is not translated into a protein;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "standard_name",
             "trans_splicing"
        ]
    },
    "5'UTR":{
        # Definition            1) region at the 5' end of a mature transcript (preceding 
        #                       the initiation codon) that is not translated into a protein;
        #                       2) region at the 5' end of an RNA virus genome (preceding the first 
        #                       initiation codon) that is not translated into a protein;
        "optional_qualifiers":[
            "allele",
             "citation",
             "db_xref",
             "experiment",
             "function",
             "gene",
             "gene_synonym",
             "inference",
             "locus_tag",
             "map",
             "note",
             "old_locus_tag",
             "standard_name",
             "trans_splicing"
        ]
    }
 }