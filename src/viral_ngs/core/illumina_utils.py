#!/usr/bin/env python3
"""
Utilities for Illumina data handling, including directory management,
run information parsing, and sample sheet processing.
"""

import csv
import logging
import os
import re
import shutil
import tempfile
import xml.etree.ElementTree

import arrow
import pandas as pd

from . import file as util_file
from . import misc as util_misc

log = logging.getLogger(__name__)


# ============================
# ***  IlluminaDirectory   ***
# ============================


class IlluminaDirectory(object):
    """A class that handles Illumina data directories"""

    def __init__(self, uri):
        self.uri = uri
        self.path = None
        self.tempDir = None
        self.runinfo = None
        self.samplesheet = None

    def __enter__(self):
        self.load()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return 0

    def load(self):
        if self.path is None:
            if "://" in self.uri:
                raise NotImplementedError("boto s3 download here uri -> tarball")
            else:
                if os.path.isdir(self.uri):
                    self.path = self.uri
                else:
                    self._extract_tarball(self.uri)
            self._fix_path()

    def _fix_path(self):
        assert self.path is not None
        # this is not the correct root-level directory
        # sometimes this points to one level up
        while True:
            if os.path.isdir(
                os.path.join(self.path, "Data", "Intensities", "BaseCalls")
            ):
                # found it! self.path is correct
                break
            else:
                subdirs = list(
                    os.path.join(self.path, x)
                    for x in os.listdir(self.path)
                    if os.path.isdir(os.path.join(self.path, x))
                )
                if len(subdirs) == 1:
                    # follow the rabbit hole
                    self.path = subdirs[0]
                else:
                    # don't know where to go now!
                    raise Exception(
                        "cannot find Data/Intensities/BaseCalls/ inside %s (%s)"
                        % (self.uri, self.path)
                    )

    def _extract_tarball(self, tarfile):
        self.tempDir = tempfile.mkdtemp(prefix="IlluminaDirectory-")
        self.path = self.tempDir
        util_file.extract_tarball(tarfile, self.tempDir)

    def close(self):
        if self.tempDir:
            shutil.rmtree(self.tempDir)
            self.tempDir = None

    def get_RunInfo(self):
        if self.runinfo is None:
            runinfo_file = os.path.join(self.path, "RunInfo.xml")
            util_file.check_paths(runinfo_file)
            self.runinfo = RunInfo(runinfo_file)
        return self.runinfo

    def get_SampleSheet(self, only_lane=None, append_run_id=None, **kwargs):
        if self.samplesheet is None:
            samplesheet_file = os.path.join(self.path, "SampleSheet.csv")
            util_file.check_paths(samplesheet_file)
            self.samplesheet = SampleSheet(
                samplesheet_file,
                only_lane=only_lane,
                append_run_id=append_run_id,
                **kwargs
            )
        return self.samplesheet

    def get_intensities_dir(self):
        return os.path.join(self.path, "Data", "Intensities")

    def get_BCLdir(self):
        return os.path.join(self.get_intensities_dir(), "BaseCalls")


# ==================
# ***  RunInfo   ***
# ==================


class RunInfo(object):
    """A class that reads the RunInfo.xml file emitted by Illumina
    MiSeq and HiSeq machines.
    """

    def __init__(self, xml_fname):
        self.fname = xml_fname
        self.root = xml.etree.ElementTree.parse(xml_fname).getroot()

    def get_fname(self):
        return self.fname

    def get_run_id(self):
        return self.root[0].attrib["Id"]

    def get_flowcell_raw(self):
        return self.root[0].find("Flowcell").text

    def get_flowcell(self):
        fc = self.get_flowcell_raw()
        # slice in the case where the ID has a prefix of zeros
        if re.match(r"^0+-", fc):
            if "-" in fc:
                # miseq often adds a bunch of leading zeros and a dash in front
                fc = "-".join(fc.split("-")[1:])
        # >=5 to avoid an exception here
        # <= 15 to limit the bytes added to each bam record
        assert len(fc) >= 5, "The flowcell ID must be five or more characters in length"
        if len(fc) > 15:
            log.warning(
                "The provided flowcell ID is longer than 15 characters. Is that correct?"
            )
        return fc

    def _get_rundate_obj(self):
        """
        Access the text of the <Date> node in the RunInfo.xml file
        and returns an arrow date object.
        """
        rundate = self.root[0].find("Date").text
        datestring_formats = [
            "YYMMDD",
            "YYYYMMDD",
            "M/D/YYYY h:mm:ss A",
            "YYYY-MM-DDTHH:mm:ssZ",
        ]
        for datestring_format in datestring_formats:
            try:
                date_parsed = arrow.get(rundate, datestring_format)
                return date_parsed
            except arrow.parser.ParserError:
                pass
        raise arrow.parser.ParserError(
            "The date string seen in RunInfo.xml ('%s') did not match known Illumina formats: %s"
            % (rundate, datestring_formats)
        )

    def get_rundate_american(self):
        return str(self._get_rundate_obj().format("MM/DD/YYYY"))

    def get_rundate_iso(self):
        return str(self._get_rundate_obj().format("YYYY-MM-DD"))

    def get_machine(self):
        return self.root[0].find("Instrument").text

    def get_read_structure(self):
        reads = []
        for x in self.root[0].find("Reads").findall("Read"):
            order = int(x.attrib["Number"])
            read = x.attrib["NumCycles"] + (
                x.attrib["IsIndexedRead"] == "Y" and "B" or "T"
            )
            reads.append((order, read))
        return "".join([r for _, r in sorted(reads)])

    def num_reads(self):
        return sum(
            1
            for x in self.root[0].find("Reads").findall("Read")
            if x.attrib["IsIndexedRead"] == "N"
        )

    def get_lane_count(self):
        layout = self.root[0].find("FlowcellLayout")
        return int(layout.attrib["LaneCount"])

    def get_surface_count(self):
        layout = self.root[0].find("FlowcellLayout")
        return int(layout.attrib["SurfaceCount"])

    def get_swath_count(self):
        layout = self.root[0].find("FlowcellLayout")
        return int(layout.attrib["SwathCount"])

    def get_tile_count(self):
        layout = self.root[0].find("FlowcellLayout")
        return int(layout.attrib["TileCount"])

    def get_section_count(self):
        layout = self.root[0].find("FlowcellLayout")
        return int(layout.attrib.get("SectionPerLane", 1))

    def tile_count(self):
        lane_count    = self.get_lane_count()
        surface_count = self.get_surface_count()
        swath_count   = self.get_swath_count()
        tile_count    = self.get_tile_count()
        section_count = self.get_section_count()

        total_tile_count = (
            lane_count * surface_count * swath_count * tile_count * section_count
        )
        return total_tile_count

    def machine_model_from_tile_count(self):
        """Return machine name and lane count based on tile count"""
        tc = self.tile_count()

        machine = None
        if tc == 2:
            log.info("Detected %s tiles, interpreting as MiSeq nano run.", tc)
            machine = {"machine": "Illumina MiSeq", "lane_count": 1}
        elif tc == 8:
            log.info("Detected %s tiles, interpreting as MiSeq micro run.", tc)
            machine = {"machine": "Illumina MiSeq", "lane_count": 1}
        elif tc == 16:
            log.info("Detected %s tiles, interpreting as iSeq run.", tc)
            machine = {"machine": "Illumina iSeq 100", "lane_count": 1}
        elif tc == 28:
            log.info("Detected %s tiles, interpreting as MiSeq run.", tc)
            machine = {"machine": "Illumina MiSeq", "lane_count": 1}
        elif tc == 38:
            log.info("Detected %s tiles, interpreting as MiSeq run.", tc)
            machine = {"machine": "Illumina MiSeq", "lane_count": 1}
        elif tc == 32:
            log.info("Detected %s tiles, interpreting as NextSeq 1000/2000 P1 run.", tc)
            machine = {"machine": "NextSeq 1000/2000", "lane_count": 1}
        elif tc == 128:
            log.info("Detected %s tiles, interpreting as HiSeq2k run.", tc)
            machine = {"machine": "Illumina HiSeq 2500", "lane_count": 2}
        elif tc == 132:
            log.info("Detected %s tiles, interpreting as NextSeq 1000/2000 P2 run.", tc)
            machine = {"machine": "NextSeq 1000/2000", "lane_count": 1}
        elif tc == 264:
            log.info("Detected %s tiles, interpreting as NextSeq 2000 P3 run.", tc)
            machine = {"machine": "NextSeq 2000", "lane_count": 2}
        elif tc == 288:
            log.info("Detected %s tiles, interpreting as NextSeq 550 (mid-output) run.", tc)
            machine = {"machine": "NextSeq 550", "lane_count": 4}
        elif tc == 624:
            log.info("Detected %s tiles, interpreting as Illumina NovaSeq 6000 run.", tc)
            machine = {"machine": "Illumina NovaSeq 6000", "lane_count": 2}
        elif tc == 768:
            log.info("Detected %s tiles, interpreting as HiSeq2500 run.", tc)
            machine = {"machine": "Illumina HiSeq 2500", "lane_count": 8}
        elif tc == 864:
            log.info("Detected %s tiles, interpreting as NextSeq 550 (high-output) run.", tc)
            machine = {"machine": "NextSeq 550", "lane_count": 4}
        elif tc == 896:
            log.info("Detected %s tiles, interpreting as HiSeq4k run.", tc)
            machine = {"machine": "Illumina HiSeq 4000", "lane_count": 8}
        elif tc == 1408:
            log.info("Detected %s tiles, interpreting as Illumina NovaSeq 6000 run.", tc)
            machine = {"machine": "Illumina NovaSeq 6000", "lane_count": 2}
        elif tc == 3744:
            log.info("Detected %s tiles, interpreting as Illumina NovaSeq 6000 run.", tc)
            machine = {"machine": "Illumina NovaSeq 6000", "lane_count": 4}
        elif tc == 6272:
            log.info("Detected %s tiles, interpreting as Illumina NovaSeq X Plus run.", tc)
            machine = {"machine": "Illumina NovaSeq X Plus", "lane_count": 8}
        else:
            log.info("Tile count: %s tiles (unknown instrument type).", tc)
        return machine

    def get_flowcell_chemistry(self):
        guessed_sequencer = self.infer_sequencer_model()
        return guessed_sequencer["chemistry"]

    def get_flowcell_lane_count(self):
        guessed_sequencer = self.infer_sequencer_model()
        try:
            return self.get_lane_count()
        except Exception:
            return guessed_sequencer.get("lane_count", None)

    def get_machine_model(self):
        guessed_sequencer = self.infer_sequencer_model()
        return guessed_sequencer["machine"]

    @classmethod
    def get_machines_for_flowcell_id(cls, fcid):
        sequencer_by_fcid = []
        for key in cls.flowcell_to_machine_model_and_chemistry:
            if re.search(key, fcid):
                sequencer_by_fcid.append(
                    cls.flowcell_to_machine_model_and_chemistry[key]
                )
        return sequencer_by_fcid

    def infer_sequencer_model(self):
        fcid = self.get_flowcell_raw()
        sequencer_by_tile_count = self.machine_model_from_tile_count()
        sequencers_by_fcid = self.get_machines_for_flowcell_id(fcid)

        if len(sequencers_by_fcid) > 1:
            raise LookupError("Multiple sequencers possible: %s", fcid)

        log.debug("self.tile_count(): %s", self.tile_count())

        if len(sequencers_by_fcid) > 0:
            if (
                sequencer_by_tile_count is not None
                and sequencers_by_fcid[0]["machine"]
                != sequencer_by_tile_count["machine"]
            ):
                log.warning(
                    "Sequencer type inferred from flowcell ID: %s does not match sequencer inferred from tile count: %s; is this a new machine type?"
                    % (
                        sequencers_by_fcid[0]["machine"],
                        sequencer_by_tile_count["machine"],
                    )
                )
            return sequencers_by_fcid[0]
        elif sequencer_by_tile_count is not None:
            log.warning(
                "Sequencer type unknown flowcell ID: %s, yet sequencer type was inferred for tile count: %s; is this a new flowcell ID pattern?"
                % (fcid, self.tile_count())
            )
            return sequencer_by_tile_count
        else:
            log.warning(
                "Tile count: %s and flowcell ID: %s are both novel; is this a new machine type?"
                % (self.tile_count(), fcid)
            )
            return {"machine": "UNKNOWN", "lane_count": self.get_lane_count()}

    flowcell_to_machine_model_and_chemistry = {
        r"[A-Z,0-9]{5}AAXX": {
            "machine"    : "Illumina Genome Analyzer IIx",
            "chemistry"  : "All",
            "lane_count" : 8,
            "note"       : "",
        },
        r"[A-Z,0-9]{5}ABXX": {
            "machine"    : "Illumina HiSeq 2000",
            "chemistry"  : "V2 Chemistry",
            "lane_count" : 8,
            "note"       : "",
        },
        r"[A-Z,0-9]{5}ACXX": {
            "machine"    : "Illumina HiSeq 2000",
            "chemistry"  : "V3 Chemistry",
            "lane_count" : 8,
            "note"       : "Also used on transient 2000E",
        },
        r"[A-Z,0-9]{5}(?:ANXX|AN\w\w)": {
            "machine"    : "Illumina HiSeq 2500",
            "chemistry"  : "V4 Chemistry",
            "lane_count" : 8,
            "note"       : "High output",
        },
        r"[A-Z,0-9]{5}(?:ADXX|AD\w\w)": {
            "machine"    : "Illumina HiSeq 2500",
            "chemistry"  : "V1 Chemistry",
            "lane_count" : 2,
            "note"       : "Rapid run",
        },
        r"[A-Z,0-9]{5}AMXX": {
            "machine"    : "Illumina HiSeq 2500",
            "chemistry"  : "V2 Chemistry (beta)",
            "lane_count" : 2,
            "note"       : "Rapid run",
        },
        r"[A-Z,0-9]{5}(?:BCXX|BC\w\w)": {
            "machine"    : "Illumina HiSeq 2500",
            "chemistry"  : "V2 Chemistry",
            "lane_count" : 2,
            "note"       : "Rapid run",
        },
        r"[A-Z,0-9]{5}AFX\w": {
            "machine"    : "NextSeq 550",
            "chemistry"  : "Mid-Output NextSeq",
            "lane_count" : 4,
            "note"       : "",
        },
        r"[A-Z,0-9]{5}AGXX": {
            "machine"    : "NextSeq 550",
            "chemistry"  : "V1 Chemistry",
            "lane_count" : 4,
            "note"       : "High-output",
        },
        r"[A-Z,0-9]{5}(?:BGXX|BG\w\w)": {
            "machine"    : "NextSeq 550",
            "chemistry"  : "V2/V2.5 Chemistry",
            "lane_count" : 4,
            "note"       : "High-output",
        },
        r"[A-Z,0-9]{5}(?:BBXX|BB\w\w)": {
            "machine"    : "Illumina HiSeq 4000",
            "chemistry"  : "Illumina HiSeq 4000",
            "lane_count" : 8,
            "note"       : "",
        },
        r"[A-Z,0-9]{5}(?:ALXX:AL\w\w)": {
            "machine"    : "HiSeq X Ten",
            "chemistry"  : "V1/V2.5 Chemistry",
            "lane_count" : 8,
            "note"       : "",
        },
        r"[A-Z,0-9]{5}(?:CCXX:CC\w\w)": {
            "machine"    : "HiSeq X Ten",
            "chemistry"  : "V2/V2.5 Chemistry",
            "lane_count" : 8,
            "note"       : "",
        },
        r"[A-Z,0-9]{5}DR\w\w": {
            "machine"    : "Illumina NovaSeq 6000",
            "chemistry"  : "V1 Chemistry",
            "lane_count" : 2,
            "note"       : "S1/SP",
        },
        r"[A-Z,0-9]{5}DM\w\w": {
            "machine"    : "Illumina NovaSeq 6000",
            "chemistry"  : "V1 Chemistry",
            "lane_count" : 2,
            "note"       : "S2",
        },
        r"[A-Z,0-9]{5}DS\w\w": {
            "machine"    : "Illumina NovaSeq 6000",
            "chemistry"  : "V1 Chemistry",
            "lane_count" : 4,
            "note"       : "S4",
        },
        r"BNS417.*": {
            "machine"    : "Illumina iSeq 100",
            "chemistry"  : "V1",
            "lane_count" : 1,
            "note"       : "AKA Firefly",
        },
        r"[0-9]{9}-\w{5}": {
            "machine"    : "Illumina MiSeq",
            "chemistry"  : "V1/V2/V3 Chemistry",
            "lane_count" : 1,
            "note"       : "",
        },
    }


# ======================
# ***  SampleSheet   ***
# ======================


class SampleSheetError(Exception):
    def __init__(self, message, fname):
        super(SampleSheetError, self).__init__(
            "Failed to load {} ({})".format(fname, message)
        )


class SampleSheet(object):
    """A class that reads an Illumina SampleSheet.csv or alternative/simplified
    tab-delimited versions as well.
    """

    def __init__(
        self,
        infile,
        use_sample_name            = True,
        only_lane                  = None,
        allow_non_unique           = False,
        append_run_id              = None,
        collapse_duplicates        = False,
        rev_comp_barcode_2         = False,
        barcode_columns_to_revcomp = None
    ):
        self.fname            = infile
        self.use_sample_name  = use_sample_name
        if only_lane is not None:
            only_lane         = str(only_lane)
        self.only_lane        = only_lane
        self.allow_non_unique = allow_non_unique
        self.append_run_id    = append_run_id
        self.rows             = []
        self._rowsOriginal    = []
        self.duplicate_rows_collapsed             = False
        self.barcodes_revcomped_relative_to_input = False
        self.barcodes_revcomped_column_names      = set()

        self._detect_and_load_sheet(infile)

        barcode_columns_to_revcomp = barcode_columns_to_revcomp or []

        if rev_comp_barcode_2 or (barcode_columns_to_revcomp and isinstance(barcode_columns_to_revcomp, (list, tuple, set, str))):
            columns_to_revcomp  = ['barcode_2'] if rev_comp_barcode_2 else []
            columns_to_revcomp += barcode_columns_to_revcomp
            self.rev_comp_barcode_values(barcode_columns_to_revcomp=columns_to_revcomp, inplace=True)

        if self.can_be_collapsed:
            if not allow_non_unique:
                raise SampleSheetError("Duplicate indices found in sample sheet", infile)
            else:
                if collapse_duplicates:
                    self.collapse_sample_index_duplicates()

    def _detect_and_load_sheet(self, infile):
        if infile.endswith((".csv", ".csv.gz")):
            with util_file.open_or_gzopen(infile, "rU") as inf:
                header = None
                miseq_skip = False
                row_num = 0
                for line_no, line in enumerate(inf):
                    if line_no == 0:
                        line = line.replace("\ufeff", "")

                    if len(line.rstrip("\r\n").strip()) == 0:
                        continue
                    csv.register_dialect(
                        "samplesheet", quoting=csv.QUOTE_MINIMAL, escapechar="\\"
                    )
                    row = next(
                        csv.reader([line.strip().rstrip("\n")], dialect="samplesheet")
                    )
                    row = [item.strip() for item in row]
                    if miseq_skip:
                        if line.startswith("[Data]"):
                            miseq_skip = False
                    elif line.startswith("["):
                        miseq_skip = True
                    elif header is None:
                        header = row
                        if all(x in header for x in ["Sample_ID", "Index"]):
                            keymapper = {
                                "Sample_ID"   : "sample",
                                "Index"       : "barcode_1",
                                "Index2"      : "barcode_2",
                                "Sample_Name" : "sample_name",
                            }
                            header = list(map(keymapper.get, header))
                        elif "Sample_ID" in header:
                            keymapper = {
                                "Sample_ID"   : "sample",
                                "index"       : "barcode_1",
                                "index2"      : "barcode_2",
                                "Sample_Name" : "sample_name",
                            }
                            header = list(map(keymapper.get, header))
                        elif "SampleID" in header:
                            keymapper = {
                                "SampleID"    : "sample",
                                "Index"       : "barcode_1",
                                "Index2"      : "barcode_2",
                                "libraryName" : "library_id_per_sample",
                                "FCID"        : "flowcell",
                                "Lane"        : "lane",
                            }
                            header = list(map(keymapper.get, header))
                        elif len(row) == 3:
                            header = ["sample", "barcode_1", "barcode_2"]
                            if "sample" not in row[0].lower():
                                row_num += 1
                                self.rows.append({
                                    "sample"    : row[0],
                                    "barcode_1" : row[1],
                                    "barcode_2" : row[2],
                                    "row_num"   : str(row_num),
                                })
                        else:
                            raise SampleSheetError("unrecognized filetype", infile)
                        for h in ("sample", "barcode_1"):
                            assert h in header
                    else:
                        row_num += 1
                        while len(row) < len(header):
                            row.append("")
                        assert len(header) == len(row)
                        row = dict((k, v) for k, v in zip(header, row) if k and v)
                        row["row_num"] = str(row_num)
                        if (
                            self.only_lane is not None
                            and row.get("lane")
                            and self.only_lane != row["lane"]
                        ):
                            continue
                        if ("sample" in row and row["sample"]) and (
                            "barcode_1" in row and row["barcode_1"]
                        ):
                            self.rows.append(row)
            if (
                self.use_sample_name
                and "sample_name" in header
                and all(row.get("sample_name") for row in self.rows)
            ):
                for row in self.rows:
                    row["library_id_per_sample"] = row["sample"]
                    row["sample"] = row["sample_name"]
            for row in self.rows:
                if "sample_name" in row:
                    del row["sample_name"]

        elif infile.endswith((".txt", ".txt.gz", ".tsv")):
            self.rows = []
            row_num = 0
            for row in util_file.read_tabfile_dict(infile):
                assert row.get("sample") and row.get("barcode_1")
                row_num += 1
                row["row_num"] = str(row_num)
                self.rows.append(row)
        else:
            raise SampleSheetError("unrecognized filetype", infile)

        if not self.rows:
            raise SampleSheetError("empty file", infile)

        for row in self.rows:
            row["library"] = row["sample"]
            if row.get("library_id_per_sample"):
                row["library"] += ".l" + row["library_id_per_sample"]
            row["run"] = row["library"]
        if len(set(row["run"] for row in self.rows)) != len(self.rows):
            if self.allow_non_unique:
                log.warning("non-unique library IDs in this lane")
                unique_count = {}
                for row in self.rows:
                    unique_count.setdefault(row["library"], 0)
                    unique_count[row["library"]] += 1
                    row["run"] += ".r" + str(unique_count[row["library"]])
            else:
                raise SampleSheetError("non-unique library IDs in this lane", infile)
        if self.append_run_id:
            for row in self.rows:
                row["run"] += "." + self.append_run_id

        for row in self.rows:
            row["sample_original"] = row["sample"]
            row["sample"]  = util_file.string_to_file_name(row["sample"])
            row["library"] = util_file.string_to_file_name(row["library"])
            row["run"]     = util_file.string_to_file_name(row["run"])

        if all(row.get("barcode_2") for row in self.rows):
            self.indexes = 2
        elif any(row.get("barcode_2") for row in self.rows):
            raise SampleSheetError(
                "inconsistent single/double barcoding in sample sheet", infile
            )
        else:
            self.indexes = 1

    @property
    def can_be_collapsed(self) -> bool:
        """Return True if duplicate barcodes exist that could be collapsed."""
        assert len(self.rows) > 0, "No sample sheet rows to collapse"
        df = pd.json_normalize(self.rows).astype(str).fillna("")
        grouping_cols = ["barcode_1"]
        if "barcode_2" in df.columns:
            grouping_cols.append("barcode_2")
        duplicated_mask = df.duplicated(subset=grouping_cols, keep=False)
        return duplicated_mask.any()

    def collapse_sample_index_duplicates(self, output_tsv=None, overwrite_instance_data=True):
        """Collapse duplicate barcode rows into single entries."""
        assert len(self.rows) > 0, "No sample sheet rows to collapse"
        hash_if_longer_than = 32

        df = pd.json_normalize(self.rows).astype(str).fillna("")
        original_cols = df.columns.tolist()

        grouping_cols = ["barcode_1"]
        has_barcode2 = "barcode_2" in df.columns
        if has_barcode2:
            grouping_cols.append("barcode_2")

        duplicated_mask = df.duplicated(subset=grouping_cols, keep=False)
        df_unique = df[~duplicated_mask].copy()
        df_duplicates = df[duplicated_mask].copy()
        grouped_dups = df_duplicates.groupby(grouping_cols, sort=False)

        collapsed_rows = []
        for group_keys, group_df in grouped_dups:
            if not isinstance(group_keys, tuple):
                group_keys = (group_keys,)

            b1 = group_keys[0]
            b2 = group_keys[1] if has_barcode2 else None

            if len(group_df) == 1:
                collapsed_rows.append(group_df.iloc[0])
            else:
                row_dict = {}
                if "barcode_1" in df.columns:
                    row_dict["barcode_1"] = b1
                if has_barcode2:
                    row_dict["barcode_2"] = b2

                if not has_barcode2:
                    row_dict["sample"] = f"{b1}"
                else:
                    if b2 and b2.strip() != "":
                        row_dict["sample"] = f"{b1}-{b2}"
                    else:
                        row_dict["sample"] = f"{b1}"

                for col in original_cols:
                    if col in ("barcode_1", "barcode_2", "sample"):
                        continue
                    col_values = group_df[col].tolist()
                    row_dict[col] = util_misc.collapse_dup_strs_to_str_or_md5(
                        col_values,
                        suffix="_muxed",
                        hash_if_longer_than=hash_if_longer_than
                    )

                row_dict["run"] = f"{row_dict['sample']}.l{row_dict['library_id_per_sample']}"
                if self.append_run_id:
                    row_dict["run"] += "." + self.append_run_id

                collapsed_rows.append(pd.Series(row_dict))

        df_collapsed = pd.DataFrame(collapsed_rows, dtype=str) if collapsed_rows else pd.DataFrame(columns=original_cols, dtype=str)
        out_df = pd.concat([df_collapsed, df_unique], ignore_index=True)
        out_df = out_df[original_cols]

        rows_collapsed = out_df.to_dict(orient="records")

        if len(rows_collapsed) < len(self.rows):
            log.info("%s: %i rows collapsed ----> %i rows", os.path.basename(self.fname), len(self.rows), len(rows_collapsed))
            self.duplicate_rows_collapsed = True
        else:
            log.info("%s: ZERO rows collapsed", os.path.basename(self.fname))

        excluded_tsv_output_columns = ('row_num', 'library', 'run', 'sample_original')
        filtered_cols = [col for col in original_cols if col not in excluded_tsv_output_columns]
        if output_tsv is not None:
            log.info("Saving collapsed sample sheet to: %s", os.path.realpath(output_tsv))
            out_df[filtered_cols].to_csv(output_tsv, sep="\t", index=False)

        if overwrite_instance_data:
            self._rowsOriginal = self.rows
            self.rows = rows_collapsed

        return rows_collapsed

    def inner_demux_mapper(self):
        """Build a DataFrame mapping barcode groups to barcode_3 values."""
        df = pd.json_normalize(self.rows).astype(str)
        has_barcode2 = ("barcode_2" in df.columns)
        columns = df.columns

        grouping_cols = ["barcode_1"]
        if has_barcode2:
            grouping_cols.append("barcode_2")

        def fill_inline_index_ids(group):
            if "Inline_Index_ID" not in group.columns:
                group["Inline_Index_ID"] = [str(i + 1) for i in range(len(group))]
            else:
                mask = (group["Inline_Index_ID"].isna()) | (group["Inline_Index_ID"] == "")
                next_ids = (str(i + 1) for i in range(len(group)))
                group.loc[mask, "Inline_Index_ID"] = [next(next_ids) for _ in range(sum(mask))]
            return group

        df = df.groupby(grouping_cols, group_keys=False, as_index=False)
        df = df[columns].apply(fill_inline_index_ids, include_groups=False)

        def build_run_string(row):
            sample_val = row.get("sample", "")
            lib_val = row.get("library_id_per_sample", "")
            run_str = f"{sample_val}.l{lib_val}"
            if self.append_run_id:
                run_str += f".{self.append_run_id}"
            return run_str

        df["run"] = df.apply(build_run_string, axis=1)

        def build_muxed_run_string(row):
            b1 = row.get("barcode_1", "")
            b2 = row.get("barcode_2", "") if has_barcode2 else None
            if b2 and b2.strip():
                muxed_sample = f"{b1}-{b2}"
            else:
                muxed_sample = b1
            lib_val = row.get("library_id_per_sample", "")
            muxed_run_str = f"{muxed_sample}.l{lib_val}"
            if self.append_run_id:
                muxed_run_str += f".{self.append_run_id}"
            return muxed_run_str

        df["muxed_run"] = df.apply(build_muxed_run_string, axis=1)

        out_cols = ["barcode_1"]
        if has_barcode2:
            out_cols.append("barcode_2")
        out_cols.extend(["barcode_3", "sample", "library_id_per_sample", "Inline_Index_ID", "run", "muxed_run"])

        existing_cols = [c for c in out_cols if c in df.columns]
        out_df = df[existing_cols].copy()
        out_df.set_index("sample", inplace=True)
        return out_df

    def make_barcodes_file(self, outFile):
        """Create input file for Picard ExtractBarcodes"""
        if self.num_indexes() == 2:
            header = ["barcode_name", "library_name", "barcode_sequence_1", "barcode_sequence_2"]
        else:
            header = ["barcode_name", "library_name", "barcode_sequence_1"]
        with open(outFile, "wt") as outf:
            outf.write("\t".join(header) + "\n")
            for row in self.rows:
                out = {
                    "barcode_sequence_1": row["barcode_1"],
                    "barcode_sequence_2": row.get("barcode_2", ""),
                    "barcode_name": row["sample"],
                    "library_name": row["library"],
                }
                outf.write("\t".join(out[h] for h in header) + "\n")

    def write_tsv(self, outFile, force=False):
        """Write sample sheet to a tab-delimited file."""
        if os.path.exists(outFile) and not force:
            raise FileExistsError(f"Output file {outFile} already exists. Use force=True to overwrite.")
        with open(outFile, "wt") as outf:
            writer = csv.DictWriter(outf, self.rows[0].keys(), delimiter="\t")
            writer.writeheader()
            for row in self.rows:
                writer.writerow(row)

    def rev_comp_barcode_values(self, barcode_columns_to_revcomp=None, inplace=True):
        """Reverse-complement barcode values in specified columns."""
        barcode_columns_to_revcomp = barcode_columns_to_revcomp or ["barcode_2"]

        if type(barcode_columns_to_revcomp) is str:
            barcode_columns_to_revcomp = [barcode_columns_to_revcomp]

        if inplace:
            self._rowsOriginal = self.rows.copy()
            for row_idx, row in enumerate(self.rows):
                for column_name in barcode_columns_to_revcomp:
                    if column_name in row:
                        try:
                            row[column_name] = util_misc.reverse_complement(row[column_name])
                        except Exception:
                            log.warning("Failed to reverse-complement barcode value on line %s of %s: '%s'",
                                       row_idx + 1, self.fname, row[column_name])
                            pass
                        self.barcodes_revcomped_relative_to_input = True
                        self.barcodes_revcomped_column_names.add(column_name)
            return self
        else:
            new_sheet_fp = util_file.mkstempfname(f'{os.path.basename(self.fname)}_rev-comped-{"-".join(barcode_columns_to_revcomp)}.txt')
            log.debug("Creating a new SampleSheet object with reverse-complemented barcodes: %s", new_sheet_fp)
            shutil.copyfile(self.fname, new_sheet_fp)
            new_ss = SampleSheet(
                new_sheet_fp,
                barcode_columns_to_revcomp=barcode_columns_to_revcomp,
                allow_non_unique=self.allow_non_unique
            )
            return new_ss

    def make_params_file(self, bamDir, outFile):
        """Create input file for Picard IlluminaBasecallsToXXX"""
        if self.num_indexes() == 2:
            header = ["OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME", "BARCODE_1", "BARCODE_2"]
        else:
            header = ["OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME", "BARCODE_1"]
        with open(outFile, "wt") as outf:
            outf.write("\t".join(header) + "\n")
            rows = self.rows + [{
                "barcode_1": "N",
                "barcode_2": "N",
                "sample": "Unmatched",
                "library": "Unmatched",
                "run": "Unmatched",
            }]
            for row in rows:
                out = {
                    "BARCODE_1": row["barcode_1"],
                    "BARCODE_2": row.get("barcode_2", ""),
                    "SAMPLE_ALIAS": row["sample"],
                    "LIBRARY_NAME": row["library"],
                }
                out["OUTPUT"] = os.path.join(bamDir, row["run"] + ".bam")
                outf.write("\t".join(out[h] for h in header) + "\n")

    def get_fname(self):
        return self.fname

    def get_rows(self):
        return self.rows

    def print_rows(self, row_indices=None):
        """Print rows, optionally only a subset."""
        rows_selected = [self.rows[i] for i in row_indices] if row_indices else self.rows
        for r in rows_selected:
            longest_key_len = len(max(r.keys(), key=len))
            print(f"Sample: {r['sample']}")
            for k, v in r.items():
                if k != 'sample':
                    print(f"\t{k:<{longest_key_len + 1}}: {v}")
            print("")

    def num_indexes(self):
        """Return 1 or 2 depending on whether pools are single or double indexed."""
        return self.indexes

    @property
    def num_samples(self):
        return len(self.rows)

    def fetch_by_index(self, idx):
        idx = str(idx)
        for row in self.rows:
            if idx == row["row_num"]:
                return row
        return None
