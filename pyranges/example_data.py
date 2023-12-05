"""Module of example data.

See Also
--------

pyranges.random : generate random PyRanges

Examples
--------

>>> pr.data.f1
Chromosome      Start      End  Name         Score  Strand
category        int64    int64  object       int64  category
------------  -------  -------  ---------  -------  ----------
chr1                3        6  interval1        0  +
chr1                5        7  interval2        0  -
chr1                8        9  interval3        0  +
PyRanges with 3 rows and 6 columns.
Contains 1 chromosomes and 2 strands.
"""
import functools
import importlib
import tempfile
from functools import cached_property
from importlib.resources import files
from pathlib import Path

import pandas as pd
import pkg_resources

import pyranges as pr

__all__ = [
    "f1",
    "f2",
    "chipseq",
    "chipseq_background",
    "aorta",
    "aorta2",
    "ensembl_gtf",
    "gencode_gtf",
    "ucsc_bed",
    "control_bam",
    "cpg",
    "exons",
    "chromsizes",
]


class ExampleData:
    """"""

    _files: dict[str, Path] = {}

    @classmethod
    @property
    def files(cls) -> dict[str, Path]:
        """Return a dict of the basenames to full paths of the example data in the project.

        Examples
        --------
        >>> bam = ExampleData.files["smaller.bam"]
        >>> bam.exists()
        True
        >>> bam == importlib.resources.files().joinpath("data/smaller.bam")
        True
        """
        if cls._files:
            return cls._files
        cls._files = {
            f.name: f
            for f in files("pyranges").joinpath("data").iterdir()
            if "__" not in f.name
        }
        return cls._files

    @staticmethod
    def _read_bed_from_string(contents):
        with tempfile.NamedTemporaryFile("w") as f:
            f.write(contents)
            f.flush()
            return pr.read_bed(f.name)

    @staticmethod
    def _read_gtf_from_string(contents):
        with tempfile.NamedTemporaryFile("w") as f:
            f.write(contents)
            f.flush()
            return pr.read_gtf(f.name)

    @cached_property
    def chipseq(self) -> "pr.PyRanges":
        contents = """chr8	28510032	28510057	U0	0	-
chr7	107153363	107153388	U0	0	-
chr5	135821802	135821827	U0	0	-
chr14	19418999	19419024	U0	0	-
chr12	106679761	106679786	U0	0	-
chr21	40099618	40099643	U0	0	+
chr8	22714402	22714427	U0	0	-
chr19	19571102	19571127	U0	0	+
chr3	140986358	140986383	U0	0	-
chr10	35419784	35419809	U0	0	-
chr4	98488749	98488774	U0	0	+
chr11	22225193	22225218	U0	0	+
chr1	38457520	38457545	U0	0	+
chr1	80668132	80668157	U0	0	-
chr2	152562484	152562509	U0	0	-
chr4	153155301	153155326	U0	0	+
chr9	120803448	120803473	U0	0	+
chr6	89296757	89296782	U0	0	-
chr1	194245558	194245583	U0	0	+
chr8	57916061	57916086	U0	0	+"""
        return self._read_bed_from_string(contents)

    @property
    def chipseq_background(self) -> "pr.PyRanges":
        contents = """chr7	20246668	20246693	U0	0	+
chr1	39036822	39036847	U0	0	+
chr19	47109000	47109025	U0	0	-
chr10	90059861	90059886	U0	0	-
chr3	55648137	55648162	U0	0	+
chr7	91135110	91135135	U0	0	+
chr13	100938475	100938500	U0	0	+
chr3	115816130	115816155	U0	0	+
chr19	43528773	43528798	U0	0	+
chr10	73781101	73781126	U0	0	+"""
        return self._read_bed_from_string(contents)

    @cached_property
    def chromsizes(self) -> "pr.PyRanges":
        contents = """chr1	0	249250621
chr2	0	243199373
chr3	0	198022430
chr4	0	191154276
chr5	0	180915260
chr6	0	171115067
chr7	0	159138663
chrX	0	155270560
chr8	0	146364022
chr9	0	141213431
chr10	0	135534747
chr11	0	135006516
chr12	0	133851895
chr13	0	115169878
chr14	0	107349540
chr15	0	102531392
chr16	0	90354753
chr17	0	81195210
chr18	0	78077248
chr20	0	63025520
chrY	0	59373566
chr19	0	59128983
chr22	0	51304566
chr21	0	48129895
chrM	0	16571"""
        return self._read_bed_from_string(contents)

    @cached_property
    def ensembl_gtf(self) -> "pr.PyRanges":
        """Example gtf file from Ensembl."""

        contents = """#!genome-build GRCh38.p10
#!genome-version GRCh38
#!genome-date 2013-12
#!genome-build-accession NCBI:GCA_000001405.25
#!genebuild-last-updated 2017-06
1	havana	gene	11869	14409	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
1	havana	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; tag "basic"; transcript_support_level "1";
1	havana	exon	11869	12227	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "1"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002234944"; exon_version "1"; tag "basic"; transcript_support_level "1";
1	havana	exon	12613	12721	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00003582793"; exon_version "1"; tag "basic"; transcript_support_level "1";
1	havana	exon	13221	14409	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002312635"; exon_version "1"; tag "basic"; transcript_support_level "1";
1	havana	exon	112700	112804	.	-	.	gene_id "ENSG00000238009"; gene_version "6"; transcript_id "ENST00000471248"; transcript_version "1"; exon_number "2"; gene_name "AL627309.1"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "AL627309.1-203"; transcript_source "havana"; transcript_biotype "lincRNA"; exon_id "ENSE00001957285"; exon_version "1"; transcript_support_level "5";
1	havana	exon	110953	111357	.	-	.	gene_id "ENSG00000238009"; gene_version "6"; transcript_id "ENST00000471248"; transcript_version "1"; exon_number "3"; gene_name "AL627309.1"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "AL627309.1-203"; transcript_source "havana"; transcript_biotype "lincRNA"; exon_id "ENSE00001879696"; exon_version "1"; transcript_support_level "5";
1	ensembl	transcript	120725	133723	.	-	.	gene_id "ENSG00000238009"; gene_version "6"; transcript_id "ENST00000610542"; transcript_version "1"; gene_name "AL627309.1"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "AL627309.1-205"; transcript_source "ensembl"; transcript_biotype "lincRNA"; tag "basic"; transcript_support_level "5";
1	ensembl	exon	133374	133723	.	-	.	gene_id "ENSG00000238009"; gene_version "6"; transcript_id "ENST00000610542"; transcript_version "1"; exon_number "1"; gene_name "AL627309.1"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "AL627309.1-205"; transcript_source "ensembl"; transcript_biotype "lincRNA"; exon_id "ENSE00003748456"; exon_version "1"; tag "basic"; transcript_support_level "5";
1	ensembl	exon	129055	129223	.	-	.	gene_id "ENSG00000238009"; gene_version "6"; transcript_id "ENST00000610542"; transcript_version "1"; exon_number "2"; gene_name "AL627309.1"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "AL627309.1-205"; transcript_source "ensembl"; transcript_biotype "lincRNA"; exon_id "ENSE00003734824"; exon_version "1"; tag "basic"; transcript_support_level "5";
1	ensembl	exon	120874	120932	.	-	.	gene_id "ENSG00000238009"; gene_version "6"; transcript_id "ENST00000610542"; transcript_version "1"; exon_number "3"; gene_name "AL627309.1"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "AL627309.1-205"; transcript_source "ensembl"; transcript_biotype "lincRNA"; exon_id "ENSE00003740919"; exon_version "1"; tag "basic"; transcript_support_level "5";"""
        return self._read_gtf_from_string(contents)

    @cached_property
    def f1(self) -> "pr.PyRanges":
        contents = """chr1	3	6	interval1	0	+
chr1	5	7	interval2	0	-
chr1	8	9	interval3	0	+"""
        return self._read_bed_from_string(contents)

    @cached_property
    def f2(self) -> "pr.PyRanges":
        contents = """chr1	1	2	a	0	+
chr1	6	7	b	0	-"""
        with tempfile.NamedTemporaryFile("w") as f:
            f.write(contents)
            f.flush()
            return pr.read_bed(f.name)

    @property
    def aorta(self) -> "pr.PyRanges":
        return pr.read_bed(self.files["aorta.bed"])

    @property
    def aorta2(self) -> "pr.PyRanges":
        return pr.read_bed(self.files["aorta2.bed"])


def get_example_path(basename) -> Path:
    full_path = pkg_resources.resource_filename("pyranges", "data/{}".format(basename))

    if full_path.endswith(".bam"):
        # hack to load index too
        pkg_resources.resource_filename("pyranges", "data/{}.bai".format(basename))

    return Path(full_path)


def aorta() -> "pr.PyRanges":
    """
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # | Chromosome   | Start     | End       | Name       | Score     | Strand       |
    >>> # | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
    >>> # |--------------+-----------+-----------+------------+-----------+--------------|
    >>> # | chr1         | 9939      | 10138     | H3K27me3   | 7         | +            |
    >>> # | chr1         | 9953      | 10152     | H3K27me3   | 5         | +            |
    >>> # | chr1         | 10024     | 10223     | H3K27me3   | 1         | +            |
    >>> # | chr1         | 10246     | 10445     | H3K27me3   | 4         | +            |
    >>> # | ...          | ...       | ...       | ...        | ...       | ...          |
    >>> # | chr1         | 9978      | 10177     | H3K27me3   | 7         | -            |
    >>> # | chr1         | 10001     | 10200     | H3K27me3   | 5         | -            |
    >>> # | chr1         | 10127     | 10326     | H3K27me3   | 1         | -            |
    >>> # | chr1         | 10241     | 10440     | H3K27me3   | 6         | -            |
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # Stranded PyRanges object has 11 rows and 6 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("aorta.bed")

    return pr.read_bed(full_path)


def aorta2() -> "pr.PyRanges":
    """
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # | Chromosome   | Start     | End       | Name       | Score     | Strand       |
    >>> # | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
    >>> # |--------------+-----------+-----------+------------+-----------+--------------|
    >>> # | chr1         | 10073     | 10272     | Input      | 1         | +            |
    >>> # | chr1         | 10280     | 10479     | Input      | 1         | +            |
    >>> # | chr1         | 16056     | 16255     | Input      | 1         | +            |
    >>> # | chr1         | 16064     | 16263     | Input      | 1         | +            |
    >>> # | ...          | ...       | ...       | ...        | ...       | ...          |
    >>> # | chr1         | 10079     | 10278     | Input      | 1         | -            |
    >>> # | chr1         | 10082     | 10281     | Input      | 1         | -            |
    >>> # | chr1         | 10149     | 10348     | Input      | 1         | -            |
    >>> # | chr1         | 19958     | 20157     | Input      | 1         | -            |
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # Stranded PyRanges object has 10 rows and 6 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("aorta2.bed")

    return pr.read_bed(full_path)


def bw() -> "pr.PyRanges":
    full_path = get_example_path("bw.bw")

    return pr.read_bigwig(full_path)


def chipseq() -> "pr.PyRanges":
    """
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # | Chromosome   | Start     | End       | Name       | Score     | Strand       |
    >>> # | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
    >>> # |--------------+-----------+-----------+------------+-----------+--------------|
    >>> # | chr1         | 212609534 | 212609559 | U0         | 0         | +            |
    >>> # | chr1         | 169887529 | 169887554 | U0         | 0         | +            |
    >>> # | chr1         | 216711011 | 216711036 | U0         | 0         | +            |
    >>> # | chr1         | 144227079 | 144227104 | U0         | 0         | +            |
    >>> # | ...          | ...       | ...       | ...        | ...       | ...          |
    >>> # | chrY         | 15224235  | 15224260  | U0         | 0         | -            |
    >>> # | chrY         | 13517892  | 13517917  | U0         | 0         | -            |
    >>> # | chrY         | 8010951   | 8010976   | U0         | 0         | -            |
    >>> # | chrY         | 7405376   | 7405401   | U0         | 0         | -            |
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # Stranded PyRanges object has 10,000 rows and 6 columns from 24 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("chipseq.bed")

    return pr.read_bed(full_path)


def chipseq_background() -> "pr.PyRanges":
    """
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # | Chromosome   | Start     | End       | Name       | Score     | Strand       |
    >>> # | (category)   | (int64)   | (int64)   | (object)   | (int64)   | (category)   |
    >>> # |--------------+-----------+-----------+------------+-----------+--------------|
    >>> # | chr1         | 39036822  | 39036847  | U0         | 0         | +            |
    >>> # | chr1         | 224145989 | 224146014 | U0         | 0         | +            |
    >>> # | chr1         | 167802964 | 167802989 | U0         | 0         | +            |
    >>> # | chr1         | 69101066  | 69101091  | U0         | 0         | +            |
    >>> # | ...          | ...       | ...       | ...        | ...       | ...          |
    >>> # | chrY         | 11936866  | 11936891  | U0         | 0         | -            |
    >>> # | chrY         | 10629111  | 10629136  | U0         | 0         | -            |
    >>> # | chrY         | 10632456  | 10632481  | U0         | 0         | -            |
    >>> # | chrY         | 11918814  | 11918839  | U0         | 0         | -            |
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # Stranded PyRanges object has 10,000 rows and 6 columns from 25 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("chipseq_background.bed")

    return pr.read_bed(full_path)


def chromsizes() -> "pr.PyRanges":
    """
    >>> # +--------------+-----------+-----------+
    >>> # | Chromosome   | Start     | End       |
    >>> # | (category)   | (int64)   | (int64)   |
    >>> # |--------------+-----------+-----------|
    >>> # | chr1         | 0         | 249250621 |
    >>> # | chr2         | 0         | 243199373 |
    >>> # | chr3         | 0         | 198022430 |
    >>> # | chr4         | 0         | 191154276 |
    >>> # | ...          | ...       | ...       |
    >>> # | chrY         | 0         | 59373566  |
    >>> # | chrX         | 0         | 155270560 |
    >>> # | chrM         | 0         | 16571     |
    >>> # | chr22        | 0         | 51304566  |
    >>> # +--------------+-----------+-----------+
    >>> # Unstranded PyRanges object has 25 rows and 3 columns from 25 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome.
    """

    full_path = get_example_path("chromsizes.bed")

    return pr.read_bed(full_path)


def control_bam() -> "pr.PyRanges":
    """
    >>> # +--------------+-----------+-----------+--------------+------------+
    >>> # | Chromosome   | Start     | End       | Strand       | Flag       |
    >>> # | (category)   | (int64)   | (int64)   | (category)   | (uint16)   |
    >>> # |--------------+-----------+-----------+--------------+------------|
    >>> # | chr1         | 887771    | 887796    | +            | 16         |
    >>> # | chr1         | 994660    | 994685    | +            | 16         |
    >>> # | chr1         | 1770383   | 1770408   | +            | 16         |
    >>> # | chr1         | 1995141   | 1995166   | +            | 16         |
    >>> # | ...          | ...       | ...       | ...          | ...        |
    >>> # | chrY         | 57402214  | 57402239  | +            | 16         |
    >>> # | chrY         | 10643526  | 10643551  | -            | 0          |
    >>> # | chrY         | 11776321  | 11776346  | -            | 0          |
    >>> # | chrY         | 20557165  | 20557190  | -            | 0          |
    >>> # +--------------+-----------+-----------+--------------+------------+
    >>> # Stranded PyRanges object has 10,000 rows and 5 columns from 25 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("control.bam")

    return pr.read_bam(full_path)


def cpg() -> "pr.PyRanges":
    """
    >>> # +--------------+-----------+-----------+-----------+
    >>> # | Chromosome   | Start     | End       | CpG       |
    >>> # | (category)   | (int64)   | (int64)   | (int64)   |
    >>> # |--------------+-----------+-----------+-----------|
    >>> # | chrX         | 64181     | 64793     | 62        |
    >>> # | chrX         | 69133     | 70029     | 100       |
    >>> # | chrX         | 148685    | 149461    | 85        |
    >>> # | chrX         | 166504    | 167721    | 96        |
    >>> # | ...          | ...       | ...       | ...       |
    >>> # | chrY         | 28555535  | 28555932  | 32        |
    >>> # | chrY         | 28773315  | 28773544  | 25        |
    >>> # | chrY         | 59213794  | 59214183  | 36        |
    >>> # | chrY         | 59349266  | 59349574  | 29        |
    >>> # +--------------+-----------+-----------+-----------+
    >>> # Unstranded PyRanges object has 1,077 rows and 4 columns from 2 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome.
    """

    full_path = get_example_path("cpg.bed")

    df = pd.read_csv(
        full_path, sep="\t", header=None, names="Chromosome Start End CpG".split()
    )

    return pr.PyRanges(df)


def ensembl_gtf() -> "pr.PyRanges":
    """
    >>> # +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
    >>> # | Chromosome   | Source     | Feature      | Start     | End       | Score      | Strand       | Frame      | gene_biotype                       | +19   |
    >>> # | (category)   | (object)   | (category)   | (int64)   | (int64)   | (object)   | (category)   | (object)   | (object)                           | ...   |
    >>> # |--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------|
    >>> # | 1            | havana     | gene         | 11868     | 14409     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
    >>> # | 1            | havana     | transcript   | 11868     | 14409     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
    >>> # | 1            | havana     | exon         | 11868     | 12227     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
    >>> # | 1            | havana     | exon         | 12612     | 12721     | .          | +            | .          | transcribed_unprocessed_pseudogene | ...   |
    >>> # | ...          | ...        | ...          | ...       | ...       | ...        | ...          | ...        | ...                                | ...   |
    >>> # | 1            | havana     | gene         | 1173055   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
    >>> # | 1            | havana     | transcript   | 1173055   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
    >>> # | 1            | havana     | exon         | 1179364   | 1179555   | .          | -            | .          | lncRNA                             | ...   |
    >>> # | 1            | havana     | exon         | 1173055   | 1176396   | .          | -            | .          | lncRNA                             | ...   |
    >>> # +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+------------------------------------+-------+
    >>> # Stranded PyRanges object has 2,446 rows and 28 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    >>> # 19 hidden columns: gene_id, gene_name, gene_source, gene_version, tag, transcript_biotype, transcript_id, transcript_name, transcript_source, transcript_support_level, ... (+ 9 more.)
    """

    full_path = get_example_path("ensembl_human.gtf.gz")

    return pr.read_gtf(full_path)


def exons() -> "pr.PyRanges":
    """
    >>> # +--------------+-----------+-----------+----------------------------------------+-----------+--------------+
    >>> # | Chromosome   | Start     | End       | Name                                   | Score     | Strand       |
    >>> # | (category)   | (int64)   | (int64)   | (object)                               | (int64)   | (category)   |
    >>> # |--------------+-----------+-----------+----------------------------------------+-----------+--------------|
    >>> # | chrX         | 135721701 | 135721963 | NR_038462_exon_0_0_chrX_135721702_f    | 0         | +            |
    >>> # | chrX         | 135574120 | 135574598 | NM_001727_exon_2_0_chrX_135574121_f    | 0         | +            |
    >>> # | chrX         | 47868945  | 47869126  | NM_205856_exon_4_0_chrX_47868946_f     | 0         | +            |
    >>> # | chrX         | 77294333  | 77294480  | NM_000052_exon_17_0_chrX_77294334_f    | 0         | +            |
    >>> # | ...          | ...       | ...       | ...                                    | ...       | ...          |
    >>> # | chrY         | 15409586  | 15409728  | NR_047633_exon_3_0_chrY_15409587_r     | 0         | -            |
    >>> # | chrY         | 15478146  | 15478273  | NR_047634_exon_18_0_chrY_15478147_r    | 0         | -            |
    >>> # | chrY         | 15360258  | 15361762  | NR_047601_exon_0_0_chrY_15360259_r     | 0         | -            |
    >>> # | chrY         | 15467254  | 15467278  | NM_001258270_exon_13_0_chrY_15467255_r | 0         | -            |
    >>> # +--------------+-----------+-----------+----------------------------------------+-----------+--------------+
    >>> # Stranded PyRanges object has 1,000 rows and 6 columns from 2 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("exons.bed")

    return pr.read_bed(full_path)


def f1() -> "pr.PyRanges":
    """
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # | Chromosome   |     Start |       End | Name       |     Score | Strand       |
    >>> # | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
    >>> # |--------------+-----------+-----------+------------+-----------+--------------|
    >>> # | chr1         |         3 |         6 | interval1  |         0 | +            |
    >>> # | chr1         |         8 |         9 | interval3  |         0 | +            |
    >>> # | chr1         |         5 |         7 | interval2  |         0 | -            |
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # Stranded PyRanges object has 3 rows and 6 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("f1.bed")

    return pr.read_bed(full_path)


def f2() -> "pr.PyRanges":
    """
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # | Chromosome   |     Start |       End | Name       |     Score | Strand       |
    >>> # | (category)   |   (int64) |   (int64) | (object)   |   (int64) | (category)   |
    >>> # |--------------+-----------+-----------+------------+-----------+--------------|
    >>> # | chr1         |         1 |         2 | a          |         0 | +            |
    >>> # | chr1         |         6 |         7 | b          |         0 | -            |
    >>> # +--------------+-----------+-----------+------------+-----------+--------------+
    >>> # Stranded PyRanges object has 2 rows and 6 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("f2.bed")

    return pr.read_bed(full_path)


def gencode_gtf() -> "pr.PyRanges":
    """
    >>> # +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+-------------------+-------+
    >>> # | Chromosome   | Source     | Feature      | Start     | End       | Score      | Strand       | Frame      | gene_id           | +15   |
    >>> # | (category)   | (object)   | (category)   | (int64)   | (int64)   | (object)   | (category)   | (object)   | (object)          | ...   |
    >>> # |--------------+------------+--------------+-----------+-----------+------------+--------------+------------+-------------------+-------|
    >>> # | chr1         | HAVANA     | gene         | 11868     | 14409     | .          | +            | .          | ENSG00000223972.5 | ...   |
    >>> # | chr1         | HAVANA     | transcript   | 11868     | 14409     | .          | +            | .          | ENSG00000223972.5 | ...   |
    >>> # | chr1         | HAVANA     | exon         | 11868     | 12227     | .          | +            | .          | ENSG00000223972.5 | ...   |
    >>> # | chr1         | HAVANA     | exon         | 12612     | 12721     | .          | +            | .          | ENSG00000223972.5 | ...   |
    >>> # | ...          | ...        | ...          | ...       | ...       | ...        | ...          | ...        | ...               | ...   |
    >>> # | chr1         | HAVANA     | exon         | 1430549   | 1430662   | .          | -            | .          | ENSG00000225285.1 | ...   |
    >>> # | chr1         | HAVANA     | transcript   | 1430663   | 1434520   | .          | -            | .          | ENSG00000225285.1 | ...   |
    >>> # | chr1         | HAVANA     | exon         | 1434177   | 1434520   | .          | -            | .          | ENSG00000225285.1 | ...   |
    >>> # | chr1         | HAVANA     | exon         | 1430663   | 1430954   | .          | -            | .          | ENSG00000225285.1 | ...   |
    >>> # +--------------+------------+--------------+-----------+-----------+------------+--------------+------------+-------------------+-------+
    >>> # Stranded PyRanges object has 4,995 rows and 24 columns from 1 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    >>> # 15 hidden columns: gene_type, gene_name, level, havana_gene, transcript_id, transcript_type, transcript_name, transcript_support_level, tag, ... (+ 6 more.)
    """

    full_path = get_example_path("gencode_human.gtf.gz")

    return pr.read_gtf(full_path)


def ucsc_bed() -> "pr.PyRanges":
    """
    >>> # +--------------+-----------+-----------+------------+------------+-----------------+--------------+---------------+-------------------+
    >>> # | Chromosome   | Start     | End       | Feature    | gene_id    | transcript_id   | Strand       | exon_number   | transcript_name   |
    >>> # | (category)   | (int64)   | (int64)   | (object)   | (object)   | (float64)       | (category)   | (float64)     | (object)          |
    >>> # |--------------+-----------+-----------+------------+------------+-----------------+--------------+---------------+-------------------|
    >>> # | chr1         | 12776117  | 12788726  | gene       | AADACL3    | nan             | +            | nan           | nan               |
    >>> # | chr1         | 169075927 | 169101957 | gene       | ATP1B1     | nan             | +            | nan           | nan               |
    >>> # | chr1         | 6845383   | 7829766   | gene       | CAMTA1     | nan             | +            | nan           | nan               |
    >>> # | chr1         | 20915589  | 20945396  | gene       | CDA        | nan             | +            | nan           | nan               |
    >>> # | ...          | ...       | ...       | ...        | ...        | ...             | ...          | ...           | ...               |
    >>> # | chrX         | 152661096 | 152663330 | exon       | PNMA6E     | 260.0           | -            | 0.0           | NM_001351293      |
    >>> # | chrX         | 152661096 | 152666808 | transcript | PNMA6E     | 260.0           | -            | nan           | NM_001351293      |
    >>> # | chrX         | 152664164 | 152664378 | exon       | PNMA6E     | 260.0           | -            | 1.0           | NM_001351293      |
    >>> # | chrX         | 152666701 | 152666808 | exon       | PNMA6E     | 260.0           | -            | 2.0           | NM_001351293      |
    >>> # +--------------+-----------+-----------+------------+------------+-----------------+--------------+---------------+-------------------+
    >>> # Stranded PyRanges object has 5,519 rows and 9 columns from 30 chromosomes.
    >>> # For printing, the PyRanges was sorted on Chromosome and Strand.
    """

    full_path = get_example_path("ucsc_human.bed.gz")

    names = "Chromosome Start End Feature gene_id transcript_id Strand exon_number transcript_name".split()
    df = pd.read_csv(full_path, sep="\t", names=names)

    return pr.PyRanges(df)
