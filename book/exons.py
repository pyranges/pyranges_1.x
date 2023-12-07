import pyranges as pr

exons = pr.example_data.exons()
cpg = pr.example_data.cpg()
cpg.join_overlaps(exons.remove_strand()).subset(lambda df: df.CpG > 25)[["CpG"]].assign(
    lambda df: df.CpG % 10, "CpGDecile"
)["chrX"].slack(500)
from piedpiper import Debug

with Debug():
    cpg.join_overlaps(exons.remove_strand()).subset(lambda df: df.CpG > 25)[["CpG"]].assign(
        lambda df: df.CpG % 10, "CpGDecile"
    )["chrX"].slack(500)
