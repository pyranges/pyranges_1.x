import pytest

from tests.helpers import assert_df_equal

from pyranges import GRanges

import pandas as pd

from io import StringIO


@pytest.fixture
def exons():

    df = pd.read_table("tests/exon.txt", sep="\t", header=None)
    df.columns = "Chromosome Start End".split() + list(df.columns[3:])

    return GRanges(df)


@pytest.fixture
def introns():

    df = pd.read_table("tests/intron.txt", sep="\t", header=None)

    print(df.head())
    print(df.shape)
    df.columns = "Chromosome Start End".split() + list(df.columns[3:])
    print(df.columns)

    return GRanges(df)


@pytest.fixture
def f1():

    df = pd.read_table("tests/f1.bed", sep="\t", header=None, names="Chromosome  Start  End  Name Score Strand".split())

    return GRanges(df)


@pytest.fixture
def f2():

    df = pd.read_table("tests/f2.bed", sep="\t", header=None, names="Chromosome  Start  End  Name Score Strand".split())

    return GRanges(df)


def assert_df_equal(df1, df2):

    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    return df1.equals(df2)
