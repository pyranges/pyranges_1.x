from typing import TYPE_CHECKING

from ruranges import complement_overlaps_numpy  # type: ignore[import]

from pyranges.core.pyranges_helpers import factorize_binary, mypy_ensure_rangeframe

if TYPE_CHECKING:
    from pyranges.range_frame.range_frame import RangeFrame


def _complement_overlaps(
    df: "RangeFrame",
    df2: "RangeFrame",
    by: list[str],
    slack: int | None = None,
) -> "RangeFrame":
    if df.empty:
        return df

    factorized, factorized2 = factorize_binary(df, df2, by)

    indices = complement_overlaps_numpy(
        chrs=factorized,
        starts=df.Start.to_numpy(),
        ends=df.End.to_numpy(),
        chrs2=factorized2,
        starts2=df2.Start.to_numpy(),
        ends2=df2.End.to_numpy(),
        slack=slack,
    )

    return mypy_ensure_rangeframe(df.take(indices))
