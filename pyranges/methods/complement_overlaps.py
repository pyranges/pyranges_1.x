from typing import TYPE_CHECKING

from pyranges.core.pyranges_helpers import ensure_rangeframe, factorize_binary

if TYPE_CHECKING:
    from pyranges.range_frame.range_frame import RangeFrame


def _complement_overlaps(
    df: "RangeFrame",
    df2: "RangeFrame",
    by: list[str],
    slack: int | None = None,
) -> "RangeFrame":
    import ruranges

    if df.empty:
        return df

    factorized, factorized2 = factorize_binary(df, df2, by)

    indices = ruranges.complement_overlaps(
        groups=factorized,  # type: ignore[arg-type]
        starts=df.Start.to_numpy(),
        ends=df.End.to_numpy(),
        groups2=factorized2,  # type: ignore[arg-type]
        starts2=df2.Start.to_numpy(),
        ends2=df2.End.to_numpy(),
        slack=slack or 0,
    )

    return ensure_rangeframe(df.take(indices))  # type: ignore[arg-type]
