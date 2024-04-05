from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import contextlib

    from pyranges import PyRanges

    with contextlib.suppress(ImportError):
        from pyrle import PyRles  # type: ignore[import-not-found]


def _to_rle(
    ranges: "PyRanges",
    value_col: str | None = None,
    *,
    strand: bool = True,
    rpm: bool = False,
    **kwargs,
) -> "PyRles":
    try:
        from pyrle import PyRles  # type: ignore[import]
        from pyrle.methods import coverage  # type: ignore[import]
    except ImportError as e:
        msg = "Using the coverage method requires that pyrle is installed."
        raise ImportError(msg) from e

    ranges = ranges.remove_strand() if not strand else ranges
    _kwargs = {
        "strand": strand,
        "value_col": value_col,
        "sparse": {"self": False},
    }  # already sparse
    kwargs.update(_kwargs)

    result = {k: coverage(v, **kwargs) for k, v in ranges.groupby(ranges.loc_columns)}

    if rpm:
        multiplier = 1e6 / len(ranges)
        result = {k: v * multiplier for k, v in result.items()}

    return PyRles(result)
