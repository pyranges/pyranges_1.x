import shutil
from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd
from tabulate import tabulate

if TYPE_CHECKING:
    import pyranges as pr

from pyranges import options


@dataclass
class AdjustedTableData:
    truncated_data: list[list[str]]
    truncated_headers: list[str]
    truncated_dtypes: list[str]
    included_columns: int


def console_width(max_total_width: int | None = None) -> int:
    """Return console width."""
    if max_total_width is not None:
        return max_total_width
    if width := options.get_option("console_width"):
        return width
    return shutil.get_terminal_size().columns


def _find_aliases_for_columns(df: "pr.RangeFrame") -> dict[str, str]:
    import random

    overlapping_cols = _get_indices_with_same_name_as_columns(df)
    return {c: "".join(random.sample(str(c.__hash__()), len(c))) for c in overlapping_cols}


def _get_indices_with_same_name_as_columns(df: "pr.RangeFrame") -> list[str]:
    if not df.index.name and not any(df.index.names):
        return []

    columns = {*df.columns}
    if df.index.name and df.index.name in columns:
        return [str(df.index.name)]

    if cols_in_both := columns.intersection(df.index.names):
        return [*cols_in_both]

    return []


def tostring(
    self: "pr.RangeFrame",
    max_col_width: int | None = None,
    max_total_width: int | None = None,
) -> str:
    """Return string representation."""
    number_index_levels = self.index.nlevels
    number_duplicated_indices = self.index.duplicated().sum()
    has_duplicated_index = number_duplicated_indices > 0

    from pyranges import options

    max_rows_to_show = options.get_option("max_rows_to_show")
    max_column_names_to_show = options.get_option("max_column_names_to_show")

    half_of_max_rows_to_show = max_rows_to_show // 2
    _self = (
        pd.concat([self.head(half_of_max_rows_to_show), self.tail(half_of_max_rows_to_show + 1)])
        if len(self) > max_rows_to_show
        else pd.DataFrame(self.copy())
    )

    replacement_names = _find_aliases_for_columns(self)
    inverted_replacement_names = {v: k for k, v in replacement_names.items()}
    _self.columns = [replacement_names.get(c, c) for c in _self.columns]
    _self = pd.DataFrame(_self).reset_index()

    truncation_marker = ["..."]
    if len(_self) > max_rows_to_show:
        head = [list(v) for _, v in _self.head(half_of_max_rows_to_show).iterrows()]
        tail = [list(v) for _, v in _self.tail(half_of_max_rows_to_show).iterrows()]
        data = [*head, truncation_marker * _self.shape[1], *tail]
    else:
        data = [list(v) for _, v in _self.iterrows()]

    adjusted_data = adjust_table_width(
        data=data,
        headers=list(_self.columns),
        dtypes=[str(t) for t in _self.dtypes],
        max_col_width=max_col_width,
        max_total_width=console_width(max_total_width),
    )
    columns_not_shown = ""
    truncated_data = adjusted_data.truncated_data
    truncated_headers = adjusted_data.truncated_headers
    truncated_headers = [inverted_replacement_names.get(c, c) for c in truncated_headers]
    truncated_dtypes = adjusted_data.truncated_dtypes
    if len(truncated_headers) != len(_self.columns):
        num_not_shown = len(_self.columns) - len(truncated_headers)
        not_shown = [
            f'"{e}"'
            for e in _self.columns[
                adjusted_data.included_columns : adjusted_data.included_columns + max_column_names_to_show
            ]
        ]
        if num_not_shown > max_column_names_to_show:
            not_shown.append("...")
        columns_not_shown = f" ({num_not_shown} columns not shown: {', '.join(not_shown)})."
        truncated_data = [row + truncation_marker for row in truncated_data]
        truncated_headers += truncation_marker
        truncated_dtypes += truncation_marker
    headers_with_dtype = [f"{h}\n{d}" for h, d in zip(truncated_headers, truncated_dtypes, strict=True)]
    contains_duplicates_string = (
        f" (with {number_duplicated_indices} index duplicates)." if has_duplicated_index else "."
    )
    class_and_shape_info = (
        f"{self.__class__.__name__} with {self.shape[0]} rows, "
        f"{_self.shape[1] - number_index_levels} columns, and {number_index_levels} index columns{contains_duplicates_string}"
    )
    truncated_df = pd.DataFrame.from_records(truncated_data, columns=truncated_headers)
    headers_with_dtype.insert(number_index_levels, "|\n|")
    truncated_df.insert(number_index_levels, "|", "|")
    return f"{tabulate(truncated_df.to_numpy(), headers_with_dtype, showindex=False)}\n{class_and_shape_info}{columns_not_shown}"


def adjust_table_width(
    data: list[list[str]],
    headers: list[str],
    dtypes: list[str],
    max_total_width: int | None = None,
    max_col_width: int | None = None,
) -> AdjustedTableData:
    """Adjust table width to fit screen or given width."""
    # Truncate individual columns to max_col_width
    truncated_headers = truncate_data([headers], max_col_width)
    truncated_dtypes = truncate_data([dtypes], max_col_width)
    truncated_data = truncate_data(data, max_col_width)

    max_total_width = int(1e5) if max_total_width is None else max_total_width
    # Calculate the cumulative width of the columns after truncation
    cumulative_width = 0
    included_columns = 0
    for index, column in enumerate(zip(*(truncated_data + truncated_dtypes + truncated_headers), strict=True)):
        column_width = max(len(str(x)) + 4 for x in column)
        if (new_width := cumulative_width + column_width) <= max_total_width:
            cumulative_width = new_width
            included_columns = index + 1
        else:
            break

    # Slice the data and headers to include only the columns that fit into the max_total_width
    final_data = [row[:included_columns] for row in truncated_data]
    final_headers = truncated_headers[0][:included_columns]
    final_dtypes = truncated_dtypes[0][:included_columns]

    return AdjustedTableData(final_data, final_headers, final_dtypes, included_columns)


def truncate_data(data: list[list[str]], max_col_width: int | None) -> list[list[str]]:
    """Truncate data to max_col_width."""
    truncated_data = []
    for row in data:
        new_row = []
        for item in row:
            item_str = str(item)
            if max_col_width is not None:
                new_row.append(
                    item_str[: max_col_width - 3] + "..."
                    if len(item_str) > max_col_width
                    else item_str[:max_col_width],
                )
            else:
                new_row.append(item_str)
        truncated_data.append(new_row)
    return truncated_data
