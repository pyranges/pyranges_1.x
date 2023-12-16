import shutil
from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd
from tabulate import tabulate

import pyranges

if TYPE_CHECKING:
    import pyranges as pr

MAX_COLUMN_NAMES_TO_SHOW = 3
MAX_ROWS_TO_SHOW = 8
HALF_OF_MAX_ROWS_TO_SHOW = MAX_ROWS_TO_SHOW // 2


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
    if width := pyranges.TOSTRING_CONSOLE_WIDTH:
        return width
    return shutil.get_terminal_size().columns


def tostring(
    self: "pr.RangeFrame",
    max_col_width: int | None = None,
    max_total_width: int | None = None,
) -> str:
    """Return string representation."""
    number_index_levels = self.index.nlevels
    number_duplicated_indices = self.index.duplicated().sum()
    has_duplicated_index = number_duplicated_indices > 0
    _self = self.reset_index()

    truncation_marker = ["..."]
    if len(_self) > MAX_ROWS_TO_SHOW:
        head = [list(v) for _, v in _self.head(HALF_OF_MAX_ROWS_TO_SHOW).iterrows()]
        tail = [list(v) for _, v in _self.tail(HALF_OF_MAX_ROWS_TO_SHOW).iterrows()]
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
    truncated_dtypes = adjusted_data.truncated_dtypes
    if len(adjusted_data.truncated_headers) != len(_self.columns):
        num_not_shown = len(_self.columns) - len(truncated_headers)
        not_shown = [
            f'"{e}"'
            for e in _self.columns[
                adjusted_data.included_columns : adjusted_data.included_columns + MAX_COLUMN_NAMES_TO_SHOW
            ]
        ]
        if num_not_shown > MAX_COLUMN_NAMES_TO_SHOW:
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
        f"{_self.__class__.__name__} with {_self.shape[0]} rows, "
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


# Define your data
data = [
    ["Alice", "Developer", "some very long text that might not fit on the screen", 30],
    ["Bob", "Manager", "another very long piece of text", 27],
]

# Define headers
headers = ["Name", "Occupation", "Description", "Age"]

# Maximum width for any individual column
max_col_width = 20

# Maximum total width for the display
max_total_width = 60  # For example, for a terminal width of 60 characters
