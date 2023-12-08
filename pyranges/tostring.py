import shutil
from dataclasses import dataclass
from typing import TYPE_CHECKING

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
    truncation_marker = ["..."]
    if len(self) >= MAX_ROWS_TO_SHOW:
        head = [list(v) for _, v in self.head(HALF_OF_MAX_ROWS_TO_SHOW).iterrows()]
        tail = [list(v) for _, v in self.tail(HALF_OF_MAX_ROWS_TO_SHOW).iterrows()]
        data = [*head, truncation_marker * self.shape[1], *(tail if len(self) > MAX_ROWS_TO_SHOW else head + tail)]
    else:
        data = [list(v) for _, v in self.iterrows()]

    adjusted_data = adjust_table_width(
        data=data,
        headers=list(self.columns),
        dtypes=[str(t) for t in self.dtypes],
        max_col_width=max_col_width,
        max_total_width=console_width(max_total_width),
    )
    columns_not_shown = "."
    truncated_data = adjusted_data.truncated_data
    truncated_headers = adjusted_data.truncated_headers
    truncated_dtypes = adjusted_data.truncated_dtypes
    if len(adjusted_data.truncated_headers) != len(self.columns):
        num_not_shown = len(self.columns) - len(truncated_headers)
        not_shown = [
            f'"{e}"'
            for e in self.columns[
                adjusted_data.included_columns: adjusted_data.included_columns + MAX_COLUMN_NAMES_TO_SHOW
            ]
        ]
        if num_not_shown > MAX_COLUMN_NAMES_TO_SHOW:
            not_shown.append("...")
        columns_not_shown = f" ({num_not_shown} columns not shown: {', '.join(not_shown)})."
        truncated_data = [row + truncation_marker for row in truncated_data]
        truncated_headers += truncation_marker
        truncated_dtypes += truncation_marker
    headers_with_dtype = [f"{h}\n{d}" for h, d in zip(truncated_headers, truncated_dtypes, strict=True)]
    class_and_shape_info = f"{self.__class__.__name__} with {self.shape[0]} rows and {self.shape[1]} columns"
    return f"{tabulate(truncated_data, headers_with_dtype)}\n{class_and_shape_info}{columns_not_shown}"


def adjust_table_width(
    data: list,
    headers: list,
    dtypes: list,
    max_total_width: int | None = None,
    max_col_width: int | None = None,
) -> AdjustedTableData:
    # Truncate individual columns to max_col_width
    truncated_headers = truncate_data([headers], max_col_width)
    truncated_dtypes = truncate_data([dtypes], max_col_width)
    truncated_data = truncate_data(data, max_col_width)

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
