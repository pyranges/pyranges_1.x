from typing import Optional, Any

from tabulate import tabulate


def adjust_table_width(
    data: list,
    headers: list,
    dtypes: Optional[list] = None,
    max_total_width: int | None = None,
    max_col_width: int | None = None,
) -> tuple[list[list[str]], list[str]]:
    # Truncate individual columns to max_col_width
    truncated_data = []
    for row in data:
        new_row = []
        for item in row:
            item_str = str(item)
            if max_col_width is not None:
                new_row.append(item_str[:max_col_width])
            else:
                new_row.append(item_str)
        truncated_data.append(new_row)

    # Adjust column widths for headers if provided
    if headers and max_col_width is not None:
        headers = [h[:max_col_width] for h in headers]

    if dtypes and max_col_width is not None:
        dtypes = [d[:max_col_width] for d in dtypes]

    # Calculate the cumulative width of the columns after truncation
    cumulative_width = 0
    included_columns = 0
    for index, item in enumerate(zip(*truncated_data)):
        column_width = max(len(str(x)) for x in item) + 2  # Add 2 for padding
        if cumulative_width + column_width <= max_total_width:
            cumulative_width += column_width
            included_columns = index + 1
        else:
            break

    # Slice the data and headers to include only the columns that fit into the max_total_width
    final_data = [row[:included_columns] for row in truncated_data]
    final_headers = headers[:included_columns] if headers else None
    final_dtypes = dtypes[:included_columns] if dtypes else None

    final_headers = (
        [f"{h}\n{t}" for h, t in zip(final_headers, final_dtypes)]
        if dtypes
        else headers
    )

    return final_data, final_headers


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
