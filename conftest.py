def pytest_configure(config) -> None:
    import pyranges as pr

    pr.options.set_option("max_rows_to_show", 8)
    pr.options.set_option("max_column_names_to_show", 3)
    pr.options.set_option("console_width", 120)
