import pytest


@pytest.fixture(autouse=True)
def set_max_console_width():
    # Assuming MAX_CONSOLE_WIDTH is a global variable in your module
    import pyranges

    pyranges.TOSTRING_CONSOLE_WIDTH = 120
