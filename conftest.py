import pytest
from hypothesis import settings


@pytest.fixture(autouse=True)
def set_max_console_width() -> None:
    # Assuming MAX_CONSOLE_WIDTH is a global variable in your module
    import pyranges

    pyranges.TOSTRING_CONSOLE_WIDTH = 120


def pytest_addoption(parser):
    parser.addoption(
        "--derandomize",
        action="store_true",
        default=False,
        help="Derandomize Hypothesis tests"
    )


def pytest_configure(config):
    if config.getoption("--derandomize"):
        settings.register_profile("derandomize", derandomize=True)
        settings.load_profile("derandomize")
