import pytest
from hypothesis import settings


@pytest.fixture(autouse=True)
def set_max_console_width() -> None:
    import pyranges

    pyranges.options.set_option("console_width", 120)
    # below: edit default dictionary so that reset_option() doesn't screw up testing
    pyranges.options.options_default["console_width"] = (120, pyranges.options.options_in_use["console_width"][1])


def pytest_addoption(parser):
    parser.addoption("--derandomize", action="store_true", default=False, help="Derandomize Hypothesis tests")


def pytest_configure(config):
    if config.getoption("--derandomize"):
        settings.register_profile("derandomize", derandomize=True)
        settings.load_profile("derandomize")
