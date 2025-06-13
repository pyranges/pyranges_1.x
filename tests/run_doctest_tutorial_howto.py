import glob

import pytest


def main():
    rst_files = ["docs/tutorial.rst", *glob.glob("docs/how_to*.rst")]
    # Build the command with all .rst files
    # Set the root directory explicitly
    args = ["--rootdir=.", *rst_files]
    # Call pytest with the list of .rst files
    return pytest.main(args)


if __name__ == "__main__":
    import sys

    sys.exit(main())
