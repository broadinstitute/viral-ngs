import argparse

from viral_ngs.core import cmd


def test_attach_main_cleans_docstring_indentation():
    def main_example():
        """
            Run an example command.

            This description is indented like normal function docstrings.
        """

    parser = argparse.ArgumentParser()
    cmd.attach_main(parser, main_example)

    assert parser.description == (
        "Run an example command.\n\n"
        "This description is indented like normal function docstrings."
    )
