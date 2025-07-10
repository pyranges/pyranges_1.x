info = """Utilities to instruct a AI coding assistant for pyranges prompts.

Get a prompt to copy-paste into an AI assistant to prime it for pyranges coding tasks:
    >>> import pyranges as pr
    >>> pr.assistant.prompt()

Make a file with pyranges documentation to upload to the AI assistant:
    >>> pr.assistant.export_docs("pr_docs.txt")"""

prompt_default = """Act as an expert bioinformatician programmer experienced in pyranges (complete documentation attached for you to learn). Next, answer my requests for code by first explaining the workflow, followed by oneliner-style code snippets, as concise as possible but elegant, preceded by the text of task as commented code. """
prompt_add_concise = """Output code that is as concise as possible but elegant. Assume pyranges is fully installed. No import statements. Use aptly named variables, no need to declare them. """


class Assistant:
    """A class to instruct a AI coding assistant for pyranges prompts."""

    def __init__(self) -> None:
        pass

    def __str__(self) -> str:
        return info

    def __repr__(self) -> str:
        return self.__str__()

    def prompt(self, *, concise: bool = False) -> str:
        """Get an example prompt to use for the AI coding assistant.

        Parameters
        ----------
        concise : bool, default False
            If True, the prompt will request concise code snippets.

        Returns
        -------
        str
            The prompt string to use with the AI coding assistant.

        """
        if concise:
            return prompt_default + prompt_add_concise
        return prompt_default

    def export_docs(self, to_file=None, *, include_df=False) -> str | None:
        """Build a single string containing all RST sources **plus** all public PyRanges docstrings.

        Parameters
        ----------
        to_file : str or pathlib.Path, optional
            Destination path.  If *None* (default) the text is returned.
        include_df : bool, default False
            Also include methods inherited unchanged from :class:`pandas.DataFrame`.

        Returns
        -------
        str | None
            The documentation blob, or None if written to a file.

        Examples
        --------
        >>> import pyranges
        >>> pyranges.assistant.export_docs(to_file="pr_docs.txt")

        """
        from pyranges.methods.export_docs import _export_docs

        return _export_docs(to_file=to_file, include_df=include_df)


assistant = Assistant()
