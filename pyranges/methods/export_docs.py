import importlib
import inspect
import pkgutil
from collections.abc import Iterable
from pathlib import Path


# --------------------------------------------------------------------------- #
# helpers                                                                     #
# --------------------------------------------------------------------------- #
def _iter_rst(project_root: Path) -> Iterable[tuple[str, str]]:
    """Yield every ``*.rst`` under docs/doc/tutorials (if present)."""
    for folder in ("docs", "doc", "tutorials"):
        base = project_root / folder
        if not base.exists():
            continue
        for path in base.rglob("*.rst"):
            yield str(path.relative_to(project_root)), path.read_text(encoding="utf-8")


def _public_from_module(module) -> Iterable[tuple[str, object]]:
    """Yield `(name, obj)` for symbols *defined in* ``module`` (no re-exports)."""
    for name, obj in vars(module).items():
        if name.startswith("_"):
            continue
        if getattr(obj, "__module__", None) == module.__name__:
            yield name, obj


def _iter_public_members(modname: str) -> Iterable[tuple[str, object]]:
    """Yield public `(name, obj)` pairs for *modname* and selected sub-packages."""
    try:
        mod = importlib.import_module(modname)
    except ModuleNotFoundError:
        return

    # 1 top-level names in the requested module
    yield from _public_from_module(mod)

    # 2  walk sub-packages, but only inside the public namespaces we care about
    if not hasattr(mod, "__path__"):
        return

    prefix = f"{modname}."
    allowed_prefixes = ("pyranges.seqs", "pyranges.orfs")

    for _, subname, _ in pkgutil.walk_packages(mod.__path__, prefix):
        if not subname.startswith(allowed_prefixes):
            continue
        try:
            subm = importlib.import_module(subname)
        except ImportError:
            continue
        yield from _public_from_module(subm)


def _format_doc(sig: str, obj) -> str:
    """Format *obj*'s docstring under a reST header."""
    doc = inspect.getdoc(obj) or ""
    return f"{'#' * 100}\n{sig}\n{'#' * len(sig)}\n{doc}\n"


def _iter_class_methods(
    cls,
    *,
    include_df: bool = False,
) -> Iterable[tuple[str, object]]:
    """Yield (name, obj) for methods defined **on cls itself**.

    If *include_df* is False we skip attributes not present in ``cls.__dict__``,
    i.e. those inherited unchanged from pandas DataFrame.
    """
    for n, o in inspect.getmembers(cls):
        if n.startswith("_"):
            continue
        if not callable(o):
            continue
        if include_df or n in cls.__dict__:  # defined or overridden here
            yield n, o


def _export_docs(
    to_file: str | Path | None = None,
    *,
    include_df: bool = False,
) -> str | None:
    pkg_root = Path(__file__).resolve().parents[1]  # …/pyranges
    pieces: list[str] = []

    # 1. RST sources
    for relpath, txt in _iter_rst(pkg_root.parent):
        pieces.append(f".. *** {relpath} ***\n\n{txt}")

    # 2. public API docstrings
    from pyranges import PyRanges, RangeFrame

    class_map = {"PyRanges": PyRanges, "RangeFrame": RangeFrame}
    public_modules = ["pyranges", "pyranges.seqs", "pyranges.orfs"]

    pieces.append("\n\n====================\nPublic API docstrings\n====================")

    for cls_name, cls in class_map.items():
        pieces.append(_format_doc(f"class {cls.__module__}.{cls_name}", cls))
        for meth_name, meth in _iter_class_methods(cls, include_df=include_df):
            pieces.append(_format_doc(f"{cls_name}.{meth_name}()", meth))

    for mod in public_modules:
        for name, obj in _iter_public_members(mod):
            if inspect.isfunction(obj):
                pieces.append(_format_doc(f"function {obj.__module__}.{name}()", obj))

    # 3. output
    full_text = "\n\n".join(pieces)

    if to_file is not None:
        path = Path(to_file)
        path.write_text(full_text, encoding="utf-8")
        return None

    return full_text
