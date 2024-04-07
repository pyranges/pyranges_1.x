import typing
from collections.abc import Callable
from functools import wraps
from typing import Any

if typing.TYPE_CHECKING:
    from pyranges.core.pyranges_main import PyRanges


def decorate_to_pyranges_method(func: Callable[..., Any]) -> Callable[..., Any]:
    """Decorate a function ext.moduleX.something.func_name to be accessed also as PyRanges.moduleX.func_name."""

    @wraps(func)
    def wrapper(self: "PyRanges", *args, **kwargs) -> Any:
        return func(self.pyranges_instance, *args, **kwargs)

    return wrapper
