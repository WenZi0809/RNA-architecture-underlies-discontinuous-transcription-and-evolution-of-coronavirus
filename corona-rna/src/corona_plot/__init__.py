from . import corona_plot

_public_functions = []
for name in dir(corona_plot):
    if not name.startswith('_') and callable(getattr(corona_plot, name)):
        globals()[name] = getattr(corona_plot, name)
        _public_functions.append(name)

__all__ = _public_functions

__version__ = "0.1.0"