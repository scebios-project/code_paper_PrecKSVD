clear classes
mod = py.importlib.import_module('GLS');
py.importlib.reload(mod);
%py.GLS.greedy_least_squares(7,7,7)
