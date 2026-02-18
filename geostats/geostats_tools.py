""" A suite of functions supporting geotatistical analysis.
    Default back-end: Python. 
    Optional back-end: Cython (user should specify). If fails, defaults to Python. 
"""

__author__ = 'A. Renmin Pretell Ductram'

def select_backend(use_cython=False):

    global _impl, __implementation__

    if use_cython:
        try:
            import geostats.geostats_tools_cy as geostats_tools
            _impl = geostats_tools
            __implementation__ = 'Cython'
            return
        except ImportError:
            print('Cython implementation failed to load. Defaulting to Python.')
            pass

    import geostats.geostats_tools_py as geostats_tools
    _impl = geostats_tools
    __implementation__ = 'Python'

select_backend(use_cython=False)

def get_E_epi_dist(*args, **kwargs):
    return _impl.get_E_epi_dist(*args, **kwargs)

def get_A_epi_dist(*args, **kwargs):
    return _impl.get_A_epi_dist(*args, **kwargs)

def get_Vs30_dist(*args, **kwargs):
    return _impl.get_Vs30_dist(*args, **kwargs)

def get_E_dist(*args, **kwargs):
    return _impl.get_E_dist(*args, **kwargs)

def get_A_dist(*args, **kwargs):
    return _impl.get_A_dist(*args, **kwargs)
