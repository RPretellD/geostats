""" A suite of functions supporting geotatistical analysis.
    Default back-end: Python. 
    Optional back-end: Cython (user should specify). If fails, defaults to Python. 
"""

__author__ = 'A. Renmin Pretell Ductram'

def select_backend(backend):

    global _impl, __implementation__
    
    backend = backend.strip().lower()
    
    if backend == 'cython':
        try:
            import geostats.geostats_tools_cy as geostats_tools
            _impl = geostats_tools
            __implementation__ = 'Cython'
            return
        except ImportError:
            print('Cython implementation failed to load. Defaulting to Python.')
            backend = 'python'

    if backend == 'python':
        import geostats.geostats_tools_py as geostats_tools
        _impl = geostats_tools
        __implementation__ = 'Python'
        return
        
    raise ValueError('backend must be "python" or "cython"')


select_backend(backend="Python")

def get_backend():
    return __implementation__

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
