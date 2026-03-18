""" A suite of functions to compute correlation matrices for geostatistical analysis.
    Default back-end: Python. 
    Optional back-end: Cython (user should specify). If fails, defaults to Python. 
"""

__author__ = 'A. Renmin Pretell Ductram'

def select_backend(backend):

    global _impl, __implementation__
    
    backend = backend.strip().lower()
    
    if backend == 'cython':
        try:
            import geostats.build_Mrho_cy as build_Mrho
            _impl = build_Mrho
            __implementation__ = 'Cython'
            return
        except ImportError:
            print('Cython implementation failed to load. Defaulting to Python.')
            backend = 'python'

    if backend == 'python':
        import geostats.build_Mrho_py as build_Mrho
        _impl = build_Mrho
        __implementation__ = 'Python'
        return
        
    raise ValueError('backend must be "python" or "cython"')


select_backend(backend="Python")

def get_backend():
    return __implementation__

def MrhoE_sam2sam(*args, **kwargs):
    return _impl.MrhoE_sam2sam(*args, **kwargs)

def MrhoEA_sam2sam(*args, **kwargs):
    return _impl.MrhoEA_sam2sam(*args, **kwargs)

def MrhoEAS_sam2sam(*args, **kwargs):
    return _impl.MrhoEAS_sam2sam(*args, **kwargs)

def MrhoE_sam2uns(*args, **kwargs):
    return _impl.MrhoE_sam2uns(*args, **kwargs)

def MrhoEA_sam2uns(*args, **kwargs):
    return _impl.MrhoEA_sam2uns(*args, **kwargs)

def MrhoEAS_sam2uns(*args, **kwargs):
    return _impl.MrhoEAS_sam2uns(*args, **kwargs)

def MrhoE_sam2grid(*args, **kwargs):
    return _impl.MrhoE_sam2grid(*args, **kwargs)

def MrhoEA_sam2grid(*args, **kwargs):
    return _impl.MrhoEA_sam2grid(*args, **kwargs)
