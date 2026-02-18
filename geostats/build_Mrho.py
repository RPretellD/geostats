""" A suite of functions to compute correlation matrices for geostatistical analysis.
    Default back-end: Python. 
    Optional back-end: Cython (user should specify). If fails, defaults to Python. 
"""

__author__ = 'A. Renmin Pretell Ductram'

def select_backend(use_cython=False):

    global _impl, __implementation__

    if use_cython:
        try:
            import geostats.build_Mrho_cy as build_Mrho
            _impl = build_Mrho
            __implementation__ = 'Cython'
            return
        except ImportError:
            print('Cython implementation failed to load. Defaulting to Python.')
            pass

    import geostats.build_Mrho_py as build_Mrho
    _impl = build_Mrho
    __implementation__ = 'Python'

select_backend(use_cython=False)

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
