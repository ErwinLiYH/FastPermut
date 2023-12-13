from . cimport pyconf

def get_clib_version():
    cdef int major, minor, patch
    pyconf.get_version(&major, &minor, &patch)
    return [major, minor, patch]