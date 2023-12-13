from . cimport pyconf

def version():
    cdef int major, minor, patch
    get_version(&major, &minor, &patch)
    return [major, minor, patch]