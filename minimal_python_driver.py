import ctypes

lib = ctypes.cdll.LoadLibrary("./build/libvlmc.so")

train_kl = lib.train_kl
train_kl.argtypes = [
    ctypes.c_char_p,
    ctypes.c_char_p,
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_float,
]
train_kl.restype = ctypes.c_char_p

lib.train_ps.argtypes = [
    ctypes.c_char_p,
    ctypes.c_char_p,
    ctypes.c_size_t,
    ctypes.c_size_t,
]
lib.train_ps.restype = ctypes.c_char_p

lib.train_kl_parameters.argtypes = [
    ctypes.c_char_p,
    ctypes.c_char_p,
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_size_t,
]
lib.train_kl_parameters.restype = ctypes.c_char_p

if __name__ == '__main__':
    print(train_kl(b"KL trained", b"ACGTACGACACACTTGT", 5, 2, 0.2))
    print(lib.train_ps(b"PS trained", b"ACGTACGACACACTTGT", 5, 2))
    print(lib.train_kl_parameters(b"Parameter trained", b"ACGTACGACACACTTGT", 5, 2, 192))
