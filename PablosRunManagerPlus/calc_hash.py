#! /usr/bin/env python
import sys
import os
import hashlib
from ctypes import *

def main(argv):
    a = str(argv[0])
    a = a.strip()
    a = a.split('/')
    a = [i for i in a if i!='']
    a = '_'.join(a[-3:])
    m = hashlib.md5()
    m.update(a)
    dig = m.hexdigest()
    b = c_uint(long(dig,base=16))
    print b.value


if __name__ == "__main__":
    main(sys.argv[1:])
