from __future__ import print_function
import sys

def eprint(*args, **kwargs):
    # No longer used and will be deprecated in the future
    tag = kwargs.pop('tag', 'qapa')
    print("[{}] {}".format(tag, *args), file=sys.stderr, **kwargs)
