from __future__ import print_function
import sys
from tqdm import tqdm

def eprint(*args, **kwargs):
    tag = kwargs.pop('tag', 'qapa')
    tqdm.write("[{}] {}".format(tag, *args), file=sys.stderr, **kwargs)
