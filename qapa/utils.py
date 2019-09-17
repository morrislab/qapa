import sys
from tqdm import tqdm

def eprint(*args, **kwargs):
    # No longer used and will be deprecated in the future
    tag = kwargs.pop('tag', 'qapa')
    tqdm.write("[{}] {}".format(tag, *args), file=sys.stderr, **kwargs)
