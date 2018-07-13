

import logging

import os.path

_logger = logging.getLogger(__name__)



def insert_adj_filename(filename, adj, ext=None):
    head, tail = os.path.split(filename)
    root, mext = os.path.splitext(tail)
    if ext is not None:
        mext = ext
    return os.path.join(head, "{}_{}{}".format(root, adj, mext))


