import time
import os
import glob
import gzip

from pathlib import Path
from collections import defaultdict

import dask.bag as db
import numpy as np


def cluster_writer(prefix, partition_index, data):
    """Takes an iterable of the form list(list) or similar structure with tuples and writes the the second value in the
    nested list to a file based on the first value in the nested list.

    Parameters
    ----------
    prefix: filename

    Returns
    -------
    None
    """
    clusters = defaultdict(list)
    for i in data:
        clusters[i[0]].append(i[1])
    for key, values in clusters.items():
        with open(f"{prefix}_cluster_{key}.{partition_index}.txt", 'w') as f:

            values = ['\t'.join(map(str, value)) for value in values]
            f.writelines(map(lambda x: f"{x}\n", values))


def store_np_matrix(matrix, path):
    np.savetxt(path, matrix, delimiter='/t')


def backup_dir(directory):
    """Backup a directory before writing to it.

    Parameters
    ----------
    directory: str

    Returns
    -------
    backup_dir: str OR None
    """
    if type(directory) is str:
        if directory.endswith('/'):
            directory = directory[:-1]
    target_path = Path(directory)
    if target_path.exists():
        t = time.localtime(target_path.stat().st_ctime)
        new_name = f"{directory}_backup_{t[0]}{t[1]}{t[2]}_{t[3]}{t[4]}{t[5]}"
        target_path.rename(new_name)
        return new_name
    return None


def combine_files(prefix):
    """Takes a prefix and finds all files of the form {prefix}[*1].[*2].txt and combines all the files with [*1] into
    one .txt file with name {prefix}[*1].txt. update 2-03-2020: now also enumerates compounds.

    Parameters
    ----------
    prefix: str

    Returns
    -------
    None
    """
    files = glob.glob(f"{prefix}*_*.txt")
    cluster_prefixes = set(['.'.join(x.split('.')[:-2]) for x in files])
    prefix_dict = {}

    for p in cluster_prefixes:
        prefix_dict[p] = glob.glob(f"{p}.*.txt")

    for k, v in prefix_dict.items():
        bag = db.read_text(v)
        with gzip.open(f'{k}.txt.gz', 'wt') as f:
            lines = [line.strip('\n') for line in bag.compute()]
            content = ''.join([f'{b}\t{a}\n' for a, b in enumerate(lines)])
            f.write(content)
        [os.remove(x) for x in v]


def zfill_enum(iterable, padding):
    """ enumerates an iterable with a number with a certain amount of padding
    Example:
    zfill_enum((a, b), 4) -> ([a, 0000], [b, 0001])

    Parameters
    ----------
    iterable: iterable(objects)
    padding: int

    Returns
    -------
    result: list(list)
    """
    result = []
    for i, item in enumerate(iterable):
        result.append([item, str(i).zfill(padding)])
    return result


def enum(iterable):
    """ enumerates an iterable with a number
    Example:
    zfill_enum((a, b), 4) -> ([a, 0], [b, 1])

    Parameters
    ----------
    iterable: iterable(objects)
    padding: int

    Returns
    -------
    result: list(list)
    """
    result = []
    for i, item in enumerate(iterable):
        result.append([item, i])
    return result


def closest_node(node, nodes):
    ''' Computes closest point to node in nodes array

    Parameters
    ----------
    node: numpy.array
    nodes: numpy.array

    Returns
    -------
    int
    '''
    deltas = nodes - node
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)


def match_node(node, nodes):
    ''' Computes closest point to node in nodes array
    #TODO make this actually different form closest_node.
    Parameters
    ----------
    node: numpy.array
    nodes: numpy.array

    Returns
    -------
    int
    '''

    deltas = nodes - node
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)


def flatten_tuple(input):
    """Takes a tuple as input and outputs a flattened tuple if needed.
    Example:
    flatten_tuple((a, b), c) -> (a, b, c)

    Parameters
    ----------
    input: tuple

    Returns
    -------
    flatten: tuple
    """
    if type(input) != tuple:
        return (input,)
    elif len(input) == 0:
        return ()
    else:
        return flatten_tuple(input[0]) + flatten_tuple(input[1:])

