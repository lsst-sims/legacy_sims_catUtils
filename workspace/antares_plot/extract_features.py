import argparse
import os
import scipy
import sqlite3

from FeatureExtraction import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, default=None)
    parser.add_argument('--prefix', type=str, default=None)
    parser.add_argument('--out_file', type=str, default=None)

    args = parser.parse_args()
    if args.data_dir is None:
        raise RuntimeError("Must specify data_dir")
    if args.prefix is None:
        raise RuntimeError("Must specify prefix")
    if args.out_file is None:
        raise RuntimeError("Must specify out_file")
