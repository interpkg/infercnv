"""
    Hua Sun
    2025-06-17 v0.33

"""

import argparse
import pandas as pd
import os
import yaml

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, default='', help='input file')
parser.add_argument('-r', '--run', default='all', help='run step')
args = parser.parse_args()


def main():
    path = os.path.dirname(__file__)
    with open(args.input, 'r') as f:
        config = yaml.safe_load(f)
    r_script = config['rscript']

    if args.run == 'all':
        print('[INFO] Extraction data ...')
        fscript = f'{path}/src/seu2mtx.v3.R'
        cmd = f'{r_script} {fscript} {args.input}'
        os.system(cmd)

    if args.run in ['all', 'cnv']:
        print('[INFO] Call CNV ...')
        fscript = f'{path}/src/infercnv.R'
        cmd = f'{r_script} {fscript} {args.input}'
        os.system(cmd)



if __name__ == '__main__':
    main()


