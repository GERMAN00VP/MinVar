import argparse
from .main import run_minvar
import glob
from multiprocessing import Pool


def run_folder(args):
    folder, min_freq, gff = args
    return run_minvar(folder, min_freq, gff)


def main():

    parser = argparse.ArgumentParser(prog="MinVar")

    parser.add_argument("--input", required=True, help="Path or pattern to IRMA folders")
    parser.add_argument("--min_freq", type=float, default=0.1)
    parser.add_argument("--gff", default=None)
    parser.add_argument("--threads", type=int, default=1)

    args = parser.parse_args()

    folders = glob.glob(args.input)

    params = [(f, args.min_freq, args.gff) for f in folders]

    with Pool(args.threads) as p:
        p.map(run_folder, params)


if __name__ == "__main__":
    main()
