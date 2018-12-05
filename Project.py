# Names: Sam Lee, Noah Kindt, Ben Novacek

import argparse
import pandas as pd

def make_arg_parser():
  parser = argparse.ArgumentParser()
  parser.add_argument("-s", "--seqs",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to sequence csv [required]")

  return parser

if __name__ == '__main__':
  parser = make_arg_parser()
  args = parser.parse_args()

  seq_file = pd.read_csv(args.seqs)

  print(seq_file)
