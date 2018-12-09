# Names: Sam Lee, Noah Kindt, Ben Novacek

import argparse
import pandas as pd
import random
import numpy as num

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
  states = num.array(["A", "C", "G", "T"])
  obs = num.empty((4,50), dtype=int)

  for i in range(4):
    for j in range(50):
      obs[i][j] = random.randint(0,2000)

  print(seq_file)
  print(obs)
