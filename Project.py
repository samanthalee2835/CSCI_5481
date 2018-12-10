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
                      help="Path to sequence tsv [required]")

  return parser

if __name__ == '__main__':
  parser = make_arg_parser()
  args = parser.parse_args()
  norm_height = 500

  seq_df = pd.read_csv(args.seqs,sep="\t")
  states = num.array(([(norm_height,0,0,0)],[(0,norm_height,0,0)],[(0,0,norm_height,0)],[(0,0,0,norm_height)]))
  obs = num.empty((50,4), dtype=(int,4))

  gb = seq_df.groupby('base_call')
  a_df = gb.get_group('A').reset_index()
  c_df = gb.get_group('C').reset_index()
  g_df = gb.get_group('G').reset_index()
  t_df = gb.get_group('T').reset_index()

  a_len = len(a_df.index)
  c_len = len(c_df.index)
  g_len = len(g_df.index)
  t_len = len(t_df.index)
  
  for i in range(50):
    a_index = random.randint(0,a_len-1)
    c_index = random.randint(0,c_len-1)
    g_index = random.randint(0,g_len-1)
    t_index = random.randint(0,t_len-1)
    obs[i][0] = (a_df['A_area'][a_index],a_df['C_area'][a_index],a_df['G_area'][a_index],a_df['T_area'][a_index])
    obs[i][1] = (c_df['A_area'][c_index],c_df['C_area'][c_index],c_df['G_area'][c_index],c_df['T_area'][c_index])
    obs[i][2] = (g_df['A_area'][g_index],g_df['C_area'][g_index],g_df['G_area'][g_index],g_df['T_area'][g_index])
    obs[i][3] = (t_df['A_area'][t_index],t_df['C_area'][t_index],t_df['G_area'][t_index],t_df['T_area'][t_index])


  transitions = num.empty((4,4), dtype=float)

  print(a_df)
#  print(obs)
