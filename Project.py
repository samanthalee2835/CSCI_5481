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

def transition_prob(df, df_len, arr, row_idx, current_base):

  current_to_a_count = 0
  current_to_c_count = 0
  current_to_t_count = 0
  current_to_g_count = 0

  for j in range(df_len):
    a_area = df['A_area'][j]
    c_area = df['C_area'][j]
    g_area = df['G_area'][j]
    t_area = df['T_area'][j]
    max_height = max(a_area, c_area, g_area, t_area)
    max_wo_current_a = max(c_area, g_area, t_area)
    max_wo_current_c = max(a_area, g_area, t_area)
    max_wo_current_g = max(a_area, c_area, t_area)
    max_wo_current_t = max(a_area, c_area, g_area)
    if (max_height == a_area):
      if (current_base == 'A'):
        if (max_wo_current_a == 0):
          current_to_a_count += 1
        elif (max_wo_current_a == c_area):
          current_to_c_count += 1
        elif (max_wo_current_a == g_area):
          current_to_g_count += 1
        elif (max_wo_current_a == t_area):
          current_to_t_count += 1
      else:
          current_to_a_count += 1
    elif (max_height == c_area):
      if (current_base == 'C'):
        if (max_wo_current_c == 0):
          current_to_c_count += 1
        elif (max_wo_current_c == a_area):
          current_to_a_count += 1
        elif (max_wo_current_c == g_area):
          current_to_g_count += 1
        elif (max_wo_current_c == t_area):
          current_to_t_count += 1
      else:
          current_to_c_count += 1
    elif (max_height == g_area):
      if (current_base == 'G'):
        if (max_wo_current_g == 0):
          current_to_g_count += 1
        elif (max_wo_current_g == a_area):
          current_to_a_count += 1
        elif (max_wo_current_g == c_area):
          current_to_c_count += 1
        elif (max_wo_current_g == t_area):
          current_to_t_count += 1
      else:
          current_to_g_count += 1
    elif (max_height == t_area):
      if (current_base == 'T'):
        if (max_wo_current_t == 0):
          current_to_t_count += 1
        elif (max_wo_current_t == a_area):
          current_to_a_count += 1
        elif (max_wo_current_t == c_area):
          current_to_c_count += 1
        elif (max_wo_current_t == g_area):
          current_to_g_count += 1
      else:
          current_to_t_count += 1

  transition_array[row_idx][0] = 1.0 * current_to_a_count / df_len
  transition_array[row_idx][1] = 1.0 * current_to_c_count / df_len
  transition_array[row_idx][2] = 1.0 * current_to_g_count / df_len
  transition_array[row_idx][3] = 1.0 * current_to_t_count / df_len

  return transition_array

if __name__ == '__main__':
  parser = make_arg_parser()
  args = parser.parse_args()
  norm_height = 500

  seq_df = pd.read_csv(args.seqs,sep="\t")
  states = num.array(([(norm_height,0,0,0)],[(0,norm_height,0,0)],[(0,0,norm_height,0)],[(0,0,0,norm_height)]))
  obs = num.empty((50,4), dtype=(int,4))
  transition_array = num.empty((4,4), dtype=float)

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

  transition_array = transition_prob(a_df, a_len, transition_array, 0, 'A')
  transition_array = transition_prob(c_df, c_len, transition_array, 1, 'C')
  transition_array = transition_prob(g_df, g_len, transition_array, 2, 'G')
  transition_array = transition_prob(t_df, t_len, transition_array, 3, 'T')

  print(transition_array)
 
