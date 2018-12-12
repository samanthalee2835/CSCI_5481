# Names: Sam Lee, Noah Kindt, Ben Novacek
# hmmlearn library is open source and available at https://github.com/hmmlearn/hmmlearn
# As written, the program will give a RuntimeWarning of division by zero, not sure where the issue is, but there is output

from __future__ import division  # The hmm library needs this imported
from hmmlearn import hmm  # Makes the Hidden Markov Model possible
import argparse  # So the tsv file(s) can be read in
import pandas as pd  # This code makes use of pandas dataframes
import random  # To gather a random sampling of heights from the tsv
import numpy as num  # Arrays are used throughout the code

def make_arg_parser():  # This iteration only expects one tsv file, with many heights
  parser = argparse.ArgumentParser()
  parser.add_argument("-s", "--seqs",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to sequence tsv [required]")

  return parser

def transition_prob(df, df_len, arr, row_idx, current_base):  # Calculates probability of bases transitioning from data
                                                              # Takes dataframe, its length, row index of transition, and current base
  current_to_a_count = 0  # Current base is generic, to allow for generality
  current_to_c_count = 0
  current_to_t_count = 0
  current_to_g_count = 0

  for j in range(df_len):
    a_area = df['A_area'][j]  # Pull the data from the appropriate columns in the dataframe and over every row
    c_area = df['C_area'][j]  # In this project, height == area
    g_area = df['G_area'][j]
    t_area = df['T_area'][j]
    max_height = max(a_area, c_area, g_area, t_area)  # Find the max height out of all heights at the current index
    max_wo_current_a = max(c_area, g_area, t_area)  # Max heights discounting the current base area/height
    max_wo_current_c = max(a_area, g_area, t_area)
    max_wo_current_g = max(a_area, c_area, t_area)
    max_wo_current_t = max(a_area, c_area, g_area)
    if (max_height == a_area):  
      if (current_base == 'A'):  # If the current base matches the max height, some extra work is needed
        if (max_wo_current_a == 0):  # If the only non zero height is the current base, current->current
          current_to_a_count += 1
        elif (max_wo_current_a == c_area):  # For any height that isn't the current base, the highest is where the base might transition to
          current_to_c_count += 1
        elif (max_wo_current_a == g_area):
          current_to_g_count += 1
        elif (max_wo_current_a == t_area):
          current_to_t_count += 1
      else:                     # If the current base isn't the max height, then definitely current->max height base.
          current_to_a_count += 1
    elif (max_height == c_area):  # Same structure as A
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
    elif (max_height == g_area):  # Same structure as A
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
    elif (max_height == t_area):  # Same structure as A
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

  arr[row_idx][0] = 1.0 * current_to_a_count / df_len  # Put probabilities in their proper place
  arr[row_idx][1] = 1.0 * current_to_c_count / df_len
  arr[row_idx][2] = 1.0 * current_to_g_count / df_len
  arr[row_idx][3] = 1.0 * current_to_t_count / df_len

  return arr

def ems_prob(obs_arr, ele_idx, arr):  # Takes the observation array, which base to use, and the emissions array to build

  ele_dict = {}  # The probabilities are created using a dictionary
  total = 0.0

  for i in range(200):  # Over the entire observation array
    tup = obs_arr[0][i]  # Pull the tuple from the array
    element = tup[ele_idx]  # Get the correct base height of the tuple
    if element in ele_dict:  # If the height is already in the dictionary, pull the count and increase it
      count = ele_dict.get(element, "none")
      count += 1
      d1 = {element: count}
      ele_dict.update(d1)
    else:  # If this is a new height, add it to the dictionary
      d1 = {element: 1}
      ele_dict.update(d1)

  for j in range(200):  # Second loop, now that the dictionary is finished.
    tup = obs_arr[0][j]  # Pull the tuple
    element = tup[ele_idx]  # Get the correct base height of the tuple
    value = ele_dict.get(element, "none")  # Get the count from the dictionary
    calc = (1.0 / (len(ele_dict) * value)) # There's a 1/(dictionary length) probability of a height appearing, with a 1/(specific height count) of that particular entry
    arr[ele_idx][j] = calc  # Store the probability accordingly
    total = total + calc  # Used to track whether probabilities sum to 1 in a row

  return arr

if __name__ == '__main__':
  parser = make_arg_parser()  # Need those arguments
  args = parser.parse_args()
  norm_height = 500  # Sets height the result will come out as.

  seq_df = pd.read_csv(args.seqs,sep="\t")  # Reads the tab delimited file into a pandas dataframe
  states = num.array(([(norm_height,0,0,0)],[(0,norm_height,0,0)],[(0,0,norm_height,0)],[(0,0,0,norm_height)]))  # Each entry represents normalized A,C,G,T respectively
  obs = num.empty((1,200), dtype=(int,4))  # There'll be 200 observations of tuples
  transition_array = num.empty((4,4), dtype=float)  # 4 states makes the transition array a 4x4

  gb = seq_df.groupby('base_call')  # The dataframe will be grouped by base call and the indicies stored in gb
  a_df = gb.get_group('A').reset_index()  # Use gb to break up the dataframe into each base call, resetting the indices for later
  c_df = gb.get_group('C').reset_index()
  g_df = gb.get_group('G').reset_index()
  t_df = gb.get_group('T').reset_index()

  a_len = len(a_df.index)  # Need the lengths of each new dataframe for bounds on randint
  c_len = len(c_df.index)
  g_len = len(g_df.index)
  t_len = len(t_df.index)

  for i in range(50):  # Pull 50 random rows from each dataframe.
    a_index = random.randint(0,a_len-1)  # Needs to be minus 1, since randint is inclusive.
    c_index = random.randint(0,c_len-1)
    g_index = random.randint(0,g_len-1)
    t_index = random.randint(0,t_len-1)
    obs[0][i] = (a_df['A_area'][a_index],a_df['C_area'][a_index],a_df['G_area'][a_index],a_df['T_area'][a_index])  # Store all A tuples in the first 50 indices
    obs[0][50+i] = (c_df['A_area'][c_index],c_df['C_area'][c_index],c_df['G_area'][c_index],c_df['T_area'][c_index])  # C's go in 50-99
    obs[0][100+i] = (g_df['A_area'][g_index],g_df['C_area'][g_index],g_df['G_area'][g_index],g_df['T_area'][g_index])  # G's go in 100-149
    obs[0][150+i] = (t_df['A_area'][t_index],t_df['C_area'][t_index],t_df['G_area'][t_index],t_df['T_area'][t_index])  # T's go in 150-199

  transition_array = transition_prob(a_df, a_len, transition_array, 0, 'A')  # Get the transition probabilities for base A
  transition_array = transition_prob(c_df, c_len, transition_array, 1, 'C')  # Get the transition probabilities for base C
  transition_array = transition_prob(g_df, g_len, transition_array, 2, 'G')  # Get the transition probabilities for base G
  transition_array = transition_prob(t_df, t_len, transition_array, 3, 'T')  # Get the transition probabilities for base T

  ems = num.empty((4,200), dtype = float)  # 4 states, 200 observations, 4x200 array
  ems = ems_prob(obs, 0, ems)  # Calculate the emission probabilities for A
  ems = ems_prob(obs, 1, ems)  # Calculate the emission probabilities for C
  ems = ems_prob(obs, 2, ems)  # Calculate the emission probabilities for G
  ems = ems_prob(obs, 3, ems)  # Calculate the emission probabilities for T
 
  model = hmm.MultinomialHMM(n_components = len(states), init_params="")  # Creates the model
  model.startprob_ = num.array([0.25,0.25,0.25,0.25])  # starting probability should be equal for A,C,G, or T
  model.transmat_ = transition_array  # Populate the model with the transition and emission probabilities
  model.emissionprob_ = ems

  lst = []  # Empty list for sequence of indicies
  for j in range(200):  # Every entry in observation needs to occur at least once in hmmlearn, so one of each index gets stored in the array
    lst.append(j)

  random.shuffle(lst)  # This step is potentially not necessary, but useful for testing
  seq = num.array([lst]).T  # The model is expecting a column array, so turn it into one
  model = model.fit(seq)  # Add the sequence to the model
  logprob, output = model.decode(seq, algorithm= "viterbi")  # Traverse the HMM, storing the output and probabilities
  print ",".join(map(lambda x: str(obs[0][int(x)]), seq.T[0]))  # Print out the sequence that was fed to the model
  print ",".join(map(lambda x: str(states[int(x)]), output))  # Print out the resulting states that fit the sequence
