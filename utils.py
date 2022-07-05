# ------------------------------------------------- LIBRARIES ----------------------------------------------------- #

import Bio.SeqIO
import numpy as np
import timeit
import matplotlib.pyplot as plt
import copy
import sys
import utils
from termcolor import colored
from Levenshtein import distance as lev




# ------------------------------------------------- FUNCTIONS ----------------------------------------------------- #

def create_dictionary(file, k = 50):
  start_time = timeit.default_timer()
  dict_frequency = {}
  for record in Bio.SeqIO.parse(file, "fasta"):
      sequence = str(record.seq)
      for i in range(len(sequence) - k + 1 ):
          key = sequence[i:i+k]

          if key not in dict_frequency:
              dict_frequency[key] = 0
          dict_frequency[key] +=1
  stop_time = timeit.default_timer()
  reading_time = stop_time - start_time
  print(f"Time required for reading the file and store it in a dictionary is {(reading_time)/60:.2f} minutes")

  return dict_frequency


def plot_frequency_occurrences(dictionary, k):
  ## -- PLOT: (x,y) = (number of occurencies, frquency) of a k-mer -- ##
  frequency = {}
  freq = list(dictionary.values())
  freq.sort()
  for el in freq:
    if el in frequency.keys():
        frequency[el] +=1 
    else:
        frequency[el] = 1

  #frequency.pop(1)
  plot_dict = {}
  for i in range(1,300):
      if i in frequency.keys():
      	
      	plot_dict[i] = frequency[i]
  
  keys = list(plot_dict.keys())
  values = list(plot_dict.values())

  results = (keys, values)

  # -- PLOT: log-scale -- ##
  x,y = results
  fig, ax = plt.subplots()
  ax.plot(x, y,  label = f'k = {k}')
  ax.set_yscale('log')
  ax.legend(loc = 'upper right')
  plt.xlim(-2)
  plt.show()
  


def filter_errors(dictionary, seq_err = 10):
    print(f'Init dimension of the dictionary: {len(dictionary)}')
    start_time = timeit.default_timer()
    for key,value in dict(dictionary).items():
        if value < seq_err:
            del dictionary[key]
    filter_time = timeit.default_timer()
    filtering_time = filter_time - start_time
    print(f"It takes: {(filtering_time):.2f} seconds to filter the dictionary")
    copied_dict = copy.deepcopy(dictionary)
    copy_time = timeit.default_timer()
    copying_time = copy_time - filter_time
    print(f"It takes: {(copying_time):.2f} seconds to copy the dictionary")
    print(f'Final dimension of the dictionary: {len(dictionary)}')
    return dictionary, copied_dict


def concatenate_kmer(dictionary):
  print("\nConcatenating k-mer:")
  concatenated_strings = []
  banned = []
  

  for key in dictionary: #real_filtered:
      string = key
      counter = 1
      
      
      if key not in banned:
          flag = False
          for key2 in dictionary: #real_filtered:
              if key != key2 and (key2 not in banned) and (key not in banned):

                  
                  ### --> the tail of the 1st kmer with the head of the 2nd <-- ###
                  key_cutted = string[counter:]
                  key2_cutted = key2[:-1]
                  
                  if key_cutted == key2_cutted:
                      flag = True
                      banned.append(key2)
                      counter+=1
                      letter = key2[-1]
                      string += key2[-1]
                      
                  else:                    
                  ### --> the head of the 1st kmer with the tail of the 2nd <-- ###
                      key_cutted = string[:-counter]
                      key2_cutted = key2[1:]
                      
                      if key_cutted == key2_cutted:
                          flag = True
                          banned.append(key2)
                          counter+=1
                          letter = key2[0]
                          string = letter + string
          
          if flag == True : 
              banned.append(key)        
          concatenated_strings.append(string)
  print("Concatenating k-mer finished.")
  return concatenated_strings



"""## Find distances"""



## distances

## GET ALL THE POSSIBLE DISTANCES ##
def all_distances(real_seqs, variant_seqs):
  distances = {}
  for real in real_seqs:
      for variant in variant_seqs:

          distance = lev(real, variant)
          distances[(real,variant)] = distance
          #print(f"Distance {real} - {variant} : {distance}")
  
  return distances

def min_distances(real_seqs, variant_seqs):
  ## GET THE MINIMUM DISTANCE for each real k-mer ##
  print("Calculate minimum distances")
  min_distances = {}
  for real in real_seqs:
      closest = ''
      min_distance = sys.maxsize
      for variant in variant_seqs:

          distance = lev(real, variant)
          if distance <= min_distance:
              min_distance = distance
              closest = variant
    
      if closest != '':
      	  if min_distance <= 10: #scelto a priori --> numero max di snp per print
          	min_distances[(real,closest)] = min_distance
      else:
          print(f"No comparison for the real k-mer: {real}")
  
  return min_distances

def print_SNP(min_distances):
  for key, val in min_distances.items():
      print('SNP')
      print(f'\tW: {key[0]}')
      print(f'\tV: {key[1]}')
      print(f'\t\tDistance: {val}\n------------------------------------------------------')

def print_colored_SNP(min_distances):
  print("Printing results:")
  for key, val in min_distances.items():

      if len(key[0]) == len(key[1]):
          diff = [i for i in range(len(key[0])) if key[0][i] != key[1][i]]
          real = ''
          variant = ''
    
          for i in range(len(key[0])):
            if i in diff:
              real += colored(key[0][i], 'red')
              variant += colored(key[1][i], 'red')
            else:
              real += key[0][i]
              variant += key[1][i]
    
    
          print('SNP')
          print(f'\tW: {real}')
          print(f'\tV: {variant}')
          print(f'\t\tDistance: {val}\n------------------------------------------------------')

