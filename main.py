# Commented out IPython magic to ensure Python compatibility.
# directory


# Import libraries

import timeit
import utils
import argparse

## -- INPUT PARAMETERS -- ##

parser = argparse.ArgumentParser()
parser.add_argument('first_file', help='directory of the file with the natural genome')
parser.add_argument('second_file', help='directory of the file with the variant genome')
parser.add_argument("-k","--kmerlength", dest = 'k', type=int, default = 50, help="choose the length of the kmer")
args = parser.parse_args()


k = args.k

## ---  READING THE FILE -- ##
first_file = args.first_file
second_file = args.second_file
files = [first_file, second_file]

print(f"Start the computation:\n\tk-mer length:{k}\n\tReal file: {first_file}\n\tVariant file: {second_file}")

#first_file = "real_files/salmonella-enterica.reads.fna"
#second_file ='real_files/salmonella-enterica-variant.reads.fna'
#files = [first_file, second_file]
#k = 50

real_sequences, real_filtered = {}, {}
variants_sequences, variant_filtered = {}, {}
real_file = True
for file in files:
  print(f"Processing file: {file}")  
  dictionary = utils.create_dictionary(file, k)
  utils.plot_frequency_occurrences(dictionary, k)

  
  correct = False
  while correct == False:
    print('Enter the value before which you wish to consider a sequencing error (suggested \'seq_err = 10\'):')
    seq_err = input()
    try:
      seq_err = int(seq_err)
      correct = True
    except ValueError:
      print(f"You have to insert a integer number. Input not valid: {seq_err}")
  if real_file == True:
    real_sequences, real_filtered = utils.filter_errors(dictionary, seq_err)
    real_file = False
  else:
    variants_sequences, variant_filtered = utils.filter_errors(dictionary, seq_err)

  print("-----------------------------------------------------------------------------------")

start_time = timeit.default_timer()

for keyST in real_sequences:
    if keyST in variants_sequences:
        del real_filtered[keyST]
        del variant_filtered[keyST]
   
            
stop_time = timeit.default_timer()
working_time = stop_time - start_time
print(f"It takes:{(working_time/60):.2f} minutes to do the filtering between the dictionaries ")

print(f'\nREAL:\nInit dim: {len(real_sequences)}')
print(f'Final dim: {len(real_filtered)}')
print(f'\nVARIANT:\nInit dim: {len(variants_sequences)}')
print(f'Final dim: {len(variant_filtered)}')

real_strings = utils.concatenate_kmer(real_filtered)
variant_strings = utils.concatenate_kmer(variant_filtered)

# calculate minimum distances
min_distances = utils.min_distances(real_strings, variant_strings)

utils.print_colored_SNP(min_distances)




