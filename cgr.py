# # (gH)   -_-  cgr.py  ;  TimeStamp (unix) : 17 Juillet 2015 vers 12:22

'''Generate a chaos game representation (CGR) as an image from a
fasta file located at PATH with visualisation
Usage: python cgr.py PATH or python cgr.py PATH &  # gh
'''

# Original author: Author: Bostjan Cigan
# http://zerocool.is-a-geek.net24
# Edited by Josuah Demangeon to suit the need of the University of Angers

import sys
import argparse
import collections
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
import math
##from Bio import SeqIO #### PAS BESOIN DU MAIN DANS CE PROJET
import operator

def main():
    parser = argparse.ArgumentParser(description =
'''Generate a chaos game representation (CGR) as an image from a fasta
file located at PATH.  Usage: python cgr.py PATH or python cgr.py PATH &''') # gh
    parser.add_argument('-f', '--file', required=True, nargs=1,
                        metavar='PATH_TO_FASTA_FILE',
                        help='Path to the FASTA formatted file')
    parser.add_argument('-k', '--k-mer', required=True, nargs=1, type=int,
                        metavar='K-MER_LENGTH',
                        help='Length of the words, size of the CGR')
    parser.add_argument('-o', '--output', nargs=1, metavar='PATH_TO_OUTPUT',
                        help='Path where the images will be saved')
    parser.add_argument('-t', '--table-length', nargs=1, type=int,
                        metavar='TABLE_LENGTH', help='''Number of
                        probabilities to send to standard output for
                        the html table''')
    parser.add_argument('-a', '--adjustment', nargs=1, type=float,
                        metavar='ADJUSTMENT_VALUE', help='''Float nuber between
                        0 and 10 to increase visibility on high k-mer sizes''')

    # Set the parameters from arguments
    args = parser.parse_args()
    k = args.k_mer[0]
    if args.adjustment:
        a = args.adjustment[0]
    else:
        a = 0

    # The CGR grid width and height
    array_size = int(math.sqrt(4**k))

    # Open the file and parse the fasta sequence
    with open(args.file[0]) as f:
        sequences = SeqIO.parse(f, 'fasta')

        # For each sequences from the fasta file
        for sequence in sequences:
            # Generate the CGR
            kmer_count = count_kmers(str(sequence.seq), k) # gh
            prob = probabilities(kmer_count, k, str(sequence.seq), a) # gh
            chaos = chaos_game_representation(prob, k)

            # Draw the figure with matplotlib

            plt.title('Chaos game representation of ' + sequence.id
                        + ' for ' + str(k) + '-mers' + "\n" ) # gh

            # Assign colors to the values (change cmap for different colors)
            im = plt.imshow(chaos,
                         # To set the plot start at 0
                         extent=[0, array_size, 0, array_size],
                         interpolation='nearest', cmap=cm.hot_r)
            # From http://stackoverflow.com/questions/11776663
            plt.colorbar(im, orientation='vertical')

            # Add black (k) lines to separate the main four regions
            plt.axhline(y=array_size / 2, color='k')
            plt.axvline(x=array_size / 2, color='k')

            # Remove the ticks (unneeded)
            plt.tick_params(
                axis='both', which='both',
                top='off', bottom='off', left='off', right='off',
                labelbottom='off', labelleft='off')

            # Change the ticks to A, T, G and C
            nucleotides = {'A': [0, 1], 'T': [1, 1],
                           'C': [0, 0], 'G': [1, 0]}
            for nuc, i in nucleotides.iteritems():
                plt.annotate(nuc, xy=(i[0], i[1]),
                             xycoords='axes fraction',
                             fontsize=16,
                             xytext=(i[0] * 30 - 20, i[1] * 30 - 23),
                             textcoords='offset points')

            # If output ("-o") is set,
            if args.output:
                # create an image in the folder
                output_path = (args.output[0] + '/'
                               + sequence.id + '_' + str(k) + '-mer'
                               + '_adj=' + str(a) + '.png')
                plt.savefig(output_path)
                # Then print the html code to display the image to std output
                print ('<img class="cgr" src="' + output_path + '"/>')
            else:
                # Else, display the figure (usefull for testing)
                plt.show()

            # If table-length ("-t") is set
            if args.table_length:
                # return an html table with the kmers and their frequency
                print ('<table class="cgr">')
                print ('  <th><td>Word</td><td>Frequency</td></th>')
                # Sort the probabilities
                sorted_prob = sorted(prob.items(), key=operator.itemgetter(1),
                                     reverse=True)

                for i in range(args.table_length[0]):
                    # For each probability of kmer
                    print ('  <tr>'
                           + '<td>' + sorted_prob[i][0] + '</td>'
                           + '<td>' + str(sorted_prob[i][1]) + '</td>'
                           + '</tr>')

                print ('</table>')


def count_kmers(sequence, k):
    '''Build a dictionnary (d) with the words (sequence of A, T, C and G)
    with the word as keys and the number of occurence as values.  The
    length of the words are given by the parameter k
    '''
    d = collections.defaultdict(int)
    for i in range(len(sequence) - (k - 1)): # for each possible word
        d[sequence[i:i + k]] += 1 # increment the number of occurence in d
    for key in d.keys():
        # remove the entries that contain "N" in their key that stand for
        # unknown residue: not represented in CGR
        if "N" in key:
            del d[key]
    return d


def probabilities(kmer_count, k, sequence, a):
    '''From a dictionnary of kmer (keys) keys and number of occurences
    (value) build a similar "probabilities" dictionnary with the
    probabilities instead of the number of occurence.  The
    probabilities are calculated from the form: P(pi) = N_pi / (N - k
    + 1) Where "N" is the length of the sequence, "k" is the length
    of the k-mer and "N_pi" is the number of occurrences of the k-mer
    "pi".
    '''
    probabilities = collections.defaultdict(float)
    for key, value in kmer_count.items():
        probability = float(value) / (len(sequence) - k + 0.5)
        probabilities[key] = ((10 - a) + a/k) ** math.log(probability, 10)
    return probabilities


def chaos_game_representation(probabilities, k):
    '''Convert the list of probabilities into an matrice of numeric values
    for which each quart of the matrice correspond to one nucleotid
    (A, T, C or G), themself subdivided in four regions with each also
    correspoonding to a nucleotidem and this as much time as the
    number of letters in each word (k).  An analogy is made between
    the sequence of the letters and rthe position in the
    matrice. Which square, sub-square, sub-sub-square give the first,
    second and third nucleotide of the sequence.
    '''
    # The width (and height if the array) is given by : sqrt(4 ^ k) 4
    # because there are 4 nucleotidsm A, T, C and G. k is the size of a
    # word because sqrt because the vblues are distributed in 2
    # dimension, and then require this much less space.
    array_size = int(math.sqrt(4**k))
    chaos = []
    for i in range(array_size):
        chaos.append([0]*array_size)

    for key, value in probabilities.items(): # for each word
        maxx = array_size
        maxy = array_size
        posx = 1 # reset the initial position for each word to 1
        posy = 1
        for char in key: # and each letter,

            # Iterate through the following rules to determin its
            # position in the matrice. The nucleotids have the following
            # coordinates:
            # A(0, 0) T(1, 0) C(0, 1) G(1, 1)
            # The rule used to distribute the words in the matrice are
            # those of the chaos game
            if char == "T":
                posx += maxx / 2
            elif char == "C":
                posy += maxy / 2
            elif char == "G":
                posx += maxx / 2
                posy += maxy / 2
            maxx /= 2
            maxy /= 2
        # Then, place the whole word in the matrice at its right calculted place
        chaos[int(posy - 1)][int(posx - 1)] = value

    return chaos


if __name__ == '__main__':
    main()#sys.argv[1])
