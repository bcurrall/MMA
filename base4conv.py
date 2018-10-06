### build MMA matrix

## imports
import numpy as np
import re
## definitions


# convert to base 4

def base_convert(i, b): # converts base 10 to other bases - usually base 4 here
    result = []
    if i == 0:
        result.insert(0, 0)
    else:
        while i > 0:
            result.insert(0, i % b)
            i = i // b
    return result


def mma_seq(i, x):
    result = []
    rlist = list(range(i, x))
    for seq in rlist:
        result.append(int("".join(map(str, base_convert(seq, 4)))))
    return result


def mma_matrix(refSeqStart, refSeqEnd, fqSeqStart, fqSeqEnd):
    fullList = []
    i = 0
    refSeq = mma_seq(refSeqStart, refSeqEnd)
    fastqSeq = mma_seq(fqSeqStart, fqSeqEnd)
    for r_seq in refSeq:
        partialList = []
        for f_seq in fastqSeq:
            partialList.append([f_seq, fastqSeq[i]])
        fullList.append(partialList)
        i += 1
    return fullList



print(mma_matrix(0, 4**2, 0, 4**2))

fullList = mma_matrix(0, 4**2, 0, 4**2)

'''
fullArray = np.empty((len(refSeq), len(fastqSeq)), int)
print("full array = ", fullArray)
'''

'''
fullList = []
for r_seq in refSeq:
    print(r_seq)
    for f_seq in fastqSeq:
        print("fastSeq[] = ", fastqSeq[r_seq])
        fullList.append(fastqSeq[r_seq])
        print("results", fullList)
'''

'''
def mma_matrix ():
    fullList = []
    i = 0
    for r_seq in refSeq:
        print("r_seq = ", r_seq)
        print("i = ", i)
        partialList = []
        for f_seq in fastqSeq:
            print(f_seq)
            print("fastSeq[] = ", fastqSeq[i])
            partialList.append([f_seq, fastqSeq[i]])
            print("partialList = ", partialList)
        fullList.append(partialList)
        i += 1
    print("fullList = ", fullList)
    return fullList
'''


'''
for n, i in enumerate(fullList):
    if i == 0:
        fullList[n] = 5
    if i == 1:
        fullList[n] = 6
    if i == 2:
        fullList[n] = 7
    if i == 3:
        fullList[n] = 8

print("fullList_conv = ", fullList)
'''

'''
a = [1, 2, 3, 4]

print(a)

for n, i in enumerate(fullList):
    print('n = ', n, 'i = ', i)
    print('fullList(n) = ', fullList[n])
    if i == 1:
        fullList[n] = 6

print('fullList_conv = ', fullList)
'''






# smith-waterman - simple code, will be slow; check out pyPaSWAS for speed up

def smith_waterman(a: str, b: str, alignment_score: float = 1, gap_cost: float = 1) -> float:
  """
  Compute the Smith-Waterman alignment score for two strings.
  See https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Algorithm
  This implementation has a fixed gap cost (i.e. extending a gap is considered
  free). In the terminology of the Wikipedia description, W_k = {c, c, c, ...}.
  This implementation also has a fixed alignment score, awarded if the relevant
  characters are equal.
  Kinda slow, especially for large (50+ char) inputs.
  """
  # H holds the alignment score at each point, computed incrementally
  H = np.zeros((len(a) + 1, len(b) + 1))
  for i in range(1, len(a) + 1):
    for j in range(1, len(b) + 1):
      # The score for substituting the letter a[i-1] for b[j-1]. Generally low
      # for mismatch, high for match.
      match = H[i-1,j-1] + (alignment_score if a[i-1] == b[j-1] else 0)
      # The scores for for introducing extra letters in one of the strings (or
      # by symmetry, deleting them from the other).
      delete = H[1:i,j].max() - gap_cost if i > 1 else 0
      insert = H[i,1:j].max() - gap_cost if j > 1 else 0
      H[i,j] = max(match, delete, insert, 0)
  # The highest score is the best local alignment.
  # For our purposes, we don't actually care _what_ the alignment was, just how
  # aligned the two strings were.
  return H.max()

print(smith_waterman("10132", "101312"))