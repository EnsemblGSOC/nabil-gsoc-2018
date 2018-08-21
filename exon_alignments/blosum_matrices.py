"""
    Obtain Blosum matrix scores for computing 
    protein similarity
"""


from Bio.SubsMat.MatrixInfo import *

                                    # maps string of blosum matrix names to blosum matrix objects
blosums = { "blosum30" : blosum30,                  
            "blosum35" : blosum35,
            "blosum40" : blosum40,
            "blosum45" : blosum45,
            "blosum50" : blosum50,
            "blosum55" : blosum55,
            "blosum60" : blosum60,
            "blosum62" : blosum62,
            "blosum65" : blosum65,
            "blosum70" : blosum70,
            "blosum75" : blosum75,
            "blosum80" : blosum80,
            "blosum85" : blosum85,
            "blosum90" : blosum90,
            "blosum95" : blosum95,
            "blosum100" : blosum100 }


def get_raw_blosum_matrix(blosum_mat):
    """
    Returns the BioPython Blosum Matrix for alignment purposes
    
    Arguments:
        blosum_mat {str} -- which blosum matrix to use
    
    Returns:
        dict -- The blosum matrix
    """


    return blosums[blosum_mat]


def get_blosum_scores(blosum_mat, gap_penalty=0.0):
    """
    Obtain the requested Blosum scoring with gap penalties
    
    Arguments:
        blosum_mat {str} -- name of the blosum matrix
    
    Keyword Arguments:
        gap_penalty {float} -- gap open penalty (default: {0.0})        
    
    Returns:
        [tuple] -- a tuple (blosum_matrix as a dictionary, max score in the matrix, min score in the matrix)
    """

    
    blosum = blosums[blosum_mat]    # obtaining the requested blosum matrix
                                    # initializing max and min scores
    blosum_scores_max = float(gap_penalty)
    blosum_scores_min = float(gap_penalty)

    out = {}                        # initializing output matrix

    for i in blosum:

        blosum_scores_max = max(blosum_scores_max, blosum[i])
        blosum_scores_min = min(blosum_scores_min, blosum[i])

        try:

            out[i[0]][i[1]] = blosum[i]

        except:             

            out[i[0]] = {}
            out[i[0]][i[1]] = blosum[i]


        try:

            out[i[1]][i[0]] = blosum[i]

        except:

            out[i[1]] = {}
            out[i[1]][i[0]] = blosum[i]

    out['-'] = {}

    for i in out:       # introducing gap open penalty

        out[i]['-'] = gap_penalty 
        out['-'][i] = gap_penalty 


                        # formating the output structure
    return (out, blosum_scores_max, blosum_scores_min)

def main():
    # example
    print(get_blosum_scores('blosum62'))

if __name__ == '__main__':
    main()