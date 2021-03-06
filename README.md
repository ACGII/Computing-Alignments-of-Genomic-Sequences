# Computing-Alignments-of-Sequences
## Coursera Fundamentals of Computer Science
### Algorithmic Thinking (Part 2)
### Rice University

Implementation of four functions,   these functions in Application 4 are used to analyze two problems involving comparison of similar sequences: 

1.  The first pair of functions will return matrices that we will use in computing the alignment of two sequences. The two types of    matrices: alignment matrices and scoring matrices. 

     A.  Alignment matrices will follow the same indexing scheme that we used for grids in "Principles of Computing". Entries in the alignment matrix will be indexed by their row and column with these integer indices starting at zero. We will model these matrices as lists of lists in Python and can access a particular entry via an expression of the form alignment_matrix[row][col]|}alignment_matrix[row][col].
     
     a.  compute_alignment_matrix(seq_x,seq_y,scoring_matrix,global_flag): Takes as input two sequences seq_y whose elements share a
common alphabet with the scoring matrix scoring_matrix. The function computes and returns the alignment matrix for seq_y as described in the Homework. If global_flag is True, each entry of the alignment matrix is computed using the method described in Question 8 of the Homework. If global_flag is False, each entry is computed using the method described in Question 12 of the Homework.
          
     B.  For scoring matrices, we take a different approach since the rows and the columns of the matrix are indexed by characters. In particular, we will represent a scoring matrix in Python as a dictionary of dictionaries. Given two characters row_char and col_char, we can access the matrix entry corresponding to this pair of characters via scoring_matrix[row_char][col_char]
     
     a.  build_scoring_matrix(alphabet,diag_score,off_diag_score,dash_score): Takes as input a set of characters alphabet and three scores diag_score, off_diag_score, and dash_score. The function returns a dictionary of dictionaries whose entries are indexed by pairs of characters in alphabet plus '-'. The score for any entry indexed by one or more dashes is dash_score. The score for the remaining diagonal entries is diag_score. Finally, the score for the remaining off-diagonal entries is off_diag_score.
            
2.  The second pair of functions will return global and local alignments of two input sequences based on a provided alignment matrix.  For the this part of Project 4, you will use the alignment matrix returned by compute_alignment_matrix to compute global and local alignments of two sequences seq_x and seq_y. The first function will implement the method ComputeAlignment discussed in Question 9 of the Homework.

     A.  compute_global_alignment(seq_x,seq_y,scoring_matrix,alignment_matrix): Takes as input two sequences seq_x and seq_y whose elements share a common alphabet with the scoring matrix scoring_matrix. This function computes a global alignment of seq_x and seq_y using the global alignment matrix alignment_matrix.The function returns a tuple of the form (score, align_x, align_y) where score is the score of the global alignment align_x and align_y. Note that align_x and align_y should have the same length and may include the padding character  ’-’.
     
     B.  The function, optimal local alignment, starting at the maximum entry of the local alignment matrix and working backwards to zero as described in Question 13 of the Homework.  compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix): Takes as input two sequences seq_x and seq_y whose elements share a common alphabet with the scoring matrix scoring_matrix. This function computes a local alignment of seq_x and seq_y using the local alignment matrix alignment_matrix.The function returns a tuple of the form (score, align_x,align_y) where score is the score of the optimal local alignment align_x and align_y. Note that align_x and align_y should have the same length and may include the padding character ’-’.
