# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 00:09:51 2016

@author: Aishwarya

@Title: Program to implement Smith-Watermann algorithm for local-pairwise alignment.

"""

#Check the validity of input
def validity(s1, s2):
    a = ''
    b = ''
    for i in s1:
        if i in ['A', 'C', 'T', 'G']:
            a = a + i
    if len(a) < len(s1):
        print 'Invalid characters found and ignored in Seq1. '
    for i in s2:
        if i in ['A', 'C', 'T', 'G']:
            b = b + i
    if len(b) < len(s2):
        print 'Invalid characters found and ignored in Seq2.' 
    return a, b

#Remove space from input
def removeSpace(s1, s2):
    a = ''
    b = ''
    for i in s1:
        if i != ' ':
            a = a + i
    for i in s2:
        if i != ' ':
            b = b + i
    return a, b

# Initilize the H matrix with all zeroes
def initialize( n1, n2 ):
    H = []
    for i in range( n1 ):
        H.append([])
        for j in range( n2 ):
            H[i].append(0)
    
    return H

#Function to caluculate score S(a, b)
def S( a, b ):
    if a == b:
        return 2 # Matching score is 2
    else:
        return -1 #Mismatch penalty is -1

# Populate the H matrix with values.        
def fill_values( H, s1, s2 ):
    n1 = len(s1)
    n2 = len(s2)
    w = -1 # Weight is -1
    dic= {} # This dictionary stores the direction of the location which gives the maximum value.
    for i in range(n1):
        for j in range(n2):
            dic[(i,j)] = ''
    for i in range( 1, n1 ):
        for j in range( 1, n2 ):
            a = H[i-1][j-1] + S(s1[i],s2[j])
            b = H[i-1][j] + w
            c = H[i][j-1] + w
            H[i][j] = max(a,b,c,0)
            if H[i][j] == a:
                dic[(i,j)] = dic[(i,j)] + 'd' # d:diagonal
            if H[i][j] == b:
                dic[(i,j)] = dic[(i,j)]+ 'u' # u:up
            if H[i][j] == c:
                dic[(i,j)] = dic[(i,j)]+ 'l' # l:left
    return dic, H    

# From H matrix find the      
def water_smith( s1, s2, dic, i, j, k ):
    r1 = []
    r2 = []

    #Check the direction of arrow and find the sequence
    
    while i > 0 and j > 0:
        """if k > (len(dic[(i, j)]) - 1):
            r1 = []
            r2 = []
            return r1, r2"""
        if dic[(i, j)][k] == 'd':
            r1.append(s1[i])
            r2.append(s2[j])
            i = i - 1
            j = j - 1
            #print r1, r2
        elif dic[(i, j)][k] == 'u':
            r1.append(s1[i])
            r2.append('-')
            i = i - 1
            #print r1, r2
        elif dic[(i, j)][k] == 'l':
            r1.append('-')
            r2.append(s2[j])
            j = j - 1
            #print r1, r2

    return r1, r2


# Function to swap two sequences, so that s1 is always the smallest sequence.
def swap( n1, n2, s1, s2 ):
    if n1 > n2:
        return n2, n1, s2, s1
    else:
        return n1, n2, s1, s2

def find_max(i, j, H):
    m = 0
    n = len(H)
    ind_i = 0
    ind_j = 0
    for a in range(i, n):
        for b in range(j, len(H[a])):
            if H[a][b] > m:
                m = H[a][b]                
                ind_i = a
                ind_j = b

    return ind_i, ind_j    
def main():
    s1 = raw_input( 'Enter the sequence s1: ' ).upper()
    s2 = raw_input( 'Enter the sequence s2: ' ).upper()
    s1, s2 = validity(s1, s2)
    #s1, s2 = removeSpace(s1, s2)
    print 'Sequences are:'
    print 'Seq1: ', s1
    print 'Seq2: ', s2
    print '***************************'
    s1 = 'o' + s1
    s2 = 'o' + s2
    length1 = len(s1)
    length2 = len(s2)
    n1, n2, s1, s2 = swap( length1, length2, s1, s2 )
    H = initialize( n1, n2 )
    #print H
    dic, H = fill_values( H, s1, s2 )
    print 'The Matrix:'
    for rows in H:
        print rows
    print dic
        
    i = 0
    j = 0
    #while ():
    i, j = find_max(i, j + 1, H)
    k = 0
    while k < 3:
        r1, r2 = water_smith( s1, s2, dic, i, j, k )
        n1, n2, r1, r2 = swap( length1, length2, r1, r2 )
        k = k + 1
        if len(r1) > 0:
            r1 = r1[::-1]
            r2 = r2[::-1]
            print '***************************'
            #print 'The optimal local alignment:'
            print 'Seq1: ', r1
            print 'Seq2: ', r2
            
    # The sequence obtained will be in reverse direction, since we are processing from last row.
    # Therefore reverse it.
    
    
main()