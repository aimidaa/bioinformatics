# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 00:22:55 2016

@author: Aishwarya

Title: ClustalW for multiple sequence alignment

"""

###################################################

import dendogram

####################################################

''' We need to align DNA sequences. 
    So invalid charachters must be eliminated.
'''
####################################################

#Check the validity of input
def validity(s1):
    a = ''
    for i in s1:
        if i in ['A', 'C', 'T', 'G']:
            a = a + i
    if len(a) < len(s1):
        print 'Invalid characters found and ignored in Seq. '
    return a

#Remove space from input
def removeSpace(s1):
    a = ''
    for i in s1:
        if i != ' ':
            a = a + i    
    return a
    
####################################################

''' Implement Needleman-Wunsch Algorithm for global alignment.
    Return the pairwise aligned sequence.
'''
####################################################

#Initialize a 2D matrix with m rows and n columns
def initialize(m, n):
    M  = []
    for i in range(m):
        M.append([])
        for j in range(n):
            M[i].append(0)
    return M

# Find scoring function
def S(a, b, match, mismatch):
    if a == b:
        return match
    else:
        return mismatch
    
# Function to implement Needleman-Wunsch    
def globalAlign(s1, s2):
    match = 1
    mismatch = -2
    gap = -1
    n = len(s1)
    m = len(s2)
    a = 0
    b = 0
    c = 0
    
    #Initialize V and D matrix
    V = initialize(m, n) #Stores the score values
    D = initialize(m, n) #Stores the direction of arrows

    #Initialize first row and column
    for i in range(1, m):
        V[i][0] = V[i-1][0] + gap
        D[i][0] = 'u'
    for i in range(1, n):
        V[0][i] = V[0][i-1] + gap
        D[0][i] = 'l'          

    for i in range(1, m):
        for j in range(1, n):
            #print '', i, j, s1[j], s2[i], S(s1[j], s2[i], match, mismatch)
            a = V[i-1][j-1] + S(s1[j], s2[i], match, mismatch)
            b = V[i-1][j] + gap
            c = V[i][j-1] + gap
            V[i][j] = max(a, b, c)
            #Save the direction
            if V[i][j] == a:
                D[i][j] = 'd'
            elif V[i][j] == b:
                D[i][j] = 'u'
            else:
                D[i][j] = 'l'
    
    #Find optimal alignment
    r1 = []
    r2 = []    
    #Start from bottom-most position
    i = m - 1
    j = n - 1
    
    while i > 0 or j > 0:
        #Check the direction of arrow and find the sequence
        if D[i][j] == 'd':
            r1.append(s1[j])
            r2.append(s2[i])
            i = i - 1
            j = j - 1
            
        elif D[i][j] == 'l':
            r1.append(s1[j])
            r2.append('-')
            j = j - 1
            
        elif D[i][j] == 'u':
            r1.append('-')
            r2.append(s2[i])
            i = i - 1
             
    r1 = r1[::-1]
    r2 = r2[::-1]
    return r1, r2

#####################################################

''' Guide is constructed using Neighbor Joining Method.
    The order in which sequence must be aligned is returned.
'''
####################################################

# Function to find the ordering of sequences
def guide_tree( matrix, n, final, li ):
    
    # Stop when length of matrix is 1
    if len(li) == 1:
        return final
    max = 0
    ind_i = 0
    ind_j = 0
    
    for i in li:
        for j in li:
            if (i, j) in matrix.keys() and matrix[ (i, j) ] > max:
                max = matrix[ (i, j) ]
                ind_i = i
                ind_j = j
    final.append( (ind_i, ind_j) )
    
    li.remove(ind_i)
    li.remove(ind_j)
    
    del matrix[(ind_i, ind_j)] 
    # Update the matrix
    for i in range( len(li) ):
        matrix[ ((ind_i, ind_j), li[i]) ] = matrix[ (li[i], ind_i) ] + matrix[ (li[i], ind_j) ] / float(2)
        matrix[ (li[i], (ind_i, ind_j)) ] = matrix[ (li[i], ind_i) ] + matrix[ (li[i], ind_j) ] / float(2)
   
    li.append((ind_i, ind_j))
    
    return guide_tree( matrix, n, final, li)

####################################################

''' Align multiple sequence with the help of profile matrix.
        
'''
####################################################

def profile_matrix( all_s ):
    
    n = len( all_s )
    m = len( all_s[0] )
    profile = initialize(5, len(all_s[0]))
    
    for seq in all_s:        
        for col in range( m ):           
            if seq[col] == 'A':
                profile[0][col] = profile[0][col] + float(1) / n
            if seq[col] == 'C':
                profile[1][col] = profile[1][col] + float(1) / n
            if seq[col] == 'T':
                profile[2][col] = profile[2][col] + float(1) / n
            if seq[col] == 'G':
                profile[3][col] = profile[3][col] + float(1) / n
            if seq[col] == '-':
                profile[4][col] = profile[4][col] + float(1) / n
                
    return profile, m

def new_sequence( profile, n ):
    
    ind_i = 0
    rows = ['A', 'C', 'T', 'G', '-']
    r = []
    for j in range( n ):
        max = 0
        for i in profile:
            #print 'i', i,
            if i[j] > max:
                max = i[j]
                ind_i = profile.index( i )
        r.append( rows[ind_i] )
            
    return r

# Multiple sequence alignment.
def msa(all_s, n):
    
    profile = []
    result = []

    r1, r2 = globalAlign( all_s[0], all_s[1])
    result.append(r1)
    result.append(r2)
    profile, m =  profile_matrix( result )
    r = new_sequence( profile, m )
    
    for i  in range(2, n):
        r.insert(0, 'o')   
        r1, r2 = globalAlign( r, all_s[i])
        li = []
        
        for j in range(len( r1 )):
            x = 0
            if j < len( r ):
                if r[j] != '-':
                    x = 1
            if r1[j] == '-' and x == 1:
                li.append(j)
        for k in range( len(result) ):
            for l in range(len(li)):
                if li[l] == '-':
                    result[k].insert( l, '-' )
                    for a in range(len(li)):
                        li[a] += 1
                    
        result.append(r2)
        profile, m =  profile_matrix( result )          
        r = new_sequence( profile, m )
    return result

####################################################

# Function to find sequence similarity to fill similarity matrix
def similarity( s1, s2 ):
    n = len( s1 )
    score = 0
    for i in range( n ):
        if s1[i] == s2[i]:
            score = score + 1
    score = float(score) / n
    return score
    
#####################################################
    
# Define the main function from where the execution begins.
def main():
    n = input( 'Enter the number of sequences: ' )
    all_s = []
    for i in range( n ):
        seq = raw_input( 'Enter the sequence: ' ).upper()
        seq = validity( seq )
        #seq = removeSpace( seq )
        all_s.append( seq )
    print 'Input Sequences are: '
    for seq in all_s:
        print seq
    aligned_s = []
    similarity_matrix = [] 
    
    for i in range( n ):
        similarity_matrix.append([])
        for j in range( n ):
            similarity_matrix[i].append(0)
    for i in range( n ):
        for j in range( n ):
            if i != j:
                s1 = 'o' + all_s[i]
                s2 = 'o' + all_s[j]
                r1, r2 = globalAlign( s1, s2 )  
                if i < j:
                    print 'Sequence ', i, j, ' : '
                    print r1
                    print r2
                aligned_s.append( r1 )
                aligned_s.append( r2 )
                score = similarity( r1, r2 )
                similarity_matrix[i][j] = round(score,2)
                
    print 'Similarity Matrix:'
    for i in range( n ):
        for j in range( n ):
            if i <= j:
                print 0.0,
            else:
                print similarity_matrix[i][j],
        print '\n'
    
    matrix = { }
    final = []
    li = []
    for i in range( n ):
        li.append(i)
        for j in range( n ):
            matrix[(i, j)] = similarity_matrix[i][j]
    tree = guide_tree( matrix, n, final, li)
    tree = tree[-1]
    tree = tree[::-1]
    print 'Guided Tree: ', tree
    order = dendogram.printDendrogram(tree, sep=3)   
    
    aligned = []
    for i in range(len(order)):
        aligned.append(all_s[int(order[i])])
        aligned[i] = 'o' + aligned[i]
        
    result = msa( aligned, n )
    
    print 'Final Alignment: '
    for r in result:
        print r
main()             
                
#####################################################