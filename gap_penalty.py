# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 23:15:44 2016

@author: Aishwarya

Title: Implement different gap penalty schemes
"""

#Check validity of input and preprocess
def validity( s ):
    seq = []
    for i in s:
        if i in ['G', 'T', 'C', 'A', '-']:
            seq.append( i )
        else:
            print 'Invalid base ignored.'
    return seq
            
#Remove space from input
def removeSpace( s ):
    a = []
    for i in s:
        if i != ' ':
            a.append(i)
    return a

def alignment( seq1, seq2 ):
    s = []
    i = 0
    j= 0
    while i < len(seq1) and j < len(seq2):
        if seq1[i] == seq2[j]:
            s.append(seq2[j])
            i += 1
            j += 1
        else:
            s.append('-')
            i += 1
    if j < len(seq2):
        while j < len(seq2):
            s.append(seq2[j])
            j += 1
    return s

def constant_penalty( seq2, penalty, match ):
    gap_count = 0
    match_count = 0
    i = 0
    while i < len( seq2 ):
        if seq2[i] == '-':
            gap_count += 1
            while i < len( seq2 ) and seq2[i] == '-':
                i += 1
            i = i - 1
        else:
            match_count += 1
        i += 1
    #print 'constant: ', match_count, gap_count
    score = (match * match_count) + (penalty * gap_count)
    return score

def linear_penalty(  seq2, penalty, match):
    gap_count = 0
    match_count = 0
    i = 0
    while i < len(seq2):
        if(seq2[i] == '-'):
            gap_count += 1
        else:
            match_count += 1
        i += 1
    #print 'linear: ', match_count, gap_count
    score = (match * match_count) + (penalty * gap_count)
    return score
    
def affine_penalty( seq2, penalty, match, open_penalty ):
    gap_count = 0
    open_count = 0
    match_count = 0
    i = 0
    
    #Iterate through the entire sequence
    while i < len(seq2):
        #When a new gap starts: 
        if seq2[i] == '-':
            gap_count += 1
            open_count += 1
            i = i + 1
            #Increment gap count to find length of gap
            while i < len( seq2 ) and seq2[i] == '-':
                gap_count += 1
                i += 1
            i -= 1
        else:
            match_count += 1
        i += 1
    #print 'affine: ', match_count, gap_count
    score = (match * match_count) + (open_penalty * open_count) + (penalty * gap_count)
    return score
    
def main():
    #seq1 = raw_input( 'Enter the first sequence' ).upper()
    seq2 = raw_input( 'Enter the  sequence: ' ).upper()
    seq = validity( seq2 )
    seq = removeSpace( seq )
    print 'Sequence: ', seq
    #seq = raw_input( 'Enter the aligned sequence: ' ).upper()
    penalty = input( 'Enter the gap penalty: ' )
    match = input( 'Enter the match score: ' )
    #seq = alignment( seq1, seq2 )
    open_penalty = input( 'Enter the opening penalty for affine scheme: ' )
    
    print 'Constant gap penalty = ', constant_penalty( seq, penalty, match )
    print 'Linear gap penalty = ', linear_penalty( seq, penalty, match )
    print 'Affine gap penalty = ', affine_penalty( seq, penalty, match, open_penalty )
    
main()