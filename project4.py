""" project 4 """

def build_scoring_matrix(alpha, diag_sc, off_diag_sc, dash_sc):
    """Creates a scoring matrix as specified"""
    scm =  [[0]*(len(alpha) + 1) for _ in range(len(alpha) + 1)]
    alph = alpha.copy()
    abc =  str('-')
    while len(alph) > 0:
        abc += alph.pop()
    abc = str(abc)
    #print 'abc',abc,'alpha',alpha,len(abc),len(alpha)
    for iii in range(len(abc)):
        for jjj in range(len(abc)):
            if iii == 0 or jjj == 0:
                scm[iii][jjj] = dash_sc
            elif iii != 0 and jjj != 0:
                if iii == jjj:
                    scm[iii][jjj] = diag_sc
                else:
                    scm[iii][jjj] = off_diag_sc
    
    dima = {}
    for iii in range(len(abc)):
        dima[abc[iii]] = {}
        dima[abc[iii]][abc[iii]]=scm[iii][iii]
    for iii in range(len(abc)):
        for jjj in range(len(abc)):
            dima[abc[iii]][abc[jjj]] = scm[iii][jjj]
                    
                
    return dima

def compute_alignment_matrix(seqx, seqy, scm, flag):
    """ Docstring """
    #print 'seqx',seqx,'seqy',seqy
    mmm = len(seqx)
    nnn = len(seqy)
    #print 'mmm',mmm,'nnn',nnn
    alm =  [[0]*(nnn + 1) for _ in range(mmm + 1)]
    #print alm
    if flag == True:
        for iii in range(1,mmm+1):
            alm[iii][0] = alm[iii - 1][0] + scm[seqx[iii-1]]['-']
            #print 'alm',alm[iii][0],'scm',scm[seqx[iii-1]]['-']
        for jjj in range(1,nnn+1):
            alm[0][jjj] = alm[0][jjj-1] + scm['-'][seqy[jjj-1]]
        
        for iii in range(1,mmm+1):
            for jjj in range(1,nnn+1):
                mval = max( alm[iii-1][jjj-1] + scm[seqx[iii - 1]][seqy[jjj - 1]],
                  alm[iii-1][jjj] + scm[seqx[iii - 1]]['-'],
                   alm[iii][jjj-1] + scm['-'][seqy[jjj - 1]] )
                #print 'mval',mval
                alm[iii][jjj] = mval
                #print 'alm-1',alm[iii-1][jjj-1],'scmseqseq',scm[seqx[iii - 1]][seqy[jjj - 1]]
    elif flag == False:
        for iii in range(1,mmm+1):
            mval = alm[iii - 1][0] + scm[seqx[iii-1]]['-']
            mval = max(mval,0)
            alm[iii][0] = mval
            #print 'alm',alm[iii][0],'scm',scm[seqx[iii-1]]['-']
        for jjj in range(1,nnn+1): 
            mval = alm[0][jjj-1] + scm['-'][seqy[jjj-1]]
            if mval < 0:
                mval = 0
            alm[0][jjj] = mval
        
        for iii in range(1,mmm+1):
            for jjj in range(1,nnn+1):
                mval = max( alm[iii-1][jjj-1] + scm[seqx[iii - 1]][seqy[jjj - 1]],
                  alm[iii-1][jjj] + scm[seqx[iii - 1]]['-'],
                   alm[iii][jjj-1] + scm['-'][seqy[jjj - 1]] )
                if mval < 0:
                    mval = 0
                #print 'mval',mval
                alm[iii][jjj] = mval
                #print 'alm-1',alm[iii-1][jjj-1],'scmseqseq',scm[seqx[iii - 1]][seqy[jjj - 1]]
    return alm
def compute_global_alignment(seqx, seqy, scm, alm):
    """ Docstring """
    iii = len(seqx)
    jjj = len(seqy)
    xprime = ''
    yprime = ''
    while iii != 0 and jjj != 0:
        if alm[iii][jjj] == alm[iii-1][jjj-1] + scm[seqx[iii-1]][seqy[jjj-1]]:
            xprime = seqx[iii-1] + xprime
            yprime = seqy[jjj-1] + yprime
            iii -= 1
            jjj -= 1
        else:
            if alm[iii][jjj] == alm[iii-1][jjj] + scm[seqx[iii-1]]['-']:
                xprime = seqx[iii-1] + xprime
                yprime = '-' + yprime
                iii -= 1
            else:
                xprime = '-' + xprime
                yprime = seqy[jjj-1] + yprime
                jjj -= 1
    while iii != 0:
        xprime = seqx[iii-1] + xprime
        yprime = '-' + yprime
        iii -= 1
    while jjj != 0:
        xprime = '-' + xprime
        yprime = seqy[jjj-1] + yprime
        jjj -= 1
    #print 'xprime',xprime
    #print 'yprime',yprime
    scor = 0
    for iii in range(len(xprime)):
        scor += scm[xprime[iii]][yprime[iii]]
    return (scor,xprime,yprime)

def compute_local_alignment(seqx, seqy, scm, alm):
    """ Docstring """
    (scc, qxx, qyy) = compute_global_alignment(seqx, seqy, scm, alm)
    #sss = compute_alignment_matrix(qx, qy,scm,True)
    #print 'sc',sc,'qx',qx,'qy',qy
    #print 'sss',sss
    return (scc, qxx, qyy)

################# TESTS ###########################
#scmat = build_scoring_matrix('ACGT', 10 ,4, -4)
#print scmat['A']['T']
#print build_scoring_matrix(set(['A', 'C', 'T', 'G']), 6, 2, -4)
#expected {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2},
#'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2},
#'-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4},
#'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2}, 
#'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}} 
#but received (Exception: TypeError) "cannot concatenate 'str' and 'set' objects" at line 4, in build_scoring_matrix
#print compute_alignment_matrix('A', 'A', {'A': {'A': 6, 'C': 2, '-': -4,'T': 2, 'G': 2},
#                      'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2},
#                      '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4}, 
#                      'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2}, 
#                      'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}}, True)
#expected [[0, -4], [-4, 6]] 
#but received [[0]]
#print compute_alignment_matrix('ATG', 'ACG', {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2}, 
#                                        'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2},
#                                        '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4},
#                                        'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2}, 
#                                        'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}}, True)
##expected [[0, -4, -8, -12], [-4, 6, 2, -2], [-8, 2, 8, 4], [-12, -2, 4, 14]] 
#but received [[0, -4, -8, -12], [-4, 6, -2, -6], [-8, -2, 8, 0], [-12, -6, 0, 14]]
#print compute_alignment_matrix('ACTACT', 'AGCTA', {'A': {'A': 2, 'C': 1, '-': 0, 'T': 1, 'G': 1},
#                                        'C': {'A': 1, 'C': 2, '-': 0, 'T': 1, 'G': 1},
#                                           '-': {'A': 0, 'C': 0, '-': 0, 'T': 0, 'G': 0},
#                                          'T': {'A': 1, 'C': 1, '-': 0, 'T': 2, 'G': 1}, 
#                                         'G': {'A': 1, 'C': 1, '-': 0, 'T': 1, 'G': 2}}, True)
#expected [[0, 0, 0, 0, 0, 0], [0, 2, 2, 2, 2, 2], [0, 2, 3, 4, 4, 4],
#[0, 2, 3, 4, 6, 6], [0, 2, 3, 4, 6, 8], [0, 2, 3, 5, 6, 8], [0, 2, 3, 5, 7, 8]] 
#but received (Exception: IndexError) "list index out of range" at line 42, in compute_alignment_matrix
#print compute_alignment_matrix('A', 'A', {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2},
#                                    'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2},
#                                    '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4},
#                                    'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2},
#                                    'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}}, False)
#expected [[0, 0], [0, 6]] 
#but received [[0, 0], [0, 0]]
#print compute_global_alignment('A', 'A', {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2},
#                                    'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2}, 
#                                    '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4},
#                                    'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2},
#                                    'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}}, 
#                         [[0, -4], [-4, 6]])  
#returned incorrect score, expected 6 but received 0
#print compute_global_alignment('ATG', 'ACG', {'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2},
#                                       'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2}, 
#                                       '-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4},
#                                       'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2},
#                                       'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}},
#                         [[0, -4, -8, -12], [-4, 6, 2, -2], [-8, 2, 8, 4], [-12, -2, 4, 14]])
#expected ({'A': {'A': 6, 'C': 2, '-': -4, 'T': 2, 'G': 2},
#     'C': {'A': 2, 'C': 6, '-': -4, 'T': 2, 'G': 2},
#'-': {'A': -4, 'C': -4, '-': -4, 'T': -4, 'G': -4}, 
#'T': {'A': 2, 'C': 2, '-': -4, 'T': 6, 'G': 2}, 
#'G': {'A': 2, 'C': 2, '-': -4, 'T': 2, 'G': 6}}, 14, 'ATG', 'ACG', True)
#but received (Exception: KeyError) "'GTA'" at line 115, in compute_global_alignment

print compute_local_alignment('abddcdeffgh', 'aabcddefghij',
{'-': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1},
'a': {'-': -1, 'a': 2, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1},
'c': {'-': -1, 'a': -1, 'c': 2, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'b': {'-': -1, 'a': -1, 'c': -1, 'b': 2, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'e': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': 2, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'd': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': 2, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'g': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': 2, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'f': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': 2, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'i': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': 2, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'h': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': 2, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'k': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': 2, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'j': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': 2, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'm': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': 2, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'l': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': 2, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'o': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': 2, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'n': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': 2, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'q': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': 2, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'p': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': 2, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
's': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': 2, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'r': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': 2, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'u': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': 2, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
't': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': 2, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'w': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': 2, 'v': -1, 'y': -1, 'x': -1, 'z': -1}, 
'v': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': 2, 'y': -1, 'x': -1, 'z': -1}, 
'y': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': 2, 'x': -1, 'z': -1}, 
'x': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': 2, 'z': -1}, 
'z': {'-': -1, 'a': -1, 'c': -1, 'b': -1, 'e': -1, 'd': -1, 'g': -1, 'f': -1, 'i': -1, 'h': -1, 'k': -1, 'j': -1, 'm': -1, 'l': -1, 'o': -1, 'n': -1, 'q': -1, 'p': -1, 's': -1, 'r': -1, 'u': -1, 't': -1, 'w': -1, 'v': -1, 'y': -1, 'x': -1, 'z': 2}}, 
[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 1, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 3, 3, 5, 4, 3, 2, 1, 0, 0, 0], [0, 0, 0, 2, 2, 5, 7, 6, 5, 4, 3, 2, 1], [0, 0, 0, 1, 4, 4, 6, 6, 5, 4, 3, 2, 1],
[0, 0, 0, 0, 3, 6, 6, 5, 5, 4, 3, 2, 1], [0, 0, 0, 0, 2, 5, 5, 8, 7, 6, 5, 4, 3], [0, 0, 0, 0, 1, 4, 4, 7, 10, 9, 8, 7, 6], [0, 0, 0, 0, 0, 3, 3, 6, 9, 9, 8, 7, 6], [0, 0, 0, 0, 0, 2, 2, 5, 8, 11, 10, 9, 8], [0, 0, 0, 0, 0, 1, 1, 4, 7, 10, 13, 12, 11]])  
#returned incorrect score, expected 13 but received 10