import multiprocessing
import time
import itertools
import numpy as np
from datetime import datetime
from Bio import SeqIO

def matrix_building(ref, query,rlen,qlen,match_score=2,mismatch_cost=1, gap_cost=1):
    #H = np.zeros((qlen + 1,rlen + 1), np.int)
    #T = np.zeros((qlen, rlen), np.int)
    H=[ [ 0 for i in range(rlen+1) ] for j in range(qlen+1) ]
    T=[ [ 0 for i in range(rlen) ] for j in range(qlen) ]
    max_score=0
    max_i=0
    max_j=0
    for i, j in itertools.product(range(1,qlen+1), range(1,rlen+1)):
        match = H[i - 1][j - 1] + (match_score if query[i - 1] == ref[j - 1] else - mismatch_cost)
        delete = H[i][j-1] - gap_cost
        insert = H[i-1][ j] - gap_cost
        max_val=max(match, delete, insert, 0)
        if(max_val>=max_score):
            max_score=max_val
            max_i=i
            max_j=j
        H[i][j] = max_val
        if(max_val==0):
            T[i-1][j-1]=-1
        elif(max_val==insert):
            T[i-1][j-1]=2
        elif(max_val==delete):
            T[i-1][j-1]=3
        elif(max_val==match):
            if query[i - 1] == ref[j - 1]:
                T[i-1][j-1]=1
            else:
                T[i-1][j-1]=10

    return H, T, max_i, max_j,max_score

def tracing(T,i,j,prev_i,prev_j, match_, mismatch_,insert_, del_):
    while(i>=0 and j>=0):
        op=T[i][j]
        if(op==-1):
            return prev_i, prev_j,match_,mismatch_, insert_, del_
        prev_i=i
        prev_j=j
        if(op==1):
            match_=match_+1
            i=i-1
            j=j-1
        elif(op==10):
            mismatch_=mismatch_+1
            i=i-1
            j=j-1
        elif(op==2):
            insert_=insert_+1
            i=i-1
            j=j
        elif(op==3):
            del_=del_+1
            i=i
            j=j-1
    return prev_i,prev_j,match_,mismatch_, insert_, del_



'''
#a, b = 'GGTTGACTA', 'TGTTACGG'
a, b = 'ACACACTA','AGCACACA'
print(b) 
start, end = smith_waterman(a, b)
print(a[start:end])     # GTTGAC
'''
#print(matrix_building('CCGTACTA','CAGACCTA'))

def smith_waterman(ref,query,rlen,qlen):
    #now = datetime.now().time()
    #print("matrix building started at: ",now)
    H, T, max_i, max_j,max_score=matrix_building(ref, query,rlen,qlen)
    #now = datetime.now().time()
    #print("matrix building finished at: ", now)
    start_q,start_ref,num_match, num_mismatch, num_insert,num_del=tracing(T,max_i-1,max_j-1,0,0,0,0,0,0)
    now = datetime.now().time()
    #print("tracing finished at: ",now)
    return start_q,max_i-1,start_ref,max_j-1,num_match,num_mismatch,num_insert,num_del

def func(x):
	ref=x[0]
	con=x[1]
	r=ref
	q=con
	ref_len=len(r)
	con_len=len(con)
	start_q,end_q,start_ref,end_ref,num_match,num_mismatch,num_insert,num_del=smith_waterman(r,q,ref_len,con_len)
	ret_val=[start_q,end_q,start_ref,end_ref,num_match,num_mismatch,num_insert,num_del,ref.id]

	return ret_val

if __name__ == '__main__':
	
	references = list(SeqIO.parse('references.fasta','fasta'))
	contigs = list(SeqIO.parse('mega_assembly/final.contigs.fa','fasta'))
	num_ref=len(references)
	contig_scores={}
	goodcontig_count=0
	num_contig=0
	now=datetime.now().time()
	print(now)
	for con in contigs:
		con_len=len(con)
		inputs=[]
		for i in range(num_ref):
			temp=[references[i],con]
			inputs.append(temp)		
		pool = multiprocessing.Pool()
		pool = multiprocessing.Pool(processes=num_ref)
		outputs = pool.map(func, inputs)
		quality_max=-1
		refname_max=''
		for res in outputs:
			m=res[4]  ## number of match
			mm=res[5]+res[6]+res[7] ## number of mismatch=number of substitution + number of insertion + number of deletion
			alignment_len=m+mm
			match_ratio=m/con_len
			quality=match_ratio
			if(quality>quality_max):
				quality_max=quality
				refname_max=res[8]
				contig_scores[con.id]=[refname_max,m,res[5],res[6],res[7],res[0],res[1],res[2],res[3]]
	print(contig_scores)
	now=datetime.now().time()
	print(now)


