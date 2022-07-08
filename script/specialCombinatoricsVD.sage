'''
   This algorithm checks if a given adjoint orbits of simple Lie groups admitting 
   a canonical special locally homogeneous compatible almost-complex structure.
   For the details, see ....
   The Lie algebra of the orbit is computed using the theorems in the paper "Equivalence classes of Vogan diagrams", https://www.sciencedirect.com/science/article/pii/S0021869303007269

   Copyright (C) 2022  Alice Gatti

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''


import numpy as np
import time
import argparse

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument(
		"--type",
		default='A',
		help="Lie type of the Lie algebra",
		type=str)
	parser.add_argument(
		"--rank",
		default=10,
		help="Rank of the Lie algebra",
		type=int)
	parser.add_argument(
		"--s",
		default=[1],
		help="Indices of painted nodes (starting from 1)",
		type=int,nargs='+')
	parser.add_argument(
		"--cores",
		default='1',
		help="Number of cores to be used",
		type=int)

	
	args = parser.parse_args()

	lieType = args.type
	rank = args.rank
	S = args.s
	cores = args.cores

	# Delta function
	def delta(x,y):
		try:
		    x
		    y
		except IndexError:
		    return 0
		if type(y)!=list:
		    y=[y]
		if x in y:
		    return 1
		else:
		    return 0

	# Functions that compute the coefficients for the different types of Lie algebras

	# An
	def tau_an(x,s,n):
		s_ext=[0]+s+[n+1]
		j=np.argwhere(np.array(s)==x)[0][0]+1
		s1=sum(s_ext[j+2*k]-s_ext[j+2*k-1] for k in range(1,ceil((len(s)-j)/2)+1))
		s2=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		t1=sum(s_ext[j+2*k+1]-s_ext[j+2*k] for k in range(1,ceil((len(s)-(j+1))/2)+1))
		t2=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))
		return s1+s2-(t1+t2)

	@parallel(cores)
	def xi_coeff_an(x,s,rank):
		s_ext=[0]+s+[rank+1]
		j=np.argwhere(np.array(s)==x)[0][0]+1
		return s_ext[j-1]-s_ext[j+1]+2*(tau_an(x,s,rank))

	# Bn
	def tau_bn(x,s,n):
		s_ext=[0]+s+[n+1]
		j=np.argwhere(np.array(s)==x)[0][0]+1
		if x!=rank and s[-1]!=x:
		    s1_an=sum(s_ext[j+2*k]-s_ext[j+2*k-1] for k in range(1,ceil((len(s)-j)/2)+1))
		    s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		    t2_an=sum(s_ext[j+2*k+1]-s_ext[j+2*k] for k in range(1,ceil((len(s)-(j+1))/2)+1))
		    t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))
		    s1=2*s1_an+t1_an+(s_ext[j]-s_ext[j-1]-1)-delta(s_ext[-1],s_ext[j+2*ceil((len(s)-j)/2)])
		    s2=s2_an
		    t1=t1_an
		    t2=2*t2_an+(s_ext[j+1]-s_ext[j]-1)+s2_an-delta(s_ext[-1],s_ext[j+2*ceil((len(s)-(j+1))/2)+1])
		    return s1+s2-(t1+t2)
		elif x!=rank and s[-1]==x:
		    s1_an=sum(s_ext[j+2*k]-s_ext[j+2*k-1] for k in range(1,ceil((len(s)-j)/2)+1))
		    s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		    t2_an=sum(s_ext[j+2*k+1]-s_ext[j+2*k] for k in range(1,ceil((len(s)-(j+1))/2)+1))
		    t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))
		    s1=2*s1_an+t1_an+(s_ext[j]-s_ext[j-1]-1)
		    s2=s2_an
		    t1=t1_an
		    t2=t2_an+s2_an
		    return s1+s2-(t1+t2)
		else:
		    s1_an=sum(s_ext[j+2*k]-s_ext[j+2*k-1] for k in range(1,ceil((len(s)-j)/2)+1))
		    t2_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))
		    s1=s1_an+sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(floor((j+1)/2)))-1-delta(s_ext[-1],s_ext[j+2*ceil((len(s)-j)/2)])
		    t2=t2_an
		    return 2*(s1-t2)

	@parallel(cores)
	def xi_coeff_bn(x,s,rank):
		s_ext=[0]+s+[rank+1]
		j=np.argwhere(np.array(s)==x)[0][0]+1
		if x==s[-1] and x!=rank:
		    return -2*rank+s_ext[j]+s_ext[j-1]+2*(tau_bn(x,s,rank))
		elif x!=s[-1] and x!=rank:
		    return s_ext[j-1]-s_ext[j+1]+2*(tau_bn(x,s,rank))
		else:
		    return -2*rank+2*s_ext[j-1]+2*(tau_bn(x,s,rank))

	# Cn
	def tau_cn(x,s,rank):
		s_ext=[0]+s+[rank+1]
		j=np.argwhere(np.array(s)==x)[0][0]+1
		if x!=rank-1:
		    if x!=s[-1]:
		        if x!=rank and rank not in s_ext:
		            s1_an=sum(s_ext[j+2*k]-s_ext[j+2*k-1] for k in range(1,ceil((len(s)-j)/2)+1))
		            s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		            t2_an=sum(s_ext[j+2*k+1]-s_ext[j+2*k] for k in range(1,ceil((len(s)-(j+1))/2)+1))
		            t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))

		            s1_bn=2*s1_an+t1_an+(s_ext[j]-s_ext[j-1]-1)-delta(s_ext[-1],s_ext[j+2*ceil((len(s)-j)/2)])
		            s2_bn=s2_an
		            t1_bn=t1_an
		            t2_bn=2*t2_an+(s_ext[j+1]-s_ext[j]-1)+s2_an-delta(s_ext[-1],s_ext[j+2*ceil((len(s)-(j+1))/2)+1])

		            s1=s1_bn+2-delta(s_ext[-1],s_ext[j+2*ceil((len(s)-j)/2)])
		            s2=s2_bn
		            t1=t1_bn
		            t2=t2_bn+2-delta(s_ext[-1],s_ext[j+2*ceil((len(s)-(j+1))/2)+1])
		            return s1+s2-(t1+t2)

		        elif x!=rank and rank in s_ext:
		            s1_an=sum(s_ext[j+2*k]-s_ext[j+2*k-1] for k in range(1,ceil((len(s)-j)/2)+1))
		            s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		            t2_an=sum(s_ext[j+2*k+1]-s_ext[j+2*k] for k in range(1,ceil((len(s)-(j+1))/2)+1))
		            t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))
		            s1=s1_an+t2_an+(s_ext[j+1]-s_ext[j])+s2_an-delta(s_ext[-2],s_ext[j+2*ceil((len(s)-(j+1))/2)+1])
		            s2=s2_an
		            t1=t1_an
		            t2=t2_an+s1_an+(s_ext[j]-s_ext[j-1])+t1_an-delta(s_ext[-2],s_ext[j+2*ceil((len(s)-(j+1))/2)+1])
		            return s1+s2-(t1+t2)
		    else:
		        if x==rank:
		            s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		            t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))
		            delta_n=s2_an 
		            delta_nm=s2_an+t1_an 
		            return -delta_nm+2*delta_n
		        else:
		            s1_an=sum(s_ext[j+2*k]-s_ext[j+2*k-1] for k in range(1,ceil((len(s)-j)/2)+1))
		            s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		            t2_an=sum(s_ext[j+2*k+1]-s_ext[j+2*k] for k in range(1,ceil((len(s)-(j+1))/2)+1))
		            t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))

		            s1_bn=2*s1_an+t1_an+(s_ext[j]-s_ext[j-1]-1)
		            s2_bn=s2_an
		            t1_bn=t1_an
		            t2_bn=t2_an+s2_an

		            s1=s1_bn+2-delta(s_ext[-1],s_ext[j+2*ceil((len(s)-j)/2)])
		            s2=s2_bn
		            t1=t1_bn
		            t2=t2_bn

		            return s1+s2-(t1+t2)

		else:
		    if s[-1]==x: 
		        t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))
		        s1=2+t1_an+(s_ext[j]-s_ext[j-1]-1)
		        t1=t1_an
		        return s1-t1

		    elif s[-1]==rank: 
		        s1_an=sum(s_ext[j+2*k]-s_ext[j+2*k-1] for k in range(1,ceil((len(s)-j)/2)+1))
		        s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		        t2_an=sum(s_ext[j+2*k+1]-s_ext[j+2*k] for k in range(1,ceil((len(s)-(j+1))/2)+1))
		        t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))
		        rem=s2_an-t1_an-(s_ext[j]-s_ext[j-1]-1)
		        s1=s2_an+1
		        t1=t1_an
		        return rem+s1-t1

	@parallel(cores)    
	def xi_coeff_cn(x,s,rank):
		s_ext=[0]+s+[rank+1]
		j=np.argwhere(np.array(s)==x)[0][0]+1
		if x!=rank-1:
		    if x!=s[-1] or x==rank:
		        return s_ext[j-1]-s_ext[j+1]+2*(tau_cn(x,s,rank))
		    else:
		        return -2-(s_ext[j]-s_ext[j-1]-1)-2*(rank-s_ext[j])+2*(tau_cn(x,s,rank))
		else:
		    return -3-s_ext[j]+s_ext[j-1]+2*(tau_cn(x,s,rank))

	# Dn
	def tau_dn(x,s,rank):
		s_ext=[0]+s+[rank+1]
		j=np.argwhere(np.array(s)==x)[0][0]+1
		if x not in range(rank-2,rank+1):
		    if x!=s[-1] and len(set(s).intersection(set([rank-1,rank])))==0:
		        if s_ext[-2]==rank:
		            s_ext[-2]=rank-1
		        s1_an=sum(s_ext[j+2*k]-s_ext[j+2*k-1] for k in range(1,ceil((len(s)-j)/2)+1))
		        s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		        t2_an=sum(s_ext[j+2*k+1]-s_ext[j+2*k] for k in range(1,ceil((len(s)-(j+1))/2)+1))
		        t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))

		        s1=s1_an+delta(s_ext[-1],s_ext[j+2*ceil((len(s)-j)/2)])+t1_an+s1_an-3*delta(s_ext[-1],s_ext[j+2*ceil((len(s)-j)/2)])+(s_ext[j]-s_ext[j-1]-1)
		        s2=s2_an
		        t1=t1_an
		        t2=t2_an+delta(s_ext[-1],s_ext[j+2*ceil((len(s)-(j+1))/2)+1])+s2_an+t2_an-3*delta(s_ext[-1],s_ext[j+2*ceil((len(s)-(j+1))/2)+1])+(s_ext[j+1]-s_ext[j]-1)

		        return s1+s2-(t1+t2)

		    elif x!=s[-1] and len(set(s).intersection(set([rank-1,rank])))==1:
		        if s_ext[-2]==rank:
		            s_ext[-2]=rank-1
		        s1_an=sum(s_ext[j+2*k]-s_ext[j+2*k-1] for k in range(1,ceil((len(s)-j)/2)+1))
		        s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		        t2_an=sum(s_ext[j+2*k+1]-s_ext[j+2*k] for k in range(1,ceil((len(s)-(j+1))/2)+1))
		        t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))

		        s1=s1_an+s2_an+t2_an+(s_ext[j+1]-s_ext[j])
		        s2=s2_an
		        t1=t1_an
		        t2=t2_an+t1_an+s1_an+(s_ext[j]-s_ext[j-1])-delta(s_ext[-3],x)
		        return s1+s2-(t1+t2)

		    elif x!=s[-1] and len(set(s).intersection(set([rank-1,rank])))==2: 
		        s1_an=sum(s_ext[j+2*k]-s_ext[j+2*k-1] for k in range(1,ceil((len(s)-j)/2)+1))
		        s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		        t2_an=sum(s_ext[j+2*k+1]-s_ext[j+2*k] for k in range(1,ceil((len(s)-(j+1))/2)+1))
		        t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))
		        s1=t1_an+(s_ext[j]-s_ext[j-1]-1)+s1_an+s1_an-3+2*delta(rank,s_ext[j+2*ceil((len(s)-j)/2)])
		        s2=s2_an
		        t1=t1_an
		        t2=s2_an+t2_an+(s_ext[j+1]-s_ext[j]-1)+2*delta(rank,s_ext[j+2*ceil((len(s)-j-1)/2)+1])-3+t2_an
		        return s1+s2-(t1+t2)

		    elif x==s[-1]:
		        s1_an=sum(s_ext[j+2*k]-s_ext[j+2*k-1] for k in range(1,ceil((len(s)-j)/2)+1))
		        s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		        t2_an=sum(s_ext[j+2*k+1]-s_ext[j+2*k] for k in range(1,ceil((len(s)-(j+1))/2)+1))
		        t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))

		        s1=t1_an+(s_ext[j]-s_ext[j-1]-1)
		        s2=s2_an
		        t1=t1_an
		        t2=s2_an

		        return s1+s2-(t1+t2)
		else:
		    if x==rank-2:
		        if len(set([rank-1,rank]).intersection(set(s)))==0: 
		            return s_ext[j]-s_ext[j-1]-1
		        elif len(set([rank-1,rank]).intersection(set(s)))==1: 
		            s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		            t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))
		            return 1+2*s2_an-2*t1_an-(s_ext[j]-s_ext[j-1]-1)   
		        else:
		            s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		            return s_ext[j]-s_ext[j-1]+1
		    else:          
		        if len(set([rank-1,rank]).intersection(set(s)))==1:
		            s_ext[j]=rank-1
		            s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		            t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))
		            return 2*(s2_an-t1_an)
		        else:
		            if s_ext[j]==rank:
		                j=j-1
		            s2_an=sum(s_ext[j-2*k+1]-s_ext[j-2*k] for k in range(1,floor(j/2)+1))
		            t1_an=sum(s_ext[j-2*k]-s_ext[j-2*k-1] for k in range(1,floor((j+1)/2)))
		            return (rank-2-s_ext[j-1])

	@parallel(cores)        
	def xi_coeff_dn(x,s,rank):
		s_ext=[0]+s+[rank+1]
		j=np.argwhere(np.array(s)==x)[0][0]+1
		if x not in range(rank-2,rank+1):
		    if x!=s[-1]:
		        if s_ext[-2]==rank-1:
		            s_ext[-2]=rank
		        return s_ext[j-1]-s_ext[j+1]+2*(tau_dn(x,s,rank))
		    elif x==s[-1]:
		        return -2*rank+s_ext[j]+s_ext[j-1]+1+2*(tau_dn(x,s,rank))
		elif x==rank-1 or x==rank:
		    if rank-1 in s and s_ext[j]==rank:
		        j=j-1
		    return -2-(2-delta(rank,s)*delta(rank-1,s))*(rank-2-s_ext[j-1])+2*(tau_dn(x,s,rank))
		else:
		    return -rank+s_ext[j-1]+1-(1-delta(rank-1,s))-(1-delta(rank,s))+2*(tau_dn(x,s,rank))

	# Special diagram
	def is_special(rank,s,lieType):
		if lieType=='A':
		    coeff=xi_coeff_an([(l,s,rank) for l in s])
		elif lieType=='B':
		    coeff=xi_coeff_bn([(l,s,rank) for l in s])
		elif lieType=='C':
		    coeff=xi_coeff_cn([(l,s,rank) for l in s])
		elif lieType=='D':
		    coeff=xi_coeff_dn([(l,s,rank) for l in s])
		coeff=list(coeff)
		signs=[]
		for i in range(len(s)):
		    signs.append(sign(list(coeff)[i][-1]))
		if len(np.unique(signs))==1:
		    if signs[0]==0:
		    	return 'Symplectic Calabi-Yau'
		    elif signs[0]==1:
		    	return 'Symplectic Fano'
		    else:
		    	return 'Symplectic general type'
		else:
		    return 'No'
		    
	P=[S[i]-1 for i in range(len(S))]
	isSpecial=1
	
	if lieType in ['A','B','C','D']:
		print('')
		print(lieType+str(rank))
		if lieType=='A':
			d=rank*(rank+2)
		elif lieType=='B':
			d=rank*(2*rank+1)
		elif lieType=='C':
			d=rank*(2*rank+1)
		elif lieType=='D':
			d=rank*(2*rank-1)
		print('Dimension:',d)
		print('')
		if rank<25:
			print(DynkinDiagram([lieType,rank]))
			print('')
		print('Painted nodes:',S)
		print('')
		t0=time.time()
		result=is_special(rank,S,lieType)
		print('Is the diagram special? ',result)
		t1=time.time()-t0
		if result!='No':
			isSpecial=0
	else:
		print('')
		print(lieType+str(rank))
		W=WeylGroup([lieType,rank],implementation='permutation')
		positiveRoots=W.positive_roots()	# Positive roots
		C=CartanMatrix([lieType,rank])	# Cartan matrix of the Lie algebra   
		if lieType=='F':	# Scalar product induced on the roots
			B=(1/36)*matrix([[4,-2,0,0],[-2,4,-2,0],[0,-2,2,-1],[0,0,-1,2]])
		else: 
			B=matrix(QQ,gap('BilinearFormMat(RootSystem(SimpleLieAlgebra("'+lieType+'",'+str(rank)+',Rationals)))'))	
		print("Dimension:",2*len(positiveRoots)+rank)
		print(' ')
		print(DynkinDiagram([lieType,rank]))
		print('')
		print('Painted nodes:',s)
		print('')
		P=[S[i]-1 for i in range(len(S))]
		t0=time.time()
		compactroots=[root for root in positiveRoots if sum(root[k] for k in P)%2==0]	# Compact roots
		noncompactroots=[root for root in positiveRoots if sum(root[k] for k in P)%2!=0]	# Non-compact roots
		epsilon={root: (1 if root in noncompactroots else -1) for root in positiveRoots}	# Compactness coefficient
		eta=-2*sum(epsilon[alpha]*alpha for alpha in positiveRoots)	# Eta vector
		phi0=(C)*(eta-2*sum(root for root in positiveRoots if all(root[k]==0 for k in P)))	# phi0 vector
		isSpecial=1
		if all(phi0[k]==0 for k in range(len(phi0))):
			print('Is the diagram special? Symplectic Calabi-Yau')	# If phi0 is 0, then the orbit is symplectic Calabi-Yau
			isSpecial=0
		elif [sgn(phi0[k]) for k in range(len(phi0))]==[-1 if k in P else 0 for k in range(len(phi0))]:	# If phi0 has negative entries in correspondence of posnoncomproot, then the orbit
			print('Is the diagram special? Symplectic general type')
			isSpecial=0
		elif [sgn(phi0[k]) for k in range(len(phi0))]==[1 if k in P else 0 for k in range(len(phi0))]:	# If phi0 has positive entries in correspondence of posnoncomproot, then the orbit
			print('Is the diagram special? Symplectic Fano')
			isSpecial=0
		else:
			print('Is the diagram special? No')
		t1=time.time()-t0
		
	# Compute the Lie algebra of the orbit and the stabilizer
	if isSpecial==0:
		# Compute the Lie algebra of the orbit
		if lieType=='A':
			equivclass=sum((-1)^(len(P)-s)*(P[s-1]+1) for s in range(1,len(P)+1))
			if equivclass<=((rank+1)/2).floor():
				eqclass=equivclass
			else:
				eqclass=rank+1-equivclass
			print('Lie algebra: su('+str(eqclass)+','+str(rank+1-eqclass)+')')
		elif lieType=='B':
			eqclass=sum((-1)^(len(P)-s)*(P[s-1]+1) for s in range(1,len(P)+1))
			print('Lie algebra: so('+str(2*eqclass)+','+str(2*rank-2*eqclass+1)+')')
		elif lieType=='C':
			if rank-1 in P:
				print('Lie algebra: sp('+str(rank)+',R)')
			else:
				N=sum((-1)^(len(P)-s)*(P[s-1]+1) for s in range(1,len(P)+1))
				if N<= rank/2:
					eqclass=N
				else:
					eqclass=rank-N
				print('Lie algebra: sp('+str(eqclass)+','+str(rank-eqclass)+')')
		elif lieType=='D':	
			if (rank-2 in P and rank-1 not in P) or (rank-2 not in P and rank-1 in P) :
				print('Lie algebra: so*('+str(2*rank)+')')
			elif Set([rank-2,rank-1]) in P:
				N=sum((-1)^(len(P)-s)*(P[s-1]+1) for s in range(1,len(P)-1))
				if N<=rank/2:
					eqclass=N-1
				else:
					eqclass=rank-N-1
				print('Lie algebra: so('+str(2*eqclass)+','+str(2*rank-2*eqclass)+')')
			else:
				N=sum((-1)^(len(P)-s)*(P[s-1]+1) for s in range(1,len(P)+1))
				if N<=rank/2:
					eqclass=N
				else:
					eqclass=rank-N
				print('Lie algebra: so('+str(2*eqclass)+','+str(2*rank-2*eqclass)+')')
		elif lieType=='G':
			print('Lie algebra: g2(2)')
		elif lieType=='F':
			if Set(P).intersection(Set([0,1]))!=Set([]):
				print('Lie algebra: f4(4)')
			else:
				print('Lie algebra: f4(-20)')
		elif lieType=='E' and rank==6:
			II=[j for j in P if j<=3 and j!=1]
			JJ=[j for j in P if j>3]
			if 1 in P:
				s=1
			else:
				s=0
			if II!=[] or JJ!=[]:
				if II!=[] and 0 not in II:
					I=sum((-1)^(len(II)-a-1)*(II[a]) for a in range(len(II)))
				elif 0 in II:
					I=(-1)^(len(II)-1)+sum((-1)^(len(II)-a)*(II[a]) for a in range(len(II)))
				else:
					I=0
				J=sum((-1)^(len(JJ)-a-1)*(JJ[a]) for a in range(len(JJ)))
				print(I,J)
				if P==[0] or P==[5] or P==[2,4] or P==[0,3,4] or P==[0,1] or P==[1,2] or P==[1,4] or P==[1,5] or P==[1,3,5] or (len(JJ)!=1 and J==2-I and (I+s)%2==1) or (len(JJ)!=1 and J==4-I and (I+s)%2==0) or (len(JJ)!=1 and J==1-I) or (len(JJ)==1 and ((J==4+I and (I+s)%2==1) or J==1+I)):
					print('Lie algebra: e6(-14)')
				else:
					print('Lie algebra: e6(2)')				
			else: 
				print('Lie algebra: e6(2)')
		elif lieType=='E' and rank==7:
			II=[j for j in P if j<=3 and j!=1]
			JJ=[j for j in P if j>3]
			if 1 in P:
				s=1
			else:
				s=0
			if II!=[] or JJ!=[]:
				if II!=[] and 0 not in II:
					I=sum((-1)^(len(II)-a-1)*(II[a]) for a in range(len(II)))
				elif 0 in II:
					I=(-1)^(len(II)-1)+sum((-1)^(len(II)-a)*(II[a]) for a in range(len(II)))
				else:
					I=0
				J=sum((-1)^(len(JJ)-a-1)*(JJ[a]) for a in range(len(JJ)))
				if P==[0] or P==[3] or P==[5] or P==[3,5] or P==[3,4,6] or P==[1,4] or P==[1,6] or P==[1,2,4] or P==[0,1,3,4] or (len(JJ)!=1 and (((J==1-I or J==3-I) and (I+s)%2==1) or ((J==2-I or J==4-I) and (I+s)%2==0) )) or (len(JJ)==1 and (((J==1+I or J==2+I or J==3+I or J==5+I) and (I+s)%2==0) or (J==4+I and (I+s)%2==1))):
					print('Lie algebra: e7(-5)')
				elif P==[6] or P==[2,4] or P==[0,3,4] or P==[0,1] or P==[1,2] or P==[1,5] or P==[1,3,5] or P==[1,3,4,6] or P==[1,3,4,5,6] or (len(JJ)!=1 and ((J==1-I and (I+s)%2==0) or (J==2-I and (I+s)%2==1))) or (len(JJ)==1 and ((J==1+I or J==2+I or J==5+I) and (I+s)%2==1)):
					print('Lie algebra: e7(-25)')	
				else:
					print('Lie algebra: e7(7)')			
			else: 
				print('Lie algebra: e7(7)')
		elif lieType=='E' and rank==8:
			II=[j for j in P if j<=3 and j!=1]
			JJ=[j for j in P if j>3]
			if 1 in P:
				s=1
			else:
				s=0
			if II!=[] or JJ!=[]:
				if II!=[] and 0 not in II:
					I=sum((-1)^(len(II)-a-1)*(II[a]) for a in range(len(II)))
				elif 0 in II:
					I=(-1)^(len(II)-1)+sum((-1)^(len(II)-a)*(II[a]) for a in range(len(II)))
				else:
					I=0
				J=sum((-1)^(len(JJ)-a-1)*(JJ[a]) for a in range(len(JJ)))
				if P==[7] or P==[2] or P==[3] or P==[0,2] or P==[1,2] or P==[1,5] or P==[1,6] or (len(JJ)!=1 and ((J==1-I or J==5-I) and (I+s)%2==0)) or (len(JJ)!=1 and (J==3-I and (I+s)%2==1)) or (len(JJ)!=1 and (J==2-I or J==6-I)) or (len(JJ)==1 and ((J==1+I or J==5+I) and (I+s)%2==1)) or (len(JJ)==1 and ((J==3+I and (I+s)%2==0) or J==2+I or J==6+I)):
					print('Lie algebra: e8(-24)')
				else:
					print('Lie algebra: e8(8)')					
			else: 
				print('Lie algebra: e8(8)')
		# Compute the Lie algebra of the stabilizer
		if len(P)<rank:
			st=str(CartanType([lieType,rank]).subtype([i+1 for i in range(rank) if i not in P])).translate({ord('['):None,ord("'"):None,ord(','):None,ord(']'):None,ord(' '):None})
			stab=''
			for i in range(len(st)):
				if st[i] in ['A','B','C','D','E']:
					if st[i]=='A':
						stab=stab+'su('
						for j in range(i+1,len(st)):
							if st[j]=='r' or st[j]=='x' or j==len(st)-1:
								if j==len(st)-1:
									stab=stab+str(int(st[i+1:len(st)])+1)+') x '
									break
								else:
									stop=j
									stab=stab+str(int(st[i+1:stop])+1)+') x ' 
									break
					if st[i]=='B':
						stab=stab+'so('
						for j in range(i,len(st)):
							if st[j]=='r' or st[j]=='x' or j==len(st)-1:
								if j==len(st)-1:
									stab=stab+str(int(st[i+1:len(st)])*2+1)+') x '
									break
								else:
									stop=j
									stab=stab+str(int(st[i+1:stop])*2+1)+') x ' 
									break
					if st[i]=='C':
						stab=stab+'sp('
						for j in range(i,len(st)):
							if st[j]=='r' or st[j]=='x' or j==len(st)-1:
								if j==len(st)-1:
									stab=stab+str(int(st[i+1:len(st)]))+') x '
									break
								else:
									stop=j
									stab=stab+str(int(st[i+1:stop]))+') x ' 
									break
					if st[i]=='D':
						stab=stab+'so('
						for j in range(i,len(st)):
							if st[j]=='r' or st[j]=='x' or j==len(st)-1:
								if j==len(st)-1:
									stab=stab+str(int(st[i+1:len(st)])*2)+') x '
									break
								else:
									stop=j
									stab=stab+str(int(st[i+1:stop])*2)+') x ' 
									break
					if st[i]=='E':
						if st[i+1]=='6':
							stab=stab+'e6 x '
						if st[i+1]=='7':
							stab=stab+'e7 x '
			print('Stabilizer: '+stab[:-2]+' x R'+str(len(P))+'\n')
		else:
			print('Stabilizer: R'+str(len(P))+'\n')
		
		print('Time:',t1)
