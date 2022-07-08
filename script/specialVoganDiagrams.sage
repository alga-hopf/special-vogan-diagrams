import argparse

'''
   This algorithm provides adjoint orbits of simple Lie groups admitting 
   a canonical special locally homogeneous compatible almost-complex structure.
   For the details, see https://link.springer.com/article/10.1007/s00209-022-02995-9.
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
		
	args = parser.parse_args()

	lieType = args.type
	rank = args.rank

	print(DynkinDiagram([lieType,rank]))
	W=WeylGroup([lieType,rank],implementation='permutation')
	positiveRoots=W.positive_roots()	# Positive roots
	C=CartanMatrix([lieType,rank])	# Cartan matrix of the Lie algebra   
	if lieType=='F':	# Scalar product induced on the roots
		B=(1/36)*matrix([[4,-2,0,0],[-2,4,-2,0],[0,-2,2,-1],[0,0,-1,2]])
	else: 
		B=matrix(QQ,gap('BilinearFormMat(RootSystem(SimpleLieAlgebra("'+lieType+'",'+str(rank)+',Rationals)))'))	
	print("Dimension:",2*len(positiveRoots)+rank)
	print('  ')
	for P in [q for q in Combinations(range(rank)) if q!=[]]:
		compactroots=[root for root in positiveRoots if sum(root[k] for k in P)%2==0]	# Compact roots
		noncompactroots=[root for root in positiveRoots if sum(root[k] for k in P)%2!=0]	# Non-compact roots
		epsilon={root: (1 if root in noncompactroots else -1) for root in positiveRoots}	# Compactness coefficient
		eta=-2*sum(epsilon[alpha]*alpha for alpha in positiveRoots)	# Eta vector
		phi0=(C)*(eta-2*sum(root for root in positiveRoots if all(root[k]==0 for k in P)))	# phi0 vector
		isSpecial=1
		if all(phi0[k]==0 for k in range(len(phi0))):	# If phi0 is 0, then the orbit is symplectic Calabi-Yau
			print([var('t'+str(k)) if k in P else 0 for k in range(len(phi0))],' for all ti>0     Non-compact simple roots:',P,'  ','symplectic Calabi-Yau')
			print("Dimension V:",(2*sum(all(root[k]==0 for k in P) for root in positiveRoots)+rank),"    Dimension G/V:",2*(len(positiveRoots)-sum(all(root[k]==0 for k in P) for root in positiveRoots)))
			print("Hermitian scalar curvature:",0)
			isSpecial=0
		elif [sgn(phi0[k]) for k in range(len(phi0))]==[-1 if k in P else 0 for k in range(len(phi0))]:	# If phi0 has negative entries in correspondence of posnoncomproot, then the orbit
			Omega0=[root for root in Set(positiveRoots).difference(Set([root for root in positiveRoots if all(root[k]==0 for k in P)]))]	# is sympl. general type
			Omega0nc=[root for root in Set(Omega0).intersection(Set(noncompactroots))]
			typeOrbit=''
			if len(P)>1:
				typeOrbit='symplectic general type'
			else:
				for alpha in Omega0nc:
					for beta in Omega0nc:
						if alpha+beta in Omega0:
							typeOrbit='symplectic general type'
							break
				if typeOrbit!='symplectic general type':
					typeOrbit='general type'
			print(-phi0/(gcd(phi0))	,'   ','Non-compact simple roots:',P,'   ',typeOrbit)
			print("Dimension V:",(2*(len(positiveRoots)-len(Omega0))+rank),"    Dimension G/V:",2*len(Omega0))
			print("Hermitian scalar curvature:",4*gcd(phi0)*sum(sum(epsilon[alpha]*((alpha*B*beta)/(sum(phi0[z]*beta[z]*B[z,z] for z in range(rank)))) for alpha in Omega0) for beta in Omega0))
			isSpecial=0
		elif [sgn(phi0[k]) for k in range(len(phi0))]==[1 if k in P else 0 for k in range(len(phi0))]:	# If phi0 has positive entries in correspondence of posnoncomproot, then the orbit
			T=Combinations([i for i in range(rank) if i not in P])					# is symplectic Fano e there may be more orbits with this diagram
			for i in range(T.cardinality()):
				S=T[i]+P
				phi=(C)*(eta-2*sum(root for root in positiveRoots if all(root[k]==0 for k in S)))
				if [sgn(phi[k]) for k in range(len(phi))]==[1 if k in S else 0 for k in range(len(phi))]:
					print(phi*(1/gcd(phi)),'   ','Non-compact simple roots:',P,'   ','S:',S,'   ','symplectic Fano')
					Omega=[root for root in Set(positiveRoots).difference(Set([root for root in positiveRoots if all(root[k]==0 for k in S)]))]
					print("Dimension V:",(2*(len(positiveRoots)-len(Omega))+rank),"    Dimension G/V:",2*len(Omega))
					print("Hermitian scalar curvature:",-4*gcd(phi)*sum(sum(epsilon[alpha]*((alpha*B*beta)/(sum(phi[z]*beta[z]*B[z,z] for z in range(rank)))) for alpha in Omega) for beta in Omega))
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
						
		# Compute the Lie algebra of the orbit in the case symplectic Calabi-Yau and symplectic general type
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

	print('Time:',cputime())
