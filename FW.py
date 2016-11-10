import numpy as np
import pandas as pd
from scipy import optimize 
from scipy.optimize import minimize_scalar
import scipy.linalg as la
import numpy.linalg as la2
import scipy.integrate as integrate 


def estimatetime(t0,xa,ca):
	ta = t0*(1+0.15*(xa/ca)**4)
	return ta

def estimateZ(alpha,xa,ca,t0,ya):
	Z = 0
	for i in range(len(xa)):
		Z += integrate.quad(lambda x: estimatetime(t0[i],x,ca[i]),0,xa[i]+alpha*(ya[i]-xa[i]))[0]
	return Z

def linearsearch(xa,ca,t0,ya):
	alpha = minimize_scalar(estimateZ, args=(xa, ca, t0,ya), bounds = (0,1), method = 'Bounded')
	return alpha.x



# main functions
#### Step 1: Network Representation and Data Structure
## Define the link-node matrix 
LinkNode = pd.read_csv("linknode2.csv", header = None)
LinkNode = LinkNode.as_matrix()
#print LinkNode
#print LinkNode.shape()
#print LinkNode

## Import Demand matrix (Q)
Q = pd.read_csv("Q.csv", header = None)
Q = Q.as_matrix()
#print Q

## create travel time vector (ta)
n = 76 # number of total links
k = 24 # number of total nodes


## create initial link flow vector (X)
#X = np.zeros(n)


## creat link flow matrix for iterations (Y)
s = (n,k) # 76 links, 24 nodes/origins
#Y = np.zeros(s) # each entry represents the flow on link a from origin i 

## import the travel time coeff. estimation matrix (Coeff)
coeff = pd.read_csv("coeff2.csv", header = None)
coeff = coeff.as_matrix()
#print coeff
#print coeff.shape




### Step 2: Shortest Path Searching (Solve LP)

##Initialization

t0 = coeff [:,0] # free flow travel time from 1st column of *Coeff*
ca = coeff[:,1] # capacity for each link 


#ya = np.sum(Y,axis = 1)  # column sums of Y: total flow on link a from all origins
#Z = np.dot(np.transpose(t0),ya) # k*k matrix
#LEF = np.dot(np.transpose(LinkNode),Y) # a k*k matrix, each column is for each origin i; within each column, there's one row for origin constraint, one row for destination constraint, and k-2 rows for zero constraint

origq = np.sum(Q, axis = 1) # row sums of Q: total flow from origin i 
destq = -Q 
s2 = (k,k)
RHT = np.zeros(s2)# each row represents an origin, each column represents the flow on node k (with origin i )
RHT = -Q
np.fill_diagonal(RHT, origq) 
#print RHT


c0_0 = np.transpose(t0)
c_0 = np.tile(c0_0,k)

A0 = np.transpose(LinkNode) # Construct block matrix for A_eq
A1 = [A0]*k
A = la.block_diag(*A1) 
#print A

b0 = np.transpose(RHT)
b = np.ravel(b0, order = 'F')[:,np.newaxis] # construct long b_eq

ybounds = (0, None)
result = optimize.linprog(
	c_0, A_eq = A, b_eq = b, bounds = (ybounds), options = {"disp":True, "maxiter":2000,"bland":True} 
	)
#print result
#print len(result['x'])
result = np.reshape(result['x'],(k,n))
#print result
xa = np.sum(result, axis = 0) # intialization xa
#print xa
ta = estimatetime(t0,xa,ca)


###############
step = 0
tanorm = 1000000

iteration = []
Z = []

while (tanorm>7.6): # allow each link has 0.1 diff. in ta on average
	### Update 
	print "step ", step
	iteration.append(step)
	ta_old = ta
	

	### direction
	c0 = np.transpose(ta)
	c = np.tile(c0,k)
	result = optimize.linprog(
	c, A_eq = A, b_eq = b, bounds = (ybounds), options = {"disp":True, "maxiter":2000,"bland":True} 
	)
	#resultz = result['fun'] #print objective value
	#print "z is",resultz
	resultx = np.reshape(result['x'],(k,n))
	ya = np.sum(resultx, axis = 0) # yn

	### move
	alpha = linearsearch(xa,ca,t0,ya)
	print "alpha is", alpha
	xa = (1-alpha)*xa + alpha * ya 
	print "xa is ",xa

### Update 
	ta = estimatetime(t0,xa,ca)
	tanorm = la2.norm(ta-ta_old)
	z = np.dot(np.transpose(xa),ta)
	Z.append(z)
	print "ta is ", ta
	print "norm of ta is ", tanorm
	step +=1 



import matplotlib.pyplot as plt
plt.plot(iteration,Z,'ro')
plt.xlabel('iteration number')
plt.ylabel('Z(x)')
plt.show()


np.savetxt("ta", ta, delimiter = ",")
np.savetxt("xa", xa, delimiter = ",")

