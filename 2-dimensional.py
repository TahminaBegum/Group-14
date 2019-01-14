from ariadne import *
import random

print(dir())

def Intialization():
	print("1. x*x+y*y-1")
	print("2. 4*x*x+y*y-1")

def function_for_execution(function_number,x):
	if function_number == 1:
		f=x[0]*x[0]+x[1]*x[1]-1
	elif function_number == 2:
		f=4*x[0]*x[0]+x[1]*x[1]-1
	return f

def Get_Random_Number():
	#Select initial estimator
	#r=2.0
	r = random.randint(1,9)
	#r1 = r+2
	rand = FloatDPBounds(r)
	return rand

def steepest_Descent_Step(xk,f):
	dpr=DoublePrecision()
	fxk=FloatDPApproximation(f(xk))
	zero=FloatDPApproximation(0)		#zero for checking
	print("checking whether ", fxk, " < ", zero)
	print (type(fxk<zero))
	print (type(decide(fxk<zero)))
	print("fxk:",fxk)
	if (decide(fxk < zero)):
		print("It is: Box containing root found")
		return 0
	elif (fxk > zero):
		print("It is not")
		#ff=derivative(f)
		ff0=derivative(f,0)		#gradient
		ff1=derivative(f,1)
		ff=[ff0,ff1]
		print("Partial derivatives of f: ", ff)		
		#sk=-1*ff(xk)			#search direction
		#print("sk ",sk)	
		#alpha = Get_Random_Number()	#initialise random variable alpha
		#print("alpha ",alpha)
		#expr=xk+alpha*sk
		#print("expr ",expr)
		#G=f(expr)
		#print("G ",G)
		#g=derivative(G)
		#print("g ", g)		
		#h = g @ sk			# dot product between vectors g and sk
		alpha = 0.1			
		#alpha = 1DNewton(h)		# apply 1D newton method to h to find alpha
		#xk1=xk+alpha*sk			# find next point
		return xk1

	print("WHY THIS NOT CRASHING?????")
	#print(steepest_Descent_Step(xk,f))
	#while(steepest_Descent_Step(xk,f) > 0 ):

def main():
	#initialization with Some functions
	Intialization()
	function_number = int(input("Enter your function(1 to 2):"))
	print(function_number)
    
	#To get the value of Epsilon
	#E=UpperInterval({cast_exact(.00001):0})
	#X=cast_singleton(E)
	#Ep=FloatDPApproximation(X)

	#Find our Selected function

	dpr=DoublePrecision()
	argument_size=2
	x=[EffectiveScalarMultivariateFunction.coordinate(argument_size,index)
	for index in
		range(0,argument_size)]
	
	f = function_for_execution(function_number,x)
	print("my_function: ",f)
	xk=FloatDPApproximationVector([1,1],dpr)

	steepest_Descent_Step(xk,f)				#while steepest descent gives positive value, repeat

main()
