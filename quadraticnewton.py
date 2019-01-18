from ariadne import *
import random

def Intialization():
	print("1. a*x*x+b*y*y+c*x*y+d")

def function_for_execution(function_number,x,a,b,c,d):
	if function_number == 1:
		f = a*x[0]*x[0]+b*x[1]*x[1]+c*x[0]*x[1]+d
	return f

#def function_for_show(function_number,a,b,c,d):
#    if function_number == 1:
#        f = a*x*x + b*x*x + c*x*y+ d
#    return f

def Get_Random_Number():
	#Select initial estimator
	#r=2.0
	r = random.randint(1,9)
	#r1 = r+2
	rand = FloatDPBounds(r)
	return rand


	zero=FloatDPApproximation(0)		#zero for checking
	print("checking whether ", fxk, " < ", zero)
	ff=derivative(f,xk)
	print("Partial derivatives of f: ", ff)
	print (type(fxk<zero))
	print (type(decide(fxk<zero)))
	print("fxk:",fxk)
	if (decide(fxk < zero)):
		print("It is: Box containing root found")
		return 0
	elif (fxk > zero):
		print("It is not")

def steepest_Descent_Step(x0, y0, i_a, i_b, i_c, i_d):

    x = x0
    y = y0  
    a = i_a
    b = i_b
    c = i_c
    d = i_d
    while True:
        fk = a*x*x + b*y*y +c*x*y+d
        if(decide(fk<0)):
            print("Root found")
            break
        else:
            prevxk = x
            prevyk = y
            ffk1 = 2*a*x +c*y
            print("ffk1", ffk1)
            ffk2 = 2*b*y +c*x
            print("ffk2", ffk2)
            sk1 = (c*c - 4*a*b)*(2*b*ffk1 -c*ffk2)
            print("sk1", sk1)
            sk2 = (c*c - 4*a*b)*(2*a*ffk2 - c*ffk1)
            print("sk2", sk2)
            alphak = (sk1*sk1 + sk2*sk2)/(2*sk1*sk1+ 2*sk2*sk2+ 2*c*sk1*sk2) 
            print("alpha", alphak)

            xk = x + alphak * sk1
            yk = y + alphak * sk2
            print("x" ,xk, "y", yk)
            currentx = xk
            currenty = yk
            print("previous x is", prevxk," and previous y is", prevyk, "current x is", xk, "and current y is",yk)
            
            #if(decide((prevxk != x) and (prevyk != y))):
             #   prevxk = currentx
             #   prevyk = currenty
             #   currentxk = xk
              #  currentyk = yk
               # print("1.previous x and y" , prevxk, prevyk, "current x and y", currentxk, currentyk)
            #else:
             #   print("2.previous x and y" , prevxk, prevyk, "current x and y", currentx, currenty)



            if(decide((xk <= 0) and (yk <= 0))):
                break           
            
    

def main():

    #initialization with Some functions
    Intialization()
    function_number = int(input("Enter your function choice:"))
    print(function_number)

    #Find our Selected function
    #x = FloatMPBounds.identity()
    #y = FloatMPBounds.identity()
    
    
	#To get the value of Epsilon
	#E=UpperInterval({cast_exact(.00001):0})
	#X=cast_singleton(E)
	#Ep=FloatDPApproximation(X)

    #x=EffectiveScalarUnivariateFunction.identity()
    #y=EffectiveScalarUnivariateFunction.identity()
    #ff = function_for_execution(function_number,x,y,a,b,c,d)
    
    dpr=DoublePrecision()
    argument_size=2
    x=[EffectiveScalarMultivariateFunction.coordinate(argument_size,index)
	for index in
		range(0,argument_size)]
   
	#Find our Selected function
    input_a = int(input("Enter the first value of function:"))
    input_b = int(input("Enter the second value of function:")) 
    input_c = int(input("Enter the third value of function:"))
    input_d = int(input("Enter the constant value of function:")) 
    i_a = FloatMPBounds(input_a)
    i_b = FloatMPBounds(input_b)
    i_c = FloatMPBounds(input_c)
    i_d = FloatMPBounds(input_d)
    input_x = int(input("Enter the value of x0:"))
    input_y = int(input("Enter the value of y0:"))
    ff = function_for_execution(function_number,x,input_a,input_b,input_c,input_d)
    #fx = function_for_show(function_number,a,b,c,d)
    print("my function", ff)
    steepest_Descent_Step(input_x, input_y, i_a, i_b, i_c, i_d)				#while steepest descent gives positive value, repeat

main()
