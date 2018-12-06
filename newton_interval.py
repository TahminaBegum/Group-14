from ariadne import *
import random
#include <floatmp.hpp>

def Intialization():
    print("1. x*x-10")
    print("2. x*x*x-2")
    print("3. 6*x*x+4*x-50")
    print("4. x*x-200")

def function_for_execution(function_number,x):
    #f=['x*x-2','x*x*x-10','4*x*x-50','x*x-200']
    if function_number==1:
        f= x*x-10
    elif function_number==2:
        f=x*x*x-2
    elif function_number==3:
        f=6*x*x+4*x-50
    else:
        f=x*x-200
    
    return f

def Get_Random_Number():
    #Select initial estimator
    #r=2.0
    r = random.randint(1,9)
    #r1 = r+2
    rand = FloatDPBounds(r)
    return rand

def newton_step(f, qx):
    new_qx = qx - f(qx)/derivative(f,qx)

    #if(new_qx < qx):
        #print("Initial interval : ",new_qx , qx  )
        #print(new_ax)
    #elif(qx < new_qx):
        #print("Initial interval : ",qx , new_qx  )
       # print(new_ax)
    return new_qx

def random_newton_step(f):
    qx = Get_Random_Number()
    new_qx = qx - f(qx)/derivative(f,qx)

    #if(new_qx < qx):
        #print("Initial interval : ", new_qx , qx )
        #print(new_ax)
    #else:
        #print("Initial interval : ", qx , new_qx  )
        #print(new_ax)
    return qx ,new_qx

def mean_of_interval(qx, new_qx):
    mean = (qx + new_qx)/2
    return mean

#def inm_formula(f, qx, new_qx):

    


def interval_newton_method(f, qx, new_qx , Ep):

    if(qx < new_qx):
        lower = FloatDPBounds(qx)
        higher = FloatDPBounds(new_qx)
        #print(lower, higher)
    elif(qx > new_qx):   
        lower = FloatDPBounds(new_qx)
        higher = FloatDPBounds(qx)
        #print(lower, higher)
    #mean_of_interval(lower, higher) != 0.01
    i=0
    delta = derivative(f, qx) 
    mx = mean_of_interval(lower, higher)
    N_lower = mx - ( f( mx ) / derivative( f, lower ) )
    N_higher = mx - ( f( mx ) / derivative( f, higher ) )
    
    #while(delta!=0):
       # i=i+1
       # interval_newton_method(f, N_lower, N_higher , Ep)



   # if((N_lower > lower) and (N_higher < higher )):
  #      IX=UpperInterval({cast_exact(N_lower):cast_exact(N_higher)})
   # elif((N_lower > lower) and (N_higher > higher )):
   #     IX=UpperInterval({cast_exact(N_lower):cast_exact(higher)})
   # elif((N_lower < lower) and (N_higher < higher )):
   #     IX=UpperInterval({cast_exact(lower):cast_exact(N_higher)})
    #elif((N_lower < lower) and (N_higher > higher )):
    #    IX=UpperInterval({cast_exact(lower):cast_exact(higher)})

   # X=cast_singleton(IX)
   # inm=FloatDPApproximation(X)
    #print("Print X_end=",X_end)
   # return inm

    



    #print(lower, higher)
    #while(i!=5):
        #mx = mean_of_interval(lower, higher)
        #N_lower = mx - ( f( mx ) / derivative( f, lower ) )
        #N_higher = mx - ( f( mx ) / derivative( f, higher ) )
        #print(lower, higher)
        #print(N_lower, N_higher)
        #if((N_lower >= lower) and (N_higher < higher)):
            #mx = mean_of_interval(lower, N_higher)
            #N_lower = lower 
            #print("New bounds(1)", N_lower, N_higher)
            #return lower, N_higher
            #i=i+1
            
            #if((lower < N_lower) and (N_higher < higher)):
                #mx = mean_of_interval(N_lower, N_higher)
                #N_lower = mx - ( f( mx ) / derivative( f, N_lower ) )
                #print("New bounds(2)", N_lower, N_higher)
                #i=i+1
            #elif((N_lower > lower) and (N_higher > higher)):
                #mx = mean_of_interval(lower, higher)
                #N_lower = mx - ( f( mx ) / derivative( f, lower ) ) 
                #N_higher = mx - ( f( mx ) / derivative( f, higher ) )
                #print("New bounds(3)", lower, higher)
                #i=i+1
            #elif((lower < N_lower) and (N_higher >= higher)):
                #mx = mean_of_interval(N_lower, higher)
                #N_higher = mx - ( f( mx ) / derivative( f, N_higher ) )
                #print("New bounds(4)", N_lower, higher)
                #i=i+1



        #elif((lower < N_lower) and (N_higher < higher)):
            #mx = mean_of_interval(N_lower, N_higher)
            #N_lower = mx - ( f( mx ) / derivative( f, N_lower ) )
            #print("New bounds(2)", N_lower, N_higher)
            #i=i+1
        #elif((N_lower > lower) and (N_higher > higher)):
            #mx = mean_of_interval(lower, higher)
            #N_lower = mx - ( f( mx ) / derivative( f, lower ) ) 
            #N_higher = mx - ( f( mx ) / derivative( f, higher ) )
            #print("New bounds(3)", lower, higher)
            #i=i+1
        #elif((lower < N_lower) and (N_higher >= higher)):
            #mx = mean_of_interval(N_lower, higher)
            #N_higher = mx - ( f( mx ) / derivative( f, N_higher ) )
            #print("New bounds(4)", N_lower, higher)
            #i=i+1
        #if(decide(abs(delta)<=Ep)):
            #break
    #print("Root in interval", lower, higher )   
    return     

def contractor(f, qx, new_qx , Ep):

    interval_newton_method(f,qx,new_qx,Ep)





def main():
    
    #initialization with Some functions
    Intialization()
    function_number = int(input("Enter your function(1 to 4):"))
    print(function_number)    
    
    #To get the value of Epsilon
    E=UpperInterval({cast_exact(.00001):0})
    X=cast_singleton(E)
    Ep=FloatDPApproximation(X)

    #Find our Selected function
    x = EffectiveScalarUnivariateFunction.identity()
    f = function_for_execution(function_number,x)
    print("my_function: ",f)
    
    #initialize function
    print("1.Random initialization of function")
    print("2.Initialize interval by hand for 1 value")
    print("3.Initialize interval entirely by hand")
    choice = int(input("Enter your choice:"))

    IX=UpperInterval({cast_exact(1.4):2})
    X=cast_singleton(IX)
    ax=FloatDPApproximation(X)

    if choice == 1:
        random_newton_step(f)
        qx, new_qx = random_newton_step(f)
        contractor(f, qx, new_qx, Ep)

    if choice == 2:
        initialization_number = int(input("Enter the first interval:"))
        #print(initialization_number)
       # q = cast_singleton(initialization_number)
        qx = FloatDPBounds(initialization_number)
        newton_step(f, qx)
        new_qx = newton_step(f, qx)
        contractor(f, qx, new_qx, Ep)

    if choice == 3:
        first = int(input("Enter the first interval:"))
        second = int(input("Enter the second interval:")) 
        fx = FloatDPBounds(first)
        sx = FloatDPBounds(second)
        contractor(f, fx,sx, Ep)
       
    




main()
