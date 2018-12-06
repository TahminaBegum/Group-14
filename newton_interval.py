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
        lower = qx
        higher = new_qx
        #print(lower, higher)
    elif(qx > new_qx):   
        lower = new_qx
        higher = qx
        #print(lower, higher)
    #mean_of_interval(lower, higher) != 0.01
    
    #IX=UpperInterval({cast_exact(lower):cast_exact(higher)})
    #X = cast_singleton(IX)

    #mx = X.value()
    #N = mx - ( f( mx ) / derivative( f, X ) )

    i=0
    #delta = derivative(f, qx) 
    print("lower",lower ,"higher",higher)
    mx = mean_of_interval(lower, higher)
    N_lower = mx - ( f( mx ) / derivative( f, lower ) )
    N_higher = mx - ( f( mx ) / derivative( f, higher ) )
    print("Nlower",N_lower ,"Nhigher",N_higher)

    while(i!=5):
        mx = mean_of_interval(lower, higher)
        N_lower = mx - ( f( mx ) / derivative( f, lower ) )
        N_higher = mx - ( f( mx ) / derivative( f, higher ) )
        print(lower, higher)
        print(N_lower, N_higher)
        if((N_lower > lower) and (N_higher < higher)):
            mx = mean_of_interval(lower, N_higher)
            #lower = mx - ( f( mx ) / derivative( f, lower ) )
            higher = mx - ( f( mx ) / derivative( f, N_higher ) )
            #higher = N_higher
            print("New bounds(1)", lower, N_higher)
            lower = mx - ( f( mx ) / derivative( f, lower ) )
            i=i+1
        elif((lower < N_lower) and (N_higher < higher)):
            mx = mean_of_interval(N_lower, N_higher)
            N_lower = mx - ( f( mx ) / derivative( f, N_lower ) )
            print("New bounds(2)", N_lower, N_higher)
            i=i+1
        elif((N_lower > lower) and (N_higher > higher)):
            mx = mean_of_interval(lower, higher)
            N_lower = mx - ( f( mx ) / derivative( f, lower ) ) 
            N_higher = mx - ( f( mx ) / derivative( f, higher ) )
            print("New bounds(3)", lower, higher)
            i=i+1
        elif((lower < N_lower) and (N_higher >= higher)):
            mx = mean_of_interval(N_lower, higher)
            N_higher = mx - ( f( mx ) / derivative( f, N_higher ) )
            print("New bounds(4)", N_lower, higher)
            i=i+1
        #if(decide(abs(delta)<=Ep)):
            #break
    print("Root in interval", lower, higher )   












    #if((N_lower > lower)):

       # if((N_higher > higher)):
          #  N_lower = lower
          #  higher = higher
          #  print(2)
       # if((N_higher < higher)):
         #   N_lower = lower
         #   higher = N_higher
        #    print(1)
    #elif((N_higher > higher)):
        #if((N_lower > lower)):
           # N_lower = lower
            #N_higher = higher
            #print(3)
        #if((N_higher < higher)):
           # N_lower = lower
            #higher = N_higher
           # print(4)
       

    #print(lower, higher)
    #if (refines(N_lower,lower) and refines( N_higher, higher)):
        #print("Verified root in",N_lower,"within",lower)
    #elif (inconsistent(N_lower,lower) and inconsistent(N_higher, higher)):
        #print("No root in",lower, N_higher)
    #else:
            #print("Any root in",lower,"also lies in",N_lower,"so in common refinement",refinement(lower,N_lower))

    #for (nl,nh) in range(N_lower, N_higher):
        #mx = mean_of_interval(N_lower, N_higher)
        #nl = mx - ( f( mx ) / derivative( fterval, N_lower ) )
        #nh = mx - ( f( mx ) / derivative( f, N_higher ) )
        #if((nl < N_lower) and (nh < N_higher)):
            #N_higher = nh
            #print("in for in",N_lower, n_higher)
            #print(IX)
            #break
        #break

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
        interval_newton_method(f, qx, new_qx, Ep)

    if choice == 2:
        initialization_number = int(input("Enter the first interval:"))
        #print(initialization_number)
       # q = cast_singleton(initialization_number)
        qx = FloatMPBounds(initialization_number)
        newton_step(f, qx)
        new_qx = newton_step(f, qx)
        interval_newton_method(f, qx, new_qx, Ep)

    if choice == 3:
        first = int(input("Enter the first interval:"))
        second = int(input("Enter the second interval:")) 
        fx = FloatMPBounds(first)
        sx = FloatMPBounds(second)
        interval_newton_method(f, fx,sx, Ep)
       
    




main()
