from ariadne import *
import random


def Intialization():
    print("1. x*x-2")
    print("2. x*x*x-2")
    print("3. 6*x*x+4*x-50")
    print("4. x*x-200")

def function_for_execution(function_number,x):
    #f=['x*x-2','x*x*x-10','4*x*x-50','x*x-200']
    if function_number==1:
        f= x*x-2
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
    r=random.randint(1,4)
    r1=r+2
    IX=UpperInterval({cast_exact(r):cast_exact(r1)})
    X=cast_singleton(IX)
    return X


def Get_Estimator(f,X):
    k=2         # I can't apply float value
                # If lowerbound of x is negative then it suffer inf problem
    
    x_mid=X.value()
    new_X=x_mid - f(x_mid)/derivative(f,X)
  
    if inconsistent(new_X,X):
        new_X=x_mid - k*(f(x_mid)/derivative(f,X))
        print("shift")
        return Get_Estimator(f,new_X)
   
    else:
        return new_X



if __name__=='__main__':
    
    #initialization with Some functions
    Intialization()
    function_number=int(input("Enter your function(1 to 4):"))
  
    #Find our Selected function
    x = EffectiveScalarUnivariateFunction.identity()
    f=function_for_execution(function_number,x)
    print("my_function: ",f)
    
    X=Get_Random_Number()
    print("The random number",X)

    print("Root is in the range=",Get_Estimator(f,X))
    
   


