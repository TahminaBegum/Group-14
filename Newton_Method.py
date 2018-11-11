from ariadne import *

def Intialization():
    print("1. x*x-2")
    print("2. x*x*x-10")
    print("3. 6*x*x+4*x-50")
    print("4. x*x-200")

def function_for_execution(function_number,x):
    #f=['x*x-2','x*x*x-10','4*x*x-50','x*x-200']
    if function_number==1:
        f= x*x-2
    elif function_number==2:
        f=x*x*x-10
    elif function_number==3:
        f=6*x*x+4*x-50
    else:
        f=x*x-200
    
    return f


def Newton_method_iteration(f,ax):
    i=0
    while i!=1000:
        h = f(ax)/derivative(f, ax)
        print("h=",h)
        ax=ax-h
        i=i+1
    print("Root in",ax)
    print("iteration",i-1)


def Newton_method_Approximation(f,ax,Ep):
    while True:
        h = f(ax)/derivative(f, ax)
        print("Ep : ",Ep)
        print("h = ",abs(h))
        ax=ax-h
        if(abs(h)<=Ep): #I can't compare the Ariadne data type
            break
    print("Root in Approximation: ",ax)




def main():
    
    #initialization with Some functions
    Intialization()
    function_number=int(input("Enter your function(1 to 4):"))
    print(function_number)
    
    #To get the value of Epsilon
    E=UpperInterval({cast_exact(.0001):0})
    X=cast_singleton(E)
    Ep=FloatDPApproximation(X)
    
    #print our selection menu
    print("1.Iteration_based")
    print("2.Approximation_based")
    choice=int(input("Please Select(1/2): "))
    
    #Find our Selected function
    x = EffectiveScalarUnivariateFunction.identity()
    f=function_for_execution(function_number,x)
    print("my_function: ",f)
    
    #Select initial estimator
    IX=UpperInterval({cast_exact(1.4):2})
    X=cast_singleton(IX)
    ax=FloatDPApproximation(X)
    print("Initial Estimator",ax)
    
    #function which we selected
    if choice==1:
        Newton_method_iteration(f,ax)
    else:
        Newton_method_Approximation(f,ax,Ep)

main()
