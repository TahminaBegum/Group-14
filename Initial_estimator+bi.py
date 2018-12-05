from ariadne import *

def Intialization():
    print("1. x*x-2")
    print("2. (x*x-x)*(x*x+x)*(x*x-x)")
    print("3. (x-9)*(x*x-1)*(x+1)")
    print("4. (x-2)*(x-4)*(x-1)")


def function_for_execution(function_number,x):
    if function_number==1:
        f= x*x-2
    elif function_number==2:
        f=x*x*x-2
    elif function_number==3:
        f=(x-9)*(x*x-1)*(x+1)
    else:
        f=(x-2)*(x-4)*(x+1)
    
    return f



def Make_Interval(f,X_Input):
    k=2
    ax_FDA=FloatDPApproximation(X_Input)
    bx_FDA=ax_FDA - k*(f(ax_FDA)/derivative(f,ax_FDA))
    ax_Float=float(str(ax_FDA))       #This part, I also confuse because It reduce a lots of data
    bx_Float=float(str(bx_FDA))
    #ax_Float=FloatDPValue(ax_FDA)
    #bx_Float=FloatDPValue(bx_FDA)
    if ax_Float<bx_Float:
        IX=UpperInterval({cast_exact(ax_Float):cast_exact(bx_Float)})
    else:
        IX=UpperInterval({cast_exact(bx_Float):cast_exact(ax_Float)})
    
    X=cast_singleton(IX)
    X_end=FloatDPApproximation(X)
    #print("Print X_end=",X_end)
    return X_end

def Make_Interval1(f,Input_X):      #I didn't use any NM step
    IX=UpperInterval({cast_exact(0.0):Input_X})
    X=cast_singleton(IX)
    return FloatDPApproximation(X)

def Newton_method_Approximation_step(f,ax,Ep,i):
    while True:
        i=i+1
        print("ax:=",ax)
        h=f(ax)/derivative(f, ax)
        ax=ax-h
        if(decide(abs(h)<=Ep)):
            break
    print("Root in Approximation: ",ax)
    return i

def Newton_method_Approximation_step_bi(f,ax,bx,Ep,i):
    while True:
        i=i+1
        
        print("ax=",ax," and bx=",bx)
        
        h_a=f(ax)/derivative(f, ax)
        ax=ax-h_a
        if(decide(abs(h_a)<=Ep)):
            print("Root in Approximation: ",ax)
            return i
        
        h_b=f(bx)/derivative(f, bx)
        bx=bx-h_b
        if(decide(abs(h_b)<=Ep)):
            print("Root in Approximation: ",bx)
            return i


def Get_Estimator(f,X,i):
 
    fX=Evaluate(f,X)
    
    if decide(fX>0):
        Sign=1
    else:
        Sign=-1
    k_step=1
    for k in range(1,50000):
        i=i+1
        k_step=k_step*2
        print("i in GE=",i)
        X_new=X - k_step*f(X)/derivative(f,X)
        fX_new=Evaluate(f,X_new)
        
        if decide(fX_new>0):
            Sign_new=1
        else:
            Sign_new=-1
        print("Sign in GE :=",Sign," and ",Sign_new)
        print("Sign in GE :=",X," and ",X_new)
        if not(Sign==Sign_new):
           return X,X_new,i
        X=X_new
            
    print("limit need to Increass")

    

def Evaluate(f,X):
    return f(X)



if __name__=='__main__':
    
    #initialization with Some functions
    Intialization()
    function_number=int(input("Enter your function(1 to 4):"))
  
    #Find our Selected function
    x = EffectiveScalarUnivariateFunction.identity()
    f=function_for_execution(function_number,x)
    print("my_function: ",f)
    
    #To get the value of Epsilon
    E=UpperInterval({cast_exact(.00001):0})
    X=cast_singleton(E)
    Ep=FloatDPApproximation(X)
    
    Input_X=int(input("Please input a number: "))
    A=Make_Interval(f,Input_X)
    #A=Make_Interval1(f,Input_X)
    
    # print("The Initial interval", A)
    #print("The type of X_in=",type(A))
    i=0
    A_c,B,step=Get_Estimator(f,A,i)
    print("The total steps in estimator:= ",step)
    print("i=",i)
    print(" Root is in the range= ",A," to ",B)

    step_NM=Newton_method_Approximation_step(f,B,Ep,i)
    print("The total steps in NM estimator:= ",step_NM+step)

    step_E=Newton_method_Approximation_step(f,A,Ep,i)
    print("The total steps in without initial estimator:= ",step_E)

    step_E_BI=Newton_method_Approximation_step_bi(f,A_c,B,Ep,i)
    print("The total steps in with initial estimator with bi:= ",step_E_BI+step)
