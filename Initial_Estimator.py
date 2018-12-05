from ariadne import *

def Intialization():
    print("1. x*x-2")
    print("2. x*x*x-2")
    print("3. 6*x*x+4*x-50")
    print("4. (x-2)*(x-4)*(x-1)")


def function_for_execution(function_number,x):
    if function_number==1:
        f= x*x-2
    elif function_number==2:
        f=x*x*x-2
    elif function_number==3:
        f=6*x*x+4*x-50
    else:
        f=(x-2)*(x-4)*(x-1)
    
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


def Get_Estimator(f,X):
 
    #print("X = ",X)
    fX=Evaluate(f,X)
    #print("fX = ",fX)
    #print("type of the fX = ",type(fX))
    
    if decide(fX>0):
        Sign=1
    else:
        Sign=-1

            #print("Sign1 = ",Sign)

    for k in range(1,5):
       
        X_new=X - 2*k*f(X)/derivative(f,X)
        fX_new=Evaluate(f,X_new)
        
        #print("X_new = ",X_new)
        #print("K = ",k)
        #print("fX_new = ",fX_new)
        if decide(fX_new>0):
            Sign_new=1
        else:
            Sign_new=-1

        #print("Sign_new = ",Sign_new)
        if not(Sign==Sign_new):
           return X_new
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
    
    Input_X=int(input("Please input a number: "))
    A=Make_Interval(f,Input_X)
    #A=Make_Interval1(f,Input_X)
    
    # print("The Initial interval", A)
    #print("The type of X_in=",type(A))
    
    B=Get_Estimator(f,A)
    
    print(" Root is in the range= ",A," to ",B)
