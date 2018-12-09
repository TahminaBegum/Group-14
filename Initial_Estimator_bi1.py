from ariadne import *

def Intialization():
    print("1. x*x-2")
    print("2. (x*x-x)*(x*x+x)*(x*x-x)")
    print("3. (x-9)*(x*x-1)*(x+1)+1")
    print("4. (x-2)*(x-4)*(x-1)")

def function_for_execution(function_number,x):
    """Input:
        function_number -- choice of function assigned to number,
        x               -- function variable is x
        Output:
        f               -- function chosen with variable x"""
    
    if function_number==1:
        f= x*x-2
    elif function_number==2:
        f=x*x*x-2
    elif function_number==3:
        f=(x-9)*(x*x-1)*(x+1)+1
    else:
        f=(x-2)*(x-4)*(x+1)
    
    return f



def Make_Interval(f,X_Input):
    """ Make an interval [ax, bx] from single input variable of function f.
        First bound set ax = X_Input, bound bx by performing a single Newton Step (k=2): bx = ax - k*(f(ax)/f'(b)).
        Depending on direction of Newton step, the bounds are set in the interval accordingly.
        FloatDPApproximation is used to get an average value of the constructed interval, which is set as output
        
        Input:
        f               -- function chosen with variable x
        X_Input         -- initial approximation of root by user
        Output:
        X_end           -- single value, average of interval [ax, bx] or [bx, ax]
        """
    k=2
    ax_FDA=FloatDPApproximation(X_Input)
    bx_FDA=ax_FDA - k*(f(ax_FDA)/derivative(f,ax_FDA))
    ax_Float=float(str(ax_FDA))
    bx_Float=float(str(bx_FDA))
    if ax_Float<bx_Float:
        IX=UpperInterval({cast_exact(ax_Float):cast_exact(bx_Float)})
    else:
        IX=UpperInterval({cast_exact(bx_Float):cast_exact(ax_Float)})
    
    X=cast_singleton(IX)
    X_end=FloatDPApproximation(X)
    return X_end

def Make_Interval3(X1,X1_new):
    #X=FloatXP(X1)      #if anyone have idea, please fix it
    #X_new=FloatXPValue(X1_new)
    X=float(str(X1))
    X_new=float(str(X1_new))
    if X<X_new:
        IX=UpperInterval({cast_exact(X):cast_exact(X_new)})
    else:
        IX=UpperInterval({cast_exact(X_new):cast_exact(X)})
    
    X=cast_singleton(IX)
    X_end=FloatDPBounds(X)
    return X_end


def Newton_method_Approximation_step(f,ax,Ep,i):
    """ Counting newton step's till the approximation is sufficiently close to root from single point ax (<= Epsilon)
        Input:
        f           -- function chosen with variable x
        ax          -- input approximation
        Ep          -- Epsilon (= .00001)
        i           -- counter of newton step's
        Output:
        i           -- count of newton step's
        """
    while True:
        i=i+1
        #print("ax:=",ax)
        h=f(ax)/derivative(f, ax)
        ax=ax-h
        if(decide(abs(h)<=Ep)):
            break
#print("Root in Approximation: ",ax)
    return i

def Newton_method_Approximation_step_bi(f,ax,bx,Ep,i):
    """ Counting newton step's till the approximation is sufficiently close to root from both directions of the interval [ax,bx] (<= Epsilon)
        Input:
        f           -- function chosen with variable x
        ax          -- input approximation
        Ep          -- Epsilon (= .00001)
        i           -- counter of newton step's
        Output:
        i           -- count of newton step's
        """
    while True:
        i=i+1
        
        #print("ax=",ax," and bx=",bx)
        
        h_a=f(ax)/derivative(f, ax)
        ax=ax-h_a
        if(decide(abs(h_a)<=Ep)):
            #print("Root in Approximation: ",ax)
            return i
        
        h_b=f(bx)/derivative(f, bx)
        bx=bx-h_b
        if(decide(abs(h_b)<=Ep)):
            #print("Root in Approximation: ",bx)
            return i


def Get_Estimator(f,X,i):
    """ Single point taken from initial interval is evaluted by its sign, function performs interval Newton method looking for sign difference. If no sign difference is found, k (Newton step size) is doubled. Sign difference is proof of a root (intermediate value theorem).
        Function is terminated as soon as sign difference is found (limited by k=50000)
        
        Input:
        f               -- function chosen with variable x
        X               -- single point taken from inital interval
        i               -- counter of newton steps performed
        Output:
        X_new           -- upperbound of interval with proof of root
        i               -- counter of newton steps performed
        
        """
 
    fX=f(X)
    
    if decide(fX>0):
        Sign=1
    else:
        Sign=-1
    k_step=1
    h=fX/derivative(f,X)
    
    for k in range(1,50000):
        i=i+1
        k_step=k_step*2
        X_new=X - k_step*h
        fX_new=f(X_new)
        if decide(fX_new>0):
            Sign_new=1
        else:
            Sign_new=-1
        
        if not(Sign==Sign_new):
           return X_new,i

    print("limit need to Increase")






if __name__=='__main__':
    """ The main function makes use of user iteraction to define function and first approximation of root.
        A           -- average point of the interval made around the approximation input
        B           -- point that has a sign difference with A
        Interval [A,B] has therefor proof of an existing root.
        step_E      -- step counter of without initial estimator from A
        step_NM     -- step counter of Newton Method estimator evaluated from B
        step_E_BI   -- step counter of Newton Method estimator from both directions
        
        """
    
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
    
    Input_X=float(input("Input a approximation to root: "))
    
    A=Make_Interval(f,Input_X)
    
    # print("The Initial interval", A)
    #print("The type of X_in=",type(A))
    counter=0
    B,step=Get_Estimator(f,A,counter)
    #r1=Get_Estimator1(f,A,i)
    print(" Root is in the range= ",A," and",B)
    
    step_E=Newton_method_Approximation_step(f,A,Ep,counter)
    print("The total steps in without initial estimator:= ",step_E)

    step_NM=Newton_method_Approximation_step(f,B,Ep,counter)
    print("The total steps in NM estimator:= ",step_NM+step)
    
    step_E_BI=Newton_method_Approximation_step_bi(f,A,B,Ep,counter)
    print("The total steps in with initial estimator with bi:= ",step_E_BI+step)
   


