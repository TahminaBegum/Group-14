from ariadne import *
import random


def Intialization():
    print("1. x*x-2")
    print("2. x*x*x-2")
    print("3. 6*x*x+4*x-50")
    print("4. x*x-200")

def Intialization_Estimator():
    print("1. Select the Random Boundary(one bound will select in the range (1 to 9))")
    print("2. Both bound initialize by yourself")
    print("3. You want to provide one value")
    print("If you wouldn't select any option,then option 1 will select")


def function_for_execution(function_number,x):
    if function_number==1:
        f= x*x-2
    elif function_number==2:
        f=x*x*x-2
    elif function_number==3:
        f=6*x*x+4*x-50
    else:
        f=x*x-200
    
    return f

def Create_Bound(F_bound,S_bound):
    if F_bound<S_bound:
        IX=UpperInterval({cast_exact(F_bound):cast_exact(S_bound)})
    else:
        IX=UpperInterval({cast_exact(S_bound):cast_exact(F_bound)})
    X=cast_singleton(IX)
    return X


def Create_SecondBound(One_bound):
    k=2
    ax_FDA=FloatDPApproximation(One_bound)
    bx_FDA=ax_FDA - k*(f(ax_FDA)/derivative(f,ax_FDA))
    ax_Float=float(str(ax_FDA))     #cast_exact{} need integer or float as paremeter so I convert it
    bx_Float=float(str(bx_FDA))
    #ax_Float=FloatDPValue(ax_FDA)
    #bx_Float=FloatDPValue(bx_FDA)
    if ax_Float<bx_Float:
        IX=UpperInterval({cast_exact(ax_Float):cast_exact(bx_Float)})
    else:
        IX=UpperInterval({cast_exact(bx_Float):cast_exact(ax_Float)})

    X=cast_singleton(IX)
    return X

def Get_Random_Number():
    #Select initial estimator
    #r=2.0
    r=random.randint(1,9)
    return Create_SecondBound(r)

def Shift(f,X):
    k=2
    X_value=X.value()
    shift_X=X_value - k*f(X_value)/derivative(f,X)
    return shift_X

def Get_Estimator(f,X):
 
    print("X=",X)
    fX=Contractor(X)
    
    if definitely(fX<0):
        print("both negative then the shift is not working properly")
        print("The value of negative FX",fX)
        X_s=Shift(f,X)
        print("the x_s",X_s)
    
    if (possibly(fX>=0) & definitely(fX<0)): #indetermine
        print("fX=",fX)
        print("first X return")
        return X

    X_value=X.value()
    new_X=X_value - f(X_value)/derivative(f,X)
    if inconsistent(new_X,X):       #no overlap then shift
        new_X=Shift(f,X)
        print("new_X inconsistance")
        return Get_Estimator(f,new_X)
    
    else:
        Function_Contractor=Contractor(new_X)
        if (possibly(Function_Contractor<0)):
            print("new_x possibly")
            return new_X
        else:
            new_X=Shift(f,X)
            print("new_x last")
            return Get_Estimator(f,new_X)


def Contractor(X):
    return f(X)




if __name__=='__main__':
    
    #initialization with Some functions
    Intialization()
    function_number=int(input("Enter your function(1 to 4):"))
  
    #Find our Selected function
    x = EffectiveScalarUnivariateFunction.identity()
    f=function_for_execution(function_number,x)
    print("my_function: ",f)
    Intialization_Estimator()
    Estimator_Choice=int(input("Enter your Choice(1 to 3):"))
    if Estimator_Choice==2:
        F_bound=int(input("Enter your First bound:"))
        S_bound=int(input("Enter your Second bound:"))
        X=Create_Bound(F_bound,S_bound)
    #print("The total bound is:=",X)
    
    elif Estimator_Choice==3:
        One_bound=int(input("Enter your bound:"))
        X=Create_SecondBound(One_bound)
    else:
        X=Get_Random_Number()
        print("The random number",X)
    root_range=Get_Estimator(f,X)

    print("Root is in the range=",root_range)
    
   


