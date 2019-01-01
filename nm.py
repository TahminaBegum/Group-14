import matplotlib.pyplot as plt
from ariadne import *


def print_function():
    
    # Print all the function in the screen
    
    print("1. x*x-2")
    print("2. x*x*x-2")
    print("3. (x-9)*(x*x-1)*(x+1)+6")
    print("4. (x-2)*(x-4)*(x-1)+7")


def intialization_method():
    
    # Print all the method in the screen
    
    print("1. General Newton Method")
    print("2. General Newton Method with Initial estimator by using sign")
    print("3. General Newton Method with Initial estimator by using contractor")
    print("4. General Newton Method with Initial estimator by using sign in bidirection way")
    print("5. General Newton Method with Initial estimator by using contractor in bidirection way")


def select_function(function_number, x):
    
    
    """
        Input:
        function_number -- choice of function assigned to number,
        x               -- function variable is x
        Output:
        f               -- function chosen with variable x
        """
    
    
    
    if function_number == 1:
        f = x * x - 2
    elif function_number == 2:
        f = x * x * x - 2
    elif function_number == 3:
        f = (x - 9) * (x * x - 1) * (x + 1) + 6
    else:
        f = (x - 2) * (x - 4) * (x - 1) + 7
    
    return f




def get_label(function_number):
    
    # use to display label in the graph
    
    if function_number == 1:
        return "x*x-2"
    elif function_number == 2:
        return "x*x*x-2"
    elif function_number == 3:
        return "(x-9)*(x*x-1)*(x+1)+6"
    else:
        return "(x-2)*(x-4)*(x-1)+7"




def make_interval(bnd1, bnd2):
    
    """
        Make an interval [bnd1, bnd2] for the contractor method.
        Input:
        bnd1 and bnd2               -- two saperate single values
        
        Output:
        x          -- single value with interval [bnd1, bnd2] or [bnd2, bnd1]
        """
    
    bnd1 = float(str(bnd1))
    bnd2 = float(str(bnd2))
    if bnd1 < bnd2:
        ix = UpperInterval({cast_exact(bnd1): cast_exact(bnd2)})
    else:
        ix = UpperInterval({cast_exact(bnd2): cast_exact(bnd1)})
    
    x = cast_singleton(ix)
    return x




def get_estimator_contractor(f, x, step, k1):
    
    """
        Single point taken as initial input (x) and then evaluate it,(f(x)).
        By the evaluation, we try to find the sign of the function at the input point,
        The sign of input point is stored in the variable sign1.
        We try to find a point x_new, when evaluated (f(x_new) and sign stored in sign2) has a different sign compared to sign1.
        x_new is calculated by x_new=x-k*h with h=f(x)/f'(x). Initially, k is any integer. If sign1 and sign2 are same then replace the value of x_new by x-k*h,
        where the value of k is doubled.
        The iteration will break if sign1 and sign2 are not same or max k is reached.
        In our implementation, we set the loop for k from 1 to 50000.
        -------------------------------------------------------------------
        Input:
        f               -- function chosen with variable x
        x               -- single point taken from inital input
        step            -- counter of newton steps performed, initially step=0
        k               ---k is a parameter, which helps to find the next bound quickly
        ----------------------------------------------------------------------
        Output:
        x               -- initial input
        x_new           -- upperbound of interval with proof of root
        step            -- total number of iteration need to find perfect x_new
        
        """
    
    
    h = f(x) / derivative(f,x)
    k_step = k1
    for j in range(1, 500):
        step = step + 1
        x_new = x - k_step * h
        fr =f(make_interval(x, x_new))
        k_step = k_step * 2
        if not (definitely((fr) >= 0) | definitely((fr) <= 0)):
            return x, x_new, step

    print("Limit need to Increass")




def newton_method(f, x, Ep, step, rootdisplay):
    
    """
        Counting newton step's till the approximation is sufficiently close to true root from single point bnd1 (<= Epsilon)
        Input:
        f           -- function chosen with variable x
        x           -- input variable
        Ep          -- Epsilon (= .00001)
        step        -- counter of newton step's, initially step=0
        rootdisplay -- if we want to see the true root
        Output:
        step        -- return  the total number of newton step's are required to reach to closer to true root
        """
    
    while True:
        step = step + 1
        fx=f(x)
        dfx= derivative(f, x)
        if decide(dfx==0):
            print("Zero derivative. No solution found.")
            return 0
        h = fx / dfx
        x = x - h
        if (decide(abs(h) <= Ep)):
            break
    if rootdisplay:
        print("Root in Approximation: ", x)
    return step




def newton_method_bi(f, bnd1, bnd2, Ep, step, rootdisplay):
    
    """
        Counting newton step's till the approximation is sufficiently close to root from both directions of the interval [bnd1,bnd2] (<= Epsilon)
        
        ----------------------------------------------------------------
        
        Input:
        f               -- function chosen with variable x
        bnd1 & bnd2     -- input variables
        Ep              -- Epsilon (= .00001)
        step            -- counter of newton step's, initially step=0
        rootdisplay     -- root display in terminal
        
        ----------------------------------------------------------------
        Output:
        step            -- return  the total number of newton step's are required to reach to closer to true root
        """
    
    while True:
        step = step + 1
        h_bnd1 = f(bnd1) / derivative(f, bnd1)
        bnd1 = bnd1 - h_bnd1
        if (decide(abs(h_bnd1) <= Ep)):
            if rootdisplay:
                print("Root in Approximation: ", bnd1)
            return step
        
        h_bnd2 = f(bnd2) / derivative(f, bnd2)
        bnd2 = bnd2 - h_bnd2
        if (decide(abs(h_bnd2) <= Ep)):
            if rootdisplay:
                print("Root in Approximation: ", bnd2)
            return step


def get_estimator_sign_second(f, x, step,k):
    
    """
      underconstruction
        """
    
    fx = f(x)
    sign1 = decide(fx > 0) and 1 or -1
    k_step = k
    deg=2
    dfx=differential(f,x,deg)
    print(dfx,type(dfx));
    print("f''(x):", dfx[(2,)])
    #h = fx / derivative(f, x)
    h = fx / FloatDPApproximation(dfx[(2,)])
    print("the value of dfx in second order = ",dfx)
    xp=x
    for j in range(1, 100):
        step = step + 1
        x_new = x - k_step * h
        k_step = k_step * 2         # make the k double in each iteration
        print("The k step = ",k_step,"j=",j)
        print("xp h=",xp,h)
        fx_new = f(x_new)
        sign2 = decide(fx_new > 0) and 1 or -1
        print("sign1 = ",sign1,"sign2 = ",sign2,"x new = ",x_new)
        if not (sign1 == sign2):
            x_new=(xp+x_new)/2
            return xp,x_new, step
        else:
            xp=x_new
    print("limit need to Increase")




def get_estimator_sign(f, x, step,k):
    
    """
        Single point taken as initial input (x) and then evaluate it,(f(x)).
        By the evaluation, we try to find the sign of the function at the input point,
        The sign of input point is stored in the variable sign1.
        We try to find a point x_new, when evaluated (f(x_new) and sign stored in sign2) has a different sign compared to sign1.
        x_new is calculated by x_new=x-k*h with h=f(x)/f'(x). Initially, k=2. If sign1 and sign2 are same then replace the value of x_new by x-k*h,
        where the value of k is doubled.
        The iteration will break if sign1 and sign2 are not same or max k is reached.
        In our implementation, we set the loop for k from 1 to 50000.
        
        ------------------------------------------------------------------
        
        Input:
        f               -- function chosen with variable x
        x               -- single point taken from inital input
        step            -- counter of newton steps performed, initially step=0
        k               ---k is a parameter, which helps to find the next bound quickly
        
        ----------------------------------------------------------------------
        
        Output:
        x               -- the initial input
        x_new           -- upperbound of interval with proof of root
        step            -- total number of iteration need to find perfect x_new
        """
    
    fx = f(x)
    sign1 = decide(fx > 0) and 1 or -1
    k_step = k
    h = fx / derivative(f, x)
    
    for j in range(1, 100):
        step = step + 1
        x_new = x - k_step * h
        k_step = k_step * 2         # make the k double in each iteration
        print("The k step = ",k_step,"j=",j)
        print("x h=",x,h)
        fx_new = f(x_new)
        sign2 = decide(fx_new > 0) and 1 or -1
        print("sign1 = ",sign1,"sign2 = ",sign2,"x new = ",x_new)
        if not (sign1 == sign2):
            return x,x_new, step


    print("limit need to Increase")


def askbool(message):
    
    # check bool type of question
    
    a = input(message)
    if a is not "":
        if a.lower() in ["y", "yes", "t", "tr", "true", "1"]:
            return True
    return False


if __name__ == '__main__':
    
    """
        Input:
        input_x                --The user define input x(can int or float)
        --------------------------------------------------------------------------------------------------------
        Ep                     --Epsilon is a smallest positive integer (type floatDPApproximation)
        f                      --User selected function
        x                      --independent variable
        bound1                 --It's a first bound for the initial interval, which is claculated from Input_x
        bound2                 --It's the second bound, which is the return from the initial_estimator function
        counter                --it's a variable, it use for the counting the iteration.
        -----------------------------------------------------------------------------------------------------------
        step_initial_estimator               -- The total steps need to find initial estimator
        step_newton_method                   -- The total steps with general Newton method
        step_newton_method_estimator         -- The total steps with Newton method with using initial estimator as pre-processor
        step_newton_method_bidirectional     -- The total steps with Newton method bidirectional with using initial estimator as pre-processor
        ---------------------------------------------------------------------------------------------------------------
        Output:
        Display the The total steps to find the solution of function with general Newton method.
        Display the The total steps to find the solution of function by using Newton method and initial estimator
        Display the The total steps to find the solution of function by using Newton method bidirectional and initial estimator
        """
    
    generel_nm = []
    input_range = []
    initialestimator_nm = []
    initialestimator_nm_con = []
    bi_nm = []
    bi_nm_con = []
    all_method = []
    
    # temp=(1/3)
    # print("temp ",temp*2.999999)
    # To get the value of Epsilon
    e = UpperInterval({cast_exact(.00000001): 0})
    Ep = FloatDPApproximation(cast_singleton(e))
    
    print("")
    print("Welcome...")
    print("Testing Section active")
    print("")
    
    # initialization with Some functions
    print_function()
    function_number = int(input("Select an equation (1 to 4):"))
    
    # Find our Selected function
    x = EffectiveScalarUnivariateFunction.identity()
    f = select_function(function_number, x)
    print("Function: ", f)
    
    intialization_method()
    
    total_method = int(input("How many method you want to use (1 to 5)?:"))
    
    for m in range(total_method):
        all_method.append(int(input("Method:")))

    input_type_single = askbool("Are you want to single input? [default: range input] ")

    if not (input_type_single):
        input_r= int(input("Input range of the approximation root X (-X,+X): "))
        k = float(input("The value of k:"))
        interval = 2
        for a in range(-input_r, input_r):
            input_x = interval * a
            input_range.append(input_x)
            x = FloatDPApproximation(input_x)
            counter = 0
            rootdisplay = 0
            
            if not (decide(f(x) == 0)):
                bound1, bound2, step = get_estimator_sign(f, x, counter, k)
                bound1_con, bound2_con, step_con = get_estimator_contractor(f, x, counter, k)
                for j in range(len(all_method)):
                    if all_method[j] == 1:
                        step_general_nm = newton_method(f, x, Ep, counter, rootdisplay)
                        generel_nm.append(step_general_nm)
                for j in range(len(all_method)):
                    if all_method[j] == 2:
                        step_nm = newton_method(f, bound2, Ep, counter, rootdisplay)
                        initialestimator_nm.append(step + step_nm)
                for j in range(len(all_method)):
                    if all_method[j] == 3:
                        step_nm_con = newton_method(f, bound2_con, Ep, counter, rootdisplay)
                        initialestimator_nm_con.append(step_con + step_nm_con)
                for j in range(len(all_method)):
                    if all_method[j] == 4:
                        step_ge_bi = newton_method_bi(f, bound1, bound2, Ep, counter, rootdisplay)
                        bi_nm.append(step+step_ge_bi)
                for j in range(len(all_method)):
                    if all_method[j] == 5:
                        step_ge_bi_con = newton_method_bi(f, bound1_con, bound2_con, Ep, counter, rootdisplay)
                        bi_nm_con.append(step_con+step_ge_bi_con)

        titlelabels = "Function: {}  k={}".format(get_label(function_number), k)
        plt.title(titlelabels)
        plt.ylabel("Number of Steps")
        plt.xlabel("Input")
        if not (len(initialestimator_nm) == 0):
            plt.plot(input_range, initialestimator_nm, label='NM_IE')
        if not (len(generel_nm) == 0):
            plt.plot(input_range, generel_nm, label='General_NM')
        if not (len(initialestimator_nm_con) == 0):
            plt.plot(input_range, initialestimator_nm_con, label='NM_IE_CON')
        if not (len(bi_nm) == 0):
            plt.plot(input_range, bi_nm, label='BI')
        if not (len(bi_nm_con) == 0):
            plt.plot(input_range, bi_nm_con, label='BI_CON')
        plt.legend()
        plt.show()

    else:
        
        input_x = float(input("Input a approximation to root: "))
        k = float(input("The value of k:"))
        x = FloatDPApproximation(input_x)
        counter = 0
        rootdisplay = 1
        if not (decide(f(x) == 0)):
            bound1, bound2, step = get_estimator_sign(f, x, counter, k)
            bound1_con, bound2_con, step_con = get_estimator_contractor(f, x, counter, k)
            
            for j in range(len(all_method)):
                if all_method[j] == 1:
                    step_general_nm = newton_method(f, x, Ep, counter, rootdisplay)
                    print("The total steps in without initial estimator(", input_x, "):= ", step_general_nm)
            for j in range(len(all_method)):
                if all_method[j] == 2:
                    step_nm = newton_method(f, bound2, Ep, counter, rootdisplay)
                    print("The total steps in with initial estimator sign method(", bound2, "):= ", step_nm + step)
            for j in range(len(all_method)):
                if all_method[j] == 3:
                    step_nm_con = newton_method(f, bound2_con, Ep, counter, rootdisplay)
                    print("The total steps in with initial estimator contractor method(", bound2_con, "):= ",
                          step_nm_con + step_con)
            for j in range(len(all_method)):
                if all_method[j] == 4:
                    step_ge_bi = newton_method_bi(f, bound1, bound2, Ep, counter, rootdisplay)
                    print("The total steps in bi  with initial estimator sign method(", bound1, bound2, "):= ",
                              step_ge_bi + step)
            for j in range(len(all_method)):
                if all_method[j] == 5:
                    step_ge_bi_con = newton_method_bi(f, bound1_con, bound2_con, Ep, counter, rootdisplay)
                    print("The total steps in bi  with initial estimator contractor method(",bound1_con,bound2_con, "):= ",
                          step_ge_bi_con + step_con)
