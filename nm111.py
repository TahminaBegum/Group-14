#
#                                 Master Project
#
#    ---- Checking for Collisions: Validating Optimality/Feasibility ----
#
#                                   Group 14
#                          Supervisor: Pieter Collins
#
#                      University Maastricht - januari 2019
#
#

import matplotlib.pyplot as plt
from ariadne import *
from statistics import mean
import numpy as np
import pandas as pd
import time


def print_function():
    
    # Print all the function in the screen
    #a=np.random.rand(1,9)
    #print("a",a)

    print("1. x*x-2")
    print("2. x*x*x-2")
    print("3. (x-9)*(x*x-1)*(x+1)+6")
    print("4. (x-2)*(x-4)*(x-1)+7")


def intialization_method():
    
    # Print all the method in the screen
    
    print("1. General Newton Method (General_NM)")
    print("2. General Newton Method with Initial estimator by using sign (NM_IE_SIGN)")
    print("3. General Newton Method with Initial estimator by using contractor (NM_IE_CON)")
    print("4. General Newton Method with Initial estimator by using sign method and second derivatives (NM_IE_SIGN_SECOND-DERIVATIVE)")
    print("5. General Newton Method with Initial estimator by using contractor method and second derivatives (NM_IE_CON_SECOND-DERIVATIVE)")


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
        #f=x*x*x+5
    #f=x*x*x+3
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
        bnd1 and bnd2               -- two separate single values
        Output:
        x          -- single value with interval [bnd1, bnd2] or [bnd2, bnd1]
        """
    
    bnd1 = str(bnd1)
    bnd2 = str(bnd2)
    dpr=DoublePrecision()
    if bnd1 < bnd2:
        x=FloatDPBounds(Decimal(bnd1),Decimal(bnd2),dpr)
    else:
        x=FloatDPBounds(Decimal(bnd2),Decimal(bnd1),dpr)
    
    #print("interval: ",x)
    return x

def get_estimator_contractor(f, x, step, k1):
    
    """
        Contracting method is done by evalutating the sign of interval,
        the input x is the starting point of the newton step to x_new.
        The fr interval is constructed by [x, x_new], and the function is
        evaluated for all values in this interval. If the interval is strictly
        positive or strictly negative no root has been found, so the program will
        proceed to the next iteration from initial x with value of k doubled.
        The program will break if there is a sign difference found in the interval.
        The output is the contracted interval with proof of root. In our implementation,
        we set the loop for maximum 100 iterations.
        -------------------------------------------------------------------
        Input:
        f               -- function chosen with variable x
        x_init          -- single point taken from inital input
        step            -- counter of newton steps performed, initially step=0
        k               -- k is a parameter, which helps to find the next bound quickly
        ----------------------------------------------------------------------
        Output:
        x_prev          -- lowerbound of the interval
        x_new           -- upperbound of interval with proof of root
        step            -- total number of iteration need to find perfect x_new
        """
    fx=f(x)
    dfx=derivative(f,x)
    
    if decide(dfx==0):
        print("Zero derivative. No solution found in NM_IE_CON method.For x=",x)
        return x,x,0

    h = fx /dfx
    k_step = k1
    xp=x
    for j in range(1, 100):
        step = step + 1
        x_new = x - k_step * h
        #print("xp x_new=",xp,x_new)
        fr =f(make_interval(xp, x_new))
        k_step = k_step * 2
        if not (definitely((fr) >= 0) | definitely((fr) <= 0)):
            return xp, x_new, step
        else:
            xp=x_new

    print("Limit need to Increass")



def newton_method(f, x, Ep, step, rootdisplay):
    
    """
        The newton method is used for finding the exact root from the interval
        constructed by previous defined estimator program (contractor, sign).
        The newton method will iterate in this interval until it is sufficiently
        close to the true root from single bound of interval. This is achieved
        when f(x)/f'(x) is smaller or equal to epsilon. The step to achieve this
        are counted and outputed, to be added to the steps counted in the estimator program.
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
            print("Zero derivative. No solution found in General_NM method.For x=",x)
            return 0
        h = fx / dfx
        x = x - h
        if (decide(abs(h) <= Ep)):
            break
    if rootdisplay:
        print("Root in Approximation: ", x)
    return step


def get_estimator_sign_second(f, x, step_initial,k):
    
    fx = f(x)
    sign1 = decide(fx > 0) and 1 or -1
    k_step = k
    deg=2
    dfx=differential(f,x,deg)
    if decide(dfx[(2,)]==0):
        print("Zero derivative. No solution found in NM_IE_SIGN_SECOND-DERIVATIVE method.For x=",x)
        return x,x,0
    h = FloatDPApproximation(dfx[(1,)])/ (2*FloatDPApproximation(dfx[(2,)]))
    #h = fx/ (2*FloatDPApproximation(dfx[(2,)]))
    #print("x fx dfx h ",x,fx,dfx[(2,)], h)
    xp=x
    step = step_initial
    for j in range(1,10):
        #print("j",j)
        step = step + 1
        x_new = x - k_step * h
        k_step = k_step * 2         # make the k double in each iteration
        fx_new = f(x_new)
        sign2 = decide(fx_new > 0) and 1 or -1
        if not (sign1 == sign2):
            return xp,x_new, step
        else:
            xp=x_new

    #print("j sign1 sign2",j,sign1,sign2)
    h = fx/ (2*FloatDPApproximation(dfx[(2,)]))
    #print("x fx dfx h ",x,fx,dfx[(2,)], h)
    xp=x
    step = step_initial
    for j in range(1,10):
        #print("j",j)
        step = step + 1
        x_new = x - k_step * h
        k_step = k_step * 2         # make the k double in each iteration
        fx_new = f(x_new)
        sign2 = decide(fx_new > 0) and 1 or -1
        if not (sign1 == sign2):
            return xp,x_new, step
        else:
            xp=x_new

    print("limit need to Increase")



def get_estimator_sign_second_con(f,x, step_initial,k):
    
    fx = f(x)
    sign1 = decide(fx > 0) and 1 or -1
    k_step = k
    deg=2
    dfx=differential(f,x,deg)
    if decide(dfx[(2,)]==0):
        print("Zero derivative. No solution found in NM_IE_CON_SECOND-DERIVATIVE method.For x=",x)
        return x,x,0
    
    h = FloatDPApproximation(dfx[(1,)]) / (2*FloatDPApproximation(dfx[(2,)]))
    xp=x
    step = step_initial
    for j in range(1,10):
        step = step + 1
        x_new = x - k_step * h
        k_step = k_step * 2         # make the k double in each iteration
        fx_new = f(make_interval(xp, x_new))
        sign2 = decide(fx_new > 0) and 1 or -1
        if not (sign1 == sign2):
            return xp,x_new, step
        else:
            xp=x_new

    h = fx/ (2*FloatDPApproximation(dfx[(2,)]))
    xp=x
    step = step_initial
    for j in range(1,10):
        step = step + 1
        x_new = x - k_step * h
        k_step = k_step * 2         # make the k double in each iteration
        fx_new = f(make_interval(xp, x_new))
        sign2 = decide(fx_new > 0) and 1 or -1
        if not (sign1 == sign2):
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
    dfx=derivative(f, x)
    #print("the dfx",dfx)
    if decide(dfx==0):
        print("Zero derivative. No solution found in NM_IE_SIGN method.For x=",x)
        return x,x,0
    h = fx /dfx
    xp=x
    for j in range(1, 100):
        step = step + 1
        x_new = x - k_step * h
        k_step = k_step * 2         # make the k double in each iteration
        fx_new = f(x_new)
        sign2 = decide(fx_new > 0) and 1 or -1
        if not (sign1 == sign2):
            return xp,x_new, step
        else:
            xp=x_new

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
        function_number        -- User defines function
        total_method           -- Total methods used
        input_type_single      -- single input by user
        input_r                -- Interval range
        k                      -- value of k
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
        The program is a user interacting interface, so the user can input its own preferences.
        First the user is asked to choose a one dimensional function that the program will use to validate a root.
        Next the program will display the available methods for solving the rootfinding problem, and the user can
        choose how many it may use. After choosing a selection of the methods the program will ask for a single input
        or input range. A single input will initialize the program from this point, an interval will initialize for
        every integer value in the interval. When the user choose the interval option, the interval needs to be specified.
        Next, the value of k is asked to the user. The results are shown for the methods the user choose,
        number of iterations and root. Also for interval input, a plot will be made where the count of loop iterations is displayed.
        """
    
    
    generel_nm = []
    input_range = []
    initialestimator_nm = []
    initialestimator_nm_con = []
    second_nm = []
    second_nm_con = []
    all_method = []
    t_start=0
    t_end=0
    k_list=[]
    generel_nm_time = []
    initialestimator_nm_time = []
    initialestimator_nm_con_time = []
    second_nm_time = []
    second_nm_con_time = []
    mean_step_ge=[]
    mean_step_ie=[]
    mean_step_ie_con=[]
    mean_step_ie_second=[]
    mean_step_ie_second_con=[]
    mean_step_ge_time=[]
    mean_step_ie_time=[]
    mean_step_ie_con_time=[]
    mean_step_ie_second_time=[]
    mean_step_ie_second_con_time=[]
    
    #temp=FloatDPBounds(Decimal(1.2),DoublePrecision())
    #print("temp ",temp.upper()-temp.lower())
    # To get the value of Epsilon
    pr = DoublePrecision()
    Ep = FloatDP.eps(pr)
    #print("type ep",type(Ep),Ep)
    Ep = FloatDPApproximation(Ep)
    
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

    input_type_single = askbool("Are you want to use single input? [default: range input] ")

    if not (input_type_single):
        input_r= int(input("Input range of the approximation root X (-X,+X): "))
        #k = float(input("The value of k:"))
        kr=0
        for k1 in range(kr, kr+5):
            k=k1+0.7
            k_list.append(k)
            interval = 1
            for a in range(-input_r, input_r+1):
                input_x = interval * a
                input_range.append(input_x)
                x = FloatDPApproximation(input_x)
                counter = 0
                rootdisplay = 0
            
                if not (decide(f(x) == 0)):
                    for j in range(len(all_method)):
                        if all_method[j] == 1:
                            t_start=time.process_time()
                            step_general_nm = newton_method(f, x, Ep, counter, rootdisplay)
                            t_end=time.process_time()
                            if step_general_nm==0:
                                generel_nm_time.append(0)
                            else:
                                generel_nm_time.append(t_end-t_start)
                                generel_nm.append(step_general_nm)
                    for j in range(len(all_method)):
                        if all_method[j] == 2:
                            t_start=time.process_time()
                            bound1, bound2, step = get_estimator_sign(f, x, counter, k)
                            if not decide(bound1==bound2):
                                step_nm = newton_method(f, bound1, Ep, counter, rootdisplay)
                                t_end=time.process_time()
                                initialestimator_nm_time.append(t_end-t_start)
                                initialestimator_nm.append(step + step_nm)
                            else:
                                initialestimator_nm_time.append(0)
                                initialestimator_nm.append(0)
                    for j in range(len(all_method)):
                        if all_method[j] == 3:
                            t_start=time.process_time()
                            bound1_con, bound2_con, step_con = get_estimator_contractor(f, x, counter, k)
                            if not decide(bound1_con==bound2_con):
                                step_nm_con = newton_method(f, bound2_con, Ep, counter, rootdisplay)
                                t_end=time.process_time()
                                initialestimator_nm_con_time.append(t_end-t_start)
                                initialestimator_nm_con.append(step_con + step_nm_con)
                            else:
                                initialestimator_nm_con_time.append(0)
                                initialestimator_nm_con.append(0)
                    for j in range(len(all_method)):
                        if all_method[j] == 4:
                            t_start=time.process_time()
                            bound1_second, bound2_second, step_second= get_estimator_sign_second(f, x, counter, k)
                            if not decide(bound1_second==bound2_second):
                                step_ge_sign_second= newton_method(f,bound2_second, Ep, counter, rootdisplay)
                                t_end=time.process_time()
                                second_nm_time.append(t_end-t_start)
                                second_nm.append(step_second+step_ge_sign_second)
                            else:
                                second_nm_time.append(0)
                                second_nm.append(0)
                    for j in range(len(all_method)):
                        if all_method[j] == 5:
                            t_start=time.process_time()
                            bound1_second_con, bound2_second_con, step_second_con= get_estimator_sign_second_con(f, x, counter, k)
                            if not decide(bound1_second_con==bound2_second_con):
                                step_ge_second_con = newton_method(f,bound2_second_con, Ep, counter, rootdisplay)
                                t_end=time.process_time()
                                second_nm_con_time.append(t_end-t_start)
                                second_nm_con.append(step_second_con+step_ge_second_con)
                            else:
                                second_nm_con_time.append(0)
                                second_nm_con.append(0)


        
            mean_step_ge.append(np.array(generel_nm)[np.nonzero(np.array(generel_nm))].mean())
            mean_step_ie.append(np.array(initialestimator_nm)[np.nonzero(np.array(initialestimator_nm))].mean())
            mean_step_ie_con.append(np.array(initialestimator_nm_con)[np.nonzero(np.array(initialestimator_nm_con))].mean())
            mean_step_ie_second.append(np.array(second_nm)[np.nonzero(np.array(second_nm))].mean())
            mean_step_ie_second_con.append(np.array(second_nm_con)[np.nonzero(np.array(second_nm_con))].mean())
            #print("time:=",np.array(generel_nm_time)[np.nonzero(np.array(generel_nm_time))].mean())
            mean_step_ge_time.append(np.array(generel_nm_time)[np.nonzero(np.array(generel_nm_time))].mean())
            mean_step_ie_time.append(np.array(initialestimator_nm_time)[np.nonzero(np.array(initialestimator_nm_time))].mean())
            mean_step_ie_con_time.append(np.array(initialestimator_nm_con_time)[np.nonzero(np.array(initialestimator_nm_con_time))].mean())
            mean_step_ie_second_time.append(np.array(second_nm_time)[np.nonzero(np.array(second_nm_time))].mean())
            mean_step_ie_second_con_time.append(np.array(second_nm_con_time)[np.nonzero(np.array(second_nm_con_time))].mean())
           
           #print(k_list)
        N=5
        ind=np.arange(N)
        width=0.12
        fig=plt.figure()
        r1=plt.bar(ind,mean_step_ge,width,color='r')
        r2=plt.bar(ind+width,mean_step_ie,width,color='g')
        r3=plt.bar(ind+width*2,mean_step_ie_con,width,color='b')
        r4=plt.bar(ind+width*3,mean_step_ie_second,width,color='y')
        r5=plt.bar(ind+width*4,mean_step_ie_second_con,width,color='m')
        plt.xticks([r+width*2 for r in range(N)],k_list)
        titlelabels = "Function: {}".format(get_label(function_number))
        plt.title(titlelabels)
        plt.ylabel("Number of Average Steps")
        plt.legend((r1[0],r2[0],r3[0],r4[0],r5[0]),('General_NM','NM_IE_SIGN','NM_IE_CON','NM_IE_SIGN_SECOND-DERIVATIVE','NM_IE_CON_SECOND-DERIVATIVE'))
        plt.show()

#print("time mean:=",mean_step_ge_time)
        r1=plt.bar(ind,mean_step_ge_time,width,color='r')
        r2=plt.bar(ind+width,mean_step_ie_time,width,color='g')
        r3=plt.bar(ind+width*2,mean_step_ie_con_time,width,color='b')
        r4=plt.bar(ind+width*3,mean_step_ie_second_time,width,color='y')
        r5=plt.bar(ind+width*4,mean_step_ie_second_con_time,width,color='m')
        plt.xticks([r+width*2 for r in range(N)],k_list)
        titlelabels = "Function: {}".format(get_label(function_number))
        plt.title(titlelabels)
        plt.ylabel("Number of Average time(nanoseconds)")
        plt.legend((r1[0],r2[0],r3[0],r4[0],r5[0]),('General_NM','NM_IE_SIGN','NM_IE_CON','NM_IE_SIGN_SECOND-DERIVATIVE','NM_IE_CON_SECOND-DERIVATIVE'))
        plt.show()

    else:
        
        input_x = float(input("Input a approximation to root: "))
        k = float(input("The value of k:"))
        x = FloatDPApproximation(input_x)
        counter = 0
        rootdisplay = 1
        
        if not (decide(f(x) == 0)):
            for j in range(len(all_method)):
                if all_method[j] == 1:
                    t_start=time.process_time()
                    step_general_nm = newton_method(f, x, Ep, counter, rootdisplay)
                    print("The total steps in General_NM method(", input_x, "):= ", step_general_nm)
                    t_end=time.process_time()
                    print("The total TIME in General_NM method(", input_x, "):= ",t_end-t_start)
        
            for j in range(len(all_method)):
                if all_method[j] == 2:
                    t_start=time.process_time()
                    bound1, bound2, step = get_estimator_sign(f, x, counter, k)
                    if not(decide(bound1==bound2)):
                        step_nm = newton_method(f, bound2, Ep, counter, rootdisplay)
                        print("The total steps in NM_IE_SIGN method(", bound2, "):= ", step_nm + step)
                        t_end=time.process_time()
                        print("The total TIME in NM_IE_SIGN method(", bound2, "):= ",t_end-t_start)
            
            for j in range(len(all_method)):
                if all_method[j] == 3:
                    t_start=time.process_time()
                    bound1_con, bound2_con, step_con = get_estimator_contractor(f, x, counter, k)
                    if not(decide(bound1_con==bound2_con)):
                        step_nm_con = newton_method(f, bound2_con, Ep, counter, rootdisplay)
                        t_end=time.process_time()
                        print("The total steps in NM_IE_CON method(", bound2_con, "):= ",step_nm_con + step_con)
                        print("The total TIME in NM_IE_CON method(", bound2_con, "):= ",t_end-t_start)
                            
            for j in range(len(all_method)):
                if all_method[j] == 4:
                    t_start=time.process_time()
                    bound1_second, bound2_second, step_second= get_estimator_sign_second(f, x, counter, k)
                    if not(decide(bound1_second==bound2_second)):
                        step_ge_sign_second= newton_method(f,bound2_second, Ep, counter, rootdisplay)
                        t_end=time.process_time()
                        print("The total steps in NM_IE_SIGN_SECOND-DERIVATIVE(", bound2_second, "):= ",step_ge_sign_second + step_second)
                        print("The total time in NM_IE_SIGN_SECOND-DERIVATIVE(", bound2_second, "):= ",t_end-t_start)
            
            for j in range(len(all_method)):
                if all_method[j] == 5:
                    t_start=time.process_time()
                    bound1_second_con, bound2_second_con, step_second_con= get_estimator_sign_second_con(f, x, counter, k)
                    if not(decide(bound1_second_con==bound2_second_con)):
                        step_ge_second_con = newton_method(f,bound2_second_con, Ep, counter, rootdisplay)
                        t_end=time.process_time()
                        print("The total steps in NM_IE_CON_SECOND-DERIVATIVE(",bound2_second_con, "):= ",
                          step_ge_second_con + step_second_con)
                        print("The total time in NM_IE_CON_SECOND-DERIVATIVE(", bound2_second_con, "):= ",t_end-t_start)
