#
# ---Created by Group 14--------
#----Designed and Implemantated By Tahmina Begum(i6185670)
#-----Add both of ur name and id as reviewing
#

import matplotlib.pyplot as plt
from ariadne import *
from statistics import mean
import numpy as np
import pandas as pd
import time


def print_function():
    # Print all the function in the screen

    print("1. x^2-2")
    print("2. x^3-2")
    print("3. f=pow((x-3),2)*(x+1)-2")
    print("4. f=(x-1)*pow((x+2),3)*pow((x-2),2)-1")


def intialization_method():
    # Print all the method in the screen

    print("1. General Newton Method")
    print("2. General Newton Method in Initial estimator by using sign")
    print("3. General Newton Method in Initial estimator by using contractor")
    print("4. Taylor Series Approx. in Initial estimator by using sign method")
    print("5. Taylor Series Approx. in Initial estimator by using contractor ")


def select_function(function_number, x):
    """
        Input:
        function_number -- choice of function assigned to number,
        x               -- function variable is x
        Output:
        f               -- function chosen with variable x
        """

    if function_number == 1:
        f = pow(x, 2)
    elif function_number == 2:
        f = pow(x, 3) - 2
    elif function_number == 3:
        f = pow((x - 3), 2) * (x + 1) - 2
    else:
        f = (x - 1) * pow((x + 2), 3) * pow((x - 2), 2) - 1

    return f


def get_label(function_number):
    # use to display label in the graph

    if function_number == 1:
        return "x*x-2"
    elif function_number == 2:
        return "x*x*x-2"
    elif function_number == 3:
        return "pow((x-3),2)*(x+1)-2"
    else:
        return "f=(x-1)*pow((x+2),3)*pow((x-2),2)-1"


def make_interval(bnd1, bnd2):
    """
        Make an interval [bnd1, bnd2] for the contractor method.
        Input:
        bnd1 and bnd2               -- two separate single values
        Output:
        x          -- single value with interval [bnd1, bnd2] or [bnd2, bnd1]
        """

    dpr = DoublePrecision()
    if bnd1 < bnd2:
        x = FloatDPBounds(Decimal(str(bnd1)), Decimal(str(bnd2)), dpr)
    else:
        x = FloatDPBounds(Decimal(str(bnd2)), Decimal(str(bnd1)), dpr)

    # print("interval: ",x)
    return x


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
        fx = f(x)
        dfx = derivative(f, x)
        if decide(dfx == 0):
            print("Zero derivative. No solution found in NM.For x=", x)
            return 0
        h = fx / dfx
        x = x - h
        if (decide(abs(h) <= Ep)):
            break
    if rootdisplay:
        print("Root in Approximation: ", x)
    return step

def get_estimator_contractor_nm(f, x, step, k1):
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
    fx = f(x)
    dfx = derivative(f, x)

    if decide(dfx == 0):
        print("Zero derivative. No solution found in IE contractor.For x=", x)
        return x, x, 0

    h = fx / dfx
    k_step = k1
    xp = x
    for j in range(1, 100):
        step = step + 1
        x_new = x - k_step * h
        fr = f(make_interval(xp, x_new))
        k_step = k_step * 2
        if not (definitely((fr) >= 0) | definitely((fr) <= 0)):
            return xp, x_new, step
        else:
            xp = x_new

    print("Limit need to Increass c1 x", x)
    return x, x, -1



def get_estimator_sign_taylor(f, x, step_initial, k):
    """
            Single point taken as initial input (x) and then evaluate it,(f(x)).
            By the evaluation, we try to find the sign of the function at the input point,
            The sign of input point is stored in the variable sign1.
            We try to find a point x_new, when evaluated (f(x_new) and sign stored in sign2) has a different sign compared to sign1.
            x_new is calculated by x_new=x-k*h with h=(-fd1 + sng(fd1) * sq) / fd2. Initially, k=2.
            If sign1 and sign2 are same then replace the value of x_new by x-k*h,
            where the value of k is doubled.

            The iteration will break if sign1 and sign2 are not same or max k is reached.
            In our implementation, we set the loop for k from 1 to 100.
            ------------------------------------------------------------------
            Input:
            f               -- function chosen with variable x
            x               -- single point taken from inital input
            step_initial    -- counter of newton steps performed, initially step=0
            k               ---k is a parameter, which helps to find the next bound quickly
            ----------------------------------------------------------------------
            Output:
            x               -- the initial input
            x_new           -- upperbound of interval with proof of root
            step            -- total number of iteration need to find perfect x_new
            """

    k_step = k
    deg = 2
    dfx = differential(f, x, deg)
    if decide(dfx[(2,)] == 0):
        print("Zero derivative. No solution found in IE(second derivative)sign method.For x=", x)
        return x, x, 0

    fd1 = FloatDPApproximation(dfx[(1,)])
    fd2 = 2 * FloatDPApproximation(dfx[(2,)])
    fx = FloatDPApproximation(dfx[(0,)])
    sign1 = decide(fx > 0) and 1 or -1
    singr = decide((fd1) > 0) and 1 or -1
    sq = sqrt(abs(pow(fd1, 2) - (2 * fd2 * fx)))
    h = (-fd1 + singr * sq) / fd2
    xp = x
    # print("dfx x sign kstep h sq:=",dfx,x,sign,k_step,h,sq)
    step = step_initial
    for j in range(1, 100):
        step = step + 1
        x_new = x + k_step * h
        k_step = k_step * 2  # make the k double in each iteration
        fx_new = f(x_new)
        sign2 = decide(fx_new > 0) and 1 or -1
        #print("x xp xnew", x, xp, x_new)
        if not (sign1 == sign2):
            return xp, x_new, step
        else:
            xp = x_new

    print("limit need to Increase ts2", x)
    return x, x, -1


def get_estimator_contrator_taylor(f, x, step_initial, k):

    """
            Contracting method is done by evalutating the sign of interval,
            the input x is the starting point of the newton step to x_new.
            x_new is calculated by x_new=x-k*h with h=(-fd1 + sng(fd1) * sq) / fd2.
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
            x               -- single point taken from inital input
            step            -- counter of newton steps performed, initially step=0
            k               -- k is a parameter, which helps to find the next bound quickly
            ----------------------------------------------------------------------
            Output:
            x_prev          -- lowerbound of the interval
            x_new           -- upperbound of interval with proof of root
            step            -- total number of iteration need to find perfect x_new
            """

    k_step = k
    deg = 2
    dfx = differential(f, x, deg)
    if decide(dfx[(2,)] == 0):
        print("Zero derivative. No solution found in IE(second derivative)contractor method.For x=", x)
        return x, x, 0

    fd1 = FloatDPApproximation(dfx[(1,)])
    fd2 = 2 * FloatDPApproximation(dfx[(2,)])
    fx = FloatDPApproximation(dfx[(0,)])
    sign = decide((fx / fd1) > 0) and 1 or -1
    singr = decide((fd1) > 0) and 1 or -1
    sq = sqrt(abs(pow(fd1, 2) - (2 * fd2 * fx)))
    h = (-fd1 + singr * sq) / fd2
    xp = x
    step = step_initial
    for j in range(1, 100):
        step = step + 1
        x_new = x + k_step * h
        k_step = k_step * 2  # make the k double in each iteration
        fr = f(make_interval(xp, x_new))
        # print("interval",make_interval(xp, x_new)#)
        #print("x xp xnew fr sign", x, xp, x_new, fr, sign)
        if not (definitely((fr) >= 0) | definitely((fr) <= 0)):
            return xp, x_new, step
        else:
            xp = x_new

    print("limit need to Increase tc2", x)
    return x, x, -1


def get_estimator_sign_nm(f, x, step, k):
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
    dfx = derivative(f, x)
    # print("the dfx",dfx)
    if decide(dfx == 0):
        print("Zero derivative. No solution found in IE sign method.For x=", x)
        return x, x, 0
    h = fx / dfx
    xp = x
    for j in range(1, 100):
        step = step + 1
        x_new = x - k_step * h
        k_step = k_step * 2  # make the k double in each iteration
        fx_new = f(x_new)
        sign2 = decide(fx_new > 0) and 1 or -1
        # print("x xp xnew sign1 sign2",x,xp,x_new,sign1,sign2)
        if not (sign1 == sign2):
            return xp, x_new, step
        else:
            xp = x_new
    # print("x xnew sign1 sign2",x,x_new,sign1,sign2)
    print("limit need to Increase x se1", x)
    return x, x, -1


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
        
        
        
        its for test1+avg value in here to find true root.......please write
        """

    generel_nm = []
    input_range = []
    initialestimator_nm = []
    initialestimator_nm_con = []
    second_nm = []
    second_nm_con = []
    all_method = []
    t_start = 0
    t_end = 0
    generel_nm_time = []
    initialestimator_nm_time = []
    initialestimator_nm_con_time = []
    second_nm_time = []
    second_nm_con_time = []
    bar_list = []
    bars_name = []
    bar_list_time = []
    bars_name_time = []

    # To get the value of Epsilon
    pr = DoublePrecision()
    Ep = FloatDP.eps(pr)
    # print("type ep",type(Ep),Ep)
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
        input_r = int(input("Input range of the approximation root X (-X,+X): "))
        k = float(input("The value of k:"))
        interval = 1

        for a in np.arange(-input_r, input_r + 1, 0.1):
            input_x = interval * a
            # input_range.append(input_x)
            x = FloatDPApproximation(input_x)
            counter = 0
            rootdisplay = 0

            if not (decide(f(x) == 0)):
                input_range.append(input_x)
                for j in range(len(all_method)):
                    if all_method[j] == 1:
                        #t_start = time.process_time()
                        step_general_nm = newton_method(f, x, Ep, counter, rootdisplay)
                        #t_end = time.process_time()
                        if step_general_nm == 0:
                            generel_nm_time.append(0)
                        else:
                            #generel_nm_time.append(t_end - t_start)
                        generel_nm.append(step_general_nm)
                for j in range(len(all_method)):
                    if all_method[j] == 2:
                        t_start = time.process_time()
                        bound1, bound2, step = get_estimator_sign_nm(f, x, counter, k)
                        t_end = time.process_time()
                        if not decide(bound1 == bound2):
                            bound1 = (bound1 + bound2) / 2
                            step_nm = newton_method(f, bound1, Ep, counter, rootdisplay)
                            initialestimator_nm_time.append(t_end - t_start)
                            initialestimator_nm.append(step)  # i change the step+stepnm
                        else:
                            initialestimator_nm_time.append(100)
                            initialestimator_nm.append(10000)
                for j in range(len(all_method)):
                    if all_method[j] == 3:
                        t_start = time.process_time()
                        bound1_con, bound2_con, step_con = get_estimator_contractor_nm(f, x, counter, k)
                        t_end = time.process_time()
                        if not decide(bound1_con == bound2_con):
                            bound1_con = (bound1_con + bound2_con) / 2
                            step_nm_con = newton_method(f, bound1_con, Ep, counter, rootdisplay)
                            initialestimator_nm_con_time.append(t_end - t_start)
                            initialestimator_nm_con.append(step_con)
                        else:
                            initialestimator_nm_con_time.append(100)
                            initialestimator_nm_con.append(10000)
                for j in range(len(all_method)):
                    if all_method[j] == 4:
                        t_start = time.process_time()
                        bound1_second, bound2_second, step_second = get_estimator_sign_taylor(f, x, counter, k)
                        t_end = time.process_time()
                        if not decide(bound1_second == bound2_second):
                            bound1_second = (bound1_second + bound2_second) / 2
                            step_ge_sign_second = newton_method(f, bound1_second, Ep, counter, rootdisplay)

                            second_nm_time.append(t_end - t_start)
                            second_nm.append(step_second)
                        else:
                            second_nm_time.append(100)
                            second_nm.append(10000)

                for j in range(len(all_method)):
                    if all_method[j] == 5:
                        t_start = time.process_time()
                        bound1_second_con, bound2_second_con, step_second_con = get_estimator_contrator_taylor(f, x,
                                                                                                              counter,
                                                                                                              k)
                        t_end = time.process_time()
                        if not decide(bound1_second_con == bound2_second_con):
                            bound1_second_con = (bound1_second_con + bound2_second_con) / 2
                            step_ge_second_con = newton_method(f, bound1_second_con, Ep, counter, rootdisplay)

                            second_nm_con_time.append(t_end - t_start)
                            second_nm_con.append(step_second_con)
                        else:
                            second_nm_con_time.append(100)
                            second_nm_con.append(10000)

        # titlelabels = "Function: {}  k={}".format(get_label(function_number), k)
        titlelabels = "Function: {}  k={}".format(f, k)
        plt.title(titlelabels)
        plt.ylabel("Number of Steps", fontsize=20)
        plt.xlabel("Approximation Input", fontsize=20)

        #if not (len(generel_nm) == 0):
            #plt.plot(input_range, generel_nm, label='General_NM')
            #bar_list.append(mean(generel_nm))
            #bars_name.append('General_NM')
        if not (len(initialestimator_nm) == 0):
            plt.plot(input_range, initialestimator_nm, label='NM_IE')
            bar_list.append(mean(initialestimator_nm))
            bars_name.append('NM_IE')
        if not (len(initialestimator_nm_con) == 0):
            plt.plot(input_range, initialestimator_nm_con, 'r:', label='NM_IE_CON')
            bar_list.append(mean(initialestimator_nm_con))
            bars_name.append('NM_IE_CON')
        if not (len(second_nm) == 0):
            plt.plot(input_range, second_nm, label='sign+second')
            bar_list.append(mean(second_nm))
            bars_name.append('sign+second')
        if not (len(second_nm_con) == 0):
            plt.plot(input_range, second_nm_con, 'c:', label='CON+second')
            bar_list.append(mean(second_nm_con))
            bars_name.append('CON+second')
        x1, x2, y1, y2 = plt.axis()
        plt.axis((x1, x2, 0, 10))
        # plt.axis((x1,x2,0,100))
        plt.legend()
        plt.show()


     titlelabels = "Function: {}  k={}".format(get_label(function_number), k)
     plt.title(titlelabels)
     plt.ylabel("Time(nanosecond)",fontsize=20)
     plt.xlabel("Approximation Input",fontsize=20)
     if not (len(generel_nm_time) == 0):
        plt.plot(input_range, generel_nm_time, label='General_NM')
        bar_list_time.append(mean(generel_nm_time))
     if not (len(initialestimator_nm_time) == 0):
        plt.plot(input_range, initialestimator_nm_time, label='NM_IE')
        bar_list_time.append(mean(initialestimator_nm_time))
     if not (len(initialestimator_nm_con_time) == 0):
        plt.plot(input_range, initialestimator_nm_con_time, label='NM_IE_CON')
        bar_list_time.append(mean(initialestimator_nm_con_time))
     if not (len(second_nm_time) == 0):
        plt.plot(input_range, second_nm_time, label='sign+second')
        bar_list_time.append(mean(second_nm_time))
     if not (len(second_nm_con_time) == 0):
        plt.plot(input_range, second_nm_con_time, label='CON+second')
        bar_list_time.append(mean(second_nm_con_time))
     x1,x2,y1,y2 = plt.axis()
     plt.axis((x1,x2,0,0.1))
     plt.legend()
     plt.show()

    # print(bar_list)
    # x_pos = np.arange(len(bar_list))
    # plt.bar(x_pos, bar_list, color=['red', 'yellow', 'green', 'blue', 'cyan'])
    # plt.xticks(x_pos, bars_name)
    # plt.show()

    # x_pos = np.arange(len(bar_list_time))
    # plt.bar(x_pos, bar_list_time, color=['red', 'yellow', 'green', 'blue', 'cyan'])
    # plt.xticks(x_pos, bars_name)
    # plt.show()

    else:

        input_x = float(input("Input a approximation to root: "))
        k = float(input("The value of k:"))
        x = FloatDPApproximation(input_x, pr)
        # print("type of x",type(f(x)),f(x))
        counter = 0
        rootdisplay = 1

        if not (decide(f(x) == 0)):
            for j in range(len(all_method)):
                if all_method[j] == 1:
                    t_start = time.process_time()
                    step_general_nm = newton_method(f, x, Ep, counter, rootdisplay)
                    print("The total steps for true root finding without initial estimator(", input_x, "):= ", step_general_nm)
                    t_end = time.process_time()
                    print("The total TIME for true root finding without initial estimator(", input_x, "):= ", t_end - t_start)

            for j in range(len(all_method)):
                if all_method[j] == 2:
                    t_start = time.process_time()
                    bound1, bound2, step = get_estimator_sign_nm(f, x, counter, k)
                    if not (decide(bound1 == bound2)):
                        step_nm = newton_method(f, bound2, Ep, counter, rootdisplay)
                        print("The total steps for true root finding with initial estimator sign method(", bound2, "):= ", step_nm + step)
                        t_end = time.process_time()
                        print("The total TIME for true root finding with initial estimator sign method(", bound2, "):= ", t_end - t_start)

            for j in range(len(all_method)):
                if all_method[j] == 3:
                    t_start = time.process_time()
                    bound1_con, bound2_con, step_con = get_estimator_contractor_nm(f, x, counter, k)
                    if not (decide(bound1_con == bound2_con)):
                        step_nm_con = newton_method(f, bound2_con, Ep, counter, rootdisplay)
                        t_end = time.process_time()
                        print("The total steps for true root finding with initial estimator contractor method(", bound2_con, "):= ",
                              step_nm_con + step_con)
                        print("The total TIME for true root finding with initial estimator contractor method(", bound2_con, "):= ",
                              t_end - t_start)

            for j in range(len(all_method)):
                if all_method[j] == 4:
                    t_start = time.process_time()
                    bound1_second, bound2_second, step_second = get_estimator_sign_taylor(f, x, counter, k)
                    if not (decide(bound1_second == bound2_second)):
                        step_ge_sign_second = newton_method(f, bound2_second, Ep, counter, rootdisplay)
                        t_end = time.process_time()
                        print("The total steps for true root finding with initial estimator sign method(second derivatives)(",
                              bound2_second, "):= ", step_ge_sign_second + step_second)
                        print("The total time for true root finding with initial estimator sign method(second derivatives)(",
                              bound2_second, "):= ", t_end - t_start)

            for j in range(len(all_method)):
                if all_method[j] == 5:
                    t_start = time.process_time()
                    bound1_second_con, bound2_second_con, step_second_con = get_estimator_contrator_taylor(f, x, counter,
                                                                                                          k)
                    if not (decide(bound1_second_con == bound2_second_con)):
                        step_ge_second_con = newton_method(f, bound2_second_con, Ep, counter, rootdisplay)
                        t_end = time.process_time()
                        print(
                            "The total steps for true root finding contractor with initial estimator contractor method(second derivatives)(",
                            bound2_second_con, "):= ",
                            step_ge_second_con + step_second_con)
                        print("The total time for true root finding with initial estimator contractor method(second derivatives)(",
                              bound2_second_con, "):= ", t_end - t_start)
