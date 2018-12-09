from ariadne import *


def print_function():

    #Print all the function in the screen

    print("1.  x*x-2")
    print("2. (x*x-x)*(x*x+x)*(x*x-x)")
    print("3. (x-9)*(x*x-1)*(x+1)+1")
    print("4. (x-2)*(x-4)*(x-1)")


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
        f = (x - 9) * (x * x - 1) * (x + 1) + 1
    else:
        f = (x - 2) * (x - 4) * (x + 1)

    return f


def set_starting_point(f, x_input):


    """
        Make an interval [bnd1, bnd2] from single input variable (x_input) of function f.

        First bound set bnd1 = x_input, bound bnd2 by performing a single Newton Step (k=2): bnd2 = bnd1 - k*(f(bnd1)/f'(bnd1)).
        Depending on direction of Newton step, the bounds are set in the interval accordingly.
        floatDPApproximation is used to get an average value of the constructed interval, which is set as output

        Input:
                f               -- function chosen with variable x
                x_input         -- initial approximation of root by user
        Output:
                x_end           -- single value, average of interval [bnd1, bnd2] or [bnd2, bnd1]


        """


    k = 2
    bnd1_fda = FloatDPApproximation(x_input)
    bnd2_fda = bnd1_fda - k * (f(bnd1_fda) / derivative(f, bnd1_fda))
    bnd1_float = float(str(bnd1_fda))
    bnd2_float = float(str(bnd2_fda))
    if bnd1_float < bnd2_float:
        ix = UpperInterval({cast_exact(bnd1_float): cast_exact(bnd2_float)})
    else:
        ix = UpperInterval({cast_exact(bnd2_float): cast_exact(bnd1_float)})

    x = cast_singleton(ix)
    x_end = FloatDPApproximation(x)
    return x_end



def newton_method(f, x, Ep, step):

    """

        Counting newton step's till the approximation is sufficiently close to true root from single point bnd1 (<= Epsilon)
        Input:
                f           -- function chosen with variable x
                x           -- input variable
                Ep          -- Epsilon (= .00001)
                step        -- counter of newton step's, initially step=0
        Output:
                step        -- return  the total number of newton step's are required to reach to closer to true root

        """




    while True:
        step = step + 1
        # print("bnd1:=",bnd1)
        h = f(x) / derivative(f,x)
        x = x - h
        if (decide(abs(h) <= Ep)):
            break
    # print("Root in Approximation: ",bnd1)
    return step





def newton_method_bidirectional(f, bnd1, bnd2, Ep, step):


    """

        Counting newton step's till the approximation is sufficiently close to root from both directions of the interval [bnd1,bnd2] (<= Epsilon)
        Input:
                f               -- function chosen with variable x
                bnd1 & bnd2     -- input variables

                Ep              -- Epsilon (= .00001)
                step            -- counter of newton step's, initially step=0
        Output:
                step            -- return  the total number of newton step's are required to reach to closer to true root


        """




    while True:
        step = step + 1

        # print("bnd1=",bnd1," and bnd2=",bnd2)

        h_bnd1 = f(bnd1) / derivative(f, bnd1)
        bnd1 = bnd1 - h_bnd1
        if (decide(abs(h_bnd1) <= Ep)):
            # print("Root in Approximation: ",bnd1)
            return step

        h_bnd2 = f(bnd2) / derivative(f, bnd2)
        bnd2 = bnd2 - h_bnd2
        if (decide(abs(h_bnd2) <= Ep)):
            # print("Root in Approximation: ",bnd2)
            return step





def initial_estimator(f, x,step):

    """

        Single point taken as initial input (x) and then evaluate it,(f(x)).
        By the evaluation, we try to find the sign of the function at the input point,
        The sign of input point is stored in the variable sign1.
        We try to find a point x_new, when evaluated (f(x_new) and sign stored in sign2) has a different sign compared to sign1.
        x_new is calculated by x_new=x-k*h with h=f(x)/f'(x). Initially, k=2. If sign1 and sign2 are same then replace the value of x_new by x-k*h,
        where the value of k is doubled.
        The iteration will break if sign1 and sign2 are not same or max k is reached.
        In our implementation, we set the loop for k from 1 to 50000.

        Input:
                f               -- function chosen with variable x
                x               -- single point taken from inital input
                step            -- counter of newton steps performed, initially step=0
        Output:
                x_new           -- upperbound of interval with proof of root
                step            -- total number of iteration need to find perfect x_new

        """





    fx = f(x)

    if decide(fx > 0):
        sign1 = 1
    else:
        sign1 = -1
    k = 1
    h = fx / derivative(f,x)

    for k in range(1, 50000):
        step = step + 1
        k = k * 2           #make the k double in each iteration
        x_new = x - k  * h
        fx_new = f(x_new)
        if decide(fx_new > 0):
            sign2 = 1
        else:
            sign2 = -1

        if not (sign1 == sign2):
            return x_new, step

    print("limit need to Increase")


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




    # initialization with Some functions
    print_function()
    function_number = int(input("Select a function(1 to 4):"))

    # Find our Selected function
    x = EffectiveScalarUnivariateFunction.identity()
    f = select_function(function_number, x)
    #print("function: ", f)

    # To get the value of Epsilon
    e = UpperInterval({cast_exact(.00001): 0})
    xx = cast_singleton(e)
    Ep = FloatDPApproximation(xx)

    input_x = float(input("Input a approximation to root: "))

    bound1 = set_starting_point(f, input_x)

    counter = 0
    bound2, step_initial_estimator = initial_estimator(f, bound1, counter)
    #print(" Root is in the range= ", bound1, " and", bound2)

    step_newton_method = newton_method(f, bound1, Ep, counter)
    print("The total steps in without initial estimator:= ", step_newton_method)

    step_newton_method_estimator = newton_method(f, bound2, Ep, counter)
    print("The total steps in Newton method estimator:= ",step_newton_method_estimator + step_initial_estimator)

    step_newton_method_bidirectional = newton_method_bidirectional(f, bound1, bound2, Ep, counter)
    print("The total steps in with initial estimator with Newton method bidirectional:= ",step_newton_method_bidirectional + step_initial_estimator)
