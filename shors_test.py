import sys
from projectq.ops import All, Measure, R, Rx, Ry, Rz, X, Y, Z, H, S, Sdag, T, Tdag, Swap, SqrtX, C, QFT, get_inverse
from projectq import MainEngine
from projectq.backends import Simulator, CommandPrinter, ResourceCounter
from projectq.meta import Control

from projectq.cengines import (MainEngine,
                               AutoReplacer,
                               LocalOptimizer,
                               TagRemover,
                               InstructionFilter,
                               DecompositionRuleSet)
from projectq.libs.math import (AddConstant,
                                AddConstantModN,
                                MultiplyByConstantModN)

import math
from math import gcd
from fractions import Fraction
import time
import random
import numpy as np
np.random.seed(0)
import copy

print("Shor's Factoring Algorithm - ProjectQ")

backend_id='simulator'                                      # Define backend

verbose = True

# Set the backend for execution
def set_execution_target(backend_id='simulator', num_runs=100, verbose=False, use_hardware=True):
    
    backend = Simulator()
    compiler_engines=None
    authentication_error_msg = "No supporting modules for {0} backend found.  Using the simulator instead."
      
    eng_in = MainEngine(backend=backend, engine_list=compiler_engines, verbose=verbose)
    return eng_in

 
############### Circuit Definition

# Execute Shor's Factoring Algorithm given a 'number' to factor,
# the 'base' of exponentiation, and the number of qubits required 'input_size'

def ShorsAlgorithm(eng_in, input_size, number, base, verbose=verbose):
    
    # Create count of qubits to use to represent the number to factor
    # NOTE: this should match the input_size, but is calculated here to be sure
    n = int(math.ceil(math.log(number, 2)))
    
    if verbose:
        print(f"... running Shors to factor number [ {number} ] with base={base} using num_qubits={n}")
    
    # Create an engine and allocate necessary qubits
    eng = copy.copy(eng_in)
    qureg = eng.allocate_qureg(n)

    # this will hold the 2n measurement results
    measurements = [0] * (2 * n)

    # create a single control qubit
    ctrl_qubit = eng.allocate_qubit()
  
    # the QFT algorithm 
    
    X | qureg[0]
    
    if verbose:
        print("  ", end="")
    
    # perform modular exponentiation 2*n times
    for k in range(2 * n):
    
        t = (2 * n - 1 - k)
        current_a = pow(base, t, number)
        
        # one iteration of 1-qubit QPE
        H | ctrl_qubit
        with Control(eng, ctrl_qubit):
            MultiplyByConstantModN(current_a, number) | qureg

        # perform inverse QFT --> Rotations conditioned on previous outcomes
        for i in range(k):
            if measurements[i]:
                R(-math.pi/(1 << (k - i))) | ctrl_qubit
        H | ctrl_qubit

        Measure | ctrl_qubit
        
        # execute the circuit and measure control qubit
        eng.flush()
        measurements[k] = int(ctrl_qubit)
        
        # reset the control qubit (flip if it is a 1)
        if measurements[k]:
            X | ctrl_qubit

        if verbose:
            print(f"{measurements[k]}", end="")
            sys.stdout.flush()
            
    if verbose:
        print("")
        
    # measure all the qubits in number register
    All(Measure) | qureg
        
    # turn the measured values into a number in [0,1) by summing their binary values
    ma = [(measurements[2 * n - 1 - i]*1. / (1 << (i + 1))) for i in range(2 * n)]
    y = sum(ma)

    # continued fraction expansion to get denominator (the period?)
    r = Fraction(y).limit_denominator(number - 1).denominator
    f = Fraction(y).limit_denominator(number - 1)
    
    if verbose:
        print(f"  ... y = {y}  fraction = {f.numerator} / {f.denominator}  r = {f.denominator}")

    # return the (potential) period
    return r


# Filter function, which defines the gate set for the first optimization
# (don't decompose QFTs and iQFTs to make cancellation easier)
def high_level_gates(eng, cmd):
    g = cmd.gate
    if g == QFT or get_inverse(g) == QFT or g == Swap:
        return True
    if isinstance(g, BasicMathGate):
        #return False
        return True
        print("***************** should never get here !")
        if isinstance(g, AddConstant):
            return True
        elif isinstance(g, AddConstantModN):
            return True
        return False
    return eng.next_engine.is_available(cmd)


# Choose a base at random < N / 2 without a common factor of N
def choose_random_base(N):

    # try up to 100 times to find a good base
    for guess in range(100):
        a = int(np.random.random() * (N / 2))        
        if gcd(a, N) == 1:
            return a
            
    print(f"Ooops, chose non relative prime {a}, gcd={gcd(a, N)}, giving up ...")
    return 0

# Determine factors from the period
def determine_factors(r, a, N):

    # try to determine the factors
    if r % 2 != 0:
        r *= 2
    
    apowrhalf = pow(a, r >> 1, N)
    
    f1 = gcd(apowrhalf + 1, N)
    f2 = gcd(apowrhalf - 1, N)
    
    # if this f1 and f2 are not the factors
    # and not both 1
    # and if multiplied together divide N evenly
    # --> then try multiplying them together and dividing into N to obtain the factors
    f12 = f1 * f2
    if ((not f12 == N)
       and f12 > 1
       and int(1. * N / (f12)) * f12 == N):
       
        f1, f2 = f12, int(N/(f12))
    
    return f1, f2

# Attempt to execute Shor's Algorithm to find factors, up to a max number of tries
# Returns number of failures, 0 if success
def attempt_factoring(blank_eng, input_size, number, verbose):

    max_tries = 5
    trials = 0
    failures = 0
    while trials < max_tries:
        trials += 1

        # choose a base at random
        base = choose_random_base(number)
        if base == 0: break

        #eng, qureg = ShorsAlgorithm(blank_eng, input_size, number, base)
        r = ShorsAlgorithm(blank_eng, input_size, number, base, verbose=verbose)

        # try to determine the factors
        f1, f2 = determine_factors(r, base, number)

        # Success! if these are the factors and both are greater than 1
        if (f1 * f2) == number and f1 > 1 and f2 > 1:
            if verbose:
                print(f"  ==> Factors found :-) : {f1} * {f2} = {number}")
            break

        else:
            failures += 1
            if verbose:
                print(f"  ==> Bad luck: Found {f1} and {f2} which are not the factors")
                print(f"  ... trying again ...")
    
    return failures

        
#################### Main      

# Execute program with default parameters
def run (num_shots=100, verbose=verbose):
    
    # Define the compiler engine using the backend information available
    blank_eng = set_execution_target(backend_id=backend_id, verbose=verbose, num_runs=num_shots)
    
    do_interactive_test(blank_eng=blank_eng, verbose=verbose)
    return;
   
# For interactive testing
def do_interactive_test(blank_eng, verbose):

    number = int(input('Enter the number to factor: '))
    print(f"Factoring number = {number}\n")
    
    input_size = int(math.ceil(math.log(number, 2)))
    
    # attempt to execute Shor's Algorithm to find factors
    failures = attempt_factoring(blank_eng, input_size, number, verbose)
    
    # Report how many factoring failures occurred
    if failures > 0:
        print(f"*** Factoring attempts failed {failures} times!")
            

# if main, execute method
if __name__ == '__main__': run() 
