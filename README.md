# shors-pq
Enhanced Version of Shors in ProjectQ

This is a modified version of the ProjectQ standard 'shor' example.
The modifications are primarily done to enhance the overall modularity and readability of the program.

The exercise has also identified a few issues that are described below.

Note that in this example, the simple (Slow) simulator is used.
We would like to add to the code only the bare minimum additional engine components that are essential to make the algorithm function correctly and return correct factors most of the time.


### Issues

#### The *high_level_gates()* method referenced in the documentation includes code that can never be reached.

There is a reference in the ProjectQ documentation to this method and a suggestion about returning True in order to make the algorithm execute more quickly. 
However, the code that follows that line could never have been reached.
This needs to be clarified with the ProjectQ team.
In the code snippet below the 'print' line indicates the place at which the subsequent code would never be reached, even in the unedited original source code in the ProjectQ documentation.

```
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
```

#### Results of factoring are not repeatable across executions

The randomization behavior of the ProjectQ simulator does not appear to be controllable.
Generation of the random base for exponentiation uses a random seed and generates the same base each time the algorithm is run from the beginning.
However, when the QFT algorithm is run, the results are unpredictable.
Is there a way to seed the Simulator?

#### Factoring results are not calculated effectively

In the output below, one can see that the first trial resulted in a period of 91.
The number 91 is nearly divisible by 17 and possibly could have been used to determine that
17 is a factor.  But the algorithm did not do this. It looks for perfect divisors that are r+1 and r-1.
The effect is that the program fails to factor most numbers above 256 or so.
Perhaps the algorithm could be modifed to test for 'close matches'? 

Additionally, the period determined by the QFT algorithm is often not evenly divisible, especially with larger numbers.
Is there a way to improve the quality of the QFT result?

```
Shor's Factoring Algorithm - ProjectQ
Enter the number to factor: 221
Factoring number = 221

... running Shors to factor number [ 221 ] with base=60 using num_qubits=8
  0100010101001100
  ... y = 0.197784423828125  fraction = 18 / 91  r = 91
  ==> Bad luck: Found 1 and 1 which are not the factors
  ... trying again ...
... running Shors to factor number [ 221 ] with base=79 using num_qubits=8
  0110110111011011
  ... y = 0.858245849609375  fraction = 109 / 127  r = 127
  ==> Factors found :-) : 13 * 17 = 221
*** Factoring attempts failed 1 times!
```
