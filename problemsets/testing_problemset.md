## `pytest` problem set

Generally, you will write your test functions in the same script file as the functions being tested and run `pytest` on your scripts to run the tests.

1. Using the following simple function:

    ```python
    def gc_content(seq):
        valid = {'A', 'C', 'G', 'T'}
        if not set(seq).issubset(valid):
            raise ValueError("Invalid characters in sequence")
        if len(seq) == 0:
            return 0
        return (seq.count('G') + seq.count('C')) / len(seq)
    ```
    Write `pytest` unit tests for this function that:  
    - confirm GC content of `"GCGC"` is `1.0`.  
    - confirm GC content of `"ATAT"` is `0.0`.  
    - confirm GC content of `"ATGC"` gives `0.5`.  
    - confirm that empty string returns `0`.  
    - confirm that `"ATGXB"` raises a `ValueError`. 

2. Modify your script to add tests of the `gc_content` function for the following inputs (below). Re-run `pytest` on your updated script. Do these new tests pass or fail? Confirm your modified `gc_content` function can:

    - Calculate the GC content of `"ATGNNNTAGC"` as `0.3`.
    - Calculate the GC content of the lower-cased sequence `"gattacaa"` as `0.25`. 

    Are the inputs above reasonable inputs? Modify the `gc_content` function so that `pytest` passes all tests.

3. Write a function to reverse complement a DNA sequence, then write unit tests for it that:

    - Confirm an all-lower-case input sequence is correctly reverse-complemented. 

    - Confirm an all-upper-case input sequence is correctly reverse-complemented.

    - Confirm a mixed-case input sequence is correctly reverse-complemented.

    - Confirm an input sequence containing non-`ATCGN` characters triggers an exception.

4. Write a function called `isnumeric` that determines whether an input value is a numeric type and outputs a boolean value (`True`/`False`). *Hint*: use a `try`-`except` block and the `float` function. Write unit tests for it that test the conditions below. If a test fails, modify your `isnumeric` function to meet that condition:

    - Confirm the number `1.3` returns `True`.
    - Confirm the string `'1.3'` returns `True`.
    - Confirm the string `'1e-6'` returns `True`.
    - Confirm the string `'-0.0001'` returns `True`.
    - Confirm the string `'not-a-number'` returns `False`.

5. Copy your `isnumeric` function (and only that function) from Problem 4 into a new script file and save it. Prompt an AI (such as Github Copilot, ChatGPT, etc.) to write unit tests for you, putting the test functions in a new `test_` script. Run the tests with `pytest`.

    - Does the AI write tests for an input value of `None` ? If not, ask it to. Re-run the tests.
    - Does the AI write tests for complex numeric inputs (e.g., `1+5j`)? If not, ask it to. Re-run the tests. Do its tests pass? Why?
    - Modify your `isnumeric` function to use `complex` instead of `float`. Re-run the tests. Do its tests pass now? Why? 
    - What do these tests suggest about how smart/rigorous AI-generated tests are?

     

### Challenge Questions

1. Using AI, write unit tests for the DNA sequence class you wrote for the Python 11 problem set. Write your test functions in a separate `test_` script that imports your sequence class and performs tests on its methods.
    - Look carefully at all the unit test inputs and expected outputs. Are they correct?
    - Now, we will have AI write more tests for a new method you will add to your DNA sequence class. 
        1. Copy the following method stub (which contains the `def` line and a doc string, but no executable code block) into your code. 

          ```python
          def six_frames_nt(self):
              """Outputs all six reading frames (3 forward, 3 reverse complement) as nucleotide seqences. 
              Arguments: None
              Returns: list of six strings, each string is a reading frame (forward 0, 1, 2; reverse 0, 1, 2).
          
              Examples:
              >>> sixframes_nt = Sequence("id1", "ATGCGTTAG").six_frames_nt()
              """
          ```

        2. Ask AI to write your tests. 
        3. Look carefully at the tests it produced. Are they correct (correct frames, length divisible by three, reverse-strand sequence is correctly reverse complemented, etc.)?
        4. Next, write the code to extract the six reading frames as a list of sequence strings (as described in the stub above). Test your code by running `pytest` on your test class methods.

2. Lastly, do the same for a new method called `six_frames_aa` that outputs all six reading frames as a list of six amino acid sequences: 
  1. Write a method stub (include a description of what the method does, input arguments, returned values, and examples).
  2. Ask AI to write your tests.
  3. Check the correctness of the tests!
  4. *Then* write the method code. Run `pytest` on your test script to test your class methods. 