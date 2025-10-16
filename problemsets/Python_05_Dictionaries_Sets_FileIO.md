Python 5 - Dictionaries, Sets, File IO
===================

Dictionaries
============

1.  Write a script in which you construct a dictionary of your favorite things.

> Some of my favorites:
>
> | Type | Favorite                                 |
> | ---- | ---------------------------------------- |
> | book | [Jitterbug Perfume](https://en.wikipedia.org/wiki/Jitterbug_Perfume) |
> | song | [Tom Petty - I Won't Back Down](https://www.youtube.com/watch?v=nvlTJrNJ5lA) |
> | tree | [Cedar](https://sciencing.com/cedar-trees-5432718.html) |


2. Print out your favorite book.
```python
print(fav_dict['book'])
```

3. Print out your favorite book but use a variable in the key.
``` python
fav_thing = 'book'
print(fav_dict[fav_thing])
```

4. Now print your favorite tree.

5. Add your favorite 'organism' to the dictionary. Make organism the new key of `fav_thing`
```python
fav_thing = 'organism'
print(fav_dict[fav_thing])
```

6. Use a `for` loop to print out each key and value of the dictionary.
7. Take a value from the command line for `fav_thing` and print the value of that item from the dictionary.
8. Print out all the keys to the user so that they know what to pick from. Check out `input()`. Here is a [link about input](https://www.w3schools.com/python/ref_func_input.asp)

9. Change the value of your favorite organism.

10. Get the `fav_thing` from the command line and a new value for that key. Change the value with the user inputted value. And print out a confirmation.




Sets
===================

1. Make a set using the two different syntaxes for creating a set `myset = set()` and `myset2 = {}`.  

```python
mySet = set('ATGTGGG')
mySet2 = {'ATGTGGG'}
```
-  What is the difference?
-  Does it matter which method you use?
-  How many items are in mySet and mySet2?



2. Write a script that creates 2 sets using the collections of numbers below. Find the intersection, difference, union, and symmetrical difference between these two sets.  
    - 3, 14, 15, 9, 26, 5, 35, 9
    - 60, 22, 14, 0, 9  


3. Create a set using the function `set()` and a DNA sequence, what will you get back? Try it with this sequence:

```
GATGGGATTGGGGTTTTCCCCTCCCATGTGCTCAAGACTGGCGCTAAAAGTTTTGAGCTTCTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGGCAGCCAGACTGCCTTCCGGGTCACTGCCATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATTCGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTNNGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGTGGGGTTTTCCCCTCCCATGTGCTCAAGACTGGCGCTAAAAGTTTTGAGCTTCTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGGCAGCCAGACTGCCTTCCGGGTCACTGCCATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATTCGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACX
```

Dictionaries and Sets and File I/O
=====================

1.  **Nucleotide Composition**. Write a script that:

  - determines the unique characters in this sequence

  ```
GAACTCCAAAAATGAAAACATAGTAGCAATCAAAGCATCCCACTATTTTTTGTCTCTCGTTTCATTAGCGTTGTAAATTACTGATACCCTACTATACCTCTACAAGGCCTTTGTCATCTTTTTACTCAAGTGTGAAATCATCACTTATTGTATGAAGGATGAGCTTTCCGTTCGCTAGTTTGCTGAAAAGGCCTTCTGCAATAAGCTCTCTATTATCTTTAAAAAAACCTGGTTCCTGGTCTTCCATTCTGCTAAAAGCTGTAGGGGTTTTATCACGAGATTCCCGTTGGCATTCTGACTTATTAAAAATGCTTACAGAAGAAATGGATTCTTTAAATGGTCAAATTAATACGTGGACAGATAATAATCCTTTATTAGATGAAATTACGAAGCCATACAGAAAATCTTCAACTCGTTTTTTTCATCCGCTTCTTGTACTTCTAATGTCTAGAGCATCAGTAAATGGGGATCCACCGAGTCAGCAACTATTTCAAAGGTACAAACAACTTGCCCGTGTAACAGAATTGATTCATGCTGCCAATATAATTCATATTAATATTGGAGAAGAACAAAGCAACGAACAGATTAAACTTGCAACGTTGGTTGGAGATTATTTACTCGGAAAGGCGTCTGTTGATTTAGCACATTTAGAAAACAACGCTATTACAGAAATTATGGCTTCTGTTATTGCAAACTTAGTTGAAGGGCACTTCGGAAGCCGACAAAATGGCTCTGTTGGTTTGTCAAACGAACGAACCATCCTTCTGCAATCAGCCTTTATGCCAGCAAAGGCATGTTTATGCGCAAGCATATTGAATAACTCATCACAATACATTAATGATGCGTGTTTCAATTATGGAAAATTTCTAGGCTTATCGCTGCAACTGGCCCATAAGCCTGTATCTCCTGACGCCCAAGTTTTGCAAAAGAATAATGACATTTTGAAAACATATGTTGAGAATGCCAAGAGCTCATTGTCTGTTTTCCCCGATATAGAGGCTAAGCAAGCTCTCATGGAAATCGCTAATAGTGTTTCGAAGTAATCGACAGGTATTGTATCCTGGATTAATATTAGGGTGGCTCATGCATGCTCGTGCAATCGTAACAAATATGTCTTTCTTTTACGAATTTTAACGCTTCAATATAAATCATATTTTTCCTCA
  ```

  - iterate over each unique character and count the number found in the sequence
  - store each count in a dictionary. example: `nt_comp['A']=2`
  - when you are done counting nucleotides, report the number of each unique nucleotide
  - also when you are done counting, calculate and report the GC content (G_count + C_count / total_nucleotides ).




5. Write your first FASTA parser script. This is a script that reads in a FASTA file ([Python_06.fasta](../files/Python_06.fasta)) and stores each FASTA record separately for easy access for future analysis.

Things to keep in mind:
   - open your file
   - read each line
   - is your line a header line? is it a sequence line?
   - does a single FASTA record have one line of sequence or multiple lines of sequence?

   HINTS: use file I/O, if statements and dictionaries to write your first FASTA parser. Some other useful functions and methods are find, split, string concatenation.

   At the end, your script should return the following:

   fastaDict = {
      'seq1' : 'AAGAGCAGCTCGCGCTAATGTGATAGATGGCGGTAAAGTAAATGTCCTATGGGCCACCAATTATGGTGTATGAGTGAATCTCTGGTCCGAGATTCACTGAGTAACTGCTGTACACAGTAGTAACACGTGGAGATCCCATAAGCTTCACGTGTGGTCCAATAAAACACTCCGTTGGTCAAC' ,
      'seq2' : 'GCCACAGAGCCTAGGACCCCAACCTAACCTAACCTAACCTAACCTACAGTTTGATCTTAACCATGAGGCTGAGAAGCGATGTCCTGACCGGCCTGTCCTAACCGCCCTGACCTAACCGGCTTGACCTAACCGCCCTGACCTAACCAGGCTAACCTAACCAAACCGTGAAAAAAGGAATCT' ,
      'seq3' : 'ATGAAAGTTACATAAAGACTATTCGATGCATAAATAGTTCAGTTTTGAAAACTTACATTTTGTTAAAGTCAGGTACTTGTGTATAATATCAACTAAAT' ,
      'seq4' : 'ATGCTAACCAAAGTTTCAGTTCGGACGTGTCGATGAGCGACGCTCAAAAAGGAAACAACATGCCAAATAGAAACGATCAATTCGGCGATGGAAATCAGAACAACGATCAGTTTGGAAATCAAAATAGAAATAACGGGAACGATCAGTTTAATAACATGATGCAGAATAAAGGGAATAATCAATTTAATCCAGGTAATCAGAACAGAGGT' }



## Extra: Expand on the nucleotide composition exercise 
  - get the raw file [Python_06.seq.txt](https://raw.githubusercontent.com/prog4biol/pfb2025/master/files/Python_06.seq.txt)
  - in a script, open this file
  - iterate over each line in this file (seqName\tsequence\n)
     - for each sequence:
         - calculate and store the count of each unique nucleotide character in a dictionary
         - report the name, total of each nucleotide count, and the GC content 


