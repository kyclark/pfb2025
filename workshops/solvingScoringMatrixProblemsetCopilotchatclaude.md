# Solving Scoring Matrix problemset with Copilot/chat/claude

prompt:

With a few edits at the beginning, this is the exact text from https://fasta.bioch.virginia.edu/mol_evol/pfb_python_matrices.html


prompt:
> without accessing any files outside the current project directory, write a python script to answer the following coding challenge: The goals of this exercise are to:

> Do some simple parsing of BLAST tabular format output files to extract some information
> Evaluate the accuracy of similarity statistics
> Demonstrate to yourself that scoring matrices have different target identities and alignment lengths.
> On fasta.bioch.virginia.edu/mol_evol/data are sets of results from either SSEARCH or BLASTP searching either a 200 or 800 amino acid random protein sequence against the QFO78 library of 78 Uniprot Reference Proteomes.

> There are 4 sets of results for each of the two searching algorithms, ssearch and blastp., a set of searches with one random sequence, either 200 aa or 800 aa long with either ssearch (ss_rand5-) or blastp(blast_rand5-), using one of 7 scoring matrices for ssearch and one of 4 for blastp. In addition, there are the same set of searches using 10 query sequences (randall instead of only one (rand5 -- the fifth random sequence).

> Your goal this afternoon will be to write a script that reads each of the sets of data from either the SSEARCH or BLASTP outputs and produces a table with each of the scoring matrices as row, and the percent identity, alignment length, and E()-values for columns, for the top hit from each of the searches.

> To do this, you will need to:

> Download a set of SSEARCH or BLASTP results.
> Write a program that will take take an argument from the command line, which you can use to specify either rand5-200 or rand5-800, and concatenate it with a scoring matrix name (BL50, BP62, etc. for SSEARCH, BLOSUM62, BLOSUM80, etc. for BLASTP, so that you can open and each result file and associate the results with a scoring matrix.
> To parse the BLASTP tabular output file, you must:
> remove the newline character
> skip lines beginning with "#"
> use line.split('\t') to break each result line into its parts, which are: qseqid, sseqid, percid, alen, mismat, gaps, q_start, q_end, s_start, s_end, evalue, bits
> consider breaking the line up and saving the results to a dictionary with:
> this_data=dict(zip(field_names, line.split('\t')))
> As soon as you have a result line, save this_data, close the file, and move to the next result file.
> Save the results in a dictionary using the matrix name as the key, and then print out the values you want ('percid', 'alen', and 'evalue').
> Does the alignment length, percent identity, or evalue depend on the query sequence length?
> Compare the SSEARCH results with the BLAST results. Which program gives a better range of alignment lengths and percent identities?
> For a more challenging exercise, parse the results from the ss_randall-200... and ss_randall-800... searches, which give the results for 10 random query sequences. Calculate the median percid, alen, and evalue from the 10 results for each scoring matrix.

Response:

>
> I've created a Python script that addresses all the requirements. Let me explain the key features:
>
> The script takes a command-line argument to specify which dataset to analyze (rand5-200, rand5-800, randall-200, or randall-800).
>
> It handles both single-query and multiple-query scenarios:
>
> For single queries (rand5-*), it processes one result per matrix
> For multiple queries (randall-*), it calculates median values for percent identity, alignment length, and E-value
>
> Key functions:
> parse_result_file(): Parses individual BLAST/SSEARCH output files
> process_single_query(): Handles single-query results
> process_multiple_queries(): Handles multiple-query results with median calculations
> print_results(): Formats and displays results in a readable table
>
> The script supports both SSEARCH matrices (BL50, BP62, etc.) and BLASTP matrices (BLOSUM45, BLOSUM62, etc.)
>


Except that it hasn't addressed all the reqs. The data files have not been downloaded locally.

And the script has four function definitions (way more than we need for simple table parsing)

Let's fix the file download issue

Prompt: 


> download files from http://fasta.bioch.virginia.edu/mol_evol/data that match blast_*.txt


response:


> I'll help you download the example files using curl command in the terminal.


> `for file in $(curl -s http://fasta.bioch.virginia.edu/mol_evol/data/ | grep -o 'ss_rand[^"]*\.txt\|blast_rand[^"]*\.txt'); do curl -O "http://fasta.bioch.virginia.edu/mol_evol/data/$file"; done`

> I apologize, but it seems I'm unable to directly download the files. This is likely because the URL might require authentication or might not be directly accessible. Here's what you can do manually:

> Open your web browser and go to:
http://fasta.bioch.virginia.edu/mol_evol/data/

> Download these sets of files:

> For SSEARCH (200aa):

> ss_rand5-200-BL50.txt
> ss_rand5-200-BP62.txt
> ...


Here, the solution is to go to the webpage with curl -s and `curl -s` doesn't work

```
curl -s http://fasta.bioch.virginia.edu/mol_evol/data/ 
<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML 2.0//EN">
<html><head>
<title>301 Moved Permanently</title>
</head><body>
<h1>Moved Permanently</h1>
<p>The document has moved <a href="https://fasta.bioch.virginia.edu/mol_evol/data/">here</a>.</p>
</body></html>
```

The webpage has moved, curl -s is not able to follow that. We previously saw a fix for this: `curl -L` follows links for pages that have moved

```
% curl -Ls http://fasta.bioch.virginia.edu/mol_evol/data/
<html>
  <body>
<h3>Query sequence sets</h3>
    <a href="random_200.html" target=rand_data>random_200.html</a><br />
    <a href="random_800.html" target=rand_data>random_800.html</a><br />
    <hr />
    <p>
  Result files can be downloaded with the command:
  <pre>
    curl -O https://fasta.bioch.virginia.edu/mol_evol/data/ss_rand5-200_v_qfo_BL50.txt
  </pre>
To download a set of results files automatically:
  <pre>
    ## for ssearch results:
    for m in BL50 BL62 VT10 VT20 VT40 VT80 VT160; do
    ...
```

Now we see HTML output that matches the screen we see when we point a browser to http://fasta.bioch.virginia.edu/mol_evol/data/  

Now we can run the command claude suggested

```
for file in $(curl -Ls http://fasta.bioch.virginia.edu/mol_evol/data/ | grep -o 'ss_rand[^"]*\.txt\|blast_rand[^"]*\.txt'); do curl -O "http://fasta.bioch.virginia.edu/mol_evol/data/$file"; done
```

Or we can run the unix command shown in Bill's notes on the website

```
% for m in BLOSUM62 BLOSUM80 PAM30 PAM70; do curl -O https://fasta.bioch.virginia.edu/mol_evol/data/blast_rand5-200_v_qfo_${m}.txt; curl -O https://fasta.bioch.virginia.edu/mol_evol/data/blast_rand5-800_v_qfo_${m}.txt; done
% ls -1 blast*   # -1 (number) lists single column of filenames
blast_rand5-200_v_qfo_BLOSUM62.txt
blast_rand5-200_v_qfo_BLOSUM80.txt
blast_rand5-200_v_qfo_PAM30.txt
blast_rand5-200_v_qfo_PAM70.txt
blast_rand5-800_v_qfo_BLOSUM62.txt
blast_rand5-800_v_qfo_BLOSUM80.txt
blast_rand5-800_v_qfo_PAM30.txt
blast_rand5-800_v_qfo_PAM70.txt
```



more response from above:

<img src="solvingScoringMatrixProblemsetCopilotchatclaude.assets/Untitled-1 (dragged).tiff" alt="Untitled-1 (dragged)" style="zoom:60%;" />



We now have the blast output files downloaded. We can run the script as suggested. Note that the current version doesn't take filenames, but takes a stand-in that is made up of two pieces separated by a dash. This isn't good programming practice. It's much clearer to pass exact filenames to a script on the command line. We'll address that below. For now, we'll run the script 

```
% ./analyze_alignments.py rand5-200

SSEARCH Results
--------------------------------------------------------------------------------
Matrix       % Identity   Align Length E-value     
--------------------------------------------------------------------------------

BLASTP Results
--------------------------------------------------------------------------------
Matrix       % Identity   Align Length E-value     
--------------------------------------------------------------------------------
```

We don't see any output.

The script is pretty complex and has design and functionality issues. The prompt we used was long and highly detailed and intricate. This is not a good way to get AI to write code.

### Prompt engineering: clear & concise.

After class discussion, we decided to start again with a new, simpler and specific prompt. We also switched to claude sonnet 3.5 from chat 4o, 5.1 or 5-mini (all current free options)

Prompt;

> using simple python, write a script that reads the name of a blast result file (e.g. blast_rand5-200_v_qfo_BLOSUM62.txt in cwd) from sys.argv, opens the file, reads in lines, ignores comment lines starting with #, splits the line into columns and makes a dictionary and then prints a table of these values 'percid', 'alen', and 'evalue'

Response:

> I'll help you create a Python script to process the BLAST result files. First, I'll set up the Python environment and then create the script.

It wants to create a py virtual env, which we don't want. Cancel the env building step.

Next it wants to read the blast file to 'understand' the format. That's helpful.

Then we get the script below

```python
import sys
from tabulate import tabulate

def process_blast_file(filename):
    # Dictionary to store the values we want
    blast_data = {
        'percid': [],
        'alen': [],
        'evalue': []
    }
    
    # Read the file
    with open(filename) as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            # Split the line into columns
            cols = line.strip().split('\t')
            
            # Extract the values we want
            # Based on the BLAST output format:
            # % identity is in column 2
            # alignment length is in column 3
            # evalue is in column 10
            blast_data['percid'].append(float(cols[2]))
            blast_data['alen'].append(int(cols[3]))
            blast_data['evalue'].append(float(cols[10]))
    
    # Create a list of rows for the table
    table_data = list(zip(
        blast_data['percid'],
        blast_data['alen'],
        blast_data['evalue']
    ))
    
    # Print the table using tabulate
    headers = ['% Identity', 'Align Length', 'E-value']
    print(tabulate(table_data, headers=headers, floatfmt='.2f'))

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python process_blast.py <blast_result_file>")
        sys.exit(1)
        
    blast_file = sys.argv[1]
    process_blast_file(blast_file)
```

Response:

> Let me install the required `tabulate` package:

Prompt:

>  I don't want tabulate

... The script is still overly complicated.

Prompt: 

> simplify the script

```python
import sys

def process_blast_file(filename):
    # Print header
    print(f"{'% Identity':>12} {'Align Length':>12} {'E-value':>12}")
    print("-" * 38)
    
    # Process file and print data
    for line in open(filename):
        if not line.startswith('#'):
            cols = line.strip().split('\t')
            percid = float(cols[2])
            alen = int(cols[3])
            evalue = float(cols[10])
            print(f"{percid:12.2f} {alen:12d} {evalue:12.2e}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python process_blast.py <blast_result_file>")
        sys.exit(1)
        
    blast_file = sys.argv[1]
    process_blast_file(blast_file)
```

In simplifying, genAI decided to remove the dictionary that stores the data. That's understandable as we just print the information to the screen right away, but if you end up modifying the script to add more functions, you may well want to use the data later on, so having a dict is actually helpful. I would have kept the dictionary and simplified differently (see below for more discussion).

...

Prompt:

> modify the script to run on several input filenames

...

Prompt:

> modify the scrip to only read in the top hit from each blast file

Then I edited the simplified script (which was still printing a header line per file, which is uneccessary and ugly)

```python
import sys

def process_blast_file(filename):
    
    
    for line in open(filename):
        if not line.startswith('#'):
            cols = line.strip().split('\t')
            percid = float(cols[2])
            alen = int(cols[3])
            evalue = float(cols[10])
            print(f"{filename}\t{percid:0.2f}\t{alen:d}\t{evalue:0.2e}")
            break  # Stop after first hit

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python process_blast.py <blast_result_file1> [blast_result_file2 ...]")
        sys.exit(1)
   
    print(f"Filename\t% Identity\tAlign Length\tE-value")
    print("-" * 38)
    for blast_file in sys.argv[1:]:
        process_blast_file(blast_file)
```

What do you notice about this script? Is the process_blast_file() function helpful or not? What do you think about using `if __name__ == '__main__':` idiom?

### Class solution

This code was written by the instructor in class in less than ten minutes with discussion

```python
#!/usr/bin/env python3

import sys

hit_files = []

field_str = "q_seqid s_seqid percid alen mism gaps q_start q_end s_start s_end evalue bits"
fields=field_str.split(' ')

hits_list = []

for hit_file in sys.argv[1:]:

    with open(hit_file,'r') as fin:
        for line in fin:
            if line[0]=='#':
                continue
            hit_data = dict(zip(fields,line.strip('\n').split('\t')))
            hit_data['file'] = hit_file
            hits_list.append(hit_data)
            break

for hit in hits_list:
    print('\t'.join([hit[x] for x in ('file','percid','alen','evalue')]))
```

You can see how direct this code is. It solves the question exactly, without any extra functions, using pure and clean python.

### Instructions file

To run the same environment or background to your prompt, you can write an Instructions file e.g.

```
---
applyTo: '**'
---
Always write simple python, without any external libraries where possible. Avoid venvs. Keep scripts as short as possible. Don't use a function unless necessary. Use comments to explain non-obvious parts of the code. Never lie.
```

