# Really Useful Tools and Bioinformatics Pipelines

We're going to review ways to get help, and dive into more advanced unix programming. Then we'll see how to design and construct bioinformatics pipelines and think about running them.

### Getting help

`pydoc`  
google 'python3 list comprehensions`  
https://docs.python.org/3/    -> Quick search  
Help is available inside python interactive shell

```python
>>>help()
```
Just like for `man`, help text appears inside a pager like `more` or `less`.  You can use the following commands:  
```
space  -> next page  
b      -> back a page  
return -> next line  
/      -> search for a string  
q      -> quits the pager
```

```
>>> help(str)
Help on class str in module builtins:

class str(object)
 |  str(object='') -> str
 |  str(bytes_or_buffer[, encoding[, errors]]) -> str
 |  
 |  Create a new string object from the given object. If encoding or
 |  errors is specified, then the object must expose a data buffer
 |  that will be decoded using the given encoding and error handler.
 |  Otherwise, returns the result of object.__str__() (if defined)
 |  or repr(object).
 |  encoding defaults to sys.getdefaultencoding().
 |  errors defaults to 'strict'.
...
|  count(...)
 |      S.count(sub[, start[, end]]) -> int
 |      
 |      Return the number of non-overlapping occurrences of substring sub in
 |      string S[start:end].  Optional arguments start and end are
 |      interpreted as in slice notation.
...
```
  and `dir()`


```  
  >>> help(str.split)

Help on method_descriptor:

split(...)
    S.split(sep=None, maxsplit=-1) -> list of strings
    
    Return a list of the words in S, using sep as the
    delimiter string.  If maxsplit is given, at most maxsplit
    splits are done. If sep is not specified or is None, any
    whitespace string is a separator and empty strings are
    removed from the result.
(END)

```


## Advanced Unix
### awk

awk is a simple unix utility for reformatting text files. An awk script would look like this

```
#!/usr/bin/env awk
BEGIN { print "File\tOwner"}   # block executed before main script
{ print $9 "\t" $3}          # main script
END { print " - DONE -" }      # block executed after main script
```

You could run it like this `awk table.awk`. Each column (whitespace-separated) in the input appears in your script as \$1, \$2, \$3 etc. A bit like sys.argv in python.

Let's ignore the BEGIN and END blocks for now.

How could you take a long file listing and print out the owner of each file?

```bash
% ls -l
total 19664
-rw-r--r--  1 simonp  staff  312 Oct 20 11:05 scope_global.py
-rw-r--r--  1 simonp  staff  201 Oct 20 11:03 scope_global.py~
-rw-r--r--  1 simonp  staff  323 Oct 20 10:40 scope_w_function.py
-rw-r--r--  1 simonp  staff  210 Oct 20 10:33 scope_w_function.py~
-rw-r--r--  1 simonp  staff    5 Oct 15 14:15 test.nt.fa
-rw-r--r--  1 simonp  staff  103 Oct 17 19:27 while.py
-rw-r--r--  1 simonp  staff  160 Oct 17 19:27 while_else.py

```

To get the value of a scalar variable in unix tools, you prefix the variable name with a dollar sign. `$1` is the value in the first column of the input, `$i` is the value stored in the variable `i`. Here are the column variables explicitly. This is not shell output. Just a picture.

```
$1         $2 $3      $4     $5  $6  $7 $8    $9
-rw-r--r--  1 simonp  staff  160 Oct 17 19:27 while_else.py
```

We want to print the file and the owner. Find the variables. The order can be whatever we want. The awk part would look like this 

`awk '{print $9 "\t"  $3 }'`

```
while_else.py      simonp
```

Note space-separated values would be coded like so

`awk  '{print $9,  $3 }'`

```
while_else.py simonp
```



Recall, how do we get the long listing? `ls -l`. Note also the first line of `ls -l` output is `total 19664`

Put these together with our friend pipe `|` and add a tab between columns

```
% ls -l | awk 'BEGIN{print "filename\towner"} $1 != "total" {print $9 "\t"  $3}'  # no commas!! why??
filename	owner
scope_global.py 	 simonp
scope_global.py~ 	 simonp
scope_w_function.py 	 simonp
scope_w_function.py~ 	 simonp
test.nt.fa 	 simonp
while.py 	 simonp
while_else.py 	 simonp
```

Printing files modified in October, on or after the 20th

```
% ls -l | awk '$7 >= 20 && $6 == "Oct" {print $9 "\t" $3}'
scope_global.py 	 simonp
scope_global.py~ 	 simonp
scope_w_function.py 	 simonp
scope_w_function.py~ 	 simonp
```

How can you print the number of records in a fastq file? Each record is 4 lines long.

```
wc -l mouse_R7385F_*fastq | awk '{print $1/4 "\t" $2}'
85534	mouse_R7385F_smpl1.fastq
688103	mouse_R7385F_smpl10.fastq
100668	mouse_R7385F_smpl2.fastq
802370	mouse_R7385F_smpl3.fastq
407677	mouse_R7385F_smpl4.fastq
```

Another example that counts how many reads align with score 60 or more (great alignment) to the first 5kb of the reference.

```
samtools view /Users/simonp/src/samtools-1.14/test/mpileup/mpileup.2.bam | awk '$4 <=5000 && $5 >=60  {print $1 "\t" $2 "\t" $5}'  | wc -l
```

#### regular expressions in awk

Much simpler than python, put the re between slashes

```
awk '/li/ { print $2 }' mail-list.txt   # match on the whole line

awk '{ if ($1 ~ /J/) print }' inventory-shipped  # match on a field (column 1)

ls -l | awk ' /smpl[0-9]\./ {print $9}'  # match a single digit 0,1,2,3 etc before .fastq
```



Use awk if you want to reorder columns in a file, do simple filtering and calculations etc. Works great in piped command lines.

### Unix aliases

Here's a way to save typing

`alias` is a unix command that goes in your ~/.profile file (bash) or ~/.zshrc (zsh). Make one with `vi` if you don't have one already. 

```bash
alias ll='ls -l'
alias lr='ls -ltrh'
```

To get these changes, `source ~/.profile` or `source ~/.zshrc` or open a new window in terminal. Now you can type `lr` instead of `ls -ltrh` 



## Workflows and approaches

### Saving time and effort.

Your coding day is time spent doing these things:

* thinking: design
* preparation, testing
* writing code
* debugging
* running code
* thinking: analysis
* more writing, thinking
* report results

Where do you spend most of your time? What can you save time on? The more you plan out coding and check your data, the faster you'll get to the important second half of this list.

* thinking: design   Lots of time!
* preparation, testing  Lots of time!
* writing code  Quick now that you've done the first two
* debugging  Quick now that you've done the first two
* running code Very quick
* thinking: analysis  Spend lots of time on this and later steps
* more writing, thinking
* report results

Assume your data is corrupted, even if it came from a good colleague. This will stress test your code before you start writing.

Check for consistent numbers of columns in your data, files that end halfway through a line are truncated or corrupted. Is a column always numbers or mixed numbers and text? Be precise about numbers. `2000-3000` is not a number. Nor is `5kb`. Do some fields have quotes or other unusual characters, accents? Do the values seem reasonable? Are values for gene lengths between 1,000 and 10,000bp for example?

Data consistency, corruption, sanity checks
  NGS data generation: illumina, pacbio  
  formats - see biopython  
  (un)compression  
filesize and md5 checksums

## Designing and Implementing a Bioinformatics Pipeline

Say you want to automate blast runs. Your first challenge is to come up with an exact step by step recipe of how you would run a single blast job from the command line. We skip several hours of reading, research, trial and error for the sake of teaching, but here are some key steps.

Install standalone unix executables (programs) for blast+ from ncbi (see help here https://www.ncbi.nlm.nih.gov/books/NBK52640/)

makeblastdb formats a database file so that you can search therein for sequences that are similar to your query. Formatting generates a series of files that end `.pXX` where X is any character. Here's an example of the output files when the input file is `EcoliO157.uniprot.fa` (see the command below)

```
EcoliO157.uniprot.fa.pin
EcoliO157.uniprot.fa.phr
EcoliO157.uniprot.fa.psq
EcoliO157.uniprot.fa.pog
EcoliO157.uniprot.fa.pot
EcoliO157.uniprot.fa.pto
EcoliO157.uniprot.fa.ptf
EcoliO157.uniprot.fa.pos
EcoliO157.uniprot.fa.pdb
```

We just need to know that blast needs these files to search with a query. The details are not important.

When we run blast, by default, the output shows up on the screen in blast report format, like so

```
% blastp -query ilvG.bacteria.prot.fa -db EcoliO157.uniprot.fa -evalue 1e-10  
BLASTP 2.10.1+


Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.
Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.
Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of
protein database search programs", Nucleic Acids Res. 25:3389-3402.


Reference for composition-based statistics: Alejandro A. Schaffer,
L. Aravind, Thomas L. Madden, Sergei Shavirin, John L. Spouge, Yuri
I. Wolf, Eugene V. Koonin, and Stephen F. Altschul (2001),
"Improving the accuracy of PSI-BLAST protein database searches with
composition-based statistics and other refinements", Nucleic Acids
Res. 29:2994-3005.



Database: EcoliO157.uniprot.fa
           4,587 sequences; 1,403,709 total letters



Query= sp|P66947|ILVG_MYCBO Probable acetolactate synthase OS=Mycobacterium
bovis (strain ATCC BAA-935 / AF2122/97) OX=233413 GN=ilvG PE=3 SV=1

Length=547
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

P0AFI0 Oxalyl-CoA decarboxylase OS=Escherichia coli (strain K12) ...  197     4e-57
P0DP90 Acetolactate synthase isozyme 2 large subunit OS=Escherich...  149     1e-39
P00893 Acetolactate synthase isozyme 3 large subunit OS=Escherich...  149     2e-39
P08142 Acetolactate synthase isozyme 1 large subunit OS=Escherich...  147     5e-39
P0AEP7 Glyoxylate carboligase OS=Escherichia coli (strain K12) OX...  139     3e-36
P0DP89 Putative acetolactate synthase isozyme 2 large subunit OS=...  112     4e-28
P07003 Pyruvate dehydrogenase [ubiquinone] OS=Escherichia coli (s...  98.6    1e-22


>P0AFI0 Oxalyl-CoA decarboxylase OS=Escherichia coli (strain K12) OX=83333 
GN=oxc PE=1 SV=1
Length=564

 Score = 197 bits (501),  Expect = 4e-57, Method: Compositional matrix adjust.
 Identities = 170/556 (31%), Positives = 275/556 (49%), Gaps = 30/556 (5%)

Query  1    MSTDTAPAQTMHAGRLIARRLKASGIDTVFTLSGGHLFSIYDGCREEGIRLIDTRHEQTA  60
            MS        MH   +I   LK + IDT++ + G  +  +    + EGIR I  RHEQ+A
Sbjct  1    MSDQLQMTDGMH---IIVEALKQNNIDTIYGVVGIPVTDMARHAQAEGIRYIGFRHEQSA  57

Query  61   AFAAEGWSKVTRVPGVAALTAGPGITNGMSAMAAAQQNQSPLVVLGGRAP--ALRWGMGS  118
             +AA     +T+ PG+    + PG  NG++A+A A  N  P++++ G +    +    G 
Sbjct  58   GYAAAASGFLTQKPGICLTVSAPGFLNGLTALANATVNGFPMIMISGSSDRAIVDLQQGD  117

Query  119  LQEIDHVPFVAPVARFAATAQSAENAGLLVDQALQAAVSAPSGVAFVDFPMD-HAFSMSS  177
             +E+D +    P A+ A      ++ G+ + +A++ +VS   G  ++D P +  A +M  
Sbjct  118  YEELDQMNAAKPYAKAAFRVNQPQDLGIALARAIRVSVSGRPGGVYLDLPANVLAATMEK  177

Query  178  DNGRPGALTELPA--GPTPA----GDALDRAAGLLSTAQRPVIMAGTNVWWGHAEAALLR  231
            D     ALT +     P+PA      ++  A  LL+ A+RP+I+ G    +  A+  L  
Sbjct  178  DE----ALTTIVKVENPSPALLPCPKSVTSAISLLAKAERPLIILGKGAAYSQADEQLRE  233

Query  232  LVEERHIPVLMNGMARGVVPADHRLAFSRARSKALGEADVALIVGVPMDFRLGFGGV-FG  290
             +E   IP L   MA+G++   H L+ + ARS AL  ADV ++VG  +++ L  G   + 
Sbjct  234  FIESAQIPFLPMSMAKGILEDTHPLSAAAARSFALANADVVMLVGARLNWLLAHGKKGWA  293

Query  291  STTQLIVADRVEPAR-EHPRPVAAGLYGDLTAT----LSALAGSGGTDHQGWIEELATAE  345
            + TQ I  D +EP   +  RP+A  + GD+ ++    L+ L  +  T    W + L   +
Sbjct  294  ADTQFIQLD-IEPQEIDSNRPIAVPVVGDIASSMQGMLAELKQNTFTTPLVWRDILNIHK  352

Query  346  TMARDLEKAELVDDRIPLHPMRVYAELAALL--ERDALVVIDAGDFGSYAGRMIDSYLPG  403
                     +L  D  PL+     + +  +L   +D  +V +  +    A  +ID Y P 
Sbjct  353  QQNAQKMHEKLSTDTQPLNYFNALSAVRDVLRENQDIYLVNEGANTLDNARNIIDMYKPR  412

Query  404  CWLDSGPFGCLGSGPGYALAAKLARPQRQVVLLQGDGAFGFSGMEWDTLVRHNVAVVSVI  463
              LD G +G +G G GYA+ A +      VV ++GD AFGFSGME +T+ R+N+ V  VI
Sbjct  413  RRLDCGTWGVMGIGMGYAIGASVTS-GSPVVAIEGDSAFGFSGMEIETICRYNLPVTIVI  471

Query  464  GNNGIWGLEKHPMEALYGYSVVA--ELRPGTRYDEVVRALGGHGELVSVPAELRPALERA  521
             NNG  G+ +     L G    +  +L    RYD+++ A  G G  V+   ELR AL   
Sbjct  472  FNNG--GIYRGDGVDLSGAGAPSPTDLLHHARYDKLMDAFRGVGYNVTTTDELRHALTTG  529

Query  522  FASGLPAVVNVLTDPS  537
              S  P ++NV+ DP+
Sbjct  530  IQSRKPTIINVVIDPA  545


>P0DP90 Acetolactate synthase isozyme 2 large subunit OS=Escherichia 
coli (strain K12) OX=83333 GN=ilvG PE=1 SV=1
Length=548
...
```

This output is hard to parse, so we add a flat `-outfmt 7` to print tab-separated output with comments like so

```
% blastp -query ilvG.bacteria.prot.fa -db EcoliO157.uniprot.fa -outfmt 7  -evalue 1e-10 
# BLASTP 2.10.1+
# Query: sp|P66947|ILVG_MYCBO Probable acetolactate synthase OS=Mycobacterium bovis (strain ATCC BAA-935 / AF2122/97) OX=233413 GN=ilvG PE=3 SV=1
# Database: EcoliO157.uniprot.fa
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 7 hits found
sp|P66947|ILVG_MYCBO	P0AFI0	30.576	556	356	14	1	537	1	545	4.20e-57	197
sp|P66947|ILVG_MYCBO	P0DP90	27.122	542	361	15	11	535	1	525	1.12e-39	149
sp|P66947|ILVG_MYCBO	P00893	26.902	539	361	12	9	522	2	532	1.52e-39	149
sp|P66947|ILVG_MYCBO	P08142	25.091	550	384	10	2	535	4	541	5.05e-39	147
sp|P66947|ILVG_MYCBO	P0AEP7	26.460	548	365	12	21	535	14	556	2.98e-36	139
sp|P66947|ILVG_MYCBO	P0DP89	28.615	325	219	6	11	326	1	321	3.59e-28	112
sp|P66947|ILVG_MYCBO	P07003	24.539	542	362	14	17	534	9	527	1.40e-22	98.6
# BLAST processed 1 queries
```

We can send the output to a file with `-out ilvG.bacteria.prot.fa.blastp.out`

```
% blastp -query ilvG.bacteria.prot.fa -db EcoliO157.uniprot.fa -outfmt 7 -out ilvG.bacteria.prot.fa.blastp.out  -evalue 1e-10 
%
```



Here's an example of what we'd run in the terminal. You can think of this as a shell script. It documents exactly what we do, but it isn't flexible: we can't run this script on different inputs or change the E-value.

```bash
#!/usr/bin/env zsh

# first you make a blast database from the fasta file of the sequences 
#you want to blast against (target or subject sequences)
makeblastdb -in EcoliO157.uniprot.fa -dbtype prot -parse_seqids

# now you can run a protein blast search with your query ilvG
blastp -query ilvG.bacteria.prot.fa -db EcoliO157.uniprot.fa -outfmt 7 -out ilvG.bacteria.prot.fa.blastp.out -evalue 1e-10

```

Make sure these commands work ok before you start making a pipeline.

Here are the steps we need our pipeline script to perform:

We need to print a usage message

We need to track when and how the script was run

We need to format the blast database if it hasn't been done already. 

We need to run blastp

We need to check blastp ran ok (unix return code)

Print summary table of hits

```python
#!/usr/bin/env python3
import subprocess
import sys
import datetime 
import os

# help message
# what's a better way to do this?
if len(sys.argv) < 4:
    print(f'Usage: {sys.argv[0]}   <query protein fasta>  <subject protein fasta>  <min E-value>')
    exit(1)
# get cmd line params
query = sys.argv[1]
db = sys.argv[2]
evalue = float(sys.argv[3]) # immediately convert to appropriate type

if not query.endswith( ('.fa','.fasta') ):
    print('Query input file needs to end with .fa or .fasta')
    exit(12)

if not db.endswith( ('.fa','.fasta') ):
    print('Subject input file needs to end with .fa or .fasta')
    exit(12)


# 2019-10-23 13:49:27.232603
now = str(datetime.datetime.now())
# cut down to 2019-10-23 13:49
now = now[0:16]

#log run command and time/date to screen
print('#' , ' '.join(sys.argv))
print('#' , 'was run on', now)

# see if we have a blast database and make one if not
# let's figure out the next line in detail -> see terminal session
if not os.path.exists(f'{db}.phr'):
    makeblastdbcmd = f'makeblastdb -in {db} -dbtype prot -parse_seqids'
    makeblastdb_run = subprocess.run(makeblastdbcmd, shell=True , stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    # Now we need to check the UNIX return code
    # always do this!
    # 0 = success
    # non-zero =failure
    if makeblastdb_run.returncode != 0:
        print("FAILED!")
        exit(5)

#generate blast output filename
blast_out = query + '.blastp.out'

# run the command
blastcmd = f'blastp -query {query} -db {db} -outfmt 7 -out {blast_out} -evalue {evalue}'

# object is returned after run command
blastcmd_run = subprocess.run(blastcmd, shell=True , stdout = subprocess.PIPE, stderr=subprocess.PIPE)

# Now we need to check the UNIX return code
# always do this!
# 0 = success
# non-zero =failure
if blastcmd_run.returncode != 0:
    print("FAILED!")
    exit(2)

# now parse results, 
homologs = {}
with open(blast_out,'r') as blast_results:
    for line in blast_results:
        line = line.rstrip()
        if line.startswith('#'): # skip comment lines
            continue
        fields = line.split('\t')
        query = fields[0]
        subject = fields[1]
        evalue = float(fields[10]) # this will be a string because 
                            # we read in from a file
                            # don't forget to convert to float
        # collect hits and their E-values into dictionary
        if query not in homologs:
            homologs[query] = [ (subject, evalue) ]
        else:
            homologs[query].append( (subject,evalue) )

print('Hit summary')
for query in sorted(homologs):
    print('Query:',query)
    for data in homologs[query]:
        subject,evalue = data
        print(f'{subject} E-value={evalue}' )

```

Here's the output from running this script

```
% ./blastPipeline.py 
Usage: ./blastPipeline.py   <query protein fasta>  <subject protein fasta>  <min E-value>
% ./blastPipeline.py ilvG.bacteria.prot.fa EcoliO157.uniprot.fa 1e-10
# ./blastPipeline.py ilvG.bacteria.prot.fa EcoliO157.uniprot.fa 1e-10
# was run on 2025-10-23 07:16
Hit summary
Query: sp|P66947|ILVG_MYCBO
P0AFI0 E-value=4.2e-57
P0DP90 E-value=1.12e-39
P00893 E-value=1.52e-39
P08142 E-value=5.05e-39
P0AEP7 E-value=2.98e-36
P0DP89 E-value=3.59e-28
P07003 E-value=1.4e-22
(base) simonp:~/Documents/programmingForBiology/PFB2024_course_lectures% 
```





## Bioinformatics How do I ...?

Here are some bare-bones guidelines to get you going.
### filtering illumina sequence data: 

cutadapt  
trimgalore   
trimmomatic

#### QC sequence data:

fastqc
multiqc (aggregate results from bioinformatics analyses across many samples into a single report)

### resequencing, variant calling

GATK, FreeBayes

### finding genes

Maker (eukaryotes),   
Prokka/prodigal (prokaryotes)

### predicting gene function

Interproscan

### Databases store large data for easy searching and retrieval

sqlite3 is the simplest. It stores your data in a single file. Portable and simple. Gets you up and running quickly.

python has a module

```python
import sqlite3
```

We won't talk about DBs more here, but they are useful for larger data projects. They use their own language: SQL = structured query language. 

### Public databases

__NCBI__ (sequences, searching)  
Landmark (model organisms)
Above are better than: nr (proteins),  nt (nucleotides)  Lots of data, uncurated, complete
Sequence Read Archive (SRA) 454, illumina, short reads

__Uniprot__ (sequences, searching)
http://www.uniprot.org
Curated, smaller, not as inclusive as nr. 
Helpful for speeding up analysis: UniRef90 (sequences clustered at 90% identity, which is approximately genus level). Much smaller than full database. 

__PDB Protein Data Bank__
For protein structures

__Genomes__
Ensembl, JGI (plants, fungi, bacteria/metagenomes),  NCBI genome
Organism data bases, beware data quality: some are excellent, some not so well resourced.

### Write web apps

See Flask python library

### Debug my script

Run your script with the debugger module `pdb` for python debugger. Not very sophisticated, but very useful. It starts an interactive debugger that's a bit like the python interactive shell, but you are inside your script. 

We were doing this
```bash
% python3 while.py
```
Now we add `-m pdb` so it becomes
```bash
% python3 -m pdb while.py
> /Users/simonp/git/pfb2017/scripts/while.py(3)<module>()
-> count = 0
(Pdb) h

Documented commands (type help <topic>):
========================================
EOF    c          d        h         list      q        rv       undisplay
a      cl         debug    help      ll        quit     s        unt      
alias  clear      disable  ignore    longlist  r        source   until    
args   commands   display  interact  n         restart  step     up       
b      condition  down     j         next      return   tbreak   w        
break  cont       enable   jump      p         retval   u        whatis   
bt     continue   exit     l         pp        run      unalias  where    

Miscellaneous help topics:
==========================
exec  pdb

(Pdb) 

```
q quits  
h gets help  

It's a good idea to make alias for python3 -m pdb  in .profile. How would we do that?

### Write bigger python coding projects? 

VS Code

### Tell if my code is slow

Even though python is much slower than C and C++, is your script running too slowly? How can you tell? 
Two things to think about

* Is debugging painfully slow? Use the smallest test data sets you can to test and debug your script
* Do you have time to get a cup of coffee while your script is running? If you come back to your script and it's still running, and you're bored, look into speeding it up. Look up profilers, parallelization, other peoples' experiences (seqanswers.com, stackoverflow.com)

Once you are a decent programmer, the speed up you'll get  (a few milliseconds) from tinkering with your script (several hours) will not be worth it.
