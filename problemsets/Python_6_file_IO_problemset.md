Python 6 - Files and IO - Problem Set
===================



1. Write a script to do the following to [Python_06.txt](https://raw.githubusercontent.com/prog4biol/pfb2025/master/files/Python_06.txt)
   - Open and read the contents.  
   - Uppercase each line
   - Print each line to the STDOUT


2. Modify the script in the previous problem to write the contents to a new file called "Python_06_uc.txt"



3. Open and print the reverse complement of each sequence in [Python_06.seq.txt](https://raw.githubusercontent.com/prog4biol/pfb2025/master/files/Python_06.seq.txt). Each line is the following format:    `seqName\tsequence\n.` Make sure to print the output in FASTA format including the sequence name and a note in the description that this is the reverse complement. Print to STDOUT and capture the output into a file with a command line redirect '>'. 
   - **Remember is is always a good idea to start with a test set for which you know the correct output.**

4. FASTQ File Parsing:
___ 
FASTQ File Overview based on wikipedia [wikipedia/FASTQ](https://en.wikipedia.org/wiki/FASTQ_format):
___  

A FASTQ file has 4 lines per sequence record:  
  1. Begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
  2. The raw sequence letters
  3. Begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
  4. The quality values for each sequence character. This line is required to contain the same number of symbols as letters in the sequence.  
  
A FASTQ file containing a single sequence will have a format like this:
```text
@SEQ_ID  
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT  
+  
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65  
```

The quality scores are denoded with ASCII characters. The byte representing quality runs from 0x21 (lowest quality; '!' in ASCII) to 0x7e (highest quality; '~' in ASCII). 

Here are the quality value characters in left-to-right increasing order of quality (ASCII):  
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

___  

For this problem open the [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file [Python_06.fastq](https://raw.githubusercontent.com/prog4biol/pfb2025/master/files/Python_06.fastq) and read each line to calculate and report:  
    - total number of lines  
    - total number of sequence IDs  
    - total number of characters  
    - total number of nucleotides  
    - average line length of all the lines  
    - average line length of the lines that contain sequences.     





5. Goal of this problem: generate a couple of gene list that are saved in files, add their contents to sets, and compare them. 

__Generate Gene Lists:__

You are going to generate gene lists and then compare the gene lists. We want to know what genes are both stem cell proliferation genes AND transcription factors. To do this you are going to use Biomart to 1) create a list of all the genes of one species 2) get genes that have been curated as being involved in stem cell proliferation ( tagged with the Gene Ontology term for stem cell proliferation ), 3) get genes that are transcription factors ( tagged with Gene ontology terms for transcription factors). Then you will write a script that using sets to compare the lists. 

_Get all genes:_

1. Go to [Ensembl Biomart](http://useast.ensembl.org/biomart/martview).
2. In dropdown box, select "Ensembl Genes 115"  (or most current version)
3. In dropdown box, select "Ferret Genes" 
4. On the left, click Attributes
5. Expand GENE:
6. Deselect "transcript stable ID", "Gene stable ID version", and "transcript stable ID version".
7. Click Results (top left)
8. Export all results to "File" "TSV" --> GO
9. Rename the file to "ferret_all_genes.tsv"

_In the same Ensembl window, follow the steps below to get genes that have been labeled with Gene Ontology term "stem cell proliferation". For extra information on stem cell proliferation, check out  [stem cell proliferation](http://purl.obolibrary.org/obo/GO_0072089)_

11. Click "Filters"
12. Under "Gene Ontology", check "Go term name" and enter "stem cell proliferation" (clear out any previous GO term names)
13. Click Results (top left)
14. Export all results to "File" "TSV" --> GO
15. Rename the file to "ferret_stemcellproliferation_genes.tsv"

_In the same Ensembl window, follow the steps below to get genes that have been labeled with Gene Ontology term "pigmentation". For extra information on pigmentation, check out [pigmentation](http://purl.obolibrary.org/obo/GO_0043473)_. Make sure that the previous Stem Cell Proliferation GO term is replaced with the Pigmintation GO term


16. Click "Filters"
17. Under "Gene Ontology", check "Go term name" and enter "pigmentation"
18. Click Results (top left)
19. Export all results to "File" "TSV" --> GO
10. Rename the file to "ferret_pigmentation_genes.tsv"


__Write a script that opens each of the three files and add the contents (which are the Gene stable IDs) to a Set. Create one Set per file.__
  - Then using set methods, indentify all the genes that **are not** cell proliferation genes.
  - And identify all genes **that are both** stem cell proliferation genes and pigment genes.
    
*Note* Make sure to NOT add the header to your set.  

__Now, let do it again with transciption factors.__

1. Go back to your Ensembl Biomart window
2. Deselect the "GO Term Name"
3. Select "GO Term Accession"
4. Enter these two accessions IDs which in most organisms will be all the transcription factors
   - GO:0006355 is "regulation of transcription, DNA-dependent”. 
   - GO:0003677 is "DNA binding"
5.  Click Results (top left)
6. Export all results to "File" "TSV" --> GO
7. Rename the file to "ferret_transcriptionFactors.tsv"

__Write a script that creates sets to find all the genes that are both transcription factors and are involved in stem cell proliferation.__  
  - Open these files: 1) the transcription factor gene list file and 2) the stem cell proliferation gene list file.
  - Add each to a Set, one Set per file
  - Identify all the genes that are transcription factors **AND** are involved in stem cell proliferation.

__Now compare gene lists on the command line with `comm` command. You might need to `sort` each file first.__
  - What genes are transcription factors involved in stem cell proliferation?
  - What genes are transcription factors involved in pigmentation?
  - What genes are involved in both stem cell proliferation and pigmentation?
  - Find out more about genes involved in stem cell proliferation and pigmentation by doing some searches on Ensembl to get the gene names and on Google to get more information about how they act in these processes.

     
Are you still committing your files as you go?

