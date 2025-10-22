# Programming For Biology 2025

__Instructors__  
Simon Prochnik  
Sofia Robb   
Jessen Bredeson  
Eric Ross  

# Workshops and Overviews

## Pipelines -- Simon Prochnik
In modern bioinformatics research, designing and managing data‑analysis pipelines is key to achieving reproducibility, scalability, and maintainability. This lecture covers pipeline engineering in computational biology: how to structure, automate, and maintain workflows that process biological data from raw inputs to interpretable results.

__Audience__   
This lecture is intended for students, researchers, or bioinformaticians who:
- perform computational analyses of biological data (for example sequence data, “omics” data, or variant data) and want to design workflows that are robust and reproducible;
- are comfortable writing code or scripts (in Python, R, or shell), and using command‑line tools, but may not yet have formal experience building full end‑to‑end pipelines or workflow systems;
- want to move beyond one‑off analyses toward reusable, well‑structured, traceable pipelines.

You do not need to be an expert software engineer, but you should be familiar with basic programming, version control, and the concept of modular analysis.

__Goals__  
By the end of this lecture, you will be able to:
- explain the key concepts behind pipeline‑based analysis: modularity, reproducibility, traceability, and automation;
- design a pipeline workflow that connects raw data through intermediate steps to final outputs, with clear dependencies and execution logic;
- implement a basic bioinformatics pipeline, execute it, monitor it, and maintain it;

__Scope__  
In this session we will focus on:
- structuring pipeline components: defining the steps, specifying inputs and outputs, sequencing tasks, tracking dependencies;
- integrating pipelines within reproducible research practices: parameterisation and modularisation;
- pitfalls in pipeline design and how to avoid them.

We will not attempt to cover every workflow engine in existence, nor delve deeply into low‑level system tuning (e.g., HPC scheduling beyond basic integration). Instead, our focus is on practical, foundational pipeline practices that you can apply directly in your computational biology research.

__Dive in__  
  - [Lecture](../lectures/bioinfPipesLectureSimon.md)
  - Problem set: None


## Package Managers  -- Jessen Bredeson
In modern computational biology workflows, managing software dependencies and isolated environments is critical for reproducibility, portability, and efficiency. This workshop covers mamba — a lightweight, high‑performance package and environment manager — and shows how to integrate it into your research projects.

__Audience__  
This workshop is intended for students and researchers in the life sciences who:
- already use or plan to use Python or R in their analyses, and 
- have a basic familiarity with command‑line tools (for example, creating virtual environments or installing packages) but may not yet have deep experience with environment managers like mamba.

You do not need to be an expert in package management or system administration.

__Goals__  
By the end of this session, you will be able to:
- explain what mamba is, and how it differs from other environment tools (such as conda).
- install mamba and create isolated environments for your analyses.
- manage package installations and environment activation using `micromamba`.
- integrate mamba into reproducible workflows (for example in containers or shared projects).
- recognise when mamba is the right tool for your workflow

__Scope__  
We will focus on:
- installing `micromamba` on your system
- creating, removing, and activating environments
- installing packages
- `micromamba` channel configuration and cache cleaning

We will not cover every advanced feature of mamba, or dive deeply into advanced packaging internals; our focus is on practical, reproducible use of mamba for research workflows. In the following sections, we will start by giving you a quick introduction to how mamba works, then walk through a live demonstration of environment setup

__Let’s get started__  
  - [Lecture](../lectures/mamba.md)
  - Problem set: None


## BioPython -- Sofia Robb
In modern bioinformatics workflows, automating the processing, analysis, and interpretation of biological sequence data is essential for reproducible, scalable research. This workshop introduces Biopython — a comprehensive open‑source Python library tailored for computational biology and bioinformatics tasks — and shows how to leverage it for practical, real‑world biological data processing.

__Audience__
This material is designed for researchers, graduate students, or advanced undergraduates who:
- work with biological sequence data (DNA, RNA, protein) or structure data (PDB) and wish to automate analyses;
- already have experience with Python programming (variables, loops, functions, modules) but may not have used domain‑specific bioinformatics libraries;
- seek to build reproducible and maintainable code, rather than ad‑hoc scripts.

You do not need to be a software engineer, nor an expert in algorithm design, but you should feel comfortable writing and running Python code.

__Goals__
By the end of this workshop, you will be able to:
- explain what Biopython is and what kinds of bioinformatics tasks it supports;
- install and import Biopython into your workflow;
- use Biopython tools to read, parse, manipulate, and write common biological data formats (e.g., FASTA, GenBank, PDB);
- apply Biopython in practical workflows (e.g., sequence operations, alignments, structure parsing, file conversions) relevant to your research;
- integrate these capabilities into reproducible code, increasing the reliability, transparency, and shareability of your analyses.

__Scope__
In this session we will focus on:
- setting up the development environment for Biopython (installation and basic import);
- core modules and functions: working with sequence objects, reading/writing file formats, accessing biological databases programmatically;
- use‑cases common in computational biology workflows (e.g., parsing FASTA/GenBank, processing PDB files, sequence manipulation and simple downstream analyses);
- best practices for writing clear, maintainable scripts that use Biopython, including commenting, modularization, and reproducibility.

We will not exhaustively cover every Biopython module, nor dive into highly advanced topics like custom extension of Biopython internals, deep algorithmic optimizations, or extensive performance tuning. The emphasis is on getting you productive and confident with practical uses of Biopython in your research context.

__To get started, follow the links below__  
  - [Lecture](BioPython/biopython.md)
  - [Workshop](BioPython)


## Pytest Testing  -- Jessen Bredeson
Software in computational biology and data‑science workflows must not only run, but run correctly, reliably, and reproducibly. This lecture covers testing — the methods and practices for systematically verifying that code does what it is supposed to do, and continues to do so as it evolves.

__Audience__
This material is intended for students and researchers who:
- are developing scripts, packages, or pipelines Python as part of their scientific work;
- have basic programming experience (defining functions, using modules, working with version control) but may not yet have adopted disciplined testing practices;
- value reproducibility, maintainability, and collaboration in code development.
You do not need to be an expert in software engineering, but familiarity with running tools on the command line, and debugging will help.

__Goals__
By the end of this lecture, you will be able to:
- explain why testing matters in scientific software (e.g., correctness, regressions, reproducibility);
- write and run automated tests for your code;
- integrate testing into your development workflow using AI
- adopt best practices for organizing, writing, and maintaining tests to support long‑term project health.

__Scope__
In this session we will focus on:
- the fundamentals of designing effective tests (what to test, how to structure tests, how to pick meaningful test cases);
- practical Python tools and frameworks to write and execute tests;
- pitfalls and common anti‑patterns

We will not cover exhaustive detail on every testing framework nor deep dive into all advanced topics (such as property‑based testing, mutation testing, or highly optimized test harnesses). Instead, our focus is practical, actionable testing guidance tailored to research‑oriented code.

__Let's take a closer look__ 
  - [Lecture](../lectures/testing.md)
  - [Problem Set](../problemsets/testing_problemset.md)


## Ethics and Responsibility in Bioinformatics -- Simon Prochnik
  - [Lecture]()


## Coding with AI -- Simon
  - [Lecture]()
  - [Workshop]()


## RNAseq -- Brian Haas
  - [Lecture](https://github.com/trinityrnaseq/CSHLProgForBio/blob/main/rnaseq_slides_PFB2023.pdf)
  - [Workshop](RNAseq)


## Protein Homology -- Bill Pearson
  - [Lecture](Sequence_homology/cshl_pfb_25a.pdf)
  - [Lecture](Sequence_homology/cshl_pfb_25b.pdf)
  - [Workshop](Sequence_homology)
