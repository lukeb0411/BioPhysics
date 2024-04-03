# BioPhysics
Biophysics research primer (Carleton group)
Malachy Guzman, Arjendu Pattanayak


This is meant to serve as a background primer for Arjendu’s biophysics group research. This contains 1) suggested papers on biophysics, transcriptomics, development, information processing, 2) general resources on biophysics ideas and techniques, and 3) general resources or links to relevant biological ideas and techniques. 

Parts of section 3 are needed/recommended for a full grasp of what is going on in the papers in section 1. The videos in section 2 are fairly self contained assuming a bit of physics background.

Overall the best way to understand our projects is to think of them as being in dynamical systems biology. In one kind of project (let's call it type A) we use dynamical systems techniques and ideas to characterize biological dynamics; these are expected to be relatively close to spatio-temporal chaos problems in physics such as with photons etc, while the dynamics are complicated, these aren’t necessarily ‘complex biological systems’. 

In the second more complex set of projects (type B), we are actually studying ‘biological systems’ that are then changing in very interesting ways in time, where the information theoretic complexity of the system is a critical feature as well as a thermodynamic change of some sort. We hope to connect these two seemingly different aspects of the system behavior. 

These terms can be puzzling so always clarify as needed. 

The critical tool in all cases is some discrete information-theoretic measure, whether K-L entropy or PI entropy or Ordinal Pattern Subpopulation Vector Spaces to understand broader changes.


Relevant Papers

Always useful are those by our collaborator Jordi Garcia-Ojalvo: Publications - Dynamical Systems Biology lab (UPF) and https://scholar.google.com/citations?user=S76yu0YAAAAJ. 


Many of our present ideas stem from this review article by Jordi, and it has tons of well organized references: Mechanistic models of cell-fate transitions from single-cell data 
An older review article by Jordi which focuses on similar things but with less state of the art techniques: Towards a statistical mechanics of cell fate decisions 
Arjendu’s article on PI entropy, our main technique right now: Permutation entropy of indexed ensembles: quantifying thermalization dynamics 
Current work on Ordinal Pattern Subpopulation Vector Spaces
Jordi recommended this article on issues with techniques in single-cell genomics, which is very relevant for us.  The specious art of single-cell genomics 
Pre reading: transcriptomics in the Biology section below



Articles with promising public data:
Drosophila data (not an article): https://www.ebi.ac.uk/gxa/experiments/E-GEOD-18068/Results 
Single-cell developmental RNA-seq in frogs: The dynamics of gene expression in vertebrate embryogenesis at single-cell resolution 
Bulk and single-cell analysis of 110 genes in chickens, thermo & differentiation Single-Cell-Based Analysis Highlights a Surge in Cell-to-Cell Molecular Variability Preceding Irreversible Commitment in a Differentiation Process 
Single-cell developmental trajectories in mouse T cells: Wishbone identifies bifurcating developmental trajectories from single-cell data 


Selected articles (to be annotated)
Thermodynamics of Biological Processes (Rob Phillips, Caltech)
Biological Maxwell's Demon
Stem Cell Differentiation as a Non-Markov Stochastic Process (may have data for us)
Dynamics of Drosophila endoderm specification (Wieschaus and Shvartsman, Princeton)
Are Biological Systems Poised at Criticality? (Bialek, Princeton)
Entropy as a Geometrical Source of Information in Biological Organizations  


A list of subfields to look out for, which can mean the same or different things:
Biophysics
Systems biology
Computational biology 
Quantitative biology
Mathematical biology

Be careful not to confuse computational biology which is computational in the physics sense with computational biology which is basically just statistics (e.g., a lot of genomics). There’s nothing wrong with this other version of comp bio, but it’s distinct from our research angle in my opinion.


Section 2: Broader background ideas biophysics 

There are plenty of good biophysics textbooks and overviews, Arjendu has some copies of these. I have one which is a big, high level overview of topics and methods across biophysics, shared in the 0-Biophys folder: Physics of Life.pdf

A lecture series and a separate video on transcription by Bill Bialek, a major figure in biophysics: 
Statistical physics of biological systems: From molecules to minds - 1 of 4 
Statistical physics of biological systems: From molecules to minds - 2 of 4 
Statistical physics of biological systems: From molecules to minds - 3 of 4 
William Bialek   4 of 4 

Separate Bialek video on biophysics of transcription:
William Bialek (Princeton): Developing Unifying Theories for Biology 


‘Popular’ books: What is Life, Into the Warm, 



Section 3: Broader background ideas biology 

You should get textbooks in biology or genetics (recommendations available) including from the library to read in more detail on all of this. Any basic overview resources are good. Take bio courses!

 Our work currently spans three fields:
Genetics, genomics, and transcriptomics (or any -omics)
 In particular we’ve been looking at transcription, one step in the process of gene expression where DNA is transcribed into RNA. Read:
https://en.wikipedia.org/wiki/Genetics 
https://en.wikipedia.org/wiki/Transcription_(biology) 
https://en.wikipedia.org/wiki/Regulation_of_gene_expression 
The central dogma of molecular biology
Translation is the next step of gene expression, where the RNA is translated into protein, which after being chemically modified carries out different jobs in a cell. Right now our work is not concerned with translation, but it’s still good to know about.
Any “-omics” word refers to studying all of the present “thing”. So transcriptomics means measuring how much of each different transcript (the RNA transcribed from the DNA) is present in a sample. 
https://en.wikipedia.org/wiki/Transcriptome 
Our present work focuses on transcriptomics data, in particular RNA-seq
https://en.wikipedia.org/wiki/RNA-Seq 
Note the differences between bulk and single-cell sequencing.

Developmental biology
How an organism goes from a cell to a highly structured body. This is a very broad field, but our current work is focused on gene expression during development, which is a bit more specific. 
Two important parts of development that our work involves are:
Forming shapes: https://en.wikipedia.org/wiki/Morphogenesis 
Cell types: https://en.wikipedia.org/wiki/Cellular_differentiation 
Control over gene expression in time and space decides how this all goes.

Biophysics Proposal details how we’re presently thinking about this connection. So far we have mainly been studying fruit flies (Drosophila) because there’s a lot of public data. But we might move to frogs, or worms, or whoever has good data. 

TYPE A problem per description above: Neuroscience
(Luke and Sara). This is more of ‘dynamical systems’ in biology take on the dynamics of brain wave frequencies (measured by EEG) as they relate to a subject performing different tasks
https://en.wikipedia.org/wiki/Behavioral_neuroscience 
https://en.wikipedia.org/wiki/Electroencephalography 
We could however imagine looking at the developmental transcriptomics of the brain, tying all these topics together.


