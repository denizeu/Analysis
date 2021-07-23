# Data Description

### Link and Data Info
[Genome sequences link] (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049&Completeness_s=complete&QualNum_i=50&Region_s=Africa&Region_s=Asia&Region_s=South%20America) This is a link to the genome sequences of SARS-CoV2 from Africa, Asia and South America. It was accessed on July 11, 2021.

### Further Data Info
The components that are included within the sequence headers consist of location, species, date, and what type of animal had the virus (human), as well as complete genomes.
After finding a link to these genomes, I downloaded it and inserted this data as a file into the data directory. 
I also created a file called ".gitignore" to make sure that git does not track the entire file and to keep from committing all of that data.

### Altered Data
For my final project analysis, I had to cut down on how many sequences I used.
Because of the long run time of certain functions, I cut down on my data by choosing the first 12 in each time period (which I specify within BioinformaticsBISC195).
I chose the first 12 sequences that were dated at 2019, the first 12 dated at 2020 and the first 12 dated at 2021.
Next, I created a new fasta file for this and stored it in "data/refined_data.fasta" to be used in my functions.
Therefore, my data is a smaller subset of the above Genome Sequences Link.
I did this in order to maintain the integrity of my code, to be able to properly run and test all my functions, to avoid extreme wait times.

### Smaller Snippet
The smaller snippet of the very large genome database was what I used for intermediate tests before I created the altered "data/refined_data.fasta".
After running a few tests with this much smaller dataset, and after an office hours conversation with Prof. Bonham, I realized it would be much smarter to create a whole new dataset (the altered data in the excerpt above this one).
However, the dataset "data/snippet.fasta" was helpful in checking if my data was being parsed properly and being stored in the proper headers for each function I ran.
It consists of 4 genome datasets.
One data set I left with its full sequence while the others I cut down.
This was to test uniqueKmers and see if functions that might rely on size of sequence would play out if they were smaller.
This data file is not representative of true data but was manipulated to make sure my functions worked properly throughout my attempts at understanding the functions I was creating with such a large data file. 
