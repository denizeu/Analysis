# Data Description

### Link and Data Info
[Genome sequences link] (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049&Completeness_s=complete&QualNum_i=50&Region_s=Africa&Region_s=Asia&Region_s=South%20America) This is a link to the genome sequences of SARS-CoV2 from Africa, Asia and South America. It was accessed on July 11, 2021.

### Further Data Info
The components that are included within the sequence headers consist of location, species, date, and what type of animal had the virus (human), as well as complete genomes.
After finding a link to these genomes, I downloaded it and inserted this data as a file into the data directory. 
I also created a file called ".gitignore" to make sure that git does not track the entire file and to keep from committing all of that data.

## Data Files Used within this Project

### Altered Data: refined_data.fasta
For my final project analysis, I had to cut down on how many sequences I used.
Because of the long run time of certain functions, I cut down on my data by choosing the first 12 in each time period (which I specify within BioinformaticsBISC195).
I chose the first 12 sequences that were dated at 2019, the first 12 dated at 2020 and the first 12 dated at 2021.
Next, I created a new fasta file for this and stored it in "data/refined_data.fasta" to be used in my functions.
Therefore, my data is a smaller subset of the above Genome Sequences Link.
I did this in order to maintain the integrity of my code, to be able to properly run and test all my functions, to avoid extreme wait times.
This data makes it easier to return the functions quicker and still gives clean data to answer each part of my analysis.

### Smaller Snippet: snippet.fasta
The smaller snippet of the very large genome database was what I used for intermediate tests before I created the altered "data/refined_data.fasta".
After running a few tests with this much smaller dataset, and after an office hours conversation with Prof. Bonham, I realized it would be much smarter to create a whole new dataset (the altered data in the excerpt above this one).
However, the dataset "data/snippet.fasta" was helpful in checking if my data was being parsed properly and being stored in the proper headers for each function I ran.
It consists of 4 genome datasets.
One data set I left with its full sequence while the others I cut down.
This was to test my uniqueKmers function and see if functions that might rely on size of sequence would play out if they were smaller.
This data file is not representative of true data but was manipulated to make sure my functions worked properly throughout my attempts at understanding the functions I was creating with such a large data file. 

### Working Data: datatry.fasta
This data was very useful in my function testing.
It consists of 3 very small datasets from each time period: 2019, 2020, 2021.
Therefore it was useful in visualizing the functions that separate sequences based on time periods.
It was especially helpful to the kmertime function because it gives a smaller set of returned vectors that we can count to verify that the next function, kmertimes (which counts the number of sequences in each) is correct.

### Another version of Snippet: Sorting Data: data_sorting.fasta
This data comes in handy for my function that sorts sequences by size (greater than 29500).
Since my refined data set that I am using for my project is extremely large (the sequences are long and there are many), this data set returns only a smaller number of sequences making it easy to see what the function does.
However, this data has a larger number of sequences in order to work for the sorting function which needs sequences greater than 29500.
This makes it much more clear how many sequences are greater than that boundary as it only gives a small number of sequences.

### Smaller Version of Snippet: datasort.fasta
This data is very similar to my snippet data. 
It is just a version of data_sorting and snippet that does not have extremely long sequences.
It takes a few sequences from only 3 headers: one version of each time period.
This is helpful for my min max function because the small data allows me to check whether the output is correct by counting on my own.

### Ex data sets: ex1.fasta & ex2.fasta
These example datasets are the ones that we have used in class.
They are very small collections of sequences and headers.
I used them within my project to show that certain functions will error with this data since it is not large enough or since it doesnt have the elements necessary to run my tests.
For example, I used it within my runtests for kmertimes and it gave an error because these example datasets do not have headers that define the year and therefore could not be sorted by time