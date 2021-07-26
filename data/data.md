# Data Description

### Link and Data Info
[Genome sequences link] (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049&Completeness_s=complete&QualNum_i=50&Region_s=Africa&Region_s=Asia&Region_s=South%20America) This is a link to the genome sequences of SARS-CoV2 from Africa, Asia and South America. It was accessed on July 11, 2021.

### Further Data Info: genomes_CoV2.fasta
The components that are included within the sequence headers consist of location, species, date, and what type of animal had the virus (human), as well as complete genomes.
I thought this was interesting to analyze because it contains a huge variety of locations and dates which would allow me for a very in-depth analysis.
After finding a link to these genomes, I downloaded it and inserted this data as a file into the data directory.
I named this file "genomes_CoV2.fasta". 
I also created a file called ".gitignore" to make sure that git does not track the entire file and to keep from committing all of that data.

## Data Files Used in this Project

### Altered Data: refined_data.fasta
This is the data I am using as a basis for my final project.
After trying to use my larger data file: "data/genomes_CoV2.fasta", I realized that it took much too long to run my sequences through that ginormous file.
Instead, I decided it would be more useful to cut down on the number of sequences and attempt to do an analysis on a smaller number of sequences.

**About this Dataset**
- I cut down on how many sequences I used by selecting a smaller portion of sequence data from my larger "data/genomes_CoV2.fasta".
- I did this by choosing the first 12 in each time period: the first 12 sequences that were dated at 2019, the first 12 dated at 2020 and the first 12 dated at 2021.
- This gives me a good number of sequences from each time period, and also is short enough that my code is runnable.
- Next, I created a new fasta file for this and stored it in "data/refined_data.fasta" so I could use these sequences within my functions.
- Therefore, my data is a smaller subset of the above Genome Sequences Link.
- I did this in order to maintain the integrity of my code, to be able to properly run and test all my functions, to avoid extreme wait times.
- This data makes it easier to return the functions quicker and still gives clean data to answer each part of my analysis.

### Smaller Snippet: snippet.fasta
The smaller snippet of the very large genome database was what I used for intermediate tests before I created the altered "data/refined_data.fasta".
After running a few tests with this much smaller dataset, and after an office hours conversation with Prof. Bonham, I realized it would be much smarter to create a whole new dataset (the altered data in the excerpt above this one).
However, the dataset "data/snippet.fasta" was helpful in checking if my data was being parsed properly and being stored in the proper headers for each function I ran.

**About this Dataset**
- It consists of 4 genome datasets.
- One data set I left with its full sequence while the others I cut down.
- Used to test my uniqueKmers function and see if functions that might rely on size of sequence would play out if they were smaller.
- It is *not* representative of true data but was manipulated to make sure my functions worked properly throughout my attempts at understanding the functions I was creating with such a large data file. 
- I stopped using this data after realizing I could better test my functions with the other files I manipulated (trial and error!).

### Working Data: datatry.fasta
This data was very useful in my function testing.

**About this Dataset**
- It consists of 3 very small datasets from each time period: 2019, 2020, 2021.
- It was useful in visualizing the functions that separate sequences based on time periods.
- It was especially helpful to the kmertime function because it gives a smaller set of returned vectors that we can count to verify that the next function, kmertimes (which counts the number of sequences in each) is correct.
- Allowed me to visualize and double check that my functions worked as I could count the number of values returned manually.

### Sorting Data: data_sorting.fasta
This data comes in handy for my function that sorts sequences by size (greater than 29800).
Since my refined data set that I am using for my project is extremely large (the sequences are long and there are many), this data set returns only a smaller number of sequences making it easy to see what the function does.

**About this Dataset**
- This data has a larger number of sequences in order to work for the sorting function which needs sequences greater than 29800.
- Makes it much more clear how many sequences are greater than that boundary as it only gives a small number of sequences.
- I basically just used this data file to check that my sortingseq() function worked before I inputed my full refined data file into the argument of the function.

### Smaller Version of Snippet: datasort.fasta
This data is very similar to my snippet data. 
It is just a version of data_sorting and snippet that does not have extremely long sequences.

**About this Dataset**
- It takes a few sequences from only 3 headers: one version of each time period.
- This is helpful for my minMax function because the small data allows me to check whether the output is correct by counting on my own.

### Example Data Sets: ex1.fasta & ex2.fasta
These example datasets are the ones that we have used in class.
They are very small collections of sequences and headers.
I used them within my project to show that certain functions will error with this data since it is not large enough or since it doesnt have the elements necessary to run my tests.

**About this Dataset**
- I used it within my runtests for kmertimes and it gave an error because these example datasets do not have headers that define the year and therefore could not be sorted by time.
- By using these data files, I realized I needed to compile better data to test my functions that relied on time period/date or on very long sequences.