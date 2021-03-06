# My Analysis Plan

### 1. My Plan
I plan to analyze by calculating sequence lengths, GC content and kmer composition of these genomes from Asia, South America and Africa. 
I will also compare kmer compositions between different geographical locations. 
Then, I will use new functions to determine whether kmer compositions of CoV2 genomes from Africa and South America are more similar to each other than to those of Asia. 
From this analysis I would like to come to a conclusion about which geographic locations *if any* of the three have a more similar kmer composition than the other. 
Another sequence feature I'd like to analyze is how the kmer composition changed over time in each particular location. 
I want to determine whether later sequences have distinct kmer compositions in comparison to earlier recorded sequences. 

### 2. Functions
- I plan to use the previously created functions: **sequence(), gc_content(), kmercount()**. 
These will provide data for the sequence lengths, gc content and kmer count of these sequences. 
These functions have already been previously generated, and I will reuse old code by inputing my new COVID data. 
It will be interesting to calculate how gc content differs between the different geographic areas & to compare differences in content. 
Similarly the sequence length function will show how different the sequences are from each other.  

- I will also create new functions. **kmertime()** will be a function that divides kmer compositions into time periods. 
I might be able to do this by separating the sequences by time period first using a sorting function. 
(I will also use our previously created kmer functions and build on them. 
Then, from there I can specify which time period I am analyzing for each by initializing arrays to store all DNA sequences within each time period. 
The output will be kmer compositions stored in an array for each time period. 
I am also thinking of using kmercount function within this function, to utilize the code that created Dictionaries. 

- I will also make a function to divide the headers by geographic location, **kmerloc()**, in order to analyze the effects of time and location on kmer composition. 
I will split the locations as we did when searching using the startswith(header) code created in an earlier assignment, in order to separate the sequences by location (since locations are stored in the headers): Africa, South America and Asia. 
Then I will store the locations in an array for each. 
I will build on previous kmer functions to create a kmer counter for each location and insert it into the tuple. 
This would also be useful if I would like to do further analyses on it.
Then I will find a way to create a function that shows how different compositions are based on a certain threshhold to see how siginificantly different the different compositions are from one another (if this is possible). 
I will be sure to utilize normalizeDNA from the function bank we have created in order to make sure the DNA sequences are being read properly.

### 3. Plotting
- I intend to plot gc content, sequence length and kmer count by location. 
I will put these data onto a scatterplot with ***x*** being location and ***y*** being the respective analyzed data. 
With a scatterplot, patterns based on location will be clear. 
The scatterplots will bring clarity of how these factors differ and how spread out these differences are across different locations. 

- Next, to determine whether there is a difference in kmer composition over time I would use a bar graph with ***x*** being the time period, grouped by early, middle and late, and ***y*** being the kmer composition of each. 
A bar graph is ideal because it will show change between periods of time. 
I will also make each bar clearly represent the different times in order for more clarity and to differentiate changes between the locations (therefore there will be 3 separate lines: one per location). 
To do this I will use labels, colors and a legend if possible.

- When determining the distances between kmers in two sequences based on time period, I might use a boxplot chart to compare which time periods had a greater distance from each other in respect to kmer distance.
This will make it  easy to separate time period as well as account for key differences in kmers.
It will show specific numerical differences of kmer distance between later periods of time that may have experienced more variation in kmer composition as the virus spread further and earlier time periods where I hypothesize there will be a lesser distance.
This will show when kmer compositions differed most from each other across time periods and shed light on how the sequences changed from each other with time.