# First Draft of Functions

## Finding the Unique Kmers
```julia 
@doc uniqueKmers
```
This is a function which takes the COVID data I compiled from NCBI and returns which kmers of a given length exist.
It takes 2 arguments: one is the sequence in which kmers appear (sequence) and the second is the length of the kmers (k).
In my function, the first key part is making sure that the function only takes bases within AGCT in order to have correct calculations.
Next, I identify the length of k and that it cannot be greater than how long the sequence is, to assure that the right number of letters is being computed.
Then I create a Dictionary to store the kmers in, so that I can easily identify them and identify their location within the sequence.
I also initialize an array "ret" which is extremely important because this is where all the unique kmers will be stored.
I push each unique kmer into the array so that the output represents how many kmer patterns of the length "k" exist within the sequences.
This is important because it can be used as a "distance metric".
This means that the sequences that are more similar are more likely to have a common ancestor.

## Distance Kmers: Comparing Unique Kmers
```julia 
@doc kmerdist
```
Kmerdist() is a function that returns how similar genomes are.
It is scored 0-1.
0 represents a kmer set that is completely identical, and 1 represents a kmer set with no similarity.
In this function, I made the argument composed of "set1" and "set2" so that two sets of kmers could be used to determine the distance between them.
Kmerdist() is composed of a return which uses a function to compute distance.
The math that is returned is this: the length of intersection of the two sets divided by the length of union.
This means that the function is comparing distances in this way: the intersection describes the set elements that are in both set1 and set2.
The union shows a set with all the elements in both sets.
Therefore, this mathematical function dividing the two and subtracting 1- that value is a simple distance fraction that shows the distance between the kmer sets. 
This is important because it can be used to see when genomes diverge and become different from one another.

## Time vs. UniqueKmers
```julia 
@doc kmertime
```
This function is my analysis of kmer compositions comparisons between time periods. 
It is important because it shows how the sequence of DNA changes as time goes on.
First, I attempted initialize 3 arrays for early, middle, and late.
Early represents 2019, middle represents 2020 and late represents 2021.
Then I began a for loop that sees if the header contains "2019", "2020" or "2021".
If so, I pushed the unique kmer count of that sequence to the time period array that it belongs to.
This is helpful to the analysis because it splits up the kmer comparisons for each time period to show how much the unique kmer compositions differ as time goes on, and the sequence can experience more changes/mutations to its code.

### Sorting Kmertime Data into a Vector
```julia 
@doc kmertimes
```
This function is created in order to make the process for making graphs easier.
It is meant to store the number of unique kmers per time period into a vector that can be used to create a bar graph to show the differences in numbers of unique kmers between 3 distinct time periods.
I used uniquekmers within this because it is illuminating to see the differences in time period grouped by number of unique kmers.

## Kmertime Bar Graph
```julia 
using Plots 
using BioinformaticsBISC195
histogram(kmertimes("data/genomes_CoV2.fasta"))
```
This is a Plots function.
It creates a bar graph with x and y
X stores time period: "early", "middle", "late".
Y stores number of unique kmers: 63, 74, 95.
I chose a bar graph because this will allow for differences in the number of unique kmers between the three periods to be clear. 
It is also a good way to visualize the actual counts through the "y" axis.

## Comparing Kmer Geo Locations: Turkey vs. Japan
```julia
@doc pairdist
```
This is the function I am trying to use to calculate pairwise kmer distance using a matrix.
My function is meant to create a matrix of sequences which lists the sequences stored in "data/refined_data.fasta".
Then, it is meant to calculate the distance between the kmers between two sequences going through the matrix.
However, unlike needleman wunch, it will only be calculated at one half of the matrix since distances between kmers will be the same no matter which way they are ordered.
This means that I will only calculate the distance between each pair a single time.
After I created the matrix of distances, I group them by early, middle and late.
This is done in order to convert this data into a boxplot which will represent the distance between kmers in early vs early, early vs. middle, etc.


## Pairwise Kmer Distance Box Plot
```julia
Plots.gr()
x= ["Turkey", "Japan"] #x-value is location: turkey or japan
y= [Tu, Ja] #y value holds the arrays of unique kmers
pie(x, y, title= "Number of Unique COVID-19 Kmers in Turkey vs. Japan")
```
I decided to use a piechart to compare Japan and Turkey's number of unique kmers because I believe this will highlight the differences in genomes.
By comparing the two side by side and as parts of a whole, it will be clear which country had much more unique kmers.
However, this is a first draft, like all my code.