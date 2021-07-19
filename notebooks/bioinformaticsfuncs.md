# First Draft of Functions

## Finding the Unique Kmers
```julia 
function uniqueKmers(sequence, k)
    for base in sequence 
        if !occursin(base, "ACGT")
            error("Invalid base $base encountered")
        end
    end
    1 <= k <= length(sequence) || error("k must be a positive integer less than the length of the sequence")
    kmers = Dict()    
    stopindex = length(sequence) - k + 1
    ret = []
    for i in 1:stopindex
        kmer= sequence[i:i+k-1]
        kmer = uppercase(kmer)
        push!(ret, kmer)
        if haskey(kmers, kmer) == true
            kmers[kmer] = kmers[kmer] + 1
        else
            kmers[kmer] = 1
        end
    end
    return Set(ret)
end
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
function kmerdist(set1, set2)
    return 1 - (length(intersect(set1, set2))/length(union(set1,set2)))
end
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
function kmertime(path)
    early= []
    middle= []
    late= [] #Initializing empty arrays for 3 time periods: early(2019), middle(2020), and late(2021)
    str= "" #Empty string used to hold the kmer patterns to be pushed into each time period
    for line in eachline(path)
        if header occursin("2019", path) #Findall occurrences of 2019 per header
            push!(early, uniqueKmers(sequence, k)) #If the header contains the date "2019", the kmer will be pushed into the "early" array, calling the uniqueKmer function to process how many unique kmers exist in the sequence
        end
        if header occursin("2020", path)
            push!(middle, uniqueKmers(sequence, k)) #If the header contains the date "2020", the kmer will be pushed into the "middle" array
        end
        if header occursin("2021", path)
            push!(late, uniqueKmers(sequence, k)) #If the header contains the date "2021", the kmer will be pushed into the "late" array
        end
    return early, middle, late
    end
end      
```
This function is my first draft of analysis of kmer compositions comparisons between time periods. 
It is important because it shows how the sequence of DNA changes as time goes on.
First, I attempted to initialize 3 arrays for early, middle, and late.
Early represents 2019, middle represents 2020 and late represents 2021.
Then I initiated a string to hold the different kmer patterns that would be found.
Next I began a for loop that sees if the header contains "2019", "2020" or "2021".
If so, I pushed the unique kmer count of that sequence to the time period array that it belongs to.
This is helpful to the analysis because it splits up the kmer comparisons for each time period to show how much the unique kmer compositions differ as time goes on, and the sequence can experience more changes/mutations to its code.
However, I need to figure out how to make the argument succinct with the ideal output as I am currently retrieving an error.

### Sorting Kmertime Data into a Vector
```julia 
function kmertimes(path)
    kmernumba= [] #array to store the unique kmers per each time period
        data = parse_fasta(path)
        for i in data[2]
             push!(kmernumba, uniqueKmers(sequence, k)) #for the data within the pos 2 of sequence data, the unique kmers are pushed to the array.
        end
        return kmernumba
    end
```
This function is created in order to make the process for making graphs easier.
It is meant to store the number of unique kmers per time period into a vector that can be taken and made into a histogram.
I use uniquekmers within this because it would be illuminating to see the differences in time period grouped by number of unique kmers.

## Kmertime Histogram
```julia 
using Plots
using BioinformaticsBISC195
histogram(kmertimes("data/genomes_CoV2.fasta")) # TODO: This throws an error (`sequence not defined`)
```
This is a Plots function.
It is supposed to create a histogram with x and y.
X is meant to be sorted by time period and Y by number of unique kmers.
I chose a histogram because this will allow for clear differences between the three periods to be seen. 
It is also a good way to visualize the actual counts through the "y" axis.

## Comparing Kmer Geo Locations: Turkey vs. Japan
```julia
function kmerloc(path)
Tu= []
Ja= []
for header in eachline(path)
    if header occursin("TUR", path)
        push!(Tu, uniqueKmers(sequence, k))
    end
    if header occursin("Japan", path)
        push!(Ja, uniqueKmers(sequence, k))
    end
        return Tu, Ja
    kmerdist(Tu, Ja)
    end
    return kmerdist
end
```
This is the function I am trying to use to calculate the unique kmers in Turkey and Japan.
My function is meant to push these unique kmers, if found in the location which contains the headers "TUR" (Turkey) or "Japan" into two arrays.
Next, I call the kmerdist function to calculate the difference between two sets of kmers, that of Japan and Turkey. 
I am trying to figure out how to properly separate the two.
I want to choose the biggest set from each and then calculate the different number of unique kmers that appear there but I am having trouble figuring out the code for this.
I am currently retrieving a "TypeError" which indicates I need to find a non-boolean way to rewrite this.
This function is important because it compares two locations and shows how different the number of unique kmers are.
This can help determine how much the sequence changes depending on the geographic location.
Japan is close to China, but Turkey is also relatively near so it will be interesting to see how different they are. 
I would have liked to try to compare them to China but it is hard for me to determine which data indicates Wuhan's data.

## Kmer Location Plot
```julia
Plots.gr()
x= ["Turkey", "Japan"] #x-value is location: turkey or japan
# TODO: Where's the code to generate `Tu` and `Ja`?
y= [Tu, Ja] #y value holds the arrays of unique kmers
pie(x, y, title= "Number of Unique COVID-19 Kmers in Turkey vs. Japan")
```
I decided to use a piechart to compare Japan and Turkey's number of unique kmers because I believe this will highlight the differences in genomes.
By comparing the two side by side and as parts of a whole, it will be clear which country had much more unique kmers.
However, this is a first draft, like all my code.