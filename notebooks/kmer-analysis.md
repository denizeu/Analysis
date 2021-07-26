# Kmer Analyses
This is my notebook for my kmer analyses.

## Finding the Unique Kmers
This is my function to find unique kmers within a sequence:
```julia 
function uniqueKmers(sequence, k)
    sequence= normalizeDNA(sequence)
    1 <= k <= length(sequence) || error("k must be a positive integer less than the length of the sequence")
    kmers = Set{String}()    
    stopindex = length(sequence) - k + 1
    for i in 1:stopindex
        kmer= sequence[i:i+k-1]
        kmer = uppercase(kmer)
        push!(kmers, kmer)
    end
    return kmers
end
```
This function which takes a sequence from the COVID data I compiled from NCBI as well as a kmer length `k` and returns the kmers of that given length that exist within the specified sequences.
- It takes 2 arguments: 
    - 1: the sequence in which kmers appear (sequence)
    - 2: the length of the kmers (k).
- The first key part is making sure that the function *only* takes bases within AGCT, and converts ambiguous bases to `N` by calling `normalizeDNA()`.
- Then I initialize `kmers` to be an empty Set of strings.
- I push each unique kmer into the empty Set of `kmers` so that the output is the many kmer patterns of the length "k" that exist within the sequences entered.
- This is important because it can be used as a "distance metric".
    - This means that the sequences that are more similar are more likely to have a common ancestor.

An example of the outputs that I got when inputing sequences and `k` values is this:

```
uniqueKmers("ATGCGATG", 4)
    Set{Any} with 5 elements:
    "TGCG"
    "ATGC"
    "GATG"
    "CGAT"
    "GCGA"
```
- This is very interesting because it shows that there is *quite a variety* of unique kmers within this single sequences.
- Although the sequence seems small, this function clearly indicates that there are many different kmer patterns which exist.
- Therefore, this function is incredibly useful as it quickly designates a method to seeing patterns that you would not automatically catch.
- With this function, it is easier to see how sequences can differ very much, even if they don't seem to at first glance.

This function was difficult for me at first, because it was not easy for me to figure out how to locate the unique kmers.
However, once I understood what the stop index must be, it became more clear to me that the way to recognize kmers would be in such a similar manner.

## Distance between Kmers
This is my function to compute the distance between two sets of kmers:
```julia 
function kmerdist(set1, set2)
    return 1 - (length(intersect(set1, set2))/length(union(set1,set2))) #length of intersection of the sets divided by length of the union
end
```
- `Kmerdist()` is a function that returns how similar genomes are.
- It is scored 0-1.
    - 0 represents a kmer set that is completely identical
    - 1 represents a kmer set with no similarity.
- In this function, I made the argument composed of "set1" and "set2" so that two sets of unique kmers could be used to determine the distance between them.
- `Kmerdist()` is composed of a returned value:
    - The length of intersection of the two sets divided by the length of the union of the two sets.
    - This means that the function is comparing distances in this way: the intersection describes the set elements that are in both set1 and set2.
    - The union shows a set with all the elements in both sets.
- Therefore, this mathematical function dividing the two and subtracting 1- that value is a simple distance fraction that shows the distance between the kmer sets. 
- This is important because it can be used to see when genomes diverge and become different from one another.

When I tested this function out with this, I could clearly understand the metrics used for this test:
```
kmerdist(uniqueKmers("AATA", 4), uniqueKmers("CGGCCCG", 4))
    1.0
```
Since these two sequences have absolutely *nothing* in common, they are scored as 1.0.

The importance of this scoring becomes obvious when comparing two sequences with *some* similarity:
```
kmerdist(uniqueKmers("GCGCAT",2), uniqueKmers("ATAT",2))
    0.8
```
These two functions share one **A** and one **T** therefore the score is close to 1.0 but not *quite* there.

Through these examples, I was able to understand the importance of this function:
- It computes quickly, an important neccessity of comparison for two large sequences.
- Simply by inputing the sets, one is able to understand how similar or different two sequences are.
This is important for analyses because it allows researchers to understand points of divergence between sequences and helps us understand how similar or different kmer sets can be.

## Time vs. UniqueKmers
This is my function used to group unique kmers by time period:
```julia 
function kmertime(headers, sequences, k=3) #whatto set k to?
    early= []
    middle= []
    late= [] #Initializing empty arrays for 3 time periods: early(2019), middle(2020), and late(2021)
    for i in 1:length(sequences)
        if occursin("2019", headers[i]) #Findall occurrences of 2019 per header
            push!(early, uniqueKmers(sequences[i], k)) #If the header contains the date "2019", the kmer will be pushed into the "early" array, calling the uniqueKmer function to process how many unique kmers exist in the sequence
        end
        if occursin("2020", headers[i])
            push!(middle, uniqueKmers(sequences[i], k)) #If the header contains the date "2020", the kmer will be pushed into the "middle" array
        end
        if occursin("2021", headers[i])
            push!(late, uniqueKmers(sequences[i], k)) #If the header contains the date "2021", the kmer will be pushed into the "late" array
        end
    end
    return (early, middle, late)
end
```
This function is my analysis of kmer compositions comparisons between time periods. 
It is important because it shows how the DNA sequence can change as time goes on.
- First, I initialized 3 arrays for early, middle, and late.
    - With this notation, early represents 2019, middle represents 2020 and late represents 2021.
- Then I began a for loop that sees if the header contains "2019", "2020" or "2021".
    - If so, I pushed the unique kmer count of that sequence to the time period array that it belongs to.
    -  This is helpful to the analysis because it splits up the kmer comparisons for each time period to show how much the unique kmer compositions differ as time goes on, and the sequence can experience more changes/mutations to its code.

I had trouble with this code because I could not figure out how to use the `for` loop.
However, I then figured out that it should go through the *length of the sequence* and after that I was able to piece together the skeleton of this function.

- After all this, I had to figure out what to use for my `k`.
- This took many trial and errors in which I chose values ranging from 7 to 20.
However, I noticed that the differences that stood out most to me was when I inputed a `k` of 3.
- With this input, for my data seemed to vary fairly clearly.
- This kmer value is fairly low, but it shows differences in unique kmers across time periods, which is why I chose to use it.

### Using my data for this function
When I used my refined data for this function, I got a very large set of values since it is evaluating a small `k` value for a large number of sequences.
- This is good, because there are many points of comparison.
- By breaking up my data into time periods it makes it easier to compare in graph form.
- I realized why we need graphs at this point.
    - When you really want to see changes with big data, using a plot can break down the large data into clear points.

### Sorting Kmertime Data into a Vector
My function for preparing kmers for plotting: 
```julia 
function kmertimes(headers, sequences, k=3)
    early, middle, late = kmertime(headers, sequences, k) #calling kmertime to separate the sequences by time period: 2019, 2020, 2021
    early_kmers = length(union(early...)) #calculating length of union of the sets within each time period
    middle_kmers = length(union(middle...))
    late_kmers = length(union(late...))
    return(early_kmers, middle_kmers, late_kmers) #returns the number of unique kmers in each time period
end
```
This function is created in order to make the process for making graphs easier.
It is also used in order for me to see how different my data for each time period is.
- It is meant to store the number of unique kmers per time period into a vector.
    - This can be used to create a bar graph to show the differences in numbers of unique kmers between 3 distinct time periods.
I used uniquekmers within this because it is illuminating to see the differences in time period grouped by number of unique kmers.

## Kmertime Bar Graph
This is my Time Period vs. Number of Unique Kmers graphing function:
```julia 
using Plots 
using BioinformaticsBISC195
bar(["early" "middle" "late"],
           [64 73 95],
           labels = ["early" "middle" "late"],
           label = "Number of Unique Kmers",
           title = "Time Period vs. Number of Unique Kmers",
           xlabel = "Time Period",
           ylabel = "Number of Unique Kmers",
           color = [:steelblue :pink :lavender],
           bg= "beige",
           legend = :topleft
```
This is a **bar graph** function using Plots!
- It creates a bar graph with `x` and `y`.
    - `X` stores time period: `early`/2019, `middle`/2020, `late`/2021.
    - `Y` stores number of unique kmers: 63, 74, 95.
- This graph clearly shows differences in the number of unique kmers between the three periods to be clear. 
- It is a perfect visualization for an `x`-variable that is categorical and a `y`-value that is numerical.
- It's a good way to visualize the actual counts through the `y` axis.
- With this function, I was able to mess around with various layout formats.

However, I found it *extremely* difficult to compile this bar graph.
It was difficult to put the notation in the correct places.
At first I included commas in the `[]`'s which led to many errors.
For example, with my labels I had written them out like this:
```
labels = ["early" "middle" "late"]
```
and also could not figure out how to differentiate between x and y labels!
However after re-reading the notation many times, I realized that the there was a difference between the word label*s* and label..
After this I was able to compile a clear graph that separated my kmer counts by time period.

I found this graph to be incredibly interesting:
- The graph showed that as time went on, the number of unique kmer counts increased.
- Although these numbers were fairly small there was interestingly almost an increase of 10 unique kmers between early and middle and then a jump of about 20 from middle to late.
- To me, this indicates that the sequences in my data changed **alot** during the last time period in order to incur such a variation in such a small `k` value.
- To me, it was all the more interesting that I manipulated such a small value for unique kmers because it indicates that there is large difference with such a small `k` pattern.

## Comparing the Distance between Unique Kmers in Sequences
This is my function to compare distances between unique kmers:
```julia
using BioinformaticsBISC195
function pairdist(path) 
    h = kmertime(parse_fasta(path)[1], parse_fasta(path)[2], 3)[1] #takes kmertime of headers and sequences and stores them with these letters
    j =  kmertime(parse_fasta(path)[1], parse_fasta(path)[2], 3)[2]
    k =  kmertime(parse_fasta(path)[1], parse_fasta(path)[2], 3)[3]
    mesh = vcat(h,j,k) #concatenates h, j, k
    ret = zeros(36, 36) #initiates empty matrix
    for i in 1:36 
        for j in 1:36
            i <= j && continue #skips upper triangle, focuses bottom half of matrix
            d = kmerdist(mesh[i], mesh[j]) #finds distance between the kmers
            ret[i, j] = d
        end
    end
    return ret #returns the matrix with the added distance between kmers of sequences
end
```
This is the function I created to calculate pairwise kmer distance using a matrix.
- My function is meant to create a matrix of sequences which takes the sequences stored in "data/refined_data.fasta".
- Then, it calculates the distance between the kmers between two sequences going through the matrix.

When  I began, I had **alot** of trouble with this function:
- This was because I was attempting to create a matrix from values I had retrieved from previous functions.
- I was unsure how to use the previous data within this function.
- My first thought was to create arrays and try to pull from each array into a matrix.
    - That turned out to be complicated and time consuming.
- I tried to create a skeleton function imitating out previous Needleman Wunch assignment as a reference.
    - However, unlike Needleman Wunch, this would only be calculated at one half of the matrix since distances between kmers will be the same no matter which way they are ordered.

- With this knowledge, I realized that I would only calculate the distance between each pair a *single* time, and so would need to include notation of this in the function.
- With some office hours guidance, I discovered the use of `i <= j && continue` this indicates that the function should retrieve the values only if `i` is less than or equal to `j` in order to only retrieve data from the lower triangle and skipping the upper triangle since the values would be equivalent.
- I also realized that I could simply concatenate my previous data into a matrix.
    - That is the reason for the: `h`, `j`, and `k`'s defined above.
    - By defining these variables I am able to call `kmertime` and `parse_fasta` on them and then use the parsed, grouped, unique kmers from this and concatenate it with `vcat`.
After this I initialized an empty matrix of size 36 which allowed me to input my concatenated values into this empty matrix.

When I used this matrix function on my data, I was given a very large matrix output with very small distance calculations.
- This function was clearly very important because it performed tedious and complex calculations of kmer distances on many sequences.
- It retrieves very informational values and allows us to see how unique and similar/different all of my data's sequences are.
- My matrix was obviously very large as it is 36x36 and so this matrix would have taken a **LOT** of time to do manually.
- Also, it can be transformed into a graph which would allow for even clearer differentiation between sequence data of genomes.


### Sorting Distances
This is my function to sort the distances from my matrix back to time periods:
```julia
function distsort(mat) #takes matrix as an argument
    early_dist = [] #3 separate one dimensional arrays initialized for each time period
    middle_dist = []
    late_dist = []
    for i in 1:12
        for j in 1:12 #goes through the early time periods (1:12 are the first 12 in my refined_data)
            push!(early_dist, mat[i,j]) 
        end
    end
    for i in 13:24 #goes thru the middle time periods
        for j in 13:24
            push!(middle_dist, mat[i,j])
        end
    end
    for i in 25:36 #goes thru the late time periods
        for j in 25:36
            push!(late_dist, mat[i,j])
        end
    end
    return (early_dist, middle_dist, late_dist)
end
```
I used this function in order to create a Tuple of vectors that I could manipulate and turn into a graph in order to visualize the difference in unique kmers between time periods. 

This function was fairly confusing.
- This is because it was hard to realize *where* to separate the data.
- It took thinking through to realize that I would indice through 36 as only half of the matrix is being used for pairing calculations.
    - Then, I realized that the smart way to do this is by going through the sections by the numbers that they fall into:
`for j in 1:12`, `for i in 13:24`, `for i in 25:36`.
- It was also difficult to realize when to use `i` or `j` as the columns vs rows were not completely visibile within the matrix.

In order to figure this function out, I relied on looking at specific points in the matrix by doing something like:
```
distsort(pairdist("data/refined_data.fasta"))[3][141]
0.16883116883116878
```
in order to find particular calculations of use.

- As you can see, calling `[3][141]` gives a clean number rather than a huge set of numbers divided by `...` in which there are invisible numbers.
- Based on this calculation you can see that there is not much difference in the sequences I analyzed.
- This is a drawback of using a small `k` value; because there is such a small number of bases within this kmer pattern, there will be many sequences that share the same pattern of bases.
- Although it does provide interesting analysis for number of unique kmers, it does appear that these distances are very close to `0.0` which represents identical sets.
- However, this is informative as it still shows that there are a variety of differences even with a `k` of 3.


## Pairwise Kmer Distance Box Plot
```julia
using Plots
using StatsPlots
early_dist= distsort(pairdist("data/refined_data.fasta"))[1]
middle_dist =distsort(pairdist("data/refined_data.fasta"))[2]
late_dist = distsort(pairdist("data/refined_data.fasta"))[3]
boxplot(["early", "middle", "late"], [[early_dist], [middle_dist], [late_dist]], title="Distances vs Time Periods BoxPlot", xlabel="Distances", ylabel="Time Periods")  
```
I decided to use a **boxplot** to compare the number of unique kmers per time period. 
Although the functionality of this is wobbly, the message is to compare the differences between each time period with one another through **stacked box plots**.
- In this way, we would be able to visualize the differences in unique kmers between each plot in particular and also altogether.
- This is important because it will highlight the differences in between time periods and will also allow to compare the differences with particular points for means, median and range values.
- With boxplots, each time period will be compared to the two other time periods.
- This would show much more variation than any of my graphs above as it would account for the differences between periods as well as altogether.
    - I hypothesize great difference between the later times and earlier times when compared with one another.
By comparing the time periods side by side *and* as parts of a whole, it will be clear which period had much more unique kmers compared to the others, and the differences in between.