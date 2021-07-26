# Sequence Calculations
This notebook contains the calculations I used for sequence analysis.

## Mean & Standard Deviation: Lengths and Counts
This is my first function:
```julia 
using Statistics
using BioinformaticsBISC195
function lengthcount(path)
    lengths= []
    counts= []
    data= parse_fasta(path)
    for i in data[2]
        push!(lengths, length(i))
        gs = count(==('G'), i)
        cs = count(==('C'), i)
        GCcontent= (gs + cs)/length(i)
        push!(counts, GCcontent)
    end
    return (mean(lengths), std(lengths), mean(counts), std(counts))
end
```
This code is calculating the *mean* and *standard deviation* of the **gc content** and **lengths** of my coronavirus genomes.
This function is important because it helps identify the different lengths and counts between sequences within a data file.

The output of this function which looks something like:
```
(13.0, 5.656854249492381, 0.46405228758169936, 0.34199935822094457)
```
As you can see, this data shows the difference in sequences within a larger dataset:
- This is a quick and easy way to see whether certain sequences are much longer than each other, or vary much from the mean number of sequences.
- Helpful in understanding how much sequences change over time, or to further analyze to see what the affects are of longer/shorter sequences.

It is calculating the means of the lengths by initializing a one dimensional array for lengths and counts.
- The `for` loop is created in which the length of each position is pushed to the length array.
- G's and C's are identified and the GC content is determined by dividing number of g's and c's by the length of each sequence.
- The GC content is then pushed to the counts array.
- The counts of GC are important because they show the relative composition of each sequence, yet another way to have a quick yet interesting analysis.

I had some trouble creating this function at first because I was not sure how to calculate the counts but I was able to work through it to come up with this function, notably this, which took some time to figure out:
```
gs = count(==('G'), i)
        cs = count(==('C'), i)
```

## Minimum and Maximum Lengths
This is my minimum and maximum function:
```julia 
function minMax(path)
    ret = parse_fasta(path)[2]
    retmin = length(minimum(ret))
    retmax = length(maximum(ret))
    return (retmin, retmax)
end
```
This calculates the minimum and maximum sequence lengths in the dataset.
1. Ret calls the parse function on the path, and with the location of `[2]`, to indicate that we are finding the length of the sequences and not of the headers.
2. The length of the minimum and the length
of the maximum of these sequences is taken.
The minimum and maximum are outputs.
This function is important because it helps us see the range of sequence data, to see how great variety is.

### Running minMax on my Refined Data

I ran the minMax function on my refined data set, with the 12 sequences from each time period. 

It worked like this:

```
minMax("data/refined_data.fasta")
(29833, 29881)
```

This is extremely interesting because it indicates that there is *some variation* in sequence length.
- The minimum value is 29833, which is 48 lower than the maximum sequence length of 29881.
- The code I created for this function was very simple and yet it returns two values that are extremely useful and important in sequence analyses.
- **However**, although these lengths are different, considering the magnitude of these lengths, 48 is not extremely large.
- Therefore we may consider that my refined data may not have enough variation and therefore is likely not very representative of COVID genomes.
- Despite this, the code is powerful as it tells us a very interesting difference in sequences with little in-between data or calculations necessary.

## Sequence Lengths
```julia 
function seqlength(path)
    lengths = []
    data = parse_fasta(path)
    for i in data[2]
         push!(lengths, length(i))
    end
    return lengths
end
```

This function returns the sequence lengths of each DNA sequence within a dataset. 
It creates an *array* in which the lengths can be stored.

- To clean up the data and separate them into header and sequences, the `parse_fasta()` function is called.
- This is important because it allows the two to be separated, and for further functions to be worked on the sequences independent of headers.
- I next initiated a `for` loop that goes through the data stored at position 2 by using `[2]`.
- This position represents the sequence data and assures that the length function is being called onto the sequences and not the header information.
- The length of each position within the data is pushed into the lengths array and then returned.
- This is important because the length of sequences can be compared for differences or used within various functions to determine important changes between genomes.

### Using the Sequence Length Function on my Data
When I ran this function for my refined data set I got a very large output:
```
seqlength("data/refined_data.fasta")
36-element Vector{Any}:
 29816
 29807
 29858
 29858
 29774
 29807
 29899
     ⋮
 29848
 29834
 29834
 29867
 29836
 29837
 29836
 ```
 Therefore, in order to properly test my function and visually see if it works, I used a smaller data set such as my `datatry.fasta` dataset.

 With this dataset as the argument, I got a much smaller output:
 ```
seqlength("data/datatry.fasta")
3-element Vector{Any}:
  3
 64
  4
```
- With this output, I could clearly see that the function worked well for my data as well since the counts matched. 
- Because of this, creating this function was difficult: it was not easy to check whether my lengths were correct.
- However, this is an interesting analysis because it allows me to see how different the lengths are as the time periods increase. 
- Based on this output, we see some examples of sequence lengths increasing, and some examples of decreasing sequence lengths as the time goes on.
- This code is useful because utilizes a path and is able to recognize the lengths of sequences by referring to these with `i`.


## Histogram of Sequence Lengths
```julia 
using Plots
data= seqlength("data/refined_data.fasta");
histogram(data, bins=10, label= "sequence lengths", xlabel= "Sequence Lengths", ylabel= "Number of Sequences")
```
Using the function for sequence lengths above and by initiating the Plots package, this creates a histogram that shows the sequence lengths as the `x`-axis and number of occurences on the `y`-axis.
This graph is a clear visualization of what the majority length of the sequences is.
- From the graph, it seems that there are significantly more occurrences within the sequence length of around 29825.
- This can be used to determine where most of the data exists.
- It also shows that, most of the genomes are above the length of 29500.
- This is important in comparing sequence sizes and showing which are more/less similar to one another.

This histogram was difficult to create because it was my first time creating them.
I had to figure out how to properly graph the dataset as well as use the seqlength function on that data.
It took me a while to realize that I could simply use: 
```data= seqlength("data/refined_data.fasta");``` 
to perform the function on my project's data set.
Although this is very obvious, it was very important in creating a succinct graph which properly calculated the sequence lengths necessary for the histogram.

## Getting Rid of 5% of Sequences
```julia 
function sortingseq(path)
    sequences = seqlength(path)
    headers = parse_fasta(path)[1]
    check = findall(x->x<29800, sequences)
    deleteat!(sequences, check)
    deleteat!(headers, check)
    return (sequences, headers)
end
```
After plotting the sequence lengths of my data, I realized that it would be very interesting to cut out the last 5% of my sequence lengths in order to have a magnified histogram.
Based on my graph, I decided that 5% of my sequences likely fall somewhere around 29778.
Therefore, I decided it would be interesting to get rid of that 5% by sorting out the sequences *below* that length value.
This function is a sorting function which cuts down the sequence length data to only show the relevant basepair lengths above 29800.
This function returns two vectors stored at positions `[1]` and `[2]`.
In `[1]` exists the sequences that are below the length 29800 and in `[2]`, exists the headers that correspond to the sequences that are below the sequence length 29800.

### Sorting my Refined Data
When I ran the function above on my data I got:
```
32-element Vector{Any}:
 29816
 29807
 29858
 29858
 29807
 29899
 29889
     ⋮
 29848
 29834
 29834
 29867
 29836
 29837
 29836
```
Based on the histogram and the function above I learned very much about my data.
- I learned that: when getting rid of around 5%, lots of the data was toward the larger end of the scale rather than at the beginning.
    - This shows that the sequences that I am analyzing fall around the larger lengths rather than smaller lengths.
    - It also shows how important differences in sequence size are in comparing genomes, as well as how helpful graphing these differences is.
- Since this function only returns a 32-element vector, it is difficult to differentiate the number of sequences that fall in each length.
- However, as you see later on, using a histogram helps to hone in on the specific sequence lengths that we are curious about.

From this function above however, we can conclude that these sequence lengths are relatively large and fit within a small-ish range of values.
This function helps us visualize these differences.

### Histogram for Sequence Lengths above 29800
```julia
data= sortingseq("data/refined_data.fasta")[1];
histogram(data, bins= 10, label= "Sorted Sequences", xlabel= "Sequence Lengths", ylabel= "Number of Sequences", legend=:topleft)
```
This histogram was more intuitive than the last one as I *finally* understood how to initalize and label my data.
In this histogram is the data that I sorted in the function above.
You can clearly see the change in distribution, as the graph moves to the right.
Making this was interesting but not as magnifying as I expected.
There seems to be a large concentration of sequence lengths around 29830.


