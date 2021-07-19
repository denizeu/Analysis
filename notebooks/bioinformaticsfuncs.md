# Bioinformatics Function Notes

## Unique Kmer Function
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

## Comparing Region and Kmer Dist
