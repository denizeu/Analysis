# Calculations

```julia 
using BioinformaticsBISC195
function read(path)
    lengths= []
    counts=[]
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
This code is calculating the mean and standard deviation of gc content and lengths of my coronavirus genomes.
It is calculating the means of the lengths by initializing a one dimensional array for lengths and counts.
Then, it identifies the data as beinf the CoV2 genomes that I have chosen to study.
Then a for loop is created in which the length of each position is pushed to the length array.
G's and C's are identified and the GC content is determined by dividing number of g's and c's by the length of each sequence.
The GC content is then pushed to the counts array.