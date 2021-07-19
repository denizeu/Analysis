# Calculations

## Mean & Standard Deviation: Lengths and Counts
```julia 
using BioinformaticsBISC195
using Statistics
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

```julia 
function minMax(path)
    ret = parse_fasta(path)[2]
    retmin = length(minimum(ret))
    retmax = length(maximum(ret))
    return (retmin, retmax)
end
```
This is my minMax function.
This function calculates the minimum sequence length in the NCBI dataset and the maximum sequence length in the NCBI dataset.
First, ret is initialized.
Ret calls the parse function, and with the location of [2], ret indicates to the code that we are finding the length of the sequences and not of the headers.
The function does the same for the maximum.
This function is important because it helps us locate which DNA sequences are the longest and shortest.
This is important for this project because it identifies differences in sequences which can be helpful in determining which analyses are best.

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
This function returns the sequence length of each DNA sequence within the dataset from NCBI. 
It creates an array in which the lengths can be stored.
Then, to clean up the data and separate them into header and sequences, the parse_fasta() function is called.
This is important because it allows the two to be separated, and for further functions to be worked on the sequences independent of headers.
Next, I initiated a for loop that goes through the data stored at position 2 by using [2].
This position represents the sequence data and assures that the length function is being called onto the sequences and not the header information.
The length of each position within the data is pushed into the lengths array and then returned.
This is important because the length of sequences can be compared for differences or used within various functions to determine important changes between genomes.

### Histogram of Sequence Lengths
```julia 
using Plots
histogram(seqlength("data/genomes_CoV2.fasta"))
```
Using the function for sequence lengths above and by initiating the Plots package, this creates a histogram that shows the sequence lengths as the x-axis and number of occurences on the y-axis.
This can be used to determine where most of the data exists.
For example, this histogram shows that most of the genomes are above 29500.
This is important in comparing sequence sizes and showing which are more/less similar to one another.

```julia 
function sorting(path)
    sequences = seqlength(path)
    headers = parse_fasta(path)[1]
    check = findall(x->x<29500, sequences)
    cpt = 0
    for i in 1:length(check)
            splice!(headers, check[i] - cpt)
            splice!(sequences, check[i] - cpt)
            cpt = cpt + 1
    end
    return (sequences, headers)
end
```
This function is a sorting function which cuts down the sequence length data to only show the relevant basepair lengths.
Since all my sequences were >25k bases, I chose to cut off all basepairs below the length 29500.
Therefore, this function will return two vectors stored at position [1] and [2].
In [1] exists the sequences that are below the length 29500 and in [2] exists the headers that correspond to the sequences that are below the sequence length 29500.