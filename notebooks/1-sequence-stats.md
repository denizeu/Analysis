# Calculations

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
This function calculates the minimum sequence length in the NCBI dataset and the maximum sequence length in the NCVI dataset.
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
Coronavirus genome lengths

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
    