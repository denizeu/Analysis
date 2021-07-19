# Review of Analysis and Code repos

**date**: 2021-07-19

## Review Notes

- See instructions in Pull Request for incorporating changes
- Comments in julia code are preceded with `#`
- Comments in markdown look like this: `<!-- this is a comment -->`
- Comments that require action begin with `TODO`, others are just for information
- You may delete comments inside other files once you've addressed them
- Please keep this file in the repository - you may add your own responses if aplicable.

## Analysis plan

- added some newlines to follow [sembr](http://sembr.org)

## Code repo

https://github.com/denizeu/BioinformaticsBISC195.jl

- Code repository currently throws an error when I try to load it
  ```
  julia> using BioinformaticsBISC195
  [ Info: Precompiling BioinformaticsBISC195 [67893044-1481-44cf-8ab6-de350e09ee4e]
  ERROR: LoadError: syntax: incomplete: "function" at /home/kevin/Repos/bisc195/  final_projects/deniz/BioinformaticsBISC195.jl/src/BioinformaticsBISC195.jl:359   requires end
  ```

- several functions were missing `end`
- A bunch of functions defined in your analysis repo should be moved to code REPO,
  and get docs and tests
- A number of functions throw `UndefVarError`s - I'm guessing you have some global variables
  hanging around in your environment when you're testing stuff.
  But your code needs to run all the way through as-written
- Fantastic set of tests, but a bunch of tests don't work or error:

  ```
  Test Summary:          | Pass  Fail  Error  Total
  BioinformaticsBISC195  |   35     3      5     43
    Using Strings        |   35     3      5     43
      normalizeDNA       |    5                   5
      composition        |   11                  11
      gc_content         |    3     1             4
      complement         |    3     1             4
      reverse_complement |    3     1             4
      parse_fasta        |    8                   8
      uniqueKmers        |                 3      3
      kmerdist           |    2                   2
      kmertime           |                 1      1
      kmertimes          |                 1      1
      kmerloc            |                    No tests
  ```


## Analysis repo

- You don't really need to re-define the functions in your analysis repo, since you can just do `using BioinformaticsBISC195` to get them
  Since you already wrote docstrings, you can do eg `@doc kmertime` to show the docstring in your notebooks.
  - This also avoids the issue of changing the function in one place and needing to update it in the other
- I can't evaluate your plots, because the code doesn't run