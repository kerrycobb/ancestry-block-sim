# ancestry-block-sim
Tool to simulate recombination of ancestry blocks

Requires the Nim compiler. https://nim-lang.org/install.html

## Examples
See `/simulations` for examples.

#### Heterozygote disadvantage with single site in panmictic population
```bash
nim c -d:release singeSiteSteppingStone.nim
./singleSitePanmictic -g=1_000 -s=100 -d=1_000 -l=1_000_000 -p=500_000 -r=1e-6 --seed=12345 -o singleSitePanmictic.out
```
Run `./singleSitePanmictic --help` for list of arguments

#### Heterozygoe disadvantage with single site in stepping stone model 
```bash
./singleSiteSteppingStone -g=1_000 -d=50 --demeSize=1000 -l=1_000_000 -p=500_000 -r=1e-6 -s=1234 -o=singleSiteSteppingStone.out
```
Run `./singleSiteSteppingStone --help` for list of arguments

## TODO:
In order of priority
- [] How exactly is recombination rate defined? Check that code is consistent with it.
- [] Calculate proportion ancestry - Kerry
- [] procs to compute block length range, std, and mean
- [] plot junctions in ggplot
- [] how to calculate expected block length and junction number
- [] procs to compute centimorgans given the haplotype length and recombination probabilities