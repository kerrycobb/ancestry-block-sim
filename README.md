# ancestry-block-sim
Tool to simulate recombination of different ancestries

Requires Nim compiler. https://nim-lang.org/install.html

## Examples
See `/simulations` for examples.

Command to compile and run one of the example simulations
```bash
nim c -d:release singeSiteSteppingStone.nim
./singleSiteSteppingStone <outfile name: string> [starting seed: int]
```


## TODO:
In order of priority
[] Calculate proportion ancestry - Kerry
[] procs to compute block length range, std, and mean
[] plot junctions in ggplot
[] how to calculate block length and junction number expectations
[] procs to compute centimorgans given the haplotype length and recombination probabilities