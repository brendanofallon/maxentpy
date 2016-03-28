# maxentpy

maxentpy is a python wrapper for MaxEntScan to calculate splice site strength.

It contains two functions. `score5` is adapted from
[MaxEntScan::score5ss](http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html)
to score 5' splice sites. `score3` is adapted from
[MaxEntScan::score3ss](http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq_acc.html)
to score 3' splice sites. They only use Maximum Entropy Model to score.


## Examples

```python
import maxentpy as maxent
scorer = maxent.SpliceScorer()

#Score splice donor, must have exactly 3 bases on exonic and 6 bases of intronic sequence
scorer.score5('cagGTAAGT')

#Score splice acceptor, must have three bases of exonic sequence and 20 bases on intronic sequence
scorer.score3('ttccaaacgaacttttgtAGgga')
```
