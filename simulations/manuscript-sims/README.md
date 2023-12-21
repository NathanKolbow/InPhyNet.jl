## Manuscript simulations

1. Monophyletic partitions
   - Ground truth constraints accuracy
   - NNI noise accuracy
   - Distance noise accuracy
   - NNI & distance noise accuracy
     - the above 4 cases are run with the function `monophyleticRobustness` in `src/simulations/robustness-fxns.jl`
   - Taxa swapping accuracy
2. Full random partitions
   - Ground truth constraints accuracy
   - NNI noise accuracy
   - Distance noise accuracy
   - NNI & distance noise accuracy
3.  Algorithm to select partitions w/ data simulated down to sequences
   - This is a pipeline: true net → true gts → sequences → est gts → subsetting algorithm → constraint network estimation (w/ model selection for reticulation count) → merged network