# Pairwise InPhyNet

### What is the pairwise algorithm?

The `inphynet(D, constraints, namelist)` function will always provide valid output so long as the networks in `constraints` are valid semi-directed, level-1 networks. The main InPhyNet algorithm attempts to utilize all of the data in `D` and `constraints` simultaneously, but when there are more than 2 constraint networks, this process occassionally fails. InPhyNet still returns a valid network because, when this process fails, a backup version of the InPhyNet algorithm is run that works as follows:

1. Instead of merging all constraints simultaneously, merge `constraints[1]` and `constraints[2]` into a new intermediate network `Nhat12` (using the portions of `D` and `namelist` that pertain only to these constraints).
2. Remove `constraints[1]` and `constraints[2]` from the list of constraints.
3. Add `Nhat12` to the list of constraints.
4. Repeat until `constraints` has only 1 network remaining: this network is the final inferred network.

The main InPhyNet algorithm is guaranteed to return a valid network when 2 constraints are provided, so this "pairwise" algorithm is thereby also guaranteed to return a valid network. This algorithm is implemented in the internal function [`InPhyNet.inphynetpairwise`](@ref).

### Example where the pairwise algorithm is needed

Here, we provide an example where a merge conflict would arise without [`inphynetpairwise`](@ref InPhyNet.inphynetpairwise). A visual depiction and deeper explanation of this example is given in the Appendices of the InPhyNet manuscript (see [here](@ref citation) for a citation containing a URL to the manuscript).

To load the data for this example, you can use the [`loadconflictexample`](@ref InPhyNet.loadconflictexample) function as follows:

```julia
using InPhyNet
loadconflictexample()
inphynet(D, constraints, namelist; refuse_pairwise=true)	# this will give an error
inphynet(D, constraints, namelist)							# this will return a network
```

### Checking whether the pairwise algorithm was used

In some cases, it may be desirable to know whether or not the pairwise algorithm was needed as a fallback. For example, failure of the main algorithm may indicate that the pairwise distance matrix `D` and the networks in `constraints` do not align very well with one another. For analyses only seeking to infer a single network, this is straightforward and can be conducted as shown in the code snippet above. For larger studies, though, this approach is wasteful as both calls to [`inphynet`](@ref InPhyNet.inphynet) will run the initial algorithm before falling back to the pairwise algorithm. In such a case, the following more efficient strategy may be preferred:

```julia
global used_pairwise = false
global Nhat
try
	global Nhat = inphynet(D, constraints, namelist; refuse_pairwise=true)
catch
	global used_pairwise = true
	global Nhat = InPhyNet.inphynetpairwise(D, constraints, namelist)
end
```

In the example we circumvent the need to redo the work of the initial InPhyNet algorithm by calling the internal function [`InPhyNet.inphynetpairwise`](@ref) directly.
