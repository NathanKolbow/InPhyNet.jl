# Constructing the full network

Finally, putting the full network together with InPhyNet is simple:

```julia
species_net = inphynet(D, snaq_networks, namelist)
```

Our example has an outgroup taxa, so we can root the network at that outgroup to obtain the network shown below:

```julia
rootatnode!(species_net, "OUTGROUP")
```

![Full network](full_net.png)

