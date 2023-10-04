# All the setup for the DEBUG scenario
include("../main.jl")

#### SETUP ####
debugnet = readTopology("(10,(#H2,(1,(2,(((9)#H1,(3,(5,(7,(8,#H1))))))#H2))))root;")
initnet, ldict, _ = _initparentaltreealgo(debugnet)

workingset = Queue{IPT}()
enqueue!(workingset, initnet)
parentaltrees = Vector{IPT}()
####

# 2 runs go fine
for i=1:2
    ipt = dequeue!(workingset)
    net = top(ipt)

    # 0. No more hybrids? This net is done!
    if length(net.hybrid) == 0
        push!(parentaltrees, ipt)
        continue
    end

    # 1. Find a hybrid node with no other hybrids below it
    node, nodeidx = _getnexthybrid(net)
    divisions = Vector{IPT}()

    # 2a. If there is already a LineageNode associated with this node,
    #     then we are ready to condition on the reticulation
    if ldict[node] !== nothing
        print("Conditioning on reticulations - ")
        divisions = _conditiononreticulation(ipt, node, ldict)
        println("done.")

    # 2b. If there is not a LineageNode associated with thisnode,
    #     we must first condition on the coalescent events
    else
        print("Conditioning on coalescences - ")
        divisions = _conditiononcoalescences(ipt, node, ldict)
        println("done.")
    end

    # 4. Repeat with the resultant networks
    for d in divisions enqueue!(workingset, d) end
end

# Next loop causes issues
ipt = dequeue!(workingset)
net = top(ipt)

node, nodeidx = _getnexthybrid(net)
divisions = Vector{IPT}()
# Done setting up

println("\n\n\n\n\n\n\n`ipt` is the next IPT to be evaluated (by conditioning on coalescences). This is where the KeyError happens")
println("`node` is the next hybrid being evaluated.")