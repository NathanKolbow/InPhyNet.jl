include("../main.jl")

# Make sure all probabilities sum to 1 and that no errors occur
ptrees, _ = getparentaltrees(readTopology("((A,((B,(C,D)))#H1),(#H1,E));"))
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")

ptrees, _ = getparentaltrees(readTopology("(((A,B),#H1), (((C,(D,#H2)))#H1,((E)#H2,F)));"))
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")

ptrees, _ = getparentaltrees("(10,(#H2,(1,(2,(((9)#H1,(7,(8,#H1))))#H2))))root;")
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")

ptrees, _ = getparentaltrees("(10,(#H2,(1,(2,(((9)#H1,(3,(6,(8,#H1)))))#H2))))root;")
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")

# This one takes a *long* time to run
ptrees, _ = getparentaltrees("(10,(#H2,(1,(2,(((9)#H1,(3,(4,(5,((6,7),(8,#H1)))))))#H2))))root;")
abs(sum([prob(t) for t in ptrees]) - 1) < 1e-12 || error("Parental tree probabilities do not sum to 1.")