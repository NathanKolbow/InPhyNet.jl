using Test, InPhyNet

t1_good = readnewick("((A,B),(C,(D,E)));")
t2_good = readnewick("((A,B),(C,(D,E)));")
t3_good = readnewick("((C,D),(A,B));")
t4_good = readnewick("(A,(B,D));")

@test are_compatible_after_merge([t1_good, t2_good, t3_good, t4_good], "D", "E")
@test are_compatible_after_merge([t1_good, t3_good, t4_good], "D", "E")
@test are_compatible_after_merge([t1_good, t2_good], "D", "E")
@test are_compatible_after_merge([t3_good, t4_good], "D", "E")
@test are_compatible_after_merge([t1_good, t3_good], "D", "E")

t1_bad = readnewick("(((A,G),C),(E,F));")
t2_bad = readnewick("(((A,B),C),(D,E));")

@test !are_compatible_after_merge(t1_bad, t2_bad, "D", "G")
@test !are_compatible_after_merge([t1_bad, t2_bad, t1_good], "D", "G")
@test !are_compatible_after_merge([t1_bad, t2_bad, t2_good], "D", "G")
@test !are_compatible_after_merge([t1_bad, t2_bad, t3_good], "D", "G")
@test !are_compatible_after_merge([t1_bad, t2_bad, t1_good, t2_good], "D", "G")
@test !are_compatible_after_merge([t1_bad, t2_bad, t2_good, t3_good], "D", "G")
@test !are_compatible_after_merge([t1_bad, t2_bad, t1_good, t2_good, t3_good, t4_good], "D", "G")
