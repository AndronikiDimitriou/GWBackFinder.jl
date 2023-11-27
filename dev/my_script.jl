using Pkg
using Revise
Pkg.activate("/home/zaldivar/Documents/Androniki/Github/GWBackFinder.jl")
using GWBackFinder
using Statistics
#%%
idx,idx28,logbins_28=GWBackFinder.binning()
Sb1_new3, Sb2_new3, Sb3_new3, Sb4_new3, Sb5_new3, Sb6_new3, Sb7_new3, Sb8_new3, Sb9_new3, Sb10_new3, Sb11_new3, Sb12_new3, Sb13_new3, Sb14_new3, Sb15_new3, Sb16_new3, Sb17_new3, Sb18_new3, Sb19_new3, Sb20_new3, Sb21_new3, Sb22_new3, Sb23_new3, Sb24_new3, Sb25_new3, Sb26_new3, Sb27_new3 = logbins_28[1], logbins_28[2], logbins_28[3], logbins_28[4], logbins_28[5], logbins_28[6], logbins_28[7], logbins_28[8], logbins_28[9], logbins_28[10], logbins_28[11], logbins_28[12], logbins_28[13], logbins_28[14], logbins_28[15], logbins_28[16], logbins_28[17], logbins_28[18], logbins_28[19], logbins_28[20], logbins_28[21], logbins_28[22], logbins_28[23], logbins_28[24], logbins_28[25], logbins_28[26], logbins_28[27]

