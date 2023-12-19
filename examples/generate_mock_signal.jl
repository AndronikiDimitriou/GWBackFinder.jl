using Pkg
using Revise
Pkg.activate("/home/zaldivar/Documents/Androniki/Github/GWBackFinder.jl")
using GWBackFinder
using PyPlot
using NPZ
using JLD2
using ProgressMeter


#%%
f = range(start=3 * 1e-5, stop=0.5, step=1e-6)
idx,idx27,logbins_27,logbins,f_filtered=GWBackFinder.binning(f)
Sb1_new3, Sb2_new3, Sb3_new3, Sb4_new3, Sb5_new3, Sb6_new3, Sb7_new3, Sb8_new3, Sb9_new3, Sb10_new3, Sb11_new3, Sb12_new3, Sb13_new3, Sb14_new3, Sb15_new3, Sb16_new3, Sb17_new3, Sb18_new3, Sb19_new3, Sb20_new3, Sb21_new3, Sb22_new3, Sb23_new3, Sb24_new3, Sb25_new3, Sb26_new3, Sb27_new3 = logbins_27[1], logbins_27[2], logbins_27[3], logbins_27[4], logbins_27[5], logbins_27[6], logbins_27[7], logbins_27[8], logbins_27[9], logbins_27[10], logbins_27[11], logbins_27[12], logbins_27[13], logbins_27[14], logbins_27[15], logbins_27[16], logbins_27[17], logbins_27[18], logbins_27[19], logbins_27[20], logbins_27[21], logbins_27[22], logbins_27[23], logbins_27[24], logbins_27[25], logbins_27[26], logbins_27[27]

z1= npzread("/home/zaldivar/Documents/Androniki/Github/GWBackFinder/src/GWBackFinder/data/thetas_for_julia.npy")
length(z1)
#z1=vcat(rand(28),3)
#length(z1)
#z1=[0.1,0.5,0.2,0.6,0.3,0.4,0.7,0.9,0.8,1.0,0.1,0.5,0.2,0.6,0.3,0.4,0.7,0.9,0.8,1.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,3]
#z1=[0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,0.525,3]
z1
z1[1500, :]

data=GWBackFinder.model(z1[10,:],f,idx,f_filtered,logbins_27,logbins, Sb1_new3, Sb2_new3, Sb3_new3, Sb4_new3, Sb5_new3, Sb6_new3, Sb7_new3, Sb8_new3, Sb9_new3, Sb10_new3, Sb11_new3, Sb12_new3, Sb13_new3, Sb14_new3, Sb15_new3, Sb16_new3, Sb17_new3, Sb18_new3, Sb19_new3, Sb20_new3, Sb21_new3, Sb22_new3, Sb23_new3, Sb24_new3, Sb25_new3, Sb26_new3)

plt.clf()
#plt.plot(freq, data[1][1])
#plt.plot(log10.(f),Data)
plt.axvline(log10(Sb14_new3))
plt.axhline(log10(10^(-14+0.525*(-6-(-14)))))

"""
plt.axvline(log10(Sb1_new3))
plt.axvline(log10(Sb2_new3))
plt.axvline(log10(Sb3_new3))
plt.axvline(log10(Sb4_new3))
plt.axvline(log10(Sb5_new3))
plt.axvline(log10(Sb6_new3))
plt.axvline(log10(Sb7_new3))
plt.axvline(log10(Sb8_new3))
plt.axvline(log10(Sb9_new3))
plt.axvline(log10(Sb10_new3))
plt.axvline(log10(Sb11_new3))
plt.axvline(log10(Sb12_new3))
plt.axvline(log10(Sb13_new3))
plt.axvline(log10(Sb14_new3))
plt.axvline(log10(Sb15_new3))
plt.axvline(log10(Sb16_new3))
plt.axvline(log10(Sb17_new3))
plt.axvline(log10(Sb18_new3))
plt.axvline(log10(Sb19_new3))
plt.axvline(log10(Sb20_new3))
plt.axvline(log10(Sb21_new3))
plt.axvline(log10(Sb22_new3))
plt.axvline(log10(Sb23_new3))
plt.axvline(log10(Sb24_new3))
plt.axvline(log10(Sb25_new3))
plt.axvline(log10(Sb26_new3))
plt.axvline(log10(Sb27_new3))
"""

plt.show()
plt.gcf()


#%%

@showprogress for i in 426422:size(z1)[1]
    GWBackFinder.write_sample(z1[i, :], i)
end

