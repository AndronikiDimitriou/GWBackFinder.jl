using Pkg
using Revise
Pkg.activate("/home/zaldivar/Documents/Androniki/Github/GWBackFinder.jl")
using GWBackFinder
using PyPlot
using NPZ
using JLD2
using ProgressMeter

### Define the frequency range
f = range(start=3 * 1e-5, stop=0.5, step=1e-6) 

### Split in 26 bins. idx26 shows in which bin the frequency belongs to and logbins_26 are the boundaries of the bins . After we coarse grain the frequencies from f = 10−3 Hz to the maximum frequency
### fmax = 0.5 Hz), i.e., we bin them in 1000 intervals of equal log-spacing. idx shows in which of the 1000 bins the frequency belongs to and logbins the boundaries of the 1000 bins.
idx,idx26,logbins_26,logbins,f_filtered=GWBackFinder.binning(f)

### Define breakpoints, i.e., read the boundaries of the 26 bins.
Sb1, Sb2, Sb3, Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12, Sb13, Sb14, Sb15, Sb16, Sb17, Sb18, Sb19, Sb20, Sb21, Sb22, Sb23, Sb24, Sb25= logbins_26[1], logbins_26[2], logbins_26[3], logbins_26[4], logbins_26[5], logbins_26[6], logbins_26[7], logbins_26[8], logbins_26[9], logbins_26[10], logbins_26[11], logbins_26[12], logbins_26[13], logbins_26[14], logbins_26[15], logbins_26[16], logbins_26[17], logbins_26[18], logbins_26[19], logbins_26[20], logbins_26[21], logbins_26[22], logbins_26[23], logbins_26[24], logbins_26[25]

### Define values for parameters randomly. 10000000 different z (z -> 27 slopes, 1 amplitude, 1 noise parameter A) 
z= npzread("/data/users/Androniki/Dani/z_new.npy")

### Look at an example and plot it 
data,freq=GWBackFinder.model(z[5,:],f,idx,f_filtered,logbins_27,logbins, Sb1, Sb2  , Sb3  , Sb4  , Sb5  , Sb6  , Sb7  , Sb8  , Sb9  , Sb10  , Sb11  , Sb12  , Sb13  , Sb14  , Sb15  , Sb16  , Sb17  , Sb18  , Sb19  , Sb20  , Sb21  , Sb22  , Sb23  , Sb24  , Sb25   )
plt.clf()
plt.plot(freq, data[1:3970])
xlabel("f")
ylabel("Ω_GW")
plt.show()
plt.gcf()
#%%

### Generate and save data
@showprogress for i in 336126:500000
    Data_total= GWBackFinder.model(z[i, :],f,idx,f_filtered,logbins_26,logbins, Sb1, Sb2, Sb3,
    Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12,
    Sb13, Sb14, Sb15, Sb16, Sb17, Sb18, Sb19, Sb20, Sb21,
    Sb22, Sb23, Sb24, Sb25)
    GWBackFinder.write_sample(Data_total,"/data/users/Androniki/Dani_new/$(i-1).jld2")
end
