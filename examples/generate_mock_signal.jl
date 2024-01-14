using Pkg
using Revise
Pkg.activate("./")
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
Sb1, Sb2, Sb3, Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12, Sb13, Sb14, Sb15, Sb16, Sb17, Sb18, Sb19, Sb20, Sb21, Sb22, Sb23, Sb24, Sb25, Sb26= logbins_26[1], logbins_26[2], logbins_26[3], logbins_26[4], logbins_26[5], logbins_26[6], logbins_26[7], logbins_26[8], logbins_26[9], logbins_26[10], logbins_26[11], logbins_26[12], logbins_26[13], logbins_26[14], logbins_26[15], logbins_26[16], logbins_26[17], logbins_26[18], logbins_26[19], logbins_26[20], logbins_26[21], logbins_26[22], logbins_26[23], logbins_26[24], logbins_26[25], logbins_26[26]

### Define values for parameters randomly. 10000000 different z (z -> 27 slopes, 1 amplitude, 1 noise parameter A) 
#z= npzread("./z_new.npy") #download dataset from Kaggle https://www.kaggle.com/andronikidimitriou/datasets generated using sbi and load parameters 
z= npzread("/data/users/Androniki/Dani/z_new.npy") #load parameters generated using sbi
#z=rand(1000000,29) # for example uniform prior

### Look at an example and plot it.
z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24, z25, z26 = -12 .+z[1000000,1:27]*(12-(-12))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
Amp =-14+z[1000000,27]*(-6-(-14)) 
A= z[1000000,28]

data,freq=GWBackFinder.model_plot(z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24, z25, z26, Amp, A, f, idx, f_filtered,logbins, Sb1, Sb2  , Sb3  , Sb4  , Sb5  , Sb6  , Sb7  , Sb8  , Sb9  , Sb10  , Sb11  , Sb12  , Sb13  , Sb14  , Sb15  , Sb16  , Sb17  , Sb18  , Sb19  , Sb20  , Sb21  , Sb22  , Sb23  , Sb24  , Sb25   ) #to run it uncomment the frequency calculation in ./GWBackFinder.jl/src/Data_generation/mock_signal.jl
plt.clf()
plt.plot(freq, data[1:3970])
xlabel("f[Hz]")
ylabel("Ω_GW")
plt.show()
plt.title("Mock signal")
plt.gcf()
#%%


### Generate and save data
@showprogress for i in 1:5

    z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24, z25, z26 = -12 .+z[i,1:27]*(12-(-12))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    Amp =-14+z[i,27]*(-6-(-14)) 
    A= z[i,28]

    Data_total= GWBackFinder.model_train_data(z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11,
    z12, z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24, z25, z26, 
    Amp, A, f, idx, f_filtered,logbins, Sb1, Sb2  , Sb3  , Sb4  , Sb5  , Sb6  ,
    Sb7  , Sb8  , Sb9  , Sb10  , Sb11  , Sb12  , Sb13  , Sb14  , Sb15  ,
     Sb16  , Sb17  , Sb18  , Sb19  , Sb20  , Sb21  , Sb22  , Sb23  , Sb24  , Sb25 )

    GWBackFinder.write_sample(Data_total,"./$(i-1).jld2")
end
