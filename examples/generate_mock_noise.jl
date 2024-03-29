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
Sb1, Sb2, Sb3, Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12, Sb13, Sb14, Sb15, Sb16, Sb17, Sb18, Sb19, Sb20, Sb21, Sb22, Sb23, Sb24, Sb25, Sb26 = logbins_26[1], logbins_26[2], logbins_26[3], logbins_26[4], logbins_26[5], logbins_26[6], logbins_26[7], logbins_26[8], logbins_26[9], logbins_26[10], logbins_26[11], logbins_26[12], logbins_26[13], logbins_26[14], logbins_26[15], logbins_26[16], logbins_26[17], logbins_26[18], logbins_26[19], logbins_26[20], logbins_26[21], logbins_26[22], logbins_26[23], logbins_26[24], logbins_26[25], logbins_26[26]

### Define values for noise parameter 
#z= npzread("./z_noise.npy") #download dataset from Kaggle https://www.kaggle.com/andronikidimitriou/datasets generated using sbi and load parameters 
z= npzread("/data/users/Androniki/Dani/z_noise.npy") #load parameters generated using sbi
#z=rand(1000000,1) # for example uniform prior


### Look at an example and plot it 
data,freq=GWBackFinder.model_noise_plot(f,z[1]) #to run it uncomment the frequency calculation in /home/zaldivar/Documents/Androniki/Github/GWBackFinder.jl/src/Data_generation/mock_noise.jl
plt.clf()
plt.plot(freq, data[1:2970])
xlabel("f[Hz]")
ylabel("Ω_GW")
plt.title("Mock noise")
plt.show()
plt.gcf()


### Generate and save data
@showprogress for i in 1:1000000
    Data_total= GWBackFinder.model_noise_train_data(f,z[i])
    #print(length(Data_total))
    GWBackFinder.write_sample(Data_total,"/data/users/Androniki/Dani_new_noise/$(i-1).jld2")
end
