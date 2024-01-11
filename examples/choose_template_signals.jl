using Pkg
using Revise
Pkg.activate("./")
#Pkg.activate("/home/zaldivar/Documents/Androniki/Github/GWBackFinder.jl")
using GWBackFinder
using PyPlot

### Define the frequency range

f = range(start=3 * 1e-5, stop=0.5, step=1e-6)

### Split in 26 bins. idx26 shows in which bin the frequency belongs to and logbins_26 are the boundaries of the bins . After we coarse grain the frequencies from f = 10−3 Hz to the maximum frequency
### fmax = 0.5 Hz), i.e., we bin them in 1000 intervals of equal log-spacing. idx shows in which of the 1000 bins the frequency belongs to and logbins the boundaries of the 1000 bins.

idx,idx26,logbins_26,logbins,f_filtered = GWBackFinder.binning(f)

### Define different different_signals

z_power_law_test         =  GWBackFinder.zPowerLaw([-9,0.6])
z_peak_test              =  GWBackFinder.zPeak([-9.,0.002,0.2])
z_wiggly_test            =  GWBackFinder.zWiggly([-10,2,0.1])
z_Broken_powerlaw_test   =  GWBackFinder.zBroken_powerlaw([-11,-1,2/3,0.002])
z_Double_peaks_test      =  GWBackFinder.zDouble_peaks([-11.,-10.,0.001,0.01,0.1,0.1])
z_three_peaks_test       =  GWBackFinder.zThree_peaks([-10.,-10.,-10.,5*10^(-4),2*10^(-3),8*10^(-3),0.1,0.1,0.1])
z_noise                  =  GWBackFinder.zNoise([15.,3.])

### Generate different signals and plot them
plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ  = GWBackFinder.different_signals(z_power_law_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ[1:3970],label="AA channel")
plt.plot(f_total_ΤΤ, Data_ΤΤ,label="TT channel")
plt.show()
plt.title("Powerlaw")
plt.legend()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ  = GWBackFinder.different_signals(z_peak_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ[1:3970],label="AA channel")
plt.plot(f_total_ΤΤ, Data_ΤΤ,label="TT channel")
plt.show()
plt.title("Peak")
plt.legend()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ = GWBackFinder.different_signals(z_wiggly_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ[1:3970],label="AA channel")
plt.plot(f_total_ΤΤ, Data_ΤΤ,label="TT channel")
plt.show()
plt.title("Wiggly")
plt.legend()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ  = GWBackFinder.different_signals(z_Broken_powerlaw_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ[1:3970],label="AA channel")
plt.plot(f_total_ΤΤ, Data_ΤΤ,label="TT channel")
plt.title("Broken powerlaw")
plt.legend()
plt.show()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ  = GWBackFinder.different_signals(z_Double_peaks_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ[1:3970],label="AA channel")
plt.plot(f_total_ΤΤ, Data_ΤΤ,label="TT channel")
plt.show()
plt.title("Double_peak")
plt.legend()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ = GWBackFinder.different_signals(z_three_peaks_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ[1:3970],label="AA channel")
plt.plot(f_total_ΤΤ, Data_ΤΤ,label="TT channel")
plt.title("Three peaks")
plt.show()
plt.legend()
plt.gcf()

