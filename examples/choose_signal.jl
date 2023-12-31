using Pkg
using Revise
Pkg.activate("/home/zaldivar/Documents/Androniki/Github/GWBackFinder.jl")
using GWBackFinder
using PyPlot

f = range(start=3 * 1e-5, stop=0.5, step=1e-6)

idx,idx27,logbins_27,logbins,f_filtered = GWBackFinder.binning(f)
Sb1_new3, Sb2_new3, Sb3_new3, Sb4_new3, Sb5_new3, Sb6_new3, Sb7_new3, Sb8_new3, Sb9_new3, Sb10_new3, Sb11_new3, Sb12_new3, Sb13_new3, Sb14_new3, Sb15_new3, Sb16_new3, Sb17_new3, Sb18_new3, Sb19_new3, Sb20_new3, Sb21_new3, Sb22_new3, Sb23_new3, Sb24_new3, Sb25_new3, Sb26_new3, Sb27_new3 = logbins_27[1], logbins_27[2], logbins_27[3], logbins_27[4], logbins_27[5], logbins_27[6], logbins_27[7], logbins_27[8], logbins_27[9], logbins_27[10], logbins_27[11], logbins_27[12], logbins_27[13], logbins_27[14], logbins_27[15], logbins_27[16], logbins_27[17], logbins_27[18], logbins_27[19], logbins_27[20], logbins_27[21], logbins_27[22], logbins_27[23], logbins_27[24], logbins_27[25], logbins_27[26], logbins_27[27]

z_power_law_test         =  GWBackFinder.zPowerLaw([-9,0.6])
z_peak_test              =  GWBackFinder.zPeak([-9.,0.002,0.2])
z_wiggly_test            =  GWBackFinder.zWiggly([-10,2,0.1])
z_Broken_powerlaw_test   =  GWBackFinder.zBroken_powerlaw([-11,-1,2/3,0.002])
z_Double_peaks_test      =  GWBackFinder.zDouble_peaks([-11.,-10.,0.001,0.01,0.1,0.1])
z_three_peaks_test       =  GWBackFinder.zThree_peaks([-10.,-10.,-10.,5*10^(-4),2*10^(-3),8*10^(-3),0.1,0.1,0.1])
z_noise                  =  GWBackFinder.zNoise([15.,3.])

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ  = GWBackFinder.different_signals(z_power_law_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ)
plt.plot(f_total_ΤΤ, Data_ΤΤ)
plt.show()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ  = GWBackFinder.different_signals(z_peak_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ)
plt.plot(f_total_ΤΤ, Data_ΤΤ)
plt.show()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ = GWBackFinder.different_signals(z_wiggly_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ)
plt.plot(f_total_ΤΤ, Data_ΤΤ)
plt.show()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ  = GWBackFinder.different_signals(z_Broken_powerlaw_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ)
plt.plot(f_total_ΤΤ, Data_ΤΤ)
plt.show()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ  = GWBackFinder.different_signals(z_Double_peaks_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ)
plt.plot(f_total_ΤΤ, Data_ΤΤ)
plt.show()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ = GWBackFinder.different_signals(z_three_peaks_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ)
plt.plot(f_total_ΤΤ, Data_ΤΤ)
plt.show()
plt.gcf()

