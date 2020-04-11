# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Julia 1.4.0
#     language: julia
#     name: julia-1.4
# ---

readdir()

using JLD2, FileIO, CSV, Plots
using FFTW

frame = CSV.read("frame.csv");
tmp1=frame[!,:KineticEnergy];
ftmp1=fft(tmp1)./50000;
ftmp1[1]=0.;
plot(abs.(ftmp1[1:1000]),xlabel="freq",ylabel="spectrum kinetic energy",legend=false)
#savefig("Kinetic_spectrum_dirac.png")
tmp2=frame[!,:PotentialEnergyE2];
ftmp2=fft(tmp2)./50000;
ftmp2[1]=0.;
#@show kmax=findmax(abs.(ftmp)[1:200])
#@show 77*2*3.1415/50.
plot(abs.(ftmp2[1:1000])./50000,xlabel="freq",ylabel="spectrum Ex energy",legend=false)
#savefig("Ex_spectrum_dirac.png")
tmp3=frame[!,:Kineticspin];
ftmp3=fft(tmp3)./50000;
ftmp3[1]=0.;
#@show kmax=findmax(abs.(ftmp)[1:200])
#@show 77*2*3.1415/50.
plot(abs.(ftmp3[1:1000]),xlabel="freq",ylabel="spectrum Spin energy",legend=false)
#savefig("Spin_spectrum_dirac.png")
tmp1=frame[!,:Momentum1];
ftmp1=fft(tmp1)./50000;
ftmp1[1]=0.;
#@show kmax=findmax(abs.(ftmp)[1:200])
#@show 77*2*3.1415/50.
plot(abs.(ftmp1[1:1000]),xlabel="freq",ylabel="Ex energy",legend=false)
#savefig("Ex_spectrum_dirac.png")
tmp1=frame[!,:PotentialEnergyE2];
ftmp1=fft(tmp1)./50000;
ftmp1[1]=0.;
#@show kmax=findmax(abs.(ftmp)[1:200])
#@show 77*2*3.1415/50.
plot(abs.(ftmp1[1:1000]),xlabel="freq",ylabel="Ey energy",legend=false)
#savefig("Ey_spectrum_dirac.png")
tmp1=frame[!,:PotentialEnergyE3];
ftmp1=fft(tmp1)./50000;
ftmp1[1]=0.;
#@show kmax=findmax(abs.(ftmp)[1:200])
#@show 77*2*3.1415/50.
plot(abs.(ftmp1[1:1000]),xlabel="freq",ylabel="Ez energy",legend=false)
#savefig("Ez_spectrum_dirac.png")
tmp1=frame[!,:PotentialEnergyB2];
ftmp1=fft(tmp1)./50000;
ftmp1[1]=0.;
#@show kmax=findmax(abs.(ftmp)[1:200])
#@show 77*2*3.1415/50.
plot(abs.(ftmp1[1:1000]),xlabel="freq",ylabel="By energy",legend=false)
#savefig("By_spectrum_dirac.png")
tmp1=frame[!,:PotentialEnergyB3];
ftmp1=fft(tmp1)./50000;
ftmp1[1]=0.;
#@show kmax=findmax(abs.(ftmp)[1:200])
#@show 77*2*3.1415/50.
plot(abs.(ftmp1[1:1000]),xlabel="freq",ylabel="Bz energy",legend=false)
#savefig("Bz_spectrum_dirac.png")
# +
#Kineticspin,Momentum1,Momentum2,Momentum3,Momentum4,
#Momentum5,Momentum6,Momentum7,PotentialEnergyE1,
#PotentialEnergyE2,PotentialEnergyE3,PotentialEnergyB2,PotentialEnergyB3
# -


