using HypothesisTests
#include("scratch.jl")
#include("scratch-time-band-physionet.jl")

phyc_raw = physionet_raw_short
phyt_raw = physionet_raw_long
phyt_dds = physionet_phd_long
phyc_dds = physionet_phd_short

phyt_tim = time_long
phyc_tim = time_short
phyt_bnd = bandpow_long
phyc_bnd = bandpow_short
 
dres_tim = time_domain_sub

drep_tim = time_domain_main
drep_dds = phd_meth
drep_raw = raw_data_acc
drep_bnd = bandpow

phyc = [phyc_dds, phyc_raw, phyc_tim, phyc_bnd]
phyt = [phyt_dds, phyt_raw, phyt_tim, phyt_bnd]
drep = [drep_dds, drep_raw, drep_tim, drep_bnd]

#Exclude trials which used an alternative electrode combination
#Note - The data was there, but the data preprocessing resulted in mistakenly
#       using the wrong channel. This resulted in basically a chance level
#       prediction.
mask = find(mapslices(minimum,hcat(bandpow[:,2],phd_meth[:,2]),2).>0.5)

function pr(a,b)
    @printf("%0.3g", pvalue(UnequalVarianceTTest(a,b)))
end
function stats(x,mask=nothing)
    if(mask != nothing)
        pr(x[1][mask,2], x[2][mask,2])
        print(" & ")
        pr(x[1][mask,2], x[3][mask,2])
        print(" & ")
        pr(x[1][mask,2], x[4][mask,2])
        println("\\\\")
    else
        pr(x[1][:,2], x[2][:,2])
        print(" & ")
        pr(x[1][:,2], x[3][:,2])
        print(" & ")
        pr(x[1][:,2], x[4][:,2])
        println("\\\\")
    end
end

stats(phyc)
stats(phyt)
stats(drep,mask)
println("Done")
