using HDF5
using JLD
using PyPlot

out = Dict{Any,Any}()
if(true)
    desc  = "This experiment takes physionet phd features and samples over number\n"
    desc *= "of trees and the number of feats on the random forest classifer\n"
    desc *= "\n"
    desc *= "The results are stored in a dict [feats, trees] => results\n"
    desc *= "result elements are [Numeric ID, Classes, Predicted, accuracy]\n"
    x = load("physionet-september-thesis-feats-labels-fixed-labels.jld");
    Feats  = x["Feats"]
    Labels = x["Labels"]

    FF = Feats[38:end]
    LL = Labels[38:end]
    println("[INFO] Running validation")
    for i=5:5:50, j=5:5:50
        println("running trees=$j, feats=$i experiement...")
        oo = cross_validate(FF, LL, feats=i, trees=j)
        r = Matrix{Any}(4,length(oo))
        for k=1:length(oo)
            r[:,k] = [37+k, LL[k], oo[k], mean(LL[k].==oo[k])]
        end
        out[(i,j)] = r
        println(mean(r[4,:]))
    end
    save("physionet-sep-super-sensitivity.jld",
         "Description", desc, "Data", out)

    println("[INFO] Stop")
end
