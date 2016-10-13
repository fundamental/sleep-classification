using HDF5
using JLD
using PyPlot

out = Dict{Any,Any}()
if(false)
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
        println(now())
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

function getSpectra(SubjectID, workDir)
    Sp = raw_spec(SubjectID)
    Cl = readcsv("$workDir/Dejunked$SubjectID-cols.csv")
    Sp[:,find(Cl)]
    #imread("$workDir/DejunkedSpectra$SubjectID.png")
end

function obtainGradient(samples::Vector{Float64})
    #Optimize for the form of y=a*x+b
    N = length(samples)
    if(N == 1)
        return samples
    end
    (a,b) = [linspace(0,1,N) ones(N)]\samples
    out = zeros(N)
    for i=1:N
        out[i] = a*(i-1)/(N-1)+b
    end
    out
end

function zoneifyGradient(img, X, Y)
    I = float(copy(img))
    xedge = find(X)
    #println("xedge = ", xedge)
    yedge = find(Y)
    prepend!(xedge, [1])
    prepend!(yedge, [1])
    push!(yedge, size(img)[1])
    push!(xedge, size(img)[2])
    #println(size(img))
    #println(xedge)
    #println(yedge)
    #println(maximum(yedge))
    for yy=1:length(yedge)-1, xx=1:length(xedge)-1
        I[yedge[yy]:yedge[yy+1],xedge[xx]:xedge[xx+1]] = obtainGradient(img[yedge[yy]:yedge[yy+1],xedge[xx]:xedge[xx+1]][:])
    end
    I
end

function viewStuff(SubjectID::String, workDir::String, amount::Float64, misc=false)
    ep  = readcsv("$workDir/edgeParameter$SubjectID.csv")
    dat = getSpectra(SubjectID, workDir)[1:1000,:]
    ll  = readcsv("$workDir/Dejunkedlabels$SubjectID.csv")

    X = find(maximum(ep,2)[:])
    Y = find(ones(size(dat)[1]))
    prepend!(X,[1])
    push!(X,size(dat,2))

    if(misc)
        mn = zoneifyGradient(dat, X, Y)
        md = zoneify(dat, X, Y, operator=median)
    else
        mn = zoneify(dat, X, Y, operator=mean)
        md = zoneify(dat, X, Y, operator=median)
    end
    out = (1.0-amount)*dat+(amount)*(mn+md)/2.0
    (ll, out)
end

if(false)
    desc  = "This experiment takes physionet phd features and samples over
    different\n"
    desc *= "amounts of denoising from 0% to 100%\n"
    desc *= "\n"
    desc *= "The results are stored in a dict denoising-fraction => results\n"
    desc *= "result elements are [Numeric ID, Classes, Predicted, accuracy]\n"
    x = load("physionet-september-thesis-feats-labels-fixed-labels.jld");
    #Feats  = x["Feats"]
    #Labels = x["Labels"]

    #FF = Feats[38:end]
    #LL = Labels[38:end]
    drange = data_range()[38:end]

    println("[INFO] Running validation")
    for i=linspace(0.0,1.0,10)
        println(now())
        println("running denoising-amount=$i experiement...")

        #Build up the data
        FF = Matrix{Float64}[]
        LL = Vector{Int}[]
        for j=drange
            if(j != "ST7221")
                tmp = viewStuff(j, "physionet/", i, false)
                push!(LL, map(Int, tmp[1][:]))
                push!(FF, tmp[2])
            end
        end

        #Validate the data
        oo = cross_validate(FF, LL, feats=20, trees=16)
        r = Matrix{Any}(4,length(oo))
        for k=1:length(oo)
            r[:,k] = [37+k, LL[k], oo[k], mean(LL[k].==oo[k])]
        end
        out[i] = r
        println(mean(r[4,:]))
    end
    save("physionet-oct-denoise-amount-sensitivity.jld",
         "Description", desc, "Data", out)

    println("[INFO] Stop")
end
if(false)
    desc  = "This experiment takes physionet phd features and samples over
    different\n"
    desc *= "amounts of denoising from 0% to 100% testing only on one subject\n"
    desc *= "\n"
    desc *= "The results are stored in a dict denoising-fraction => results\n"
    desc *= "result elements are [Numeric ID, Classes, Predicted, accuracy]\n"
    x = load("physionet-september-thesis-feats-labels-fixed-labels.jld");
    #Feats  = x["Feats"]
    #Labels = x["Labels"]

    #FF = Feats[38:end]
    #LL = Labels[38:end]
    drange = data_range()[38:end]
    test_sub = 2

    println("[INFO] Running validation")
    for i=linspace(0.0,1.0,100)
if(i<0.84)
    continue
end
        println(now())
        println("running denoising-amount=$i experiement...")

        #Build up the data
        FF = Matrix{Float64}[]
        LL = Vector{Int}[]
        for j=drange
            if(j != "ST7221")
                tmp = viewStuff(j, "physionet/", i, false)
                push!(LL, map(Int, tmp[1][:]))
                push!(FF, tmp[2])
            end
        end

        #Validate the data
        oo = validate(hcat(FF[1],FF[3:end]...), vcat(LL[1],LL[3:end]...), FF[2],
                      LL[2], feats=20, trees=16)
        r = Vector{Any}(3)
        r[1] = LL[2]
        r[2] = oo
        r[3] = mean(LL[2].==oo)
        out[i] = r
        println(r[3])
    end
    #save("physionet-oct-denoise-one-sub-amount-sensitivity.jld",
    #     "Description", desc, "Data", out)

    println("[INFO] Stop")
end
2

one_patient_mix_amount = [
0.0 0.8933518005540166
0.010101010101010102 0.8933518005540166
0.020202020202020204 0.8883656509695291
0.030303030303030304 0.8983379501385041
0.04040404040404041 0.8911357340720222
0.050505050505050504 0.896398891966759
0.06060606060606061 0.8933518005540166
0.0707070707070707 0.8958448753462603
0.08080808080808081 0.8988919667590027
0.09090909090909091 0.8955678670360111
0.10101010101010101 0.900831024930748
0.1111111111111111 0.9011080332409972
0.12121212121212122 0.8988919667590027
0.13131313131313133 0.8980609418282548
0.1414141414141414 0.9027700831024931
0.15151515151515152 0.9027700831024931
0.16161616161616163 0.900831024930748
0.1717171717171717 0.9052631578947369
0.18181818181818182 0.9066481994459834
0.1919191919191919 0.9049861495844875
0.20202020202020202 0.907202216066482
0.21212121212121213 0.9044321329639889
0.2222222222222222 0.907202216066482
0.23232323232323232 0.9083102493074793
0.24242424242424243 0.9069252077562326
0.25252525252525254 0.9091412742382271
0.26262626262626265 0.907202216066482
0.2727272727272727 0.909972299168975
0.2828282828282828 0.909972299168975
0.29292929292929293 0.9174515235457064
0.30303030303030304 0.9124653739612189
0.31313131313131315 0.9091412742382271
0.32323232323232326 0.9130193905817174
0.3333333333333333 0.9130193905817174
0.3434343434343434 0.9157894736842105
0.35353535353535354 0.9157894736842105
0.36363636363636365 0.9185595567867036
0.37373737373737376 0.9193905817174515
0.3838383838383838 0.917174515235457
0.3939393939393939 0.9174515235457064
0.40404040404040403 0.9202216066481994
0.41414141414141414 0.9193905817174515
0.42424242424242425 0.9166204986149584
0.43434343434343436 0.9191135734072022
0.4444444444444444 0.9213296398891967
0.45454545454545453 0.9202216066481994
0.46464646464646464 0.9227146814404432
0.47474747474747475 0.9224376731301939
0.48484848484848486 0.9204986149584488
0.494949494949495 0.9188365650969529
0.5050505050505051 0.9227146814404432
0.5151515151515151 0.9238227146814404
0.5252525252525253 0.9218836565096953
0.5353535353535354 0.9232686980609418
0.5454545454545454 0.9246537396121883
0.5555555555555556 0.9235457063711912
0.5656565656565656 0.9243767313019391
0.5757575757575758 0.9240997229916897
0.5858585858585859 0.9240997229916897
0.5959595959595959 0.9240997229916897
0.6060606060606061 0.9188365650969529
0.6161616161616161 0.9265927977839336
0.6262626262626263 0.9277008310249307
0.6363636363636364 0.9246537396121883
0.6464646464646465 0.9243767313019391
0.6565656565656566 0.925207756232687
0.6666666666666666 0.9268698060941828
0.6767676767676768 0.9202216066481994
0.6868686868686869 0.9263157894736842
0.696969696969697 0.925207756232687
0.7070707070707071 0.9232686980609418
0.7171717171717171 0.925207756232687
0.7272727272727273 0.9246537396121883
0.7373737373737373 0.9240997229916897
0.7474747474747475 0.9263157894736842
0.7575757575757576 0.9254847645429363
0.7676767676767676 0.9213296398891967
0.7777777777777778 0.9279778393351801
0.7878787878787878 0.9265927977839336
0.797979797979798 0.9246537396121883
0.8080808080808081 0.9191135734072022
0.8181818181818182 0.9246537396121883
0.8282828282828283 0.9263157894736842
0.8383838383838383 0.9263157894736842
0.8484848484848485 0.9229916897506926
0.8585858585858586 0.925207756232687
0.8686868686868687 0.9282548476454293
0.8787878787878788 0.9243767313019391
0.8888888888888888 0.9204986149584488
0.898989898989899 0.9191135734072022
0.9090909090909091 0.9152354570637119
0.9191919191919192 0.9199445983379502
0.9292929292929293 0.9218836565096953
0.9393939393939394 0.9155124653739612
0.9494949494949495 0.9232686980609418
0.9595959595959596 0.9132963988919668
0.9696969696969697 0.9240997229916897
0.9797979797979798 0.9196675900277008
0.98989898989899 0.9301939058171745
1.0 0.9249307479224377]

#out = load("physionet-oct-denoise-amount-sensitivity.jld")["Data"]
