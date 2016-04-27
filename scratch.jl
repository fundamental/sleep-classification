
"""
Take a sequence of M labels and return a sequence of M labels by using very
cheap interpolation
"""
function fit_labels(labels::Vector{Int}, N::Int)
    M   = length(labels)
    seq = round(Int, linspace(1, M, N))
    labels[seq]
end


if(false)
    Root = "/home/mark/current/general-sleep/"

    Feats  = Matrix{Float64}[]
    Labels = Vector{Int}[]
    for i=1:26
        println("subject ", i)
        x = h5open(string(Root, "patient$i.h5"))
        dat = map(Float64,read(x["2"]))
        F = make_features(dat)
        L = map(Int,readcsv(string(Root, "/patients/HypnogramAASM_patient$i.txt"))[2:end])[1:6:end]
        push!(Feats,  F)
        push!(Labels, L)
    end
end

if(false)
    Root = "/home/mark/code/sleep-classification/tmp/"

    Feats  = Matrix{Float64}[]
    Labels = Vector{Int}[]
    for i=1:26
        println("subject ", i)
        F = Images.imread(string(Root,"DejunkedSpectra$i.png"))[1:1000,:]'
        L = map(Int,readcsv(string(Root, "Dejunkedlabels$i.csv")))[:]
        push!(Feats,  F)
        push!(Labels, L)
    end
end

#zero mean and 2nd std is at +-1
function histNormalize(img)
    t = img-mean(img)
    map(Float64,t./(3*std(t)))
end

function doMakeThresholds(SubjectID; workingDir="tmp", doPlot=false)
    I = PyPlot.imread("$workingDir/DejunkedSpectra$SubjectID.png")

    #figure(1)
    #imshow(I, aspect="auto", interpolation="none")
    #
    #figure(2)
    #imshow(Median3(I), aspect="auto", interpolation="none")
    #
    rmed = recursiveMedian(I)
    if(doPlot)
        figure(201)
        title("Recursive Median Results")
        imshow(rmed, aspect="auto", interpolation="none")
    end

    hn = histNormalize(rmed)
    high = GeneralMedian(threshold(hn,  0.1,  Inf), 5)
    med  = GeneralMedian(threshold(hn, -0.2,  0.2), 5)
    low  = GeneralMedian(threshold(hn, -Inf, -0.1), 5)
    imwrite(high, "$workingDir/high$SubjectID.png")
    imwrite(med,  "$workingDir/med$SubjectID.png")
    imwrite(low,  "$workingDir/low$SubjectID.png")
    combined = zeros(size(high)[1],size(high)[2],3)
    combined[:,:,1] = high
    combined[:,:,2] = med
    combined[:,:,3] = low
    if(doPlot)
        figure(202)
        title("Regions")
        imshow(combined, aspect="auto", interpolation="none")
        figure(203)
        title("Histogram")
        PyPlot.plt[:hist](hn[:],128)
    end
end

if(false)
    for i=1:26
        doMakeThresholds(i,doPlot=true)
        doRectSegment(i, "low",  100, "tmp", true)
        doRectSegment(i, "med",  200, "tmp", true)
        doRectSegment(i, "high", 300, "tmp", true)
    end
end

if(true)
    for i=1:26
        @time doLikelyHoodEst(i,"tmp", true)
    end
end


#50 trees 20 feats
raw_data_acc = [1 0.628
2 0.212
3 0.307
4 0.722
5 0.414
6 0.648
7 0.62
8 0.382
9 0.617
10 0.733
11 0.742
12 0.674
13 0.611
14 0.592
15 0.558
16 0.573
17 0.377
18 0.703
19 0.4
20 0.493
21 0.409
22 0.515
23 0.644
24 0.635
25 0.46
26 0.684]


bandpow = [2 0.7829698857736241
 3 0.2333984375
 4 0.5945205479452055
 5 0.7024952015355086
 6 0.5592039800995025
 7 0.6017964071856288
 8 0.5623100303951368
 9 0.5363716038562665
 10 0.6939782823297137
 11 0.7252543940795559
 12 0.736793327154773
 13 0.639344262295082
 14 0.639344262295082
 15 0.6484962406015038
 16 0.531578947368421
 17 0.7379912663755459
 18 0.27053571428571427
 19 0.6696
 20 0.6935933147632312
 21 0.45390070921985815
 22 0.5396966993755575
 23 0.5248868778280543
 24 0.8082051282051282
 25 0.6483279395900755
 26 0.6098239110287303
 27 0.7632367632367633]


#with feats=4
cfc_one = 
[1 0.48598130841121495
 2 0.4296875
 3 0.38447488584474887
 4 0.5431861804222649
 5 0.44676616915422884
 6 0.4221556886227545
 7 0.5096251266464032
 8 0.4136722173531989
 9 0.5350444225074038
 10 0.7335800185013877
 11 0.5801668211306765
 12 0.6084860173577628
 13 0.45389344262295084
 14 0.48026315789473684
 15 0.4494736842105263
 16 0.6777292576419214
 17 0.3991071428571429
 18 0.7064
 19 0.5088207985143919
 20 0.5825734549138805
 21 0.5066904549509367
 22 0.4796380090497738
 23 0.5764102564102564
 24 0.5577130528586839
 25 0.5430954587581094
 26 0.4865134865134865]

#with feats=8
cfc_two =
[1 0.5129802699896158
 2 0.4140625
 3 0.3771689497716895
 4 0.519193857965451
 5 0.445771144278607
 6 0.4281437125748503
 7 0.5116514690982776
 8 0.4154250657318142
 9 0.5409674234945706
 10 0.7354301572617946
 11 0.5903614457831325
 12 0.6113789778206364
 13 0.45799180327868855
 14 0.48026315789473684
 15 0.4568421052631579
 16 0.6838427947598253
 17 0.4205357142857143
 18 0.6976
 19 0.5134633240482822
 20 0.574468085106383
 21 0.512042818911686
 22 0.4832579185520362
 23 0.5743589743589743
 24 0.5555555555555556
 25 0.5495829471733086
 26 0.4875124875124875]
