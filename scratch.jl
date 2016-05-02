using JLD
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

if(false)
    for i=1:26
        @time doLikelyHoodEst(i,"tmp", true)
    end
end

function viewStuff(SubjectID, workDir, doPlot=true, misc=false)
    ep  = readcsv("$workDir/edgeParameter$SubjectID.csv")
    dat = getSpectra(SubjectID, workDir)[1:1000,:]
    ll  = readcsv("$workDir/Dejunkedlabels$SubjectID.csv")
    if(doPlot)
        figure(0);PyPlot.clf();
        plot(ll)
        figure(1);PyPlot.clf();
        imshow(dat,aspect="auto",interpolation="none")
        figure(2);PyPlot.clf();
        plot(maximum(ep,2)[:])
    end

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
    iimax = zoneify(dat, X, Y, operator=maximum)
    iimin = zoneify(dat, X, Y, operator=minimum)
    out = (0.4*dat+mn+md)/2.01
    #out = dat

    if(doPlot)
        figure(3);PyPlot.clf();
        imshow(mn,aspect="auto",interpolation="none")

        figure(4);PyPlot.clf();
        imshow(md,aspect="auto",interpolation="none")

        figure(5);PyPlot.clf();
        imshow(iimax,aspect="auto",interpolation="none")

        figure(6);PyPlot.clf();
        imshow(iimin,aspect="auto",interpolation="none")
        figure(7);PyPlot.clf();
        imshow(out,aspect="auto",interpolation="none");

        figure(8);PyPlot.clf();
        imshow(mapslices(sortperm,dat,1),aspect="auto",interpolation="none");
        figure(9);PyPlot.clf();
        imshow(mapslices(sortperm,out,1),aspect="auto",interpolation="none");
    end
    (ll, dat, mn, md, out)
end

function getSpectra(SubjectID, workDir)
    Sp = raw_spec(SubjectID)
    Cl = readcsv("$workDir/Dejunked$SubjectID-cols.csv")
    Sp[:,find(Cl)]
    #imread("$workDir/DejunkedSpectra$SubjectID.png")
end

function collectThemAll()
    result = Any[]
    for i=1:26
        stuff = viewStuff(i, "tmp")
        push!(result, (stuff[5], stuff[1]))
    end
    result
end

#results = collectThemAll();
#save("patients-april.jld", "results", results)
results = load("patients-april.jld")["results"]

if(true)
    Feats  = Matrix{Float64}[]
    Labels = Vector{Int}[]
    for i=1:26
        push!(Feats, 2.0PyPlot.imread("tmp/coverhigh$i.png")+
                    (PyPlot.imread("tmp/covermed$i.png")-0.5)-
                     PyPlot.imread("tmp/coverlow$i.png"))
        #push!(Feats,  results[i][1])
        push!(Labels, round(Int, results[i][2])[:])
    end
    results = nothing
    #cross_validate(Feats, Labels)
end

#viewStuff(22, "tmp")

#Extract data from plots
if(false)
    Expected   = Vector{Int}[]
    Classified = Vector{Int}[]
    for i=1:26
        figure(i)
        push!(Expected,   round(Int,gca()[:get_lines]()[1][:get_ydata]()))
        push!(Classified, round(Int,gca()[:get_lines]()[2][:get_ydata]()))
    end
    save("patients-april-cover.jld", "expected", Expected, "classified",
    Classified)
end

cf = make_conf(vcat(Expected...), vcat(Classified...), [1,2,3,4,5])
println("total mean accuracy = ", 100sum(diag(cf))/sum(cf), "%")
for i=1:26
    cf = make_conf(Expected[i], Classified[i], [1,2,3,4,5])
    cf2 = make_conf(Expected[i], Classified[i], [1,2,4,5])
    println("sub $i mean accuracy = ", round(Int, 100sum(diag(cf))/sum(cf)), "% or ",
    round(Int, 100sum(diag(cf2))/sum(cf2)), "% (base ", round(Int,
    100mean(Expected[i].==Classified[i])), "%)")
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

 #40 tree 20 feat
phd_meth =
[1 0.7390282131661442
 2 0.555956678700361
 3 0.7433155080213903
 4 0.7376453488372093
 5 0.6068066618392469
 6 0.8169219547775346
 7 0.5963855421686747
 8 0.4929305912596401
 9 0.8301886792452831
 10 0.8344327176781002
 11 0.7075471698113207
 12 0.7524613220815752
 13 0.6041512231282431
 14 0.718052738336714
 15 0.5466562986003111
 16 0.7648562300319489
 17 0.218
 18 0.8216636744464393
 19 0.7362045760430687
 20 0.513677811550152
 21 0.40992167101827676
 22 0.5604681404421327
 23 0.8602620087336245
 24 0.7040816326530612
 25 0.3807531380753138
 26 0.7427536231884058]


 #40 tree 20 feat
phd_cover = 
[1 0.6504702194357367
 2 0.3148014440433213
 3 0.6818181818181818
 4 0.8023255813953488
 5 0.5843591600289645
 6 0.8738147337709701
 7 0.6822289156626506
 8 0.4235218508997429
 9 0.7851959361393324
 10 0.8773087071240105
 11 0.7405660377358491
 12 0.7137834036568214
 13 0.6501111934766494
 14 0.691683569979716
 15 0.6026438569206843
 16 0.8849840255591054
 17 0.316
 18 0.7779772591262717
 19 0.6480484522207268
 20 0.5949848024316109
 21 0.5313315926892951
 22 0.5838751625487646
 23 0.888646288209607
 24 0.6766091051805337
 25 0.608786610878661
 26 0.7630434782608696]

