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

if(false)
    results = load("patients-april.jld")["results"]
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

if(false)
    cf = make_conf(vcat(Expected...), vcat(Classified...), [1,2,3,4,5])
    println("total mean accuracy = ", 100sum(diag(cf))/sum(cf), "%")
    for i=1:26
        cf = make_conf(Expected[i], Classified[i], [1,2,3,4,5])
        cf2 = make_conf(Expected[i], Classified[i], [1,2,4,5])
        println("sub $i mean accuracy = ", round(Int, 100sum(diag(cf))/sum(cf)), "% or ",
        round(Int, 100sum(diag(cf2))/sum(cf2)), "% (base ", round(Int,
        100mean(Expected[i].==Classified[i])), "%)")
    end
end

oo = nothing
if(false)
    println("[INFO] Start")
    #Work on the subjects database from dreams
    Feats  = Matrix{Float64}[]
    Labels = Vector{Int}[]
    for i=1:27
        f = readcsv("/home/mark/current/feat-time-patient$i.csv")
        ll = readcsv("/home/mark/current/general-sleep/patients/HypnogramAASM_patient$i.txt")[2:end]
        #ll = readcsv("/home/mark/current/general-sleep/HypnogramAASM_subject$i.txt")[2:end]
        N = size(f,1)
        M = length(ll)
        l = map(Int,ll[round(Int, linspace(1,M,N))])
        push!(Feats, f')
        push!(Labels, l)
    end
    oo = cross_validate(Feats, Labels, feats=3, trees=100)
    println("[INFO] Stop")
end
if(false)
    println("[INFO] Start")
    #work across different DREAMS databases
    FeatsA  = Matrix{Float64}[]
    LabelsA = Vector{Int}[]
    FeatsB  = Matrix{Float64}[]
    LabelsB = Vector{Int}[]
    for i=1:27
        f = readcsv("/home/mark/current/feat-time-patient$i.csv")
        ll = readcsv("/home/mark/current/general-sleep/patients/HypnogramAASM_patient$i.txt")[2:end]
        #ll = readcsv("/home/mark/current/general-sleep/HypnogramAASM_subject$i.txt")[2:end]
        N = size(f,1)
        M = length(ll)
        l = map(Int,ll[round(Int, linspace(1,M,N))])
        push!(FeatsA,  f')
        push!(LabelsA, l)
    end
    for i=1:20
        f  = readcsv("/home/mark/current/feat-time-subject$i.csv")
        ll = readcsv("/home/mark/current/general-sleep/HypnogramAASM_subject$i.txt")[2:end]
        N = size(f,1)
        M = length(ll)
        l = map(Int,ll[round(Int, linspace(1,M,N))])
        push!(FeatsB,  f')
        push!(LabelsB, l)
    end
    println("time_train_patient_test_subject")
    oo = validate_across(FeatsA, LabelsA, FeatsB, LabelsB, feats=3, trees=100)
    println("time_train_subject_test_patient")
    oo = validate_across(FeatsB, LabelsB, FeatsA, LabelsA, feats=3, trees=100)
    println("[INFO] Stop")
end

###############################
#   DREAMS Subject Database   #
# (I think this is patients   #
#  subjects only has 20)      #
###############################

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

 #40 tree 3 feat
 time_domain_main = [
1 0.7684319833852544
2 0.572265625
3 0.5589041095890411
4 0.6794625719769674
5 0.5442786069651742
6 0.7045908183632734
7 0.2998986828774063
8 0.4794040315512708
9 0.6801579466929911
10 0.7076780758556892
11 0.6431881371640408
12 0.6432015429122468
13 0.5870901639344263
14 0.6475563909774437
15 0.5789473684210527
16 0.6279475982532751
17 0.2642857142857143
18 0.58
19 0.6889507892293407
20 0.3262411347517731
21 0.3434433541480821
22 0.5755656108597285
23 0.7671794871794871
24 0.5782092772384034
25 0.3985171455050973
26 0.7552447552447552
27 0.6076779026217228]

#100 tree 3 feat
time_domain_trees = [
1 0.7757009345794392
2 0.6416015625
3 0.5707762557077626
4 0.6833013435700576
5 0.5263681592039801
6 0.7095808383233533
7 0.3069908814589666
8 0.5004382120946538
9 0.6801579466929911
10 0.7160037002775208
11 0.6459684893419834
12 0.6451301832208293
13 0.5973360655737705
14 0.6400375939849624
15 0.5789473684210527
16 0.6296943231441048
17 0.26160714285714287
18 0.5904
19 0.6898792943361188
20 0.2968591691995947
21 0.3434433541480821
22 0.5683257918552036
23 0.7651282051282051
24 0.5836030204962244
25 0.4198331788693234
26 0.7532467532467533
27 0.5814606741573034]




 #Actually subjects

 #100 tree 3 feat
 time_domain_sub =
 [1 0.6607700312174818
 2 0.7756345177664975
 3 0.6785714285714286
 4 0.6704653371320038
 5 0.6826462128475551
 6 0.7041123370110332
 7 0.7766798418972332
 8 0.6329896907216495
 9 0.7770814682184423
 10 0.5281007751937985
 11 0.6825396825396826
 12 0.7127991675338189
 13 0.6876687668766877
 14 0.6862549800796812
 15 0.7035714285714286
 16 0.6659772492244054
 17 0.7041123370110332
 18 0.7404505386875612
 19 0.5179090029041626
 20 0.7187772925764192]


 #Cross dataset
 #feat=3 tree=100
 time_train_patient_test_subject = [
 1 0.5639958376690947
 2 0.6812182741116751
 3 0.41964285714285715
 4 0.6505223171889839
 5 0.7526366251198466
 6 0.6429287863590772
 7 0.7371541501976284
 8 0.6092783505154639
 9 0.6472694717994628
 10 0.4166666666666667
 11 0.4751984126984127
 12 0.6566077003121749
 13 0.47074707470747074
 14 0.7380478087649402
 15 0.32976190476190476
 16 0.5718717683557394
 17 0.6880641925777332
 18 0.7463271302644466
 19 0.6786060019361084
 20 0.6445414847161572]

 #feat=3 tree=100
 time_train_subject_test_patient = 
 [1 0.671858774662513
 2 0.6943359375
 3 0.5844748858447488
 4 0.7591170825335892
 5 0.6079601990049751
 6 0.48602794411177647
 7 0.3789260385005066
 8 0.45135845749342685
 9 0.736426456071076
 10 0.5124884366327475
 11 0.6626506024096386
 12 0.5978784956605593
 13 0.6618852459016393
 14 0.6541353383458647
 15 0.56
 16 0.5851528384279476
 17 0.26517857142857143
 18 0.544
 19 0.6973073351903436
 20 0.21073961499493415
 21 0.6735057983942908
 22 0.6470588235294118
 23 0.7271794871794872
 24 0.5587918015102481
 25 0.4902687673772011
 26 0.6843156843156843
 27 0.5187265917602997]

