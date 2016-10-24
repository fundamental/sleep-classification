#Assume overlap parameter is constant
thresh_overlap = 1.5
thresh_bias    = -1.0:0.1:1.0
thresh_width   = 0.2:0.1:2.0
subject = "ST7192"

function doSense()
    SubjectID  = subject
    workingDir = "physionet"
    I::Matrix{Float64} = PyPlot.imread("$workingDir/DejunkedSpectra$SubjectID.png")

    rmed = recursiveMedian(I)

    hn = histNormalize(rmed)
        
    ll = round(Int,readcsv("$workingDir/Dejunkedlabels$SubjectID.csv")[:])

for i=thresh_bias, j=thresh_width
    if(get_match(readdir("physionet-thresh/"), Regex("edgeParameter$SubjectID-$i-$j.csv")) !=
        nothing)
        println("Skipping i=$i j=$j")
        continue
    end
    print("i=$i")
    print("j=$j")
    i/5.0
    j/5.0

    l1   = i/5-j/5.0
    l2   = i/5-j/5.0*thresh_overlap
    h1   = i/5+j/5.0
    h2   = i/5+j/5.0*thresh_overlap
    high = GeneralMedian(threshold(hn,  h1,  Inf), 5)
    med  = GeneralMedian(threshold(hn,  l2,  h2), 5)
    low  = GeneralMedian(threshold(hn, -Inf, l1), 5)
    Images.save("physionet-thresh/high$SubjectID-$i-$j.png", high)
    Images.save("physionet-thresh/med$SubjectID-$i-$j.png",  med)
    Images.save("physionet-thresh/low$SubjectID-$i-$j.png",  low)
    doRectSegment(  string(subject,"-$i-$j"), "low",  100, "physionet-thresh/", true)
    doRectSegment(  string(subject,"-$i-$j"), "med",  200, "physionet-thresh/", true)
    doRectSegment(  string(subject,"-$i-$j"), "high", 300, "physionet-thresh/", true)
    doLikelyHoodEst(string(subject,"-$i-$j"),"physionet-thresh/", false,
                    label_input=ll, img_input=I)
    #viewStuff(      string(subject,"-$i-$j"), "physionet-thresh/")
end
end

doSense()

