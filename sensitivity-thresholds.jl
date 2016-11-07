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

#doSense()
#include("misc.jl")
#drange  =  data_range()[38:end]
#FF = Matrix{Float64}[]
#LL = Vector{Int}[]
#if(true)
#    for j=drange
#        if(j != "ST7221" && j != subject)
#            tmp = viewStuff(j, "physionet/", 0.85, false)
#            push!(LL, map(Int, tmp[1][:]))
#            push!(FF, tmp[2])
#        end
#    end
#end

#figure(1010101)
#model = build_forest(vcat(LL...), hcat(FF...)', 20, 40)
#PyPlot.close("all")

#Test C block gradient and coarse bands
#ttt = Any[]
ii=0
for i=thresh_bias, j=thresh_width
    ii += 1
    if(length(ttt) > ii)
        continue
    end
    push!(ttt, Any[(i,j), viewStuff(string(subject,"-$i-$j"),
            "physionet-thresh/", 0.85, true, ["etc"],
            x_thresh=0.2, y_thresh=0)])
end

class_result = Vector{Int}[]
for t=ttt
    tmp = apply_forest(model, t[2][2]')
    println("classification accuracy ", t[1], " ", mean(t[2][1].==tmp))
    push!(class_result, tmp)
end

min_img = nothing
med_img = nothing
hgh_img = nothing
for i=thresh_bias, j=thresh_width
    if(sns[find((sns[:,1] .== i)&(sns[:,2].==j))[1],3] > 0.89)
        continue
    end
    l=PyPlot.imread("physionet-thresh/coverlow$subject-$i-$j.png")
    m=PyPlot.imread("physionet-thresh/covermed$subject-$i-$j.png")
    h=PyPlot.imread("physionet-thresh/coverhigh$subject-$i-$j.png")
    if(min_img == nothing)
        min_img = l
        med_img = m
        hgh_img = h
    else
        min_img += l
        med_img += m
        hgh_img += h
    end
end
