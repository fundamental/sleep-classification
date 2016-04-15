using Images
using PyPlot

"""
Calculate the area of a given rectangle
"""
function area(rect::Vector{Int})
    abs(1+rect[2]-rect[1])*abs(1+rect[4]-rect[3])
end

"""
Identify if a rectangle is within the bounds of a given image
"""
function inbounds(rect::Vector{Int}, img::BitArray{2})
    t1 = rect[1] < 4 || rect[2] < 4 || rect[3] < 4 || rect[4] < 4
    c2 = size(img,2)-2
    t2 = rect[1] > c2 || rect[2] > c2
    c3 = size(img,1)-2
    t3 = rect[3] > c3 || rect[4] > c3

    !(t1 || t2 || t3)
end

"""
Count the positive pixels within a rectangle
"""
function rawScore(rect::Vector{Int}, img)
    score::Int = 0
    for r=rect[3]:rect[4], c=rect[1]:rect[2]
        score += img[r,c]
    end
    score
end

"""
Produce a value to represent how good a candidate rectangle is
Larger values are better
"""
function score(rect::Vector{Int},img)
    if(!inbounds(rect,img))
        return -1000
    end
    raw = rawScore(rect,img)
    raw - 2abs(raw-area(rect))
end

"""
Identify when a rectangle has a valid initialization
"""
function validInit(rect::Vector{Int},img::BitArray{2})
    inbounds(rect, img) && rawScore(rect,img) == area(rect)
end

"""
Create a random initialization which may or may not be valid
"""
function randInit(img)
    y = rand(1:size(img,1))
    x = rand(1:size(img,2))
    [x,x+1,y,y+1]
end

#rectangle [x0,x1,y0,y1]
"""
Identify if a rectangle is within an instance of a rectangle from another set of
    rectangles
"""
function within(rect::Vector{Int}, prev::Vector{Vector{Int}})
    r = rect
    for i=prev
        if(i[1]<=r[1] && i[2]>=r[2] && i[3]<=r[3] && i[4]>=r[4])
            print("%")
            return true
        end
    end
    false
end

"""
Ensure that a rectangle vector R=[a, b, c, d]
maintains the property

R.a <= R.b
R.c <= R.d
"""
function fixRect(rect::Vector{Int})
    if(rect[1]>rect[2])
        fixRect(rect[[2,1,3,4]])
    elseif(rect[3]>rect[4])
        fixRect(rect[[1,2,4,3]])
    else
        rect
    end
end

"""
R1 is the old rectangle
R2 is the new rectangle
Return the differing rectangle
Return a negative sign if it's lost area
"""
function deltaRect(r1::Vector{Int}, r2::Vector{Int})
    #println("deltaRect(",r1,",",r2,")")
    sign_ = area(r1)<area(r2) ? 1 : -1
    #println("sign = ", sign_)
    rd::Vector{Int} = abs(r1-r2)
    #println("rd=",rd)
    if(rd[1] == 1)
        (sign_, fixRect([r2[1], r1[1], r1[3],r1[4]]))
    elseif(rd[2] == 1)
        (sign_, fixRect([r2[2], r1[2], r1[3],r1[4]]))
    elseif(rd[3] == 1)
        (sign_, fixRect([r1[1], r1[2], r2[3],r2[3]]))
    elseif(rd[4] == 1)
        (sign_, fixRect([r1[1], r1[2], r2[4],r2[4]]))
    else
        nothing
    end
end

"""
Generate a score for a rectangle on a given image
    a previous rectangle score is used to speed up the process
"""
function scoreRect(img::BitArray{2}, r1::Vector{Int}, s1::Int, r2::Vector{Int})
    if(r1 == r2)
        s1
    else
        (sign_, delta) = deltaRect(r1,r2)
        s1 + sign_*score(delta,img)
    end
end

"""
Generate a candidate rectangle from a binary image and a collection of
previously generated rectangles
"""
function generateRect(img::BitArray{2}, prev::Vector{Vector{Int}})
    for iters=1:30
        #print(".")
        #get a random initialization point
        rect = Int[0,0,0,0]
        while(!validInit(rect,img))# && !within(rect,prev))
            rect = randInit(img)
            #println("rect = ", rect)
        end

        #Find best expansion
        #TODO this can be optimized as
        #the delta score can be calculated rather than the total score each time
        tmp = Int[0,0,0,0]
        oldScore = score(rect, img)
        for i=1:2000
            bestExpand::Int = 0
            bestDir::Int    = 0
            bestDelta::Int  = 0
            #oldScore = score(rect,img)
            for dir=1:4
                for delta=[(-1:1);]
                    tmp[:] = rect[:]
                    tmp[dir] += delta
                    s::Int = 0
                    #if(!within(tmp, prev))
                        s = scoreRect(img, rect, oldScore, tmp)
                    #end

                    if(s > bestExpand || (s == bestExpand && rand() > 0.5))
                        bestExpand = s
                        bestDir    = dir
                        bestDelta  = delta
                    end
                end
            end
            if(bestDelta == 0)
                break
            end
            oldScore = bestExpand
            rect[bestDir] += bestDelta
            #print("&")
        end

        if(area(rect) >= 20 && !within(rect, prev))
            return rect
        end
    end
    nothing
end

#Solve the Approximate Set Cover Problem Greedly

"""
Identify when two intervals overlap
"""
function intervalIntersect(l1::Vector{Int},l2::Vector{Int})
    #Avoiding intersection means that l1[*]>l2[*] || l1[*]<l2[*]
    !(l1[2] < l2[1] || l1[1]>l2[2])
end

"""
Identify when rectangles intersect
"""
function intersect(r1::Vector{Int}, r2::Vector{Int})
    intervalIntersect(r1[1:2],r2[1:2]) && intervalIntersect(r1[3:4],r2[3:4])
end

function calculateCover(Cover, rect, prev, damage)
    if(prev == 0)
        return 0
    end
    if(!intersect(rect, damage))
        return prev
    end
    score_ = 0
    for r=rect[3]:rect[4], c=rect[1]:rect[2]
        score_ += 1*!Cover[r,c]
    end
    score_ + 50*score_/area(rect)
end

"""
Perform rectangle segmentation on an input image

output - a collection of rectangles in the form
[[x1 x2 y1 y2]
 [x1 x2 y1 y2]
 .
 .
 .
 [x1 x2 y1 y2]]
"""
function rectSegment(Img::BitArray{2}, figNum::Int=100, doPlot::Bool=false)
    R = Vector{Vector{Int}}()
    #tic()
    for i=1:4096
        nx = generateRect(Img,R)
        if(nx == nothing)
            break
        end
        if(within(nx, R))
            print("!!!")
        end
        push!(R, nx)

        print("*")#;area(R[end]))
    end
    #toc()


    calculateInitialCover(rect)=area(rect)

    Cover = Array(Bool, size(Img,1), size(Img,2))
    Cover[:] = false
    position = Array(Int, length(R))
    position[:] = 2000
    scores   = Array(Float64, length(R))
    println("available rectangles = ", length(R))
    for i=1:length(R)
        scores[i] = calculateInitialCover(R[i])
    end

    damageRect = [1,size(Img,2),1,size(Img,1)]
    #println(damageRect)
    #TODO this can be done *MUCH* faster by only updating
    #the score of intersecting rectangles
    #Also, the cover of some element is strictly upper bound
    #by the previous value
    #Thus a zero element is zero for all future calc
    for i=1:400
        print("=")
        for j=1:length(R)
            scores[j] = calculateCover(Cover, R[j], scores[j], damageRect)
        end
        best = reverse(sortperm(scores))[1]
        if(scores[best] < 100 && sum(position.<1000) > 20)
            break
        end
        position[best] = i
        scores[best] = 0
        rect = R[best]
        damageRect = rect
        for r=rect[3]:rect[4], c=rect[1]:rect[2]
            Cover[r,c] = true
        end
    end

    if(doPlot)
        figure(figNum)
        #title("Cover for $state")
        PyPlot.clf()
        imshow(Cover, aspect="auto", interpolation="none")
        #imwrite(float64(Cover), "$workingDir/cover$state$SubjectID.png")
    end

    #figure(3)
    #plt.clf()
    I = zeros(size(Img))#copy(Img)
    I[:,:] = 0
    total_area = 0
    for i=1:length(R)
        x = 1#rand()#0.5+0.2*(rand()-0.5)
        rect = R[i]
        for r=rect[3]:rect[4], c=rect[1]:rect[2]
            I[r,c] += x
            total_area += 1
        end
    end
    real_area = sum(I.== 1)
    if(doPlot)
        figure(figNum+1)
        PyPlot.clf();
        imshow(I, aspect="auto",interpolation="none")
    end
    println("Total Area = ", total_area)
    println("Real  Area = ", real_area)

    #II = copy(Img)
    #for i=1:length(R)
    #    if(position[i] < 1000)
    #        x = 0.25#rand()#0.5+0.2*(rand()-0.5)
    #        rect = R[i]
    #        for r=rect[3]:rect[4], c=rect[1]:rect[2]
    #            II[r,c] += x
    #        end
    #    end
    #end

    result = Vector{Vector{Int}}()
    for i=1:length(R)
        if(position[i] < 1000)
            push!(result, R[i])
        end
    end
    result
end

function doRectSegment(SubjectID::Int, state::ASCIIString, figNum::Int, workingDir::ASCIIString, doPlot::Bool)
    #File = "before-labeling.png"
    #File = "real-example.png"
    ##File = "real-example2.png"
    ##File = "real-example3.png"
    #File = "DejunkedSpectra1-high.png"
    #File = "DejunkedSpectra1-low.png"
    #State = "low"
    #SubjectID = 1
    File = "$workingDir/$state$SubjectID.png"

    img = map(float, PyPlot.imread(File))
    Img = map(x->x.>0.5, img)
    result = rectSegment(Img, figNum, doPlot)

    writecsv("$workingDir/interest$SubjectID-$state.csv", hcat(result...)')
end

"""
Identify the space that a collection of rectangles cover

sz - input matrix size
R  - collection of input rectangles

outputs a binary image indicating if a rectangle covers a pixel or not
"""
function findCover(sz::Tuple{Int,Int}, R::Vector{Vector{Int}})
    I = zeros(sz)#copy(Img)
    I[:,:] = 0
    total_area = 0
    for i=1:length(R)
        x = 1#rand()#0.5+0.2*(rand()-0.5)
        rect = R[i]
        for r=rect[3]:rect[4], c=rect[1]:rect[2]
            I[r,c] = x
            total_area += 1
        end
    end
    I
end


function demo_cover()

    sz = (800,800)
    rects = Vector{Int}[]
    push!(rects,[8,20,100,500])
    push!(rects,[400,500,450,700])
    push!(rects,[600,700,450,550])
    push!(rects,[100,450,200,280])


    cover = findCover(sz, rects)

    c2    = cover+1.0randn(sz)
    c2    = c2.>0.3

    R = rectSegment(c2, 100, true)

    c3    = findCover(sz, R)

    rpow = zeros(sz[1])
    cpow = zeros(sz[2])
    for r=R
        rpower = abs(r[2]-r[1])
        cpower = abs(r[4]-r[3])
        cpow[r[1]] += cpower
        cpow[r[2]] += cpower
        rpow[r[3]] += rpower
        rpow[r[4]] += rpower
    end

    seg_likelyhood = zeros(sz)
    for i=1:sz[1],j=1:sz[2]
        seg_likelyhood[i,j] = rpow[i] + cpow[j]
    end

    figure(1)
    imshow(cover)

    figure(2)
    imshow(c2)

    figure(3)
    imshow(c3)

    figure(4)
    imshow(seg_likelyhood)

    figure(5)
    imshow(cov(cover))
    figure(6)
    imshow(cov(c2))
    figure(7)
    imshow(cov(c3))
end

function seg_likely(sz, R)
    rpow = zeros(sz[1])
    cpow = zeros(sz[2])
    for r=R
        rpower = abs(r[2]-r[1])
        cpower = abs(r[4]-r[3])
        cpow[r[1]] += cpower
        cpow[r[2]] += cpower
        rpow[r[3]] += rpower
        rpow[r[4]] += rpower
    end

    seg_likelyhood = zeros(sz)
    for i=1:sz[1],j=1:sz[2]
        seg_likelyhood[i,j] = rpow[i] + cpow[j]
    end
    seg_likelyhood
end

function seg_seg(rpow, cpow)
    sz = (length(rpow), length(cpow))
    seg_likelyhood = zeros(sz)
    for i=1:sz[1],j=1:sz[2]
        seg_likelyhood[i,j] = rpow[i] + cpow[j]
    end
    seg_likelyhood
end

#FOR SYNTHETIC DATA ONLY
function seg_true(sz, synth)
    rpow = zeros(sz[1])
    cpow = zeros(sz[2])
    for i=2:sz[1]
        if(synth[i-1,:] != synth[i,:])
            rpow[i] = 1
        end
    end
    for i=2:sz[2]
        if(synth[:,i-1] != synth[:,i])
            cpow[i] = 1
        end
    end

    seg_likelyhood = zeros(sz)
    for i=1:sz[1],j=1:sz[2]
        seg_likelyhood[i,j] = rpow[i] + cpow[j]
    end
    seg_likelyhood
end

function demo_triple_cover(seed=rand(Int64))
    println("seed = ", seed)
    srand(seed)


    sz = (800,1200)
    general_noise_level = 0.1
    sparse_noise_level  = 8.0
    sparse_noise_cnt    = 0.15


    data_raw      = make_irregular_grid(height=sz[1], width=sz[2], levels=[-1.0,0.0,1.0,2.0])
    general_noise = general_noise_level*randn(sz...)
    sparse_noise  = zeros(sz...)

    for i=1:length(sparse_noise)
        sparse_noise[i] += (sparse_noise_cnt > rand()) ?
        sparse_noise_level*randn() : 0
    end
    
    data_in = data_raw + general_noise + sparse_noise
    din_med = recursiveMedian(data_in)

    lvl_low = din_med .< -0.5
    lvl_med = (din_med .> -0.7) & (din_med .< 0.7)
    lvl_hgh = din_med .> 0.5

    Rlow = rectSegment(lvl_low, 100, false)
    Rmed = rectSegment(lvl_med, 200, false)
    Rhgh = rectSegment(lvl_hgh, 300, false)


    figure(1);PyPlot.clf();
    imshow(data_raw,aspect="auto",interpolation="none");title("underlying model")

    
    figure(2);PyPlot.clf();
    imshow(data_in,aspect="auto",interpolation="none");title("model + noise")
    
    figure(3);PyPlot.clf();
    imshow(din_med,aspect="auto",interpolation="none");title("model + noise + median")

    figure(4);PyPlot.clf();
    imshow(hcat(lvl_low, lvl_med, lvl_hgh));title("perceived levels")

    figure(5);PyPlot.clf();
    imshow(hcat(findCover(sz, Rlow), findCover(sz,Rmed), findCover(sz,Rhgh)));

    figure(6);PyPlot.clf();
    imshow(hcat(seg_likely(sz, Rlow), seg_likely(sz, Rmed), seg_likely(sz,Rhgh)));
    title("seg likelyhood")
    
    figure(7);PyPlot.clf();
    imshow(seg_likely(sz, Rlow)+seg_likely(sz, Rmed)+seg_likely(sz,Rhgh),aspect="auto",interpolation="none");

    figure(8);PyPlot.clf()
    imshow(seg_true(sz, data_raw), aspect="auto", interpolation="none")

    figure(9);PyPlot.clf()
    imshow(findCover(sz,Rhgh)-findCover(sz,Rlow),aspect="auto",interpolation="none")

    figure(10);PyPlot.clf();
    plt[:hist](data_raw[:], 128)
    figure(11);PyPlot.clf();
    plt[:hist](data_in[:], 128)
    figure(12);PyPlot.clf();
    plt[:hist](din_med[:], 128)

    Cover = Matrix{Float64}[findCover(sz, Rlow), findCover(sz,Rmed), findCover(sz,Rhgh)]
    Rects = Vector{Vector{Int}}[Rlow, Rmed, Rhgh]
    (din_med, data_raw, data_in, Cover, Rects)
end
