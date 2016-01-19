using Images
using PyPlot

function area(rect::Vector{Int})
    abs(1+rect[2]-rect[1])*abs(1+rect[4]-rect[3])
end

function inbounds(rect::Vector{Int}, img::BitArray{2})
    t1 = rect[1] < 4 || rect[2] < 4 || rect[3] < 4 || rect[4] < 4
    c2 = size(img,2)-2
    t2 = rect[1] > c2 || rect[2] > c2
    c3 = size(img,1)-2
    t3 = rect[3] > c3 || rect[4] > c3

    !(t1 || t2 || t3)
end

function rawScore(rect::Vector{Int}, img)
    score = 0
    for r=rect[3]:rect[4], c=rect[1]:rect[2]
        score += img[r,c]
    end
    score
end

function score(rect::Vector{Int},img)
    if(!inbounds(rect,img))
        return -1000
    end
    raw = rawScore(rect,img)
    raw - 2abs(raw-area(rect))
end

function validInit(rect::Vector{Int},img::BitArray{2})
    inbounds(rect, img) && rawScore(rect,img) == area(rect)
end
function randInit(img)
    y = rand(1:size(img,1))
    x = rand(1:size(img,2))
    [x,x+1,y,y+1]
end

#rectangle [x0,x1,y0,y1]
function within(rect::Vector{Int}, prev)
    r = rect
    for i=prev
        if(i[1]<=r[1] && i[2]>=r[2] && i[3]<=r[3] && i[4]>=r[4])
            print("%")
            return true
        end
    end
    false
end

function fixRect(rect::Vector{Int})
    if(rect[1]>rect[2])
        fixRect(rect[[2,1,3,4]])
    elseif(rect[3]>rect[4])
        fixRect(rect[[1,2,4,3]])
    else
        rect
    end
end

#R1 is the old rectangle
#R2 is the new rectangle
#Return the differing rectangle
#Return a negative sign if it's lost area
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

function scoreRect(img, r1::Vector{Int}, s1::Int, r2::Vector{Int})
    #println("scoreRect()")
    if(r1 == r2)
        s1
    else
        #println("Doing Math...")
        (sign_, delta) = deltaRect(r1,r2)
        #println("sign=",sign_)
        #println("delta=",delta)
        s1 + sign_*score(delta,img)
    end
end


function generateRect(img::BitArray{2}, prev::Vector{Vector{Int}})
    while true
        #print(".")
        #get a random initialization point
        rect = Int[0,0,0,0]
        while(!validInit(rect,img) && !within(rect,prev))
            rect = randInit(img)
            #println("rect = ", rect)
        end

        #Find best expansion
        #TODO this can be optimized as
        #the delta score can be calculated rather than the total score each time
        tmp = Int[0,0,0,0]
        for i=1:1000
            bestExpand::Int = 0
            bestDir::Int    = 0
            bestDelta::Int  = 0
            oldScore = score(rect,img)
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
            rect[bestDir] += bestDelta
            #print("&")
        end

        if(area(rect) >= 20 && !within(rect, prev))
            return rect
        end
    end
end

#Solve the Approximate Set Cover Problem Greedly
function intervalIntersect(l1,l2)
    #Avoiding intersection means that l1[*]>l2[*] || l1[*]<l2[*]
    !(l1[2] < l2[1] || l1[1]>l2[2])
end

function intersect(r1, r2)
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
    
function rectSegment(Img, figNum=100, doPlot=false)
    R = Vector{Vector{Int}}()
    tic()
    for i=1:1024
        nx = generateRect(Img,R)
        if(within(nx, R))
            print("!!!")
        end
        push!(R, nx)
        
        print("*")#;area(R[end]))
    end
    toc()


    calculateInitialCover(rect)=area(rect)

    Cover = Array(Bool, size(Img,1), size(Img,2))
    Cover[:] = false
    position = Array(Int, length(R))
    position[:] = 2000
    scores   = Array(Float64, length(R))
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
        print("*")
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
    figure(101)
    PyPlot.clf();
    imshow(I, aspect="auto",interpolation="none")
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

function doRectSegment(SubjectID, state, figNum, workingDir, doPlot)
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


function findCover(sz, R)
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
