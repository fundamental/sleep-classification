using Images
using PyPlot

getX(img) = size(img)[2]
getY(img) = size(img)[1]
px(img, x, y) = img[Int(max(1,min(getY(img),y))), Int(max(1,min(getX(img),x)))]
function Median3(img)
    result = copy(img)

    for x=1:getX(img), y=1:getY(img)
        opts = [
        px(img,x-1,y-1),px(img,x+0,y-1),px(img,x+1,y-1),
        px(img,x-1,y-0),px(img,x+0,y-0),px(img,x+1,y-0),
        px(img,x-1,y+1),px(img,x+0,y+1),px(img,x+1,y+1)]
        result[y,x] = median(opts)
    end
    result
end

function GeneralMedian(img, med)
    #return img
    @assert isodd(med)
    scale = (med-1)/2
    result = copy(img)

    for x=1:getX(img), y=1:getY(img)
        opts = Float64[]
        for i=-scale:scale, j=-scale:scale
            push!(opts, px(img,x+i,y+j))
        end
        result[y,x] = median(opts)
    end
    result
end


function recursiveMedian(img)
    for i=1:20
        img = Median3(img)
    end
    img
end

#zero mean and 2nd std is at +-1
function histNormalize(img)
    t = img-mean(img)
    t./(3*std(t))
end

function threshold(img, low, high)
    t = copy(img)
    t[img.>low] = 1
    t[img.<low]  = 0
    t[img.>high] = 0
    t
end
 
function labelRow(row)
    i=0
    rowLabel = zeros(length(row))
    currentState = row[1]
    for x=2:length(row)
        if(currentState != row[x])
            currentState = row[x]
            i+=1
        end
        rowLabel[x] = i
    end
    int(rowLabel)
end

function mergeEquivSet(rename, a::Int, b::Int)
    #Find a set for a/b
    rA = rename[a]
    rB = rename[b]

    #Merge two sets
    if(rA != rB)
        dest = min(rA,rB)
        src  = max(rA,rB)
        for x=keys(rename)
            if(rename[x] == src)
                rename[x] = dest
            end
        end
    end
end

function initialRename(N)
    rename = Dict{Int,Int}()
    for i=0:N
        rename[i] = i
    end
    rename
end

function findEquiv(rename, imgCol, labelCol)
    curClass = imgCol[1]
    curLabel = labelCol[1]
    for y=2:length(labelCol)
        if(curClass == imgCol[y])
            mergeEquivSet(rename, int(curLabel), int(labelCol[y]))
        else
            curClass = imgCol[y]
            curLabel = labelCol[y]
        end
    end
end

function labelClusters(img)
    #Give each row a unique labeling
    labelFirst = mapslices(labelRow, img,2)
    #apply accumulator
    for i=2:size(labelFirst)[1]
        labelFirst[i,:] += labelFirst[i-1,end]
    end

    renameTable = initialRename(int(maximum(labelFirst)))

    for i=1:size(img)[2]
        print(".")
        findEquiv(renameTable, img[:,i],labelFirst[:,i])
    end

    #Look at the image column wise and note
    #identifiers that are to be merged

    #Optimise renameTable to result in the sequence 1:N
    #Where N is the number of elements which rename n->n

    #Perform Renumbering
    for x=1:getX(img),y=1:getY(img)
        labelFirst[y,x] = renameTable[int(labelFirst[y,x])]
    end
    labelFirst

end


function doMakeThresholds(SubjectID, workingDir, doPlot)
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
        plt.hist(hn[:],128)
    end
end
