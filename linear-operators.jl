function kClassify(data::Matrix{Float64}, centers::Int; plotID::Int=-1,
    maxiter::Int=40)
    k = kmeans(data, centers)

    res = Float64[]
    for i=1:size(data)[2]
        choice = indmin(mapslices(norm, k.centers.-data[:,i],1))
        push!(res, choice)
    end

    if(plotID >= 0)
        figure(plotID)
        plot(res+0.01randn(length(res)))
        figure(plotID+1)
        plot(res[1:end-1].!=res[2:end])
    end
    res
    (k.centers, res)
end
    
function kClassify_nway(data, centers; plotID=-1, maxiter=20)
    k = kmeans(data, centers, maxiter=maxiter)
    return k.assignments
    freq = zeros(centers)
    res = Float64[]
    for i=1:size(data)[2]
        choice = k.assignments[i];#indmin(sum((k.centers.-data[:,i]).^2,1))
        push!(res, choice)
        freq[choice] += 1
    end
    per = reverse(sortperm(freq))
    for i=1:size(data)[2]
        res[i] = per[round(Int,res[i])]
    end

    #figure(plotID)
    #plot(res+0.01randn(length(res)))
    #figure(plotID+1)
    #plot(res[1:end-1].!=res[2:end])
    res
end

"""
Perform a median over a 1D vector

Arguments:

- vec - 1D vector
- med - window of the vector
"""
function linearMedian{T}(vec::Vector{T}, med::Int)
    #return img
    @assert isodd(med)
    scale = (med-1)/2
    result = copy(vec)
    N = length(vec)
    getIt(x) = Int(x<1?1:x>N?N:x)

    for x=1:length(vec)
        opts = Float64[]
        for i=-scale:scale
            push!(opts, vec[getIt(x+i)])
        end
        result[x] = median(opts)
    end
    result
end

