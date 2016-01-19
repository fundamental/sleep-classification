#Features used in Philip Low's Dissertation...

function max_frequency(spectra::Matrix{Float64})
    N       = size(spectra,2)        #Temporal observations
    feat    = zeros(size(spectra,2)) #Output sequence of features
    M       = size(spectra,1)        #Frequency bins
    weights = linspace(0,100,M)      #Range of signal in Hz
    for i=1:N
        feat[i] = weights[indmax(spectra[:,i])]
    end
    feat
end


function central_frequency(spectra::Matrix{Float64}, n::Float64=2.0)
    N       = size(spectra,2)        #Temporal observations
    feat    = zeros(size(spectra,2)) #Output sequence of features
    M       = size(spectra,1)        #Frequency bins
    weights = linspace(0,100,M)      #Range of signal in Hz
    for i=1:N
        feat[i] = sum(weights.*spectra[:,i][:].^n)/sum(spectra[:,i].^n)
    end
    feat
end

function min_frequency(spectra::Matrix{Float64})
    N       = size(spectra,2)        #Temporal observations
    feat    = zeros(size(spectra,2)) #Output sequence of features
    M       = size(spectra,1)        #Frequency bins
    weights = linspace(0,100,M)      #Range of signal in Hz
    for i=1:N
        feat[i] = weights[indmin(spectra[:,i])]
    end
    feat
end
