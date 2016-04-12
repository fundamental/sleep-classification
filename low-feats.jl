#Features used in Philip Low's Dissertation...

low_cf(i,w,s) = sum(w.*s[:,i][:].^2.0)/sum(s[:,i].^2.0)
low_hf(i,w,s) = w[indmax(s[:,i])]
low_lf(i,w,s) = w[indmin(s[:,i])]


"""
Produce 1D Feature representation for Low style features

    These features don't generalize terribly well, but they have helped motivate
    the work in the rest of this repository
"""
function low_feat(spectra::Matrix{Float64}; func=low_cf)
    N       = size(spectra,2)        #Temporal observations
    feat    = zeros(size(spectra,2)) #Output sequence of features
    M       = size(spectra,1)        #Frequency bins
    weights = linspace(0,100,M)      #Range of signal in Hz
    for i=1:N
        feat[i] = func(i,weights,spectra)
    end
    feat
end
