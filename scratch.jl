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
