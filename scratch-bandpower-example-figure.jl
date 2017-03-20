using TotalVariation

subject = 18
f = make_features(import_data_band(subject))
l = import_labels(subject)
println(size(f))
println(size(l))
println(size(f,2)*30/60/60, "hours")
println(length(l)*5/60/60, "hours")
#return
f = f[:,1:900]
l = l[1:6:end][1:900]

figure(1)
PyPlot.clf();
samples = size(f,2)
tmp = 2*[0.12, 0.04, 0.04, 0.04, 0.04, 0.04]
for i=1:6#[3,4,6]
    plot(gstv(log(f[i,:])[:], 50, tmp[i]))
    #plot(log(f[i,:])[:])
end

plot(linspace(1,size(f,2),length(l)), l+11.9)
legend(["delta", "theta", "alpha", "low beta", "high beta",
"gamma", "sleep stage"],
loc="lower right")
savefig("band-power-gstv.png")

figure(2)
PyPlot.clf();
samples = size(f,2)
tmp = 2*[0.12, 0.04, 0.04, 0.04, 0.04, 0.04]
for i=1:6#[3,4,6]
    #plot(gstv(log(f[i,:])[:], 50, tmp[i]))
    plot(log(f[i,:])[:])
end

plot(linspace(1,size(f,2),length(l)), l+11.9)
legend(["delta", "theta", "alpha", "low beta", "high beta",
"gamma", "sleep stage"],
loc="lower right")
savefig("band-power-raw.png")
