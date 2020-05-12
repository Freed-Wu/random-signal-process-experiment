#=
This program is used to simulate LFMCW.
      ,-.
     / \  `.  __..-,O
    :   \ --''_..-'.'
    |    . .-' `. '.
    :     .     .`.'
     \     `.  /  ..
      \      `.   ' .
       `,       `.   \
      ,|,`.        `-.\
     '.||  ``-...__..-`
      |  |
      |__|
      /||\
     //||\\
    // || \\
 __//__||__\\__
'--------------'
=#
# import libraries
using Plots # data visualization
using FFTW # FFT
using DSP # convolution, window, filter
using Statistics # mean
using LaTeXStrings # image output
using CSV # table output
using DataFrames # table create
using Shell # shell inject
Shell.run("mkdir tab") # directory to save tables
Shell.run("mkdir fig") # directory to save figures
Shell.run("mkdir bib") # directory to save reference
Shell.run("cp ~/.julia/packages/FFTW/qqcBj/CITATION.bib bib")

# input
B = 10e6 # band width
τ = 200e-6 # pulse retain time
η = 0.1 # duty ratio
f = 10e9 # carrier frequence
t_d_min = 100e-3 # coherent accumlation time
c = 3e8 # light speed
R₁ = 1000 # distance of object 1
v₁ = 1000 # velocity of object 1
dBSNR = 0 # ratio of noise and signal

# emit
T = τ / η # pulse repeat period
k = B / τ # LFM accelarative frequence
α = 2π * k # LFM accelarative angle
f₀ = 2B # sample frequence
t₀ = 1 / f₀ # sample time
"noise wave"
noise = ts -> randn(ComplexF64, length(ts))
histogram(
	xlabel = L"u/\mathrm{V}",
	ylabel = L"N(u)",
	label = [L"\mathrm{Re}n" L"\mathrm{Im}n"],
	[
		real(noise(0:t₀:T)),
		imag(noise(0:t₀:T))
	]
)
savefig("fig/noise_pdf.pdf")
"get half of point number of low pass"
n_lp(xs) = Int64(round(length(xs) / 4))
"low pass, in fact, DSP exist LowPass function, but I don't know usage"
lp(xs) = [zeros(n_lp(xs)); xs[range(n_lp(xs) + 1, length = length(xs) - 2n_lp(xs))]; zeros(n_lp(xs))]
"fft in fast time"
fft1(xs) = fftshift(fft(xs, 1), 1)
"ifft in fast time"
ifft1(xs) = ifft(ifftshift(xs, 1), 1)
"low pass filter"
lpf(xs) = ifft1(lp(fft1(xs)))
A₀ = db2amp(dBSNR) * rms(lpf(noise(0:t₀:T))) # signal amplitude
"signal wave"
@. s₀(t) = A₀ * exp(im * α * (t - T / 2)^2 / 2) * (abs(t % T - T / 2) < τ / 2)
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"\mathrm{Re}s_0(t)/\mathrm{V}",
	label = [L"\mathrm{Re}s_0(t)" L"|s_0(t)|"],
	0:t₀:T,
	[
		real(s₀(0:t₀:T)),
		abs.(s₀(0:t₀:T))
	]
)
savefig("fig/baseband_wave.pdf")
"get band of signal"
band(xs, a::Float64 = 1 / 1000) = findfirst(abs.(xs) .> a * maximum(abs.(xs))):findlast(abs.(xs) .> a * maximum(abs.(xs)))
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"\mathrm{Re}s_0(t)" L"|s_0(t)|"],
	(0:t₀:T)[band(s₀(0:t₀:T))],
	[
		real(s₀(0:t₀:T)[band(s₀(0:t₀:T))]),
		abs.(s₀(0:t₀:T)[band(s₀(0:t₀:T))])
	]
)
savefig("fig/baseband_wave_band.pdf")
plot(
	xlabel = L"f/\mathrm{Hz}",
	ylabel = L"U/(\mathrm{V}\times\mathrm{s})",
	label = [L"\mathrm{Cs}S_0(f)" L"S_0(f)"],
	range(-f₀, stop = f₀, length = length(0:t₀:T)),
	[
		abs.(fft1(real(s₀(0:t₀:T)))),
		abs.(fft1(s₀(0:t₀:T)))
	]
)
savefig("fig/baseband_spectrum.pdf")

# self coherent
"match filter response"
@. s⁰(t) = conj(s₀((T + τ) / 2 - t))
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"\mathrm{Re}s^0(t)" L"|s^0(t)|"],
	0:t₀:T,
	[
		real(s⁰(0:t₀:T)),
		abs.(s⁰(0:t₀:T)),
	]
)
savefig("fig/mf_wave.pdf")
"match filter"
mf(xs, t = t₀) = conv(s⁰(0:t:T), xs)[1:length(xs)]
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"\mathrm{Re}s_0^0(t)" L"|s_0^0(t)|"],
	0:t₀:T,
	[
		real(mf(s₀(0:t₀:T))),
		abs.(mf(s₀(0:t₀:T)))
	]
)
savefig("fig/self_coherent_wave.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"\mathrm{Re}s_0^0(t)" L"|s_0^0(t)|"],
	(0:t₀:T)[band(mf(s₀(0:t₀:T)), 1 / 1000)],
	[
		real(mf(s₀(0:t₀:T)))[band(mf(s₀(0:t₀:T)), 1 / 1000)],
		abs.(mf(s₀(0:t₀:T)))[band(mf(s₀(0:t₀:T)), 1 / 1000)]
	]
)
savefig("fig/self_coherent_band.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = [L"\mathrm{Re}s_0^0(t)" L"|s_0^0(t)|"],
	(0:t₀:T)[band(mf(s₀(0:t₀:T)), 1 / 1000)],
	[
		amp2db.(abs.(real(mf(s₀(0:t₀:T)))[band(mf(s₀(0:t₀:T)), 1 / 1000)])),
		amp2db.(abs.(mf(s₀(0:t₀:T)))[band(mf(s₀(0:t₀:T)), 1 / 1000)])
	]
)
savefig("fig/self_coherent_band_db.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = [L"\mathrm{Re}s_0^0(t)" L"|s_0^0(t)|"],
	(0:t₀:T)[band(mf(s₀(0:t₀:T)), 1 / 20)],
	[
		amp2db.(abs.(real(mf(s₀(0:t₀:T)))[band(mf(s₀(0:t₀:T)), 1 / 20)])),
		amp2db.(abs.(mf(s₀(0:t₀:T)))[band(mf(s₀(0:t₀:T)), 1 / 20)])
	]
)
savefig("fig/self_coherent_enlargement_db.pdf")
f⁰ = 20B # change sample frequence
t⁰ = 1 / f⁰
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = [L"\mathrm{Re}s_0^0(t)" L"|s_0^0(t)|"],
	(0:t⁰:T)[band(mf(s₀(0:t⁰:T), t⁰), 1 / 20)],
	[
		amp2db.(abs.(real(mf(s₀(0:t⁰:T), t⁰))[band(mf(s₀(0:t⁰:T), t⁰), 1 / 20)])),
		amp2db.(abs.(mf(s₀(0:t⁰:T), t⁰))[band(mf(s₀(0:t⁰:T), t⁰), 1 / 20)])
	]
)
savefig("fig/self_coherent_enlargement_hf_db.pdf")
default_windows = [rect, cosine, hanning, hamming, lanczos, triang, bartlett, blackman]
default_windows_label = ["rect" "cosine" "hanning" "hamming" "lanczos" "triang" "bartlett" "blackman"]
"add window for 1 element function"
function compare_window(
		xs,
		process = x -> x,
		windows = default_windows
		)
	ys = []
	for window in windows
		push!(ys, process(xs .* window(length(xs))))
	end
	ys
end
mf_windows = compare_window(abs.(mf(s₀(0:t⁰:T), t⁰)[band(mf(s₀(0:t⁰:T), t⁰), 1 / 20)]), x -> amp2db.(abs.(x)))
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = default_windows_label,
	mf_windows
)
savefig("fig/self_coherent_window_enlargement_hf_db.pdf")
mf_fft_windows = compare_window(abs.(mf(s₀(0:t⁰:T), t⁰)), x -> amp2db.(abs.(fft1(x))))
plot(
	xlabel = L"f/\mathrm{Hz}",
	ylabel = L"u(t)/\mathrm{dBVs}",
	label = default_windows_label,
	range(-f⁰ / 2, stop = f⁰ / 2, length = length(mf_fft_windows[1])),
	mf_fft_windows
)
savefig("fig/self_coherent_window_fft_hf_db.pdf")
"get index of local minimum"
arglocalmin(xs) = findall(diff(sign.(diff(xs))) .> 0) .+ 1
"get index of main lobe"
function main_lobe(xs)
	ys = arglocalmin(xs)
	zs = argmax(xs)
	posts = ys .>= zs
	pres = ys .<= zs
	ys[findlast(pres)]:ys[findfirst(posts)]
end
main_lobe(xs, L) = findfirst(abs.(xs[main_lobe(xs)]) .> maximum(abs.(xs[main_lobe(xs)])) - L):findlast(abs.(xs[main_lobe(xs)]) .> maximum(abs.(xs[main_lobe(xs)])) - L)
"get width of main lobe"
main_lobe_width(xs) = length(main_lobe(xs))
main_lobe_width(xs, L) = length(main_lobe(xs, L))
"get index of side lobe"
function arg_side_lobe(xs)
	ys = copy(xs)
	ys[main_lobe(ys)] .= minimum(xs)
	argmax(ys)
end
"get height of side lobe"
side_lobe_height(xs) = xs[arg_side_lobe(xs)]
mf_fft_windows_main_lobe = main_lobe_width.(mf_fft_windows) ./ length.(mf_fft_windows) .* f⁰
Shell.run("touch tab/self_coherent_window.csv")
CSV.write(
	"tab/self_coherent_window.csv",
	DataFrame(
		窗函数 = default_windows_label[1,:],
		输出脉冲4dB宽度 = main_lobe_width.(mf_windows, 4) .* t⁰,
		旁瓣高度 = side_lobe_height.(mf_fft_windows),
		主瓣宽度 = mf_fft_windows_main_lobe,
		主瓣宽度展开倍数 = mf_fft_windows_main_lobe ./ mf_fft_windows_main_lobe[1]
	)
)

# single_mtd
λ = c / f
N = 2 ^ Int64(round(ceil(log2(t_d_min / T)))) # FFT point number
prf = 1 / T
f_d = prf / N
t_d = 1 / f_d
t2r(t) = c * t / 2
f2v(f) = -1 / 2(1 / 2c + 1 / (λ * f))
r2t(R) = 2 * R / c
v2f(v) = -2 * v / λ * c / (c - v)
Δv = f2v(-f_d)
ΔR = t2r(1 / B)
v₂ = v₁ + Δv
R₂ = R₁ + ΔR
Shell.run("touch tab/parameter.csv")
CSV.write(
	"tab/parameter.csv",
	DataFrame(
		物理量 = ["不模糊距离", "不模糊速度", "距离分辨率", "速度分辨率"],
		数值 = abs.([t2r(T), f2v(prf), ΔR, Δv])
	)
)
"echo wave about distance and vector"
@. echo_wave(R = 0, v = 0) = t -> s₀(t - r2t(R)) * exp( im * 2π * v2f(v) * t)
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"\mathrm{Re}n(t)" L"\mathrm{Re}s_1(t)"],
	t₀:t₀:T,
	[
		real.(lpf(noise(t₀:t₀:T))),
		real.(echo_wave(R₁, v₁)(t₀:t₀:T))
	]
)
savefig("fig/echo_wave.pdf")
plot(
	xlabel = L"f/\mathrm{Hz}",
	ylabel = L"U(f)/(\mathrm{V}\times\mathrm{s})",
	label = [L"|N(f)|" L"|S_1(f)|"],
	range(-f₀, stop = f₀, length = length(t₀:t₀:T)),
	[
		abs.(fft1(lpf(noise(t₀:t₀:T)))),
		abs.(fft1(echo_wave(R₁, v₁)(t₀:t₀:T)))
	]
)
savefig("fig/echo_spectrum.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"|n^0(t)|" L"|s_1^0(t)|"],
	t₀:t₀:T,
	[
		abs.(mf(lpf(noise(t₀:t₀:T)))),
		abs.(mf(echo_wave(R₁, v₁)(t₀:t₀:T)))
	]
)
savefig("fig/pulse_compression_wave.pdf")
"distance door rearrange"
rearrange(xs) = reshape(xs, Int64(round(length(xs) / N)), N)
surface(
	xlabel = L"n",
	ylabel = L"t/\mathrm{s}",
	zlabel = L"u(t, n)/\mathrm{V}",
	label = [L"|n^0(t, n)|" L"|s_1^0(t, n)|"],
	1:N,
	t₀:t₀:T,
	[
		abs.(rearrange(mf(lpf(noise(t₀:t₀:N * T))))),
		abs.(rearrange(mf(echo_wave(R₁, v₁)(t₀:t₀:N * T))))
	]
)
savefig("fig/mtd_wave.pdf")
"fft in slow time"
fft2(xs) = fftshift(fft(xs, 2), 2)
surface(
	xlabel = L"f/\mathrm{Hz}",
	ylabel = L"t/\mathrm{s}",
	zlabel = L"u(t, n)/\mathrm{V}",
	label = [L"|n^0(t)|" L"|s_1^0(t)|"],
	range((f_d - prf) / 2, stop = (prf - f_d) / 2, step = f_d),
	t₀:t₀:T,
	[
		abs.(fft2(rearrange(mf(lpf(noise(t₀:t₀:N * T)))))),
		abs.(fft2(rearrange(mf(echo_wave(R₁, v₁)(t₀:t₀:N * T)))))
	]
)
savefig("fig/mtd_slow_spectrum.pdf")
"view in slow time"
view1(xs) = maximum(xs, dims = 1)[1, :]
plot(
	xlabel = L"f/\mathrm{Hz}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"|n^0(t)|" L"|s_1^0(t)|"],
	range((f_d - prf) / 2, stop = (prf - f_d) / 2, step = f_d),
	[
		view1(abs.(fft2(rearrange(mf(lpf(noise(t₀:t₀:N * T))))))),
		view1(abs.(fft2(rearrange(mf(echo_wave(R₁, v₁)(t₀:t₀:N * T))))))
	]
)
savefig("fig/mtd_slow_spectrum_slow_view.pdf")
plot(
	xlabel = L"f/\mathrm{Hz}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = [L"|n^0(t)|" L"|s_1^0(t)|"],
	range((f_d - prf) / 2, stop = (prf - f_d) / 2, step = f_d),
	[
		amp2db.(view1(abs.(fft2(rearrange(mf(lpf(noise(t₀:t₀:N * T)))))))),
		amp2db.(view1(abs.(fft2(rearrange(mf(echo_wave(R₁, v₁)(t₀:t₀:N * T)))))))
	]
)
savefig("fig/mtd_slow_spectrum_slow_view_db.pdf")
"distance door rearrange, but row and columne swap"
rearrange2(xs) = transpose(reshape(xs, N, Int64(round(length(xs) / N))))
"add window for 2 element function in dim 2"
function compare_window2(
		xs,
		process = x -> x,
		windows = default_windows
		)
	ys = []
	for window in windows
		push!(ys, process(xs .* rearrange2(repeat(window(size(xs)[2]), size(xs)[1]))))
	end
	ys
end
plot(
	xlabel = L"f/\mathrm{Hz}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = default_windows_label,
	range((f_d - prf) / 2, stop = (prf - f_d) / 2, step = f_d),
	compare_window2(rearrange(mf(echo_wave(R₁, v₁)(t₀:t₀:N * T))), x ->  amp2db.(view1(abs.(fft2(x)))))
)
savefig("fig/mtd_window_slow_spectrum_slow_view_db.pdf")
"the SNR gain"
snr(ss, ns) = maximum(abs2.(ss)) / mean(abs2.(ns))
snr_pulse_compression = snr(mf(echo_wave(R₁, v₁)(0:t₀:T)), mf(lpf(noise(0:t₀:T))))
snr_fft2 = snr(fft2(rearrange(mf(echo_wave(R₁, v₁)(t₀:t₀:N * T)))), fft2(rearrange(mf(lpf(noise(t₀:t₀:N * T)))))) / snr_pulse_compression
"bandwidth point number"
bandwidth_n(xs, a = 1 / 1000) = length(band(xs, a))
mf_bandwidth = bandwidth_n(mf(echo_wave(R₁, v₁)(t₀:t₀:T)), 1 / 10000) * t₀
fft_bandwidth = bandwidth_n(view1(abs.(fft2(rearrange(mf(echo_wave(R₁, v₁)(t₀:t₀:N * T)))))), 1 / √2) * f_d / 2
theorem = [B * τ, N, 2τ, f_d / 2]
measure = [snr_pulse_compression, snr_fft2, mf_bandwidth, fft_bandwidth]
"relative error"
@. err(x, y) = abs(x - y) / y
Shell.run("touch tab/measure.csv")
CSV.write(
	"tab/measure.csv",
	DataFrame(
		物理量 = ["脉冲压缩增益", "傅立叶变换增益", "脉冲压缩时宽", "傅立叶变换带宽"],
		理论值 = theorem,
		测量值 = measure,
		相对误差 = err(theorem, measure)
	)
)
"compare sensitivity of a parameter, xs is the range of parameter"
function compare_sensitivity(f, xs, process, postprocess = Float64)
	ys = []
	for x in xs
		ys = [ys; process(f(x))]
	end
	zs = postprocess.(reshape(ys, Int64(round(length(ys) / length(xs))), Int64(round(length(xs)))))
end
surface(
	xlabel = L"v/(\mathrm{m}/\mathrm{s})",
	ylabel = L"t/\mathrm{s}",
	zlabel = L"u(t, v)/\mathrm{dBV}",
	label = ["s(t, v)"],
	0:100:1000,
	(t₀:t₀:T)[band(mf(echo_wave(R₁, v₁)(t₀:t₀:T)), 1 / 1000)],
	compare_sensitivity(v -> mf(echo_wave(R₁, v)(t₀:t₀:T))[band(mf(echo_wave(R₁, v₁)(t₀:t₀:T)), 1 / 1000)], 0:100:1000, x -> amp2db.(abs.(x)))
)
savefig("fig/pulse_compression_wave_velocity_db.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t, v)/\mathrm{dBV}",
	label = hcat(0:100:1000 ...),
	(t₀:t₀:T)[band(mf(echo_wave(R₁, v₁)(t₀:t₀:T)), 1 / 1000)],
	compare_sensitivity(v -> mf(echo_wave(R₁, v)(t₀:t₀:T))[band(mf(echo_wave(R₁, v₁)(t₀:t₀:T)), 1 / 1000)], 0:100:1000, x -> amp2db.(abs.(x)))
)
savefig("fig/pulse_compression_wave_velocity_section_view_time_db.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t, v)/\mathrm{dBV}",
	label = hcat(0:100:1000 ...),
	(t⁰:t⁰:T)[band(mf(echo_wave(R₁, v₁)(t⁰:t⁰:T), t⁰), 1 / 20)],
	compare_sensitivity(v -> mf(echo_wave(R₁, v)(t⁰:t⁰:T), t⁰)[band(mf(echo_wave(R₁, v₁)(t⁰:t⁰:T), t⁰), 1 / 20)], 0:100:1000, x -> amp2db.(abs.(x)))
)
savefig("fig/pulse_compression_wave_velocity_section_view_time_enlargement_db.pdf")

# double_mtd
"signals and their sum"
add(xss) = [xss; [sum(xss)]]
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = [L"|s_{11}^0(t)|" L"|s_1^0(t)|" L"s_{11}^0(t) + s_1^0(t)"],
	(t⁰:t⁰:T)[band(mf(echo_wave(R₁, v₁)(t⁰:t⁰:T), t⁰), 1 / 20)],
	add(
		[
			amp2db.(abs.(mf(1000 * echo_wave(1.1R₁, 1.1v₁)(t⁰:t⁰:T), t⁰)[band(mf(echo_wave(R₁, v₁)(t⁰:t⁰:T), t⁰), 1 / 20)])),
			amp2db.(abs.(mf(echo_wave(R₁, v₁)(t⁰:t⁰:T), t⁰)[band(mf(echo_wave(R₁, v₁)(t⁰:t⁰:T), t⁰), 1 / 20)]))
		]
	)
)
savefig("fig/amplitude_distance_distinguish_db.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = [L"|s_{11}^0(t)|" L"|s_1^0(t)|" L"s_{11}^0(t) + s_1^0(t)"],
	add(
		[
			amp2db.(view1(abs.(fft2(rearrange(mf(1000 * echo_wave(1.1R₁, 1.1v₁)(t₀:t₀:N * T))))))),
			amp2db.(view1(abs.(fft2(rearrange(mf(echo_wave(R₁, v₁)(t₀:t₀:N * T)))))))
		]
	)
)
savefig("fig/amplitude_velocity_distinguish_db.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = [L"|s_2^0(t)|" L"|s_1^0(t)|" L"s_2^0(t) + s_1^0(t)"],
	(t₀:t₀:T)[band(mf(echo_wave(R₁, v₁)(t₀:t₀:T)), 1 / 1000)],
	add(
		[
			amp2db.(abs.(mf(echo_wave(R₂, v₂)(t₀:t₀:T))[band(mf(echo_wave(R₁, v₁)(t₀:t₀:T)), 1 / 1000)])),
			amp2db.(abs.(mf(echo_wave(R₁, v₁)(t₀:t₀:T))[band(mf(echo_wave(R₁, v₁)(t₀:t₀:T)), 1 / 1000)]))
		]
	)
)
savefig("fig/distance_distinguish_db.pdf")
plot(
	xlabel = L"f/\mathrm{Hz}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = [L"|s_2^0(t)|" L"|s_1^0(t)|" L"s_2^0(t) + s_1^0(t)"],
	range((f_d - prf) / 2, stop = (prf - f_d) / 2, step = f_d),
	add(
		[
			amp2db.(view1(abs.(fft2(rearrange(mf(echo_wave(R₂, v₂)(t₀:t₀:N * T))))))),
			amp2db.(view1(abs.(fft2(rearrange(mf(echo_wave(R₁, v₁)(t₀:t₀:N * T)))))))
		]
	)
)
savefig("fig/velocity_distinguish.pdf")

