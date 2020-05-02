"""
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
"""
# import libraries
using Plots # data visualization
using FFTW # FFT
using DSP # convolution, window, filter
using Statistics # mean
using Distributions # probibility
using LaTeXStrings # image output

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
noise = ts -> randn(ComplexF64, length(ts))
histogram(
	xlabel = L"u/\mathrm{V}",
	ylabel = L"N(u)",
	label = [L"\mathrm{Re}(u_n)" L"\mathrm{Im}(u_n)"],
	[
		real(noise(0:t₀:T)),
		imag(noise(0:t₀:T))
	]
)
savefig("fig/noise_pdf.pdf")
count(x -> abs(x) < 3 / √2, real(noise(0:t₀:T))) /
	length(0:t₀:T)
cdf(Normal(), 3) - cdf(Normal(), -3)
n_lpf(xs) = Int64(round(length(xs) / 4))
lpf_f(xs) = [zeros(n_lpf(xs)); xs[n_lpf(xs) + 1:3n_lpf(xs)];
	zeros(n_lpf(xs))]
lpf(xs) = ifft(ifftshift(lpf_f(fftshift(fft(xs)))))
A = db2amp(dBSNR) * rms(lpf(noise(0:t₀:T)))
@. s₀(t) = A * exp(im * α * (t - T / 2)^2 / 2) * (abs(t % T - T / 2) < τ / 2)
# signal
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"\mathrm{Re}u_{s_0}(t)/\mathrm{V}",
	label = L"\mathrm{Re}u_{s_0}(t)",
	0:t₀:T,
	real(s₀(0:t₀:T))
)
savefig("fig/baseband_wave.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"\mathrm{Re}u_{s_0}(t)/\mathrm{V}",
	label = L"\mathrm{Re}u_{s_0}(t)",
	(T - τ) / 2:t₀:(T + τ) / 2,
	real(s₀((T - τ) / 2:t₀:(T + τ) / 2))
)
savefig("fig/baseband_wave_band.pdf")
plot(
	xlabel = L"f/\mathrm{Hz}",
	ylabel = L"U/(\mathrm{V}\times\mathrm{s})",
	label = [L"U_{\mathrm{Re}s_0}(f)" L"U_{s_0}(f)"],
	range(-f₀, stop = f₀, length = length(0:t₀:T)),
	[
		abs.(fftshift(fft(real(s₀(0:t₀:T))))),
		abs.(fftshift(fft(s₀(0:t₀:T))))
	]
)
savefig("fig/baseband_spectrum.pdf")

# self coherent
@. s⁰(t) = conj(s₀((T + τ) / 2 - t)) # MF response
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u_{s^0}(t)/\mathrm{V}",
	label = L"u_{s^0}(t)",
	0:t₀:T,
	real(s⁰(0:t₀:T))
)
savefig("fig/mf_wave.pdf")
mf(xs, t = t₀) = conv(s⁰(0:t:T), xs)[1:length(xs)]
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"\mathrm{Re}u_{s_0^0}(t)" L"|u_{s_0^0}(t)|"],
	0:t₀:T,
	[
		real(mf(s₀(0:t₀:T))),
		abs.(mf(s₀(0:t₀:T)))
	]
)
savefig("fig/self_coherent_wave.pdf")
feature(xs) = findfirst(round.(abs.(xs)) .> 0):findlast(round.(abs.(xs)) .> 0)
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"\mathrm{Re}u_{s_0^0}(t)" L"|u_{s_0^0}(t)|"],
	(0:t₀:T)[feature(mf(s₀(0:t₀:T)))],
	[
		real(mf(s₀(0:t₀:T)))[feature(mf(s₀(0:t₀:T)))],
		abs.(mf(s₀(0:t₀:T)))[feature(mf(s₀(0:t₀:T)))]
	]
)
savefig("fig/self_coherent_band.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = [L"\mathrm{Re}u_{s_0^0}(t)" L"|u_{s_0^0}(t)|"],
	(0:t₀:T)[feature(mf(s₀(0:t₀:T)))],
	[
		amp2db.(abs.(real(mf(s₀(0:t₀:T)))[feature(mf(s₀(0:t₀:T)))])),
		amp2db.(abs.(mf(s₀(0:t₀:T)))[feature(mf(s₀(0:t₀:T)))])
	]
)
savefig("fig/self_coherent_db.pdf")
feature(xs, n) = feature(xs)[Int64.(round.((length(feature(xs)) + 1) / 2 - n)): Int64.(round.((length(feature(xs)) + 1) / 2 + n))]
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = [L"\mathrm{Re}u_{s_0^0}(t)" L"|u_{s_0^0}(t)|"],
	(0:t₀:T)[feature(mf(s₀(0:t₀:T)),10)],
	[
		amp2db.(abs.(real(mf(s₀(0:t₀:T)))[feature(mf(s₀(0:t₀:T)),10)])),
		amp2db.(abs.(mf(s₀(0:t₀:T)))[feature(mf(s₀(0:t₀:T)),10)])
	]
)
savefig("fig/self_coherent_enlargement_db.pdf")
f⁰ = 20B # change sample frequence
t⁰ = 1 / f⁰
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = [L"\mathrm{Re}u_{s_0^0}(t)" L"|u_{s_0^0}(t)|"],
	(0:t⁰:T)[feature(t⁰,100)],
	[
		amp2db.(abs.(real(mf(s₀(0:t⁰:T), t⁰))[feature(t⁰,100)])),
		amp2db.(abs.(mf(s₀(0:t⁰:T), t⁰))[feature(t⁰,100)])
	]
)
savefig("fig/self_coherent_band_hf_enlargement_db.pdf")
mf_feature(xs, window = rect) = abs.(mf(xs, t⁰)[feature(t⁰)]) .* window(length(feature(t⁰)))
function compare_window(
			process,
			xs,
			ts,
			windows = [rect, cosine, hanning, hamming, lanczos,
				triang, bartlett, blackman],
			)
	curves = []
	for window in windows
		push!(curves, amp2db.(process(xs, window)))
	end
	plot(
		xlabel = L"t/\mathrm{s}",
		ylabel = L"u(t)/\mathrm{dBV}",
		# label = [L"\mathrm{Re}u_{s_0^0}(t)" L"|u_{s_0^0}(t)|"],
		ts,
		curves
	)
	savefig("fig/self_coherent_window_hf_enlargement_db.pdf")
end
compare_window(mf_feature, s₀(0:t⁰:T), (0:t⁰:T)[feature(t⁰)])

# single_mtd
λ = c / f
N = 2 ^ Int64(round(ceil(log2(t_d_min / T)))) # FFT point number
prf = 1 / T
f_d = prf / N
t_d = 1 / f_d
t2r(t) = c * t / 2
f2v(f) = -1 / 2(1 / 2c + 1 / λ * f)
r2t(R) = 2 * R / c
v2f(v) = -2 * v / λ * c / (c - v)
Δv = f2v(-f_d)
ΔR = t2r(1 / B)
v₂ = v₁ + Δv
R₂ = R₁ + ΔR
@. echo_wave(R = 0, v = 0) = t -> s₀(t - r2t(R)) *
	exp( im * 2π * v2f(v) * t)
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"\mathrm{Re}u_n(t)" L"\mathrm{Re}u_{s_1}(t)"],
	0:t₀:T,
	[
		real.(lpf(noise(0:t₀:T))),
		real.(echo_wave(R₁, v₁)(0:t₀:T))
	]
)
savefig("fig/echo_wave_wave.pdf")
plot(
	xlabel = L"f/\mathrm{Hz}",
	ylabel = L"U(f)/(\mathrm{V}\times\mathrm{s})",
	label = [L"|U_n(f)|" L"|U_{s_1}(f)|"],
	range(-f₀, stop = f₀, length = length(0:t₀:T)),
	[
		abs.(lpf_f(fftshift(fft(noise(0:t₀:T))))),
		abs.(fftshift(fft(echo_wave(R₁, v₁)(0:t₀:T))))
	]
)
savefig("fig/echo_wave_spectrum.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"|u_n^0(t)|" L"|u_{s_1}^0(t)|"],
	0:t₀:T,
	[
		abs.(mf(lpf(noise(0:t₀:T)))),
		abs.(mf(echo_wave(R₁, v₁)(0:t₀:T)))
	]
)
savefig("fig/pulse_compression_wave.pdf")
snr(ss, ns) = maximum(abs2.(ss)) / mean(abs2.(ns))
snr(mf(echo_wave(R₁, v₁)(0:t₀:T)), mf(lpf(noise(0:t₀:T))))
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"|u_n^0(t)|" L"|u_{s_1}^0(t)|"],
	0:t₀:2T,
	[
		abs.(mf(lpf(noise(0:t₀:2T)))),
		abs.(mf(echo_wave(R₁, v₁)(0:t₀:2T)))
	]
)
savefig("fig/echo_wave_wave2.pdf")
rearrange(xs) = reshape(xs, Int64(round(length(xs) / N)), N)
rearrange(xs, window, t = t₀) = rearrange(xs)[feature(t), :] .* repeat(window(length(feature(t))), 1, N)
rearrange(echo_wave(R₁, v₁)(t₀:t₀:N * T))[Int64(round(r2t(R₁) / t₀)) .+ feature(t₀), :]
echo_wave(R₁, v₁)(t₀:t₀:N * T)[Int64(round(r2t(R₁) / t₀)) .+ feature(t₀)]
surface(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	zlabel = L"u(t)/\mathrm{V}",
	label = [L"|u_n^0(t)|" L"|u_{s_1}^0(t)|"],
	1:N,
	t₀:t₀:T,
	[
		abs.(rearrange(mf(lpf(noise(t₀:t₀:N * T))))),
		abs.(rearrange(mf(echo_wave(R₁, v₁)(t₀:t₀:N * T))))
	]
)
savefig("fig/mtd_wave.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"|u_n^0(t)|" L"|u_{s_1}^0(t)|"],
	1:N,
	[
		maximum(abs.(rearrange(mf(noise(t₀:t₀:N * T)))), dims = 1)[1, :],
		maximum(abs.(rearrange(mf(echo_wave(R₁, v₁)(t₀:t₀:N * T)))), dims = 1)[1, :]
	]
)
savefig("fig/mtd_wave_slow_view.pdf")
mtd(xs, window) = fftshift(fft(rearrange(mf(xs), window), 2), 2)
surface(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	zlabel = L"u(t)/\mathrm{V}",
	label = [L"|u_n^0(t)|" L"|u_{s_1}^0(t)|"],
	range((f_d - prf) / 2, stop = (prf - f_d) / 2, step = f_d),
	t₀:t₀:T,
	[
		abs.(mtd(lpf(noise(t₀:t₀:N * T)))),
		abs.(mtd(echo_wave(R₁, v₁)(t₀:t₀:N * T)))
	]
)
savefig("fig/mtd_slow_spectrum.pdf")
mtd_f(xs, window) = maximum(abs.(mtd(xs, window)), dims = 1)[1, :]
mtd_t(xs, window) = maximum(abs.(mtd(xs, window)), dims = 2)[:, 1]
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"|u_n^0(t)|" L"|u_{s_1}^0(t)|"],
	range((f_d - prf) / 2, stop = (prf - f_d) / 2, step = f_d),
	[
		mtd_f(lpf(noise(t₀:t₀:N * T))),
		mtd_f(echo_wave(R₁, v₁)(t₀:t₀:N * T))
	]
)
savefig("fig/mtd_slow_spectrum_slow_view.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = [L"|u_n^0(t)|" L"|u_{s_1}^0(t)|"],
	range((f_d - prf) / 2, stop = (prf - f_d) / 2, step = f_d),
	[
		amp2db.(mtd_f(lpf(noise(t₀:t₀:N * T)))),
		amp2db.(mtd_f(echo_wave(R₁, v₁)(t₀:t₀:N * T)))
	]
)
savefig("fig/mtd_slow_spectrum_slow_view_db.pdf")
compare_window(mtd_f, echo_wave(R₁, v₁)(t₀:t₀:N * T), range((f_d - prf) / 2, stop = (prf - f_d) / 2, step = f_d))
compare_window(mtd_t, echo_wave(R₁, v₁)(t₀:t₀:N * T), t₀:t₀:T, [rect, triang])
snr(mtd(echo_wave(R₁, v₁)(t₀:t₀:N * T)), mtd(lpf(noise(t₀:t₀:N * T))))
get_v(xs) = f2v((argmax(mtd_f(xs)) - N / 2) * f_d)
get_r(xs) = t2r(argmax(xs) * t₀ - (T + τ) / 2)
bandwidth_ratio(xs) = sum(xs .> maximum(xs) / √2)
bandwidth_ratio(mtd_f(echo_wave(R₁, v₁)(t₀:t₀:N * T))) * f_d / 2
bandwidth_ratio(abs.(echo_wave(R₁, v₁)(t₀:t₀:N * T))) * t₀

# double_mtd
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{dBV}",
	label = L"s(t)",
)
savefig("fig/amplitude_distinguish_db.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = L"s(t)",
)
savefig("fig/distance_distinguish.pdf")
plot(
	xlabel = L"t/\mathrm{s}",
	ylabel = L"u(t)/\mathrm{V}",
	label = L"s(t)",
)
savefig("fig/velocity_distinguish.pdf")
