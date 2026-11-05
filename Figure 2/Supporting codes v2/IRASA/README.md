# IRASA

Irregular-resampling auto-spectral analysis for one LFP epoch.

Wen H, Liu Z (2016) Separating fractal and oscillatory components in the
power spectrum of neurophysiological signal. *Brain Topography* 29:13–26.
doi:10.1007/s10548-015-0448-0. The steps below are section 2.4 of that
paper.

`irasa_separate` takes a time-by-channel matrix, or one row or column,
and a sampling rate in Hz.

1. Cut the epoch into 10 windows, each 90% of the epoch, evenly spaced.
2. Hanning-tapered FFT. The FFT length is twice the smallest power of
   two greater than the window.
3. Resample each window by `h` and by `1/h`, with `h` from 1.1 to 1.9
   in steps of 0.05. Upsampling is a cubic spline. Downsampling is an
   anti-alias low-pass and then a cubic spline. Both spectra use the
   same FFT length, so they share one frequency grid. On that grid a
   peak at `f` moves to `f/h` and to `f*h`.
4. Geometric mean of each pair, then the median across `h`. That median
   is the fractal spectrum of the window (equations 6 and 8).
5. Average the fractal spectra across windows. Average the original
   spectra across the same windows.
6. Subtract. The difference is the oscillatory spectrum (equation 9).
   Negative values are kept.

Frequencies above `sampling_rate / (2 * max(h))` are omitted. The
discussion section "Choosing the resampling factors" is why: a larger
`h` reads the upsampled branch past the original Nyquist.

`irasa_power_law_fit` resamples the fractal spectrum onto equal
log-frequency spacing and fits the line described after section 2.4.
`exponent` is the positive beta in `P(f) ∝ f^(-beta)`.

`irasa_simulate_fractal` builds the inverse-FFT power-law series used
in section 2.5. `test_irasa` checks the defaults, a pure power law, a
10 Hz peak, and an 8 Hz plus 16 Hz harmonic, without the experiment
files.

The authors' released toolbox (`amri_sig_fractal` 0.10,
doi:10.4231/R7WQ01R7) is a later code revision of the same method. It
uses 15 power-of-two windows, MATLAB `resample`, and it averages the
resampled spectra before the geometric mean. This folder follows the
order and the cubic-spline resampling written in section 2.4.

Run `add_repo_paths` once per MATLAB session. The IRASA folder is
already included by that script.
