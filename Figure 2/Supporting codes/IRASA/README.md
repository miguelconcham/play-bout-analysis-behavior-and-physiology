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
10 Hz peak, an 8 Hz plus 16 Hz harmonic, and that the play-session
output keeps one spectrum row per window, without the experiment files.

## Play-session pipeline

The recording loop matches the mid-LPAG play-session spectra. It uses
`Session_structure` / `Partners session` rows only. It does not detect
play bouts.

Run `add_repo_paths` once per MATLAB session, then:

1. `Estimate_play_session_lpag_irasa` loops Neuropixels folders and
   calls `GENERATE_PLAY_SESSION_LPAG_IRASA` for each recording.
2. Each Partners-session row is split into windows. The default length
   is 10 s with a 10 s step, so the raw spectrum has 0.1 Hz bins and the
   windows do not overlap. Change `window_seconds` at the top of
   `Estimate_play_session_lpag_irasa` to use another length. Set
   `window_step_seconds` when the step should differ from the window.
   Each window stores its raw Welch and Fourier power, its fractal
   spectrum, and its periodic spectrum. Those matrices keep one row per
   window. The estimate step does not average them. The saved file is
   `Data/Analysis results/psth power by frequency and behavior/play_session_lpag_irasa.mat`
   on `\\experimentfs.home.bccn-berlin.de`.
3. `Analyse_play_session_lpag_irasa` loads that file, averages the
   windows in linear power, prints one row per Partners session, and
   draws the raw, fractal, and periodic spectra. It saves each figure
   beside the `.mat` file and leaves the figures open.
   `close_figures_after_saving` at the top of the script closes them
   after saving when set to true.

`irasa_separate` remains the single-epoch function used inside each
window. Section 2.5 of the paper simulated one 8.5 s epoch (8,500
samples at 1,000 Hz). Section 2.6 used 3 s windows stepped by 1 s on
the ECoG and MEG recordings. The 10 s default here is the longer epoch
that includes 0.1 Hz.

The authors' released toolbox (`amri_sig_fractal` 0.10,
doi:10.4231/R7WQ01R7) is a later code revision of the same method. It
uses 15 power-of-two windows, MATLAB `resample`, and it averages the
resampled spectra before the geometric mean. This folder follows the
order and the cubic-spline resampling written in section 2.4.

Run `add_repo_paths` once per MATLAB session. The IRASA folder is
already included by that script.
