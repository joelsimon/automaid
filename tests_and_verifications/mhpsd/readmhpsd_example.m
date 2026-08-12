function readmhpsd_example(filename)
% READMHPSD_EXAMPLE(filename)
%
% Read and plot an .mhpsd file in MATLAB.
%
% Relies on functions in $OMNIA, in particular `readmhpsd` --
% https://github.com/joelsimon/omnia.git
%
% Author: Joel D. Simon <jdsimon@bathymetrix.com>
% Last modified: 29-Jul-2026, 9.13.0.2553342 (R2022b) Update 9 on MACI64 (geo_mac)
% (in reality: Intel MATLAB in Rosetta 2 running on an Apple silicon Mac)

clc
close all

defval('filename', '20211116T125142.0002_6194A40E.MER.STD.mhpsd')

% Read the MERMAID Hydrophone Power-Spectral-Density file.
[hdr, psd] = readmhpsd(filename);

% Display the header info (similar to a SAC header).
disp(hdr)

% Recreate the 50-95% percentile figure from automaid v3.6.0-X
figure
semilogx(psd.freq, psd.perc50, 'b', 'LineWidth', 2)
hold on
semilogx(psd.freq, psd.perc95, 'r', 'LineWidth', 2)
hold off
legend('50th percentile', '95th percentile')
title(sprintf('Network: %s Station: %s Start Time: %s', hdr.Network, hdr.Station, ...
              hdr.StartTime))

% Note that psd.freq(1) = 0, log(0) is undef, so the first point MATLAB plots is
% psd.freq(2) (and Python draws a straight line headed towards 0 to the left and
% bounds the left XLim seemingly arbitrarily...)
xlim([psd.freq(1) psd.freq(end)])
ylim([-110 -40])

% Add cosmetic finishes.
shrink([], 1, 1.5)
longticks([], 2)
xlabel('Freq. [Hz]')
ylabel('dBfs$^2$/Hz')
latimes
%savepdf(filename)
