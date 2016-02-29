% Script to test the power complementary filter bank implemented in the
% function computeComplementaryFilterbank.m.
% The filter bank is based on the paper: 
% "Complementary N-Band IIR Filterbank Based on 2-Band Complementary Filters"
% written by Alexis Favrot and Christof Faller (2010)
%
% ----------------------------------------------------------------------
%
% Author :  (c) Nils L. Westhausen (TGM @ Jade-Hochschule)
%			
% Date   :  12 Feb 2015
% Updated:  <>
%-----------------------------------------------------------------------

% set parameters
fs = 8000;
% fs = 48000;
nfft = 4096;
iOrder = 5;
vX = [zeros(round(0.5 * nfft), 1) ; 1 ; zeros(round(0.5 * nfft) - 1, 1)];
vF = linspace(1, fs / 2, nfft / 2 + 1 );
vEdgeFreq = [1415 500 2460 250 800 1905 3260 125 370 615 1005];
% vEdgeFreq = [125 250 500 1000 2000 4000 8000 16000];

% compute filter bank
[caFiltObjs, mY] = computeComplementaryFilterbank(vEdgeFreq, iOrder, fs, vX);

% plotting the transfer function of each filter and the sum
figure; 
hold on
for idx = 1:size(mY, 2)
    vTf = fft(mY(:, idx)) ./ fft(vX);
    vAbsTf = 20 * log10(abs(vTf(1:end / 2 + 1)));
    plot(vF, vAbsTf)
end
vYsum = sum(mY, 2);
vTfSum = fft(vYsum) ./ fft(vX);
vAbsTfSum = 20 * log10(abs(vTfSum(1:end / 2 + 1)));
plot(vF, vAbsTfSum, 'k')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB re. FS]')
ylim([-60 10])
grid on
hold off

% set all bands to -12 dB except of the edge bands and plot again
figure; 
hold on
mY2 = mY;
mY2(:, 2:end - 1) = mY2(:, 2:end - 1) ./ 4;
for idx = 1:size(mY2, 2)
    vTf = fft(mY2(:, idx)) ./ fft(vX);
    vAbsTf = 20 * log10(abs(vTf(1:end / 2 + 1)));
    plot(vF, vAbsTf)
end
vYsum = sum(mY2, 2);
vTfSum = fft(vYsum) ./ fft(vX);
vAbsTfSum = 20 * log10(abs(vTfSum(1:end / 2 + 1)));
plot(vF, vAbsTfSum, 'k')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB re. FS]')
ylim([-60 10])
title('-12 dB re. FS except of the edge bands')
grid on
hold off
% end of file
