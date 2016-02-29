function [caBandFiltObjs, mBands] = computeComplementaryFilterbank(vEdgeFreqs, iOrder, fs, vX)
% Function to implement a power complementary filter bank.
% The filter bank is based on the paper: 
% "Complementary N-Band IIR Filterbank Based on 2-Band Complementary Filters"
% written by Alexis Favrot and Christof Faller (2010)
%
% ----------------------------------------------------------------------
% Usage: [caBandFiltObjs, mBands] = computeComplementaryFilterbank(vEdgeFreqs, iOrder, fs, vX)
%
%   input:   ---------
%     vEdgeFreqs       vector containing the edge frequencies of the filter bank
%     iOrder           order of the butterworth filter used for the filter bank
%     fs               sampling frequency
%     vX               signal to be filtered
%
%  output:   --------- 
%     caBandFiltObjs   cell array containing the dfilt objects for each band
%     mBands           filtered signal matrix (:,bands)
%
% Author :  (c) Nils L. Westhausen (TGM @ Jade-Hochschule)
%           
% Date   :  26 Dec 2015
% Updated:  22 Feb 2016
%-----------------------------------------------------------------------

    % check for enough input arguments 
    if (nargin < 3)
        error('not enough input arguments')
    end
    % checking if order is odd
    if rem(iOrder, 2) == 0
        error(' Order must be odd')
    end
    % check edge frequencies
    if any(vEdgeFreqs > fs/2)
        error('no frequencies greater than fs/2 allowed')
    end
    if any(vEdgeFreqs < 1) 
        error('frequencies have to be greater than zero')
    end
    % sort frequencies ascending
    vEdgeFreqs = sort(vEdgeFreqs);
    % set recursive index
    iRecursivIdx = 1;
    % normalization of frequencies
    vEdgeFreqs = vEdgeFreqs ./ fs;
    % calculating number of bands
    iNumBands = length(vEdgeFreqs) + 1;
    % preallocation of the cell for the band filter objects
    caBandFiltObjs = cell(1, iNumBands);
    % generation of filter objects
    [caFiltObjs, caCompFiltObjs, caAllPassObjs] = GenFiltObjs(vEdgeFreqs, iOrder);

    caFiltObj = {};
    % recursive generation of the filter objects of each band
    [caBandFiltObjs, ~] = GenBandFiltObjs(caFiltObj, caBandFiltObjs, ...
        iRecursivIdx, caFiltObjs, caCompFiltObjs, caAllPassObjs);
    % generating filtered output if called
    if nargout > 1
        mBands = filter_bands(vX, caBandFiltObjs);
    end
end

% function to generate filter objects of each band recursively 
function [caBandFiltObjs, RecursivIdx] = GenBandFiltObjs(caFiltObj, caBandFiltObjs, RecursivIdx, caFiltObjs, caCompFiltObjs, caAllPassObjs)
    % checking for number of given filter objects
    if length(caFiltObjs) > 2
        % separating filtObjs
        FiltIdx = 2^nextpow2(length(caFiltObjs)) / 2;
        FiltIdxsLeft = 1:FiltIdx - 1;
        FiltIdxsRigth = FiltIdx + 1:length(caFiltObjs);
        FiltObj1 = [caFiltObj, caFiltObjs(FiltIdx), ...
            caAllPassObjs(FiltIdxsRigth)];
        FiltObj2 = [caFiltObj, caCompFiltObjs(FiltIdx), ...
            caAllPassObjs(FiltIdxsLeft)];
        % recursively calling the function with the serrated filtObjs
        [caBandFiltObjs, RecursivIdx] = GenBandFiltObjs(FiltObj1, ...
            caBandFiltObjs, RecursivIdx, caFiltObjs(FiltIdxsLeft), ...
            caCompFiltObjs(FiltIdxsLeft), caAllPassObjs(FiltIdxsLeft));
        [caBandFiltObjs, RecursivIdx] = GenBandFiltObjs(FiltObj2, ...
            caBandFiltObjs, RecursivIdx, caFiltObjs(FiltIdxsRigth), ...
            caCompFiltObjs(FiltIdxsRigth), caAllPassObjs(FiltIdxsRigth));
    elseif length(caFiltObjs) == 1
         % saving filtObjs to the cell
        caBandFiltObjs{RecursivIdx} = dfilt.cascade(caFiltObj{1:end}, caFiltObjs{1});
        caBandFiltObjs{RecursivIdx+1} = dfilt.cascade(caFiltObj{1:end}, ...
            caCompFiltObjs{1});
        RecursivIdx = RecursivIdx + 2;
    elseif length(caFiltObjs) == 2
        y1 = [caFiltObj, caFiltObjs(2)];
        caBandFiltObjs{RecursivIdx + 2} = dfilt.cascade(caFiltObj{1:end}, ...
            caCompFiltObjs{2}, caAllPassObjs{1});
        [caBandFiltObjs, RecursivIdx] = GenBandFiltObjs(y1, caBandFiltObjs, ...
            RecursivIdx, caFiltObjs(1), caCompFiltObjs(1), caAllPassObjs(1));
        RecursivIdx = RecursivIdx + 1;
    else
        error('something unforeseen happened')
    end
end

% function to generate filter objects (lowpass, complementary filter, allpass)
function [caFiltObjs, caCompFiltObjs, caAllPassObjs] = GenFiltObjs(vEdgeFreqs, iOrder)
    LenEdgeFreq = length(vEdgeFreqs);
    caFiltObjs = cell(1, LenEdgeFreq);
    caCompFiltObjs = cell(1, LenEdgeFreq);
    caAllPassObjs = cell(1, LenEdgeFreq);
    for ff = 1:LenEdgeFreq
        [vB, vA] = butter(iOrder ,vEdgeFreqs(ff) * 2);
        [vBp,vAp] = iirpowcomp(vB,vA);
        sos = tf2sos(vB,vA);
        sos_p = tf2sos(vBp,vAp);
        caFiltObjs{ff} = dfilt.df1sos(sos);
        caCompFiltObjs{ff} = dfilt.df1sos(sos_p);
        caAllPassObjs{ff} = dfilt.parallel(caFiltObjs{ff}, ...
            caCompFiltObjs{ff});
    end
end

% function to filter the input signal with each filter object
function mY = filter_bands(vX, caBandFiltObjs)
    iNumBands = length(caBandFiltObjs);
    mY = zeros(length(vX), iNumBands);
    for bb = 1:iNumBands
        mY(:, bb) = filter(caBandFiltObjs{bb}, vX);
    end
end
% end of file
