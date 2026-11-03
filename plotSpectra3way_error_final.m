function [f, spectrum, CI]= plotSpectra3way_error_final(data, sampleInterval, varargin)

% Function to calculate a power spectrum and its 95% confidence interval.
% Each column of data represents one segment.
%
% Inputs:
% 1) data: data matrix, with each column representing one segment
% 2) sampleInterval: sampling interval in minutes
%
% Optional inputs:
% 'Detrend': remove a linear trend (default: false)
% 'Demean': remove the mean (default: false)
% 'Hanning': apply a Hanning window (default: false)
% 'Overlap': fraction of overlap between segments (default: 0.50)
%
% Outputs:
% 1) f: frequency in cycles per day
% 2) spectrum: averaged one-sided power spectrum
% 3) CI: lower and upper 95% confidence intervals

% Read optional inputs (Had to use AI, didn't know you could do optional
% inputs)
p= inputParser;

addParameter(p,'Detrend',false,@(x) islogical(x) && isscalar(x));
addParameter(p,'Demean',false,@(x) islogical(x) && isscalar(x));
addParameter(p,'Hanning',false,@(x) islogical(x) && isscalar(x));
addParameter(p,'Overlap',0.50,@(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);

parse(p,varargin{:});

doDetrend= p.Results.Detrend;
doDemean= p.Results.Demean;
doHanning= p.Results.Hanning;
overlap= p.Results.Overlap;

% Check inputs
if ~isnumeric(data) || isempty(data) || any(~isfinite(data(:)))
    error('Data must be a nonempty numeric matrix without NaNs or Inf values.')
end

if ~isscalar(sampleInterval) || ~isfinite(sampleInterval) || sampleInterval <= 0
    error('Sample interval must be a positive number in minutes.')
end

% Get number of points and segments
[N,nSegments]= size(data);

% Sampling frequency in samples per day
Fs= 24*60/sampleInterval;

% Frequency vector in cycles per day
f= (0:floor(N/2))' * Fs/N;

% Create Hanning window using the cosine-squared equation
tw= linspace(-0.5,0.5,N)';
hanwin= cos(pi*tw).^2;

% Choose window
if doHanning
    window= hanwin;
else
    window= ones(N,1);
end

% Preallocate spectra
segmentSpectra= zeros(length(f),nSegments);

% Calculate spectrum for each segment
for k= 1:nSegments

    x= data(:,k);

    % Remove linear trend if requested
    if doDetrend
        x= detrend(x);

    % Otherwise, remove mean if requested
    elseif doDemean
        x= x - mean(x);
    end

    % Apply Hanning window if requested
    x= x .* window;

    % Calculate Fourier transform
    X= fft(x);

    % Calculate power spectrum
    P= abs(X).^2/(Fs*sum(window.^2));

    % Keep positive frequencies
    P= P(1:floor(N/2)+1);

    % Convert to one-sided spectrum
    if rem(N,2) == 0
        P(2:end-1)= 2*P(2:end-1);
    else
        P(2:end)= 2*P(2:end);
    end

    % Store segment spectrum
    segmentSpectra(:,k)= P;

end

% Average spectra across all segments
spectrum= mean(segmentSpectra,2);

% Calculate effective degrees of freedom
% This accounts approximately for correlation between overlapping segments.
step= round(N*(1-overlap));
windowEnergy= sum(window.^2);
correlationSum= 0;

for m= 1:nSegments-1

    shift= m*step;

    if shift >= N
        break
    end

    overlapCorrelation= sum(window(1:N-shift).*window(1+shift:N)) / windowEnergy;

    correlationSum= correlationSum + ...
        (1-m/nSegments)*overlapCorrelation^2;

end

nu= 2*nSegments/(1+2*correlationSum);

% Calculate 95% confidence intervals
alpha= 0.05;

CI.lower= nu*spectrum/chi2inv(1-alpha/2,nu);
CI.upper= nu*spectrum/chi2inv(alpha/2,nu);
CI.dof= nu;

end
