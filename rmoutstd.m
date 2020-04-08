function [clean_series, ids] = rmoutstd(series, std_dev)
%RMOUTSTD(SERIES, STD_DEV) Gets the variances and covariances of the idiosyncratic components of
%pair returns and the variance of the market.
% USAGE:
%   [CLEAN_SERIES, IDS] = pair_vols(SERIES,STD_DEV)
%
% INPUTS:
%   SERIES       - Series we want to remove outliers
%   STD_DEV      - Number of standard deviations from the mean after which we will remove
%   outlierts (default = 3)
%
% OUTPUTS:
%   CLEAN_SERIES - Series without outliers
%   IDS          - Id for selection of series
%
% Copyright: Tales Padilha
% tales.padilha@economics.ox.ac.uk
% Date: 26/3/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 1
        std_dev = 3;
    case 2
    otherwise
        error('Number of inputs must be 1 or 2');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing Outliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
series(series>(median(series)+(std(series)*std_dev))|...
    series<(median(series)-(std(series)*std_dev)))=NaN;
ids = ~isnan(series);
clean_series = series(ids);
end