function [var_i, var_j, cov, var_m,param] = pair_vols(pair_r, market_r, const)
%PAIR_VOLS(P_I, P_J, DATES) Gets the variances and covariances of the idiosyncratic components of
%pair returns and the variance of the market.
% USAGE:
%   [VAR_I, VAR_J, COV, VAR_M] = pair_vols(PAIR_R,MARKET_R,DATES)
%
% INPUTS:
%   PAIR_R       - TX2 matrix with returns of the stocks in the pair
%   MARKET_R     - TX1 vector with the returns of the market
%   DATES        - TX1 vector with the dates
%   Constant     
%
% OUTPUTS:
%   VAR_I        - Estimated variances for stock i
%   VAR_J        - Estimated variances for stock j
%   COV          - Estimated covariances beween i and j
%   VAR_M        - Estimated variances for the market
%
% Copyright: Tales Padilha
% tales.padilha@economics.ox.ac.uk
% Date: 22/3/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 2
        const = 0;
    case 3
    otherwise
        error('Number of inputs must be 2 or 3');
end
% Adjusting possible NaN in 'pair_r':
pair_r(isnan(pair_r)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Isolating Idiosyncratic Returns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orthogonalizing with constant:
if const == 1;
    % Stock i:
    idios_i = pair_r(:,1) - ([ones(numel(market_r),1),market_r]*...
        ([ones(numel(market_r),1),market_r]\pair_r(:,1)));
    % Stock j:
    idios_j = pair_r(:,2) - ([ones(numel(market_r),1),market_r]*...
        ([ones(numel(market_r),1),market_r]\pair_r(:,2)));
% Orthogonalizing without a constant:
elseif const == 0;
    % Stock i:
    idios_i = pair_r(:,1) - ((market_r\pair_r(:,1))*market_r);
    % Stock j:
    idios_j = pair_r(:,2) - ((market_r\pair_r(:,2))*market_r);
else 
     error('CONST parameter must be 0 or 1');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimating DCC GARCH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.Display = 'none';
options.Diagnostics = 'off';
[~, ~, Ht] = dcc([idios_i, idios_j], [], 1, 1, 1,[],[],[],[],[],[],[],options);
var_i = squeeze(Ht(1,1,:));
var_j = squeeze(Ht(2,2,:));
cov = squeeze(Ht(1,2,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimating GJR-GARCH(1,1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options2 = optimset('fminunc');
options2.Display = 'none';
options2.Diagnostics = 'off';
[~, ~, var_m] = tarch(market_r, 1, 1, 1, [], [], [], options2);
end