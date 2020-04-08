%% Market Return and Pairs Trading Profitability
% This section of the empirical study will focus on analysing the relationship between the returns
% of pairs trading and the return of the "Market" (CRSP Value-Weighted Market Return Index).

clear all
clc

%% Selecting which portfolio we will analyse:

portfolio = "DM_5";
plot_title = "5 Pairs";

%% Loading and Formatting the Data

fprintf('Reading data... ');
% Loading portfolio returs:
filename = strcat('Daily_r_',portfolio);
load(filename)   
% Loading market returns:
load Market
% Formatting data to keep all with the same sample:
[x, dates_id] = intersect(dates, daily_returns(:,1));
dates = dates(dates_id);
r_m = r_m(dates_id);
r_p = daily_returns(:,2)*100;
clear daily_returns cum_returns dates_id x dates_id filename
fprintf('Done!\n');

%% Estimating the VAR 
% We begin by estimating the VAR(10) model from Y := [R_m R_pairs]'
% NOTE: This section uses the function "varm" for VAR estimation (it uses MLE). This is available in
% the MATLAB Econometrics Toolbox.
Y = [r_m, r_p];

% Getting the matrix of innovations (e) for the model with the chosen specification:
Mdl = varm(2, 10);
Mdl.Trend = NaN;
[EstMdl, ~, ~, e] = estimate(Mdl,Y);
summarize(EstMdl)

clear EstMdl i Mdl  
mean(e)
std(e)

% Normalizing residuals:
std_e = (e - mean(e))./std(e);

%% BEKK-GARCH Estimation
% In this section, we will estimate an asymmetric BEKK GARCH (1, 1, 1) model to the fitted residuals
% from the VAR (estimated in the previous sections).
% NOTE: This section uses the function "bekk" to estimate the BEKK GARCH model. This function is
% available in the MFE Toolbox (https://www.kevinsheppard.com/MFE_Toolbox).

p = 1;
o = 1;
q = 1;

[parameters, ll, Ht, VCV, scores]  = bekk(std_e, [], p, o, q, 'Full');
[C,A,G,B] = bekk_parameter_transform(parameters,p,o,q,size(Ht,1), 3)
% OBS: Matrices of coefficients are transposed in 'bekk' function. 

%% Analysing Variances and Correlation

% Creating the vectors with variances, covariance and correlation:
sigma2_m = squeeze(Ht(1,1,:));
sigma2_p = squeeze(Ht(2,2,:));
covar = squeeze(Ht(1,2,:));
rho = covar./(sqrt(sigma2_m.*sigma2_p));

% Setting the date for the x axis:
t = dates(11:end, 1);
d = num2str(t);
dates_plot = datetime(d,'InputFormat','yyyyMMdd','Format','dd/MM/yyyy');
clear t

figure()
subplot(2,1,1);
y = sqrt(sigma2_m);
x =sqrt(sigma2_p);
d = dates_plot;
plot(d,x)
yyaxis right
plot(d,y)
title(strcat('Sigmas - ', plot_title))
subplot(2,1,2);
z =rho;
plot(d,z)
yyaxis right
plot(d,y)
title(strcat('Rho and Market Sigma - ',plot_title))

%% News Impact Surfaces:
% Creating the news shocks:
[e_m,e_p] = meshgrid(-10:0.1:10,-10:0.1:10);
unconditional = cov(e);

% Market Variance Surface:
Z1 = C(1,1) + (A(1,1)^2).*(e_m.^2) + (2*A(1,1)*A(2,1)).*e_m.*e_p + (A(2,1)^2).*(e_p.^2)...
    + (G(1,1)^2).*(e_m.^2).*(e_m<0) + (2*G(1,1)*G(2,1)).*e_m.*e_p.*(e_m<0).*(e_p<0) + (G(2,1)^2).*(e_p.^2).*(e_p<0)...
    + (B(1,1)^2)*unconditional(1,1) + (2*B(1,1)*B(2,1))*unconditional(2,1) + (B(2,1)^2)*unconditional(2,2);

% Pairs Variance Surface:
Z2 = C(2,2) + (A(1,2)^2).*(e_m.^2) + (2*A(1,2)*A(2,2)).*e_m.*e_p + (A(2,2)^2).*(e_p.^2)...
    + (G(1,2)^2).*(e_m.^2).*(e_m<0) + (2*G(1,2)*G(2,2)).*e_m.*e_p.*(e_m<0).*(e_p<0) + (G(2,2)^2).*(e_p.^2).*(e_p<0)...
    + (B(1,2)^2)*unconditional(1,1) + (2*B(1,2)*B(2,2))*unconditional(2,1) + (B(2,2)^2)*unconditional(2,2);

%Covariance :
Z3 = C(2,1) + (A(1,1)*A(1,2)).*(e_m.^2) + ((A(1,1)*A(2,2))+(A(1,2)*A(2,1))).*e_m.*e_p + (A(2,1)*A(2,2)).*(e_p.^2)...
    + (G(1,1)*G(1,2)).*(e_m.^2).*(e_m<0) + ((G(1,1)*G(2,2))+(G(1,2)*G(2,1))).*e_m.*e_p.*(e_m<0).*(e_p<0) + (G(2,1)*G(2,2)).*(e_p.^2).*(e_p<0)...
    + (B(1,1)*B(1,2))*unconditional(1,1) + ((B(1,1)*B(2,2))+(B(1,2)*B(2,1)))*unconditional(2,1) + (B(2,1)*B(2,2))*unconditional(2,2);

figure()
subplot(2,2,1)
surf(e_m, e_p, Z2)
xlabel('Market News')
ylabel('Pair News')
zlabel('Pair Volatility')
title(strcat('Pairs Volatility - ', plot_title))
subplot(2,2,2)
surf(e_m, e_p, Z2)
xlabel('Market News')
ylabel('Pair News')
zlabel('Pair Volatility')
title(strcat('Pairs Volatility - ', plot_title))
subplot(2,2,3)
surf(e_m, e_p, Z3)
xlabel('Market News')
ylabel('Pair News')
zlabel('Covariance')
title(strcat('Covariance - ', plot_title))
subplot(2,2,4)
surf(e_m, e_p, Z3)
xlabel('Market News')
ylabel('Pair News')
zlabel('Covariance')
title(strcat('Covariance - ', plot_title))