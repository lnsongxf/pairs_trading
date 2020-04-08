%% Pairs Vol Decomposition 
% We will now decompose the realized volatility of the returns of the individual pairs.
clear all
clc

%% Selecting which portfolio we will analyse:
portfolio = "5";

%% Getting the Sample of all pair returns

% Loading the data:
filename = strcat('Results_DM_',portfolio);
load(filename)  

% Creating the variable vectors that we will fill:
pr = [];
var_i = [];
var_j = [];
cov = [];
var_m = [];
m_set = [];

% Now filling these vectors using the info from "Profits":
for m=1:numel(Profits)
    y = Profits{m};
    for p=1:numel(y)
       x = y{p};
       if ~isempty(x)
           pr = [pr; x(:,1)];
           var_i = [var_i; x(:,5)];
           var_j = [var_j; x(:,6)];
           cov = [cov; x(:,7)];
           var_m = [var_m; x(:,8)];
           m_set = [m_set; ones(size(x,1),1)*m];
       end
    end
end

%% Cleaning the Sample Before Estimation
% Setting the realized volatility of the pair returns:
y = (pr*100).^2;
% Setting number of stds from zero we will remove:
std_dev = 4;
%Setting these outliers to NaN:
y(y>(median(y)+(std(y)*std_dev))|y<(median(y)-(std(y)*std_dev)))=NaN;
var_m(var_m>(median(var_m)+(std(var_m)*std_dev))|var_m<(median(var_m)-(std(var_m)*std_dev)))=NaN;
var_i(var_i>(median(var_i)+(std(var_i)*std_dev))|var_i<(median(var_i)-(std(var_i)*std_dev)))=NaN;
var_j(var_j>(median(var_j)+(std(var_j)*std_dev))|var_j<(median(var_j)-(std(var_j)*std_dev)))=NaN;
cov(cov>(median(cov)+(std(cov)*std_dev))|cov<(median(cov)-(std(cov)*std_dev)))=NaN;

% Removing all i's with any NaN:
y_clean = y(~(isnan(y)|isnan(var_m)|isnan(var_i)|isnan(var_j)|isnan(cov)));
var_m_clean = var_m(~(isnan(y)|isnan(var_m)|isnan(var_i)|isnan(var_j)|isnan(cov)));
var_i_clean = var_i(~(isnan(y)|isnan(var_m)|isnan(var_i)|isnan(var_j)|isnan(cov)));
var_j_clean = var_j(~(isnan(y)|isnan(var_m)|isnan(var_i)|isnan(var_j)|isnan(cov)));
cov_clean = cov(~(isnan(y)|isnan(var_m)|isnan(var_i)|isnan(var_j)|isnan(cov)));
m_set_clean = m_set(~(isnan(y)|isnan(var_m)|isnan(var_i)|isnan(var_j)|isnan(cov)));


%% Doing the Decomposition
% Doing the regression:
X = [var_m_clean, var_i_clean, var_j_clean, cov_clean];
[B_full, TSTAT, ~, ~, ~, R2] = olsnw(y_clean,X,0)

%% Checking Stability

% How many m's (months) we want to include in each estimation?
m_sets = 1;

% Estimating the rolling betas and their std errors:
ms = unique(m_set_clean);
betas = NaN(numel(ms(m_sets:end)),4);
betas_se = NaN(numel(ms(m_sets:end)),4);
i = 1;
for m=ms(m_sets:end)'
    id = (m_set_clean>(m-m_sets))&(m_set_clean<=m);
    X = [var_m_clean(id), var_i_clean(id), var_j_clean(id), cov_clean(id)];
    [B, ~, ~, VCVNW] = olsnw(y_clean(id),X,0);
    betas(i,:) = B'; 
    betas_se(i,:) = sqrt(diag(VCVNW))';
    i=i+1;
end

% Getting the dates:
d = num2str(dates(t_month_first(ms(m_sets:end))));
dates_plot = datetime(d,'InputFormat','yyyyMMdd','Format','dd/MM/yyyy');

% Doing the Plots:
figure()
titles = {'Market Beta', 'Idiosyncratic i Beta', 'Idiosyncratic j Beta', 'Covariance Beta'};
for i = 1:4
    subplot(2,2,i)
    [b, ids] = rmoutstd(betas(:,i),4);
    %CI_up = B_full(i) + 1.96 * betas_se(:,i); 
    %CI_down = B_full(i) - 1.96 * betas_se(:,i);
    %plot(dates_plot(ids), [b, CI_up(ids), CI_down(ids)])
    plot(dates_plot(ids), b)
    title(titles{i})
end

%% Checking Specific Windows
% What time window we are going to use for estimation?
tw = [20100104, 20180601];

% Estimating decomposition for this widow:
m_inf = find(t_month_first==find(dates==tw(1)))-12;
m_sup = find(t_month_first==find(dates==tw(2)))-12;

% Creating the variable vectors that we will fill:
pr = [];
var_i = [];
var_j = [];
cov = [];
var_m = [];
m_set = [];

% Now filling these vectors using the info from "Profits":
for m=m_inf:m_sup
    y = Profits{m};
    for p=1:numel(y)
       x = y{p};
       if ~isempty(x)
           pr = [pr; x(:,1)];
           var_i = [var_i; x(:,5)];
           var_j = [var_j; x(:,6)];
           cov = [cov; x(:,7)];
           var_m = [var_m; x(:,8)];
           m_set = [m_set; ones(size(x,1),1)*m];
       end
    end
end

% Setting the realized volatility of the pair returns:
y = (pr*100).^2;
% Setting number of stds from zero we will remove:
std_dev = 4;
%Setting these outliers to NaN:
y(y>(median(y)+(std(y)*std_dev))|y<(median(y)-(std(y)*std_dev)))=NaN;
var_m(var_m>(median(var_m)+(std(var_m)*std_dev))|var_m<(median(var_m)-(std(var_m)*std_dev)))=NaN;
var_i(var_i>(median(var_i)+(std(var_i)*std_dev))|var_i<(median(var_i)-(std(var_i)*std_dev)))=NaN;
var_j(var_j>(median(var_j)+(std(var_j)*std_dev))|var_j<(median(var_j)-(std(var_j)*std_dev)))=NaN;
cov(cov>(median(cov)+(std(cov)*std_dev))|cov<(median(cov)-(std(cov)*std_dev)))=NaN;
% Removing all i's with any NaN:
y_clean = y(~(isnan(y)|isnan(var_m)|isnan(var_i)|isnan(var_j)|isnan(cov)));
var_m_clean = var_m(~(isnan(y)|isnan(var_m)|isnan(var_i)|isnan(var_j)|isnan(cov)));
var_i_clean = var_i(~(isnan(y)|isnan(var_m)|isnan(var_i)|isnan(var_j)|isnan(cov)));
var_j_clean = var_j(~(isnan(y)|isnan(var_m)|isnan(var_i)|isnan(var_j)|isnan(cov)));
cov_clean = cov(~(isnan(y)|isnan(var_m)|isnan(var_i)|isnan(var_j)|isnan(cov)));
m_set_clean = m_set(~(isnan(y)|isnan(var_m)|isnan(var_i)|isnan(var_j)|isnan(cov)));

% Doing the decomposition:
X = [var_m_clean, var_i_clean, var_j_clean, cov_clean];
[B, TSTAT, ~, ~, ~, R2] = olsnw(y_clean,X,0)
