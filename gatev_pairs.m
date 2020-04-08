%% Pairs Trading - Gatev 2006
% We will use the clean data (see Data_Cleaning.m for info on how to clean the data) to identify and
% implement Pairs Trading following Gatev (2006).
clear all
clc

%% Setting the size of the portfolio
% Before we proceed to identify the pairs, we need to define how many pairs our portfolio of pairs
% will have. This will be defined by the variable 'P'. For example, setting P=20:
P=40;

%% Matching Pairs
fprintf( 'Reading data... ' ) ;
load('Clean_Data.mat')
load('Market.mat')
fprintf('Done!\n');

% Lets now build a matrix with the normalized prices for each training set and get the index of the 
% pairs:
Pairs_id = [];
Pairs_Sigmas = [];
for m = 13:(numel(t_month_first)-5)
    if mod(m-12, 10) == 0
    fprintf('Processing %d of %d training sets...\n', m-12, numel(t_month_first)-17);
    end
    % Restraining the training set:
    t_training_months = t_month_first(m-12:m);
    training_date_low = dates(t_training_months(1)-1); 
    training_date_high = dates(t_training_months(end)-1);
  
    % Getting the indices for the dates of the training set:
    training_date_low_id = dates_LT(training_date_low - minDate + 1);
    training_date_high_id = dates_LT(training_date_high - minDate + 1);
    
    % Now getting the indices of the stocks:
    s_id = Stocks_id{m-12};
    
    % Building the matrix of normalized prices for this sample:
    norm_prices = prices(training_date_low_id:training_date_high_id, s_id)./prices(training_date_low_id, s_id);

    % Using these values to build the matrix of SSDs:
    SSDs = nan(numel(s_id) , numel(s_id));
    for i = 1:numel(s_id)  
        Dif = norm_prices - norm_prices(:, i);
        SSDs(i,:) = sum(Dif.^2);
    end
    SSDs = triu(SSDs);
    SSDs(SSDs==0)= NaN;    
    
    % Getting the top P pairs:
    i_id = nan(P, 1);
    j_id = nan(P, 1);
    Sigma_ij = nan(P, 1);
    for k=1:P 
       [A, i] = min(SSDs);
       [b, j] = min(A);
       i_id(k) = i(j);
       j_id(k) = j;
       SSDs(i(j), :) = NaN; 
       SSDs(:, i(j)) = NaN; 
       SSDs(j, :) = NaN; 
       SSDs(:, j) = NaN; 
       Sigma_ij(k) = std(norm_prices(:,i(j))-norm_prices(:,j));
    end
    Pairs_id = [Pairs_id;{[i_id,j_id]}];
    Pairs_Sigmas = [Pairs_Sigmas; {Sigma_ij}];
end
fprintf('Done!\n');
clear A b Dif i i_id j j_id k m norm_prices s_id Sigma_ij SSDs t_training_months training_date_low training_date_low_id training_date_high training_date_high_id 

%% Trading Pairs

% Quick adjustment for the index forward:
t_month_first = [t_month_first; numel(dates)+1];

% Trading the pairs:
Profits = [];
for m = 13:(numel(t_month_first)-6)
    if mod(m-12, 10) == 0
    fprintf('Processing %d of %d test sets...\n', m-12, numel(t_month_first)-18);
    end
    
    % First getting the index of the stocks that will be considered in this test set:
    s_id = Stocks_id{m-12};
    i_id = Pairs_id{m-12}(:,1);
    j_id = Pairs_id{m-12}(:,2);
    % Getting the standard deviation of the pairs:
    Sigma_ij = Pairs_Sigmas{m-12};
    
    % Restraining the test set:
    t_test_months = t_month_first(m:m+6);
    test_date_low = dates(t_test_months(1)-1); 
    test_date_high = dates(t_test_months(end)-1);
    % Getting the indices for the dates of the training set:
    test_date_low_id = dates_LT(test_date_low - minDate + 1);
    test_date_high_id = dates_LT(test_date_high - minDate + 1);
    
    % Getting the prices from the test set (required for trading):
    test_prices = prices(test_date_low_id:test_date_high_id, s_id);
    
    % Building the matrix of normalized prices for this test set:
    norm_prices = prices(test_date_low_id:test_date_high_id, s_id)./prices(test_date_low_id, s_id);
    
    % Let's trade these Pairs!!!
    PR = []; 
    paramet = [];
    for pair=1:P
        % Setting initial positions to zero:
        position_i = 0;
        position_j = 0;
        pr_ij = [];
        % Looping for all days in the test set:
        for t=1:size(norm_prices,1)
            % First calculating the spread (normalized) between the two stocks of the pair:
            spread = norm_prices(t,i_id(pair)) - norm_prices(t,j_id(pair));   
            % If we are not trading the pair in 't':
            if position_i == 0
                % Following the 2sd rule of Gatev (2006):
                if spread < -2*Sigma_ij(pair) || spread > 2*Sigma_ij(pair)  
                    if norm_prices(t,i_id(pair)) > norm_prices(t,j_id(pair))
                        position_i = -1;
                        position_j = 1;
                    else
                        position_i = 1;
                        position_j = -1;                    
                    end
                    pair_t0 = t;
                    spread_0 = spread;
                end
           % If we are currently trading the pair:
            else
                if isnan(test_prices(t,i_id(pair))) || isnan(test_prices(t,j_id(pair)))
                    position_i = 0;
                    position_j = 0;
                else
                    % Calculating the pair daily profit: 
                    r_i = log(test_prices(t,i_id(pair)))-log(test_prices(t-1,i_id(pair)));
                    r_j = log(test_prices(t,j_id(pair)))-log(test_prices(t-1,j_id(pair)));
                    profit_ij = position_i*r_i + position_j*r_j;
                    % Adjusting the position for the next day:
                    position_i = position_i*(1+r_i);
                    position_j = position_j*(1+r_j);
                    % Saving the profit, the number of days since the pair opened and the date:
                    pr_ij =[pr_ij; profit_ij, t-pair_t0, t-1+test_date_low_id, t];
                    if sign(spread) ~= sign(spread_0)
                        position_i = 0;
                        position_j = 0;
                    end
                end
            end
        end
        
        % Let's now estimate the vols used for decomposition!
        if ~isempty(pr_ij)
            % First selectiong prices for traning and test sample:
            vol_date_low = dates(t_month_first(m-12)-1); 
            vol_date_low_id = dates_LT(vol_date_low - minDate + 1);
            % Now the series of returns:
            vol_prices = prices(vol_date_low_id:test_date_high_id, s_id);
            ij_rets = diff(log([vol_prices(:,i_id(pair)), vol_prices(:,j_id(pair))]))*100;
            m_rets = r_m(vol_date_low_id+1:test_date_high_id);
            % Using the "pair_vols" to get the vols:
            [var_i, var_j, cov, var_m] = pair_vols(ij_rets, m_rets);
            % Add vol to "profits":
            n_test = size(vol_prices,1)-size(test_prices,1);
            var_i = var_i(n_test:end);
            var_j = var_j(n_test:end);
            cov =  cov(n_test:end);
            var_m = var_m(n_test:end);
            pr_ij = [pr_ij, var_i(pr_ij(:,4)), var_j(pr_ij(:,4)), cov(pr_ij(:,4)), var_m(pr_ij(:,4))];
        end
        % Saving all info into "PR":
        PR = [PR; {pr_ij}];

    end
    
    Profits = [Profits; {PR}];
    
end
clear i_id j_id norm_prices pair pair_t0 position_i position_j PR r_i r_j pr_ij profit_ij s_id ...
    Sigma_ij spread spread_0 t t_test_months test_date_high test_date_high_id test_date_low ...
    test_date_low_id test_prices vol_date_low vol_date_low_id vol_prices ij_rets m_rets n_test ...
    var_i var_j cov var_m pr_ij PR 

save('Results_DM_40', '-v7.3')
fprintf('Done!\n');
