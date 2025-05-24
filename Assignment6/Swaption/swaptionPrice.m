function [price,tree] = swaptionPrice(t_omega, strike, dates, discounts, deltaT)

% Hull-White model parameters
a = 10/100;
sigma = 0.8/100;
ttm = 10;

% Tree parameters
N_steps_in_1y = round(1/deltaT);
N_year_p1 = round(1/deltaT)+1;
N_steps = N_steps_in_1y*ttm;

[tree,l_max,N_max] = treeHW(a, sigma, t_omega, deltaT); l_max = round(l_max);

years = dates(1)+calyears(2):calyears(1):dates(1)+calyears(t_omega);
disc_year = interpolation(discounts, dates, years);
fwd_year_intra = disc_year(2:end)./disc_year(1:end-1);
fwd_year = disc_year(end)./disc_year(1:end-1); % 0 to T

sigmaHJM = @(s,T) sigma/a * (1 - exp(-a * (T - s)));

for i=1:8
    % yearly
    x = [tree{N_year_p1+N_steps_in_1y*i}.x];
    sigma_prev = sigmaHJM(i+1,i+1); % is zero
    sigma_next = sigmaHJM(i+1,10);
    
    integranda = @(t) sigmaHJM(t,i+1).^2 - sigmaHJM(t,10).^2;
    intSigma = quadgk(integranda, 0, i+1);

    B_y = fwd_year(i)*exp(-x/sigma*(sigma_next-sigma_prev)-0.5*intSigma);
    for j=1:length(B_y)
        tree{N_year_p1+N_steps_in_1y*i}(j).B_y = B_y(j);
    end

    x = [tree{N_year_p1+N_steps_in_1y*i}.x];
    sigma_prev = sigmaHJM(i+1,i+1); % is zero
    sigma_next = sigmaHJM(i+1,i+2);
    
    integranda = @(t) sigmaHJM(t,i+1).^2 - sigmaHJM(t,i+2).^2;
    intSigma = quadgk(integranda, 0, i+1);

    B_BPV = fwd_year_intra(i)*exp(-x/sigma*(sigma_next-sigma_prev)-0.5*intSigma);
    for j=1:length(B_BPV)
        tree{N_year_p1+N_steps_in_1y*i}(j).B_BPV = B_BPV(j);
    end
end

intrinsic_value = zeros(ttm-1,N_max);

for i=1:8
    B_BPV = zeros(ttm-i,N_max);

    % disp(size(B_BPV))

    for j=i:8
        % disp(size([tree{N_year_p1+N_steps_in_1y*j}.B_BPV]))
        B_BPV(j,:) = [tree{N_year_p1+N_steps_in_1y*j}.B_BPV];
    end

    B_y = [tree{N_year_p1+N_steps_in_1y*i}.B_y];
    BPV = sum(B_BPV, 1);

    for j=1:length(BPV)
        tree{N_year_p1+N_steps_in_1y*i}(j).BPV = BPV(j);
        tree{N_year_p1+N_steps_in_1y*i}(j).S = (1-B_y(j))/BPV(j);
    end
    
    S = [tree{N_year_p1+N_steps_in_1y*i}.S];
    intrinsic_value(i+1,:) = BPV.*max(S-strike,0);
end


% codice loro
x = [tree{end}.x];

node_dates = (datenum(dates(1)) + deltaT*(0:N_steps)*365)';
node_dates = datetime(node_dates, 'ConvertFrom', 'datenum');

discounts_node = interpolation(discounts, dates, node_dates);
discounts_node(1) = 1;

fwd_discount_nodes = discounts_node(2:end)./discounts_node(1:end-1);
sigma_star = (sigma/a) * sqrt(deltaT - 2 *( (1 - exp(-a*deltaT)) / a ) + (1 - exp(-2*a*deltaT)) / (2*a) );

% initialize
fwdDF_present = zeros(N_max, N_steps);

for i=1:N_steps
    sigmaNext = sigmaHJM(0, deltaT);
    integranda = @(t) sigmaHJM(0, t).^2 - sigmaHJM(0,t+deltaT).^2;
    int_sigma = quadgk(integranda, 0, i);

    fwdDF_present(:,i) = fwd_discount_nodes(i)*exp(-x * sigmaNext - 0.5*int_sigma);
end

value = zeros(N_max, (ttm)*N_steps_in_1y+1);

sigma_hat = sigma * sqrt((1-exp(-2*a*deltaT))/(2*a));
D_x = sigma_hat * sqrt(3);
l_min = -l_max;
mu_hat = 1- exp(-a*deltaT);
reset_dates = (datenum(dates(1)) + (1:ttm)*365)';
reset_dates = datetime(reset_dates, 'ConvertFrom', 'datenum');

intrinsic_value = intrinsic_value';

for i = (ttm)*N_steps_in_1y:-1:1

    % compute the continuation value
    for j = 1:N_max

        if j == 1
            value(j,i) = C_contvalue(l_max, mu_hat, sigma, sigma_star, a, deltaT, fwdDF_present(:,i), value(:,i+1), D_x, x);

        elseif j == N_max
            value(j,i) = B_contvalue(l_min, mu_hat, sigma, sigma_star, a, deltaT, fwdDF_present(:,i), value(:,i+1), D_x, x);

        else
            value(j,i) = A_contvalue(l_max - j + 1, j, mu_hat, sigma, sigma_star, a, deltaT, fwdDF_present(:,i), value(:,i+1), D_x, x);            
        end

    end
    

   if find(reset_dates == node_dates(i))
    
        % find the index of the reset date
        index = find(reset_dates == node_dates(i));

        % not sure
        if index == ttm
            index = index-1;
        end

        continuation_value = max(intrinsic_value(:,index), value(:,i));
        
        value(:,i) = continuation_value;
    
    else
        continue
    end

end

price = value(l_max+1,1);

end