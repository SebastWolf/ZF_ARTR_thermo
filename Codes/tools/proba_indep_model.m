function proba = proba_indep_model(config, h_field)

% firing_rates = mean(Rates_A,2);
% config = Rates_A;
%h_field = log(firing_rates./(1-firing_rates));
aa = 1+exp(h_field);
Z_ind = prod(aa);
for t = 1:size(config,2);
    proba(t) = exp(sum(h_field.*config(:,t)))/Z_ind;
end
end

