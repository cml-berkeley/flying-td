function htc_func = get_htc(T_slider,T_disk,h_in_mod,h_in_orig,p,p_orig,Ce,k1_phon,k2_phon,k3_phon,b_phon,k_bulk,lambda0_bulk,sigma,gamma,Pr,max_htc)
% T_slider & T_disk is in degC
% h_in is in nm
% p is pressure above ambient in MPa

dT = T_slider-T_disk;
dT(dT < 4) = 4;
T_disk(T_disk < 25) = 25;
h_in_mod(h_in_mod<0) = 0;
h_in_orig(h_in_orig<0.01) = 0.01;


% %% Phonon parameters
% k1_phon = 0.99;
% k2_phon = -0.83;
% k3_phon = -1.99;
% b_phon = 11.4;
% 
% %% Air parameters
% k_bulk = 0.0261; 
% lambda0_bulk = 67.1e-9; 
% sigma = 0.6; 
% gamma = 1.4015; 
% Pr = 0.71; 

%%
b_air = (2-sigma)/sigma*2*gamma/(gamma+1)/Pr;

p0 = 101325*1e-6; 
T0 = 273.15+25;
T = (T_slider + T_disk)/2;
lambda_bulk = lambda0_bulk*p0./(p_orig+p0).*(T+273.15)./T0; 

%% Evaluate htc_fun due to phonon using h_in_orig
htc_func_1 = exp( k1_phon*(log(T_disk+273.15)-log(273.15+25)) + k2_phon*(log(dT)-log(400)) + k3_phon*log(h_in_orig) + b_phon );

%% Evaluate htc_fun due to air using h_in_mod
lambda_eff = zeros(size(h_in_mod));
lambda_eff(h_in_mod*1e-9<lambda_bulk) = 3/4.*h_in_mod(h_in_mod*1e-9<lambda_bulk)*1e-9 ...
    - h_in_mod(h_in_mod*1e-9<lambda_bulk).*1e-9/2.*log(h_in_mod(h_in_mod*1e-9<lambda_bulk).*1e-9./lambda_bulk(h_in_mod*1e-9<lambda_bulk));
lambda_eff(h_in_mod*1e-9>=lambda_bulk) = lambda_bulk(h_in_mod*1e-9>=lambda_bulk) ...
    - lambda_bulk(h_in_mod*1e-9>=lambda_bulk).^2./4./(h_in_mod(h_in_mod*1e-9>=lambda_bulk).*1e-9);

k_eff = k_bulk./lambda_bulk.*lambda_eff; 

h_eff = h_in_mod.*1e-9 + 2.*b_air.*lambda_eff;
htc_func_2 = k_eff./h_eff.*Ce;          

%% Total htc
htc_func = htc_func_1 + htc_func_2;

htc_func(htc_func > max_htc) = max_htc;
htc_func(htc_func < 100) = 100;

end