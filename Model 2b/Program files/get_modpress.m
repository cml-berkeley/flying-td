function [pres_out,C_out] = get_modpress(pres_in,FH_in)

% sigma = 0.5
load('CML_data_mod_FH_sigma_0_5.mat','CMLdata_C', 'CMLdata_FH', 'CMLdata_modFH') 

C_out = interp1(CMLdata_FH,CMLdata_C,FH_in);
pres_out = C_out.*pres_in;

end