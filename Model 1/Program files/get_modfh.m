function FH_out = get_modfh(FH_in)

% sigma = 0.5
load('CML_data_mod_FH_sigma_0_5.mat','CMLdata_C', 'CMLdata_FH', 'CMLdata_modFH') 

FH_out = interp1(CMLdata_FH,CMLdata_modFH,FH_in);

end