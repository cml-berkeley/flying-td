function htc_3 = map_CML_ANSYS_htc(T_slider,T_disk,FH_ANSYS_mod,FH_ANSYS_orig,PressLoad_ANSYS,PressLoad_ANSYS_orig,Ce_ANSYS,k1_phon,k2_phon,k3_phon,b_phon,k_bulk,lambda0_bulk,sigma,gamma,Pr,max_htc)

%% Get FilmCoeff from FH and PressLoad
htc_3 = get_htc(T_slider,T_disk,FH_ANSYS_mod,FH_ANSYS_orig,PressLoad_ANSYS,PressLoad_ANSYS_orig,Ce_ANSYS,k1_phon,k2_phon,k3_phon,b_phon,k_bulk,lambda0_bulk,sigma,gamma,Pr,max_htc);

fid2 = fopen('Convection_coef_slider_out.dat', 'w');
fprintf(fid2,'%15.5f\t%15.5f\t%15.5f\n',htc_3');
fclose(fid2);


end