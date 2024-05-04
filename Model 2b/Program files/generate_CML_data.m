ii =4; % Select index of "usergeom_peak" for which you intend to save the CML Air data

CMLdata_usergeom_peak(1) = usergeom_peak(ii);
CMLdata_FH01{1} = FH01{ii};    
CMLdata_pressload01{1} = pressload01{ii};  
CMLdata_imf01{1} = imf01{ii};  
CMLdata_FH_save1{1} = FH_save1{ii};

save('CML_data_mat.mat','CMLdata_usergeom_peak','CMLdata_FH01','CMLdata_pressload01','CMLdata_imf01','CMLdata_FH_save1')