function [htc_3] = map_CML_ANSYS_avg(T_slider,T_disk,surf_ele_num,k1_phon,k2_phon,k3_phon,b_phon,k_bulk,lambda0_bulk,sigma,gamma,Pr,max_htc)

PressLoad01 = importdata('PressLoad01.dat');
FH01 = importdata('FH01.dat');

x = PressLoad01(1,2:end)';
y = PressLoad01(2:end,1);
nx = size(x,1); 
ny = size(y,1); 

PressLoad_orig = PressLoad01(2:end,2:end);
[PressLoad,Ce] = get_modpress(PressLoad01(2:end,2:end),FH01(2:end,2:end)); % Make this change for SBL code: get mod pres from pres

FH_orig = FH01(2:end,2:end)-4.594410468582577e-10;
FH_mod = get_modfh(FH01(2:end,2:end)); % Make this change for SBL code: get mod FH from FH

imf01 = importdata('imf01.dat');
IMF = -imf01*9.81/1e6;

elex = importdata('ele_x.dat');
eley = importdata('ele_y.dat');
Vx = elex(:);
Vy = eley(:);

%% Map using averaging
temp_P = zeros(size(Vx));
temp_P_orig = zeros(size(Vx));
temp_IMF = zeros(size(Vx));
temp_Ce = zeros(size(Vx));
temp_FH_mod = zeros(size(Vx));
temp_FH_orig = zeros(size(Vx));

for ii = 1:length(Vx)
    diffValues_x = x - Vx(ii);
    diffValues_x(diffValues_x > 0) = -inf;
    [~, index_x] = max(diffValues_x);

    diffValues_y = y - Vy(ii);
    diffValues_y(diffValues_y > 0) = -inf;
    [~, index_y] = max(diffValues_y);
     
    if index_x<nx
        index_xp1 = index_x+1;
    else
        index_xp1 = index_x-1;
    end

    if index_y<ny
        index_yp1 = index_y+1;
    else
        index_yp1 = index_y-1;
    end
    
    x1 = x(index_x);
    x2 = x(index_xp1);
    y1 = y(index_y);
    y2 = y(index_yp1);
           
    V = [1 x1 y1 x1*y1; 1 x2 y1 x2*y1; 1 x2 y2 x2*y2; 1 x1 y2 x1*y2];
    
    p1 = PressLoad(index_y,index_x);
    p2 = PressLoad(index_y,index_xp1);
    p3 = PressLoad(index_yp1,index_xp1);
    p4 = PressLoad(index_yp1,index_x);
    c = V\[p1; p2; p3; p4];
    temp_P(ii) = c(1) + c(2)*Vx(ii) + c(3)*Vy(ii) + c(4)*Vx(ii)*Vy(ii);
    
    p1_orig = PressLoad_orig(index_y,index_x);
    p2_orig = PressLoad_orig(index_y,index_xp1);
    p3_orig = PressLoad_orig(index_yp1,index_xp1);
    p4_orig = PressLoad_orig(index_yp1,index_x);
    c = V\[p1_orig; p2_orig; p3_orig; p4_orig];
    temp_P_orig(ii) = c(1) + c(2)*Vx(ii) + c(3)*Vy(ii) + c(4)*Vx(ii)*Vy(ii);
    
    ce1 = Ce(index_y,index_x);
    ce2 = Ce(index_y,index_xp1);
    ce3 = Ce(index_yp1,index_xp1);
    ce4 = Ce(index_yp1,index_x);
    c = V\[ce1; ce2; ce3; ce4];
    temp_Ce(ii) = c(1) + c(2)*Vx(ii) + c(3)*Vy(ii) + c(4)*Vx(ii)*Vy(ii);

    IMF1 = IMF(index_y,index_x);
    IMF2 = IMF(index_y,index_xp1);
    IMF3 = IMF(index_yp1,index_xp1);
    IMF4 = IMF(index_yp1,index_x);   
    c = V\[IMF1; IMF2; IMF3; IMF4];
    temp_IMF(ii) = c(1) + c(2)*Vx(ii) + c(3)*Vy(ii) + c(4)*Vx(ii)*Vy(ii);
    
    FH1 = FH_mod(index_y,index_x);
    FH2 = FH_mod(index_y,index_xp1);
    FH3 = FH_mod(index_yp1,index_xp1);
    FH4 = FH_mod(index_yp1,index_x);
    c = V\[FH1; FH2; FH3; FH4];
    temp_FH_mod(ii) = (c(1) + c(2)*Vx(ii) + c(3)*Vy(ii) + c(4)*Vx(ii)*Vy(ii))*1e9;
    
    FH1_orig = FH_orig(index_y,index_x);
    FH2_orig = FH_orig(index_y,index_xp1);
    FH3_orig = FH_orig(index_yp1,index_xp1);
    FH4_orig = FH_orig(index_yp1,index_x);
    c = V\[FH1_orig; FH2_orig; FH3_orig; FH4_orig];
    temp_FH_orig(ii) = (c(1) + c(2)*Vx(ii) + c(3)*Vy(ii) + c(4)*Vx(ii)*Vy(ii))*1e9;
end

%%
PressLoad_ANSYS = reshape(temp_P,[surf_ele_num,3]);

fid0 = fopen('Press_ANSYS.dat', 'w');
fprintf(fid0,'%15.5f\t%15.5f\t%15.5f\n',PressLoad_ANSYS');
fclose(fid0);

PressLoad_ANSYS_orig = reshape(temp_P_orig,[surf_ele_num,3]);
fid0 = fopen('Press_ANSYS_orig.dat', 'w');
fprintf(fid0,'%15.5f\t%15.5f\t%15.5f\n',PressLoad_ANSYS_orig');
fclose(fid0);

Ce_ANSYS = reshape(temp_Ce,[surf_ele_num,3]);

IMF_ANSYS = reshape(temp_IMF,[surf_ele_num,3]);

fid0 = fopen('vdW_ANSYS.dat', 'w');
fprintf(fid0,'%15.5f\t%15.5f\t%15.5f\n',IMF_ANSYS');
fclose(fid0);

FH_ANSYS_mod = reshape(temp_FH_mod,[surf_ele_num,3]);

fid1 = fopen('FH_ANSYS_mod.dat', 'w');
fprintf(fid1,'%15.5f\t%15.5f\t%15.5f\n',FH_ANSYS_mod');
fclose(fid1);

FH_ANSYS_orig = reshape(temp_FH_orig,[surf_ele_num,3]);

fid1 = fopen('FH_ANSYS_orig.dat', 'w');
fprintf(fid1,'%15.5f\t%15.5f\t%15.5f\n',FH_ANSYS_orig');
fclose(fid1);

%% Get FilmCoeff from FH and PressLoad
htc_3 = get_htc(T_slider,T_disk,FH_ANSYS_mod,FH_ANSYS_orig,PressLoad_ANSYS,PressLoad_ANSYS_orig,Ce_ANSYS,k1_phon,k2_phon,k3_phon,b_phon,k_bulk,lambda0_bulk,sigma,gamma,Pr,max_htc);

fid2 = fopen('Convection_coef_slider_out.dat', 'w');
fprintf(fid2,'%15.5f\t%15.5f\t%15.5f\n',htc_3');
fclose(fid2);

end