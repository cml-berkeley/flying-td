tic 
clc
clear
%% Input Parameters
TFCPower = 102;
RPM = 5400; % disk RPM
rr = 27e-3; % radial position of slider
mu = 0.7; % Friction coefficient

% Phonon parameters
% We model the phonon conduction htc using the following formula:
% log(htc_phon) = k1_phon*(log(T_disk+273.15)-log(273.15+25)) + k2_phon*(log(dT)-log(400)) + k3_phon*log(h_in_orig) + b_phon
% T_disk = disk temperature in deg C
% dT = slider-disk temperature difference in deg C
% h = spacing in nm
k1_phon = 0.99;
k2_phon = -0.83;
k3_phon = -1.99;
b_phon = 11.4;

% Air parameters
k_bulk = 0.0261; % Air thermal conductivity
lambda0_bulk = 67.1e-9; % Air Mean Free Path 
sigma = 0.6; % Thermal Accomodation Coefficient 
gamma = 1.4015; % Ratio of Specific Heat Capacity for Air
Pr = 0.71; % Prandtl Number for Air

max_htc = 1.5e7; % Interface thermal conductance for the HDI (applied as an upper limit for the HDI htc)

% Slider FEM parameters
surf_ele_num = 1505; % number of surface elements on ABS surface of slider ANSYS model
nx = 401;
ny = 401;

%% ROUGHNESS PARAMETERS
R = 0.020e-6; % asperity radius of curvature
eta = 5000e12; % asperity areal density
sigma_z = 0.5e-9; % std deviation of surface heights

sigma_a = sqrt(sigma_z^2-3.717e-4/(eta^2*R^2)); % std deviation of apserity height
alpha = 2412.74*(sigma_z*R*eta)^2;
ys = 4*sigma_z/sqrt(pi*alpha);

%%
load('CML_data_mat.mat')
max_iteration = 10;
tol_FH = 0.05; % tol for convergence = diff in FH between succesive iterations of 0.05 nm
tol_conv = 1e-2;  % tol for convergence = rel error between input and output htc of 1e-3

U = RPM*2*pi/60*rr; 

for i=1:length(TFCPower)
    
fid_a=fopen('MOAI.txt','w');
fprintf(fid_a,'%8.2f',TFCPower(i));
fclose(fid_a);
disp(['TFC Power: ',num2str(TFCPower(i))])

error_save_FH = zeros(max_iteration,1);
error_save_conv = zeros(max_iteration,1);
error_save_conv_ideal = zeros(max_iteration,1);
FH_save1{i}=zeros(120,1);
    
main_iteration_flag = 1;
main_iteration_count = 1;

while main_iteration_flag
    
    fid_b=fopen('main_iteration_count.dat','w');
    fprintf(fid_b,'%15.5f',main_iteration_count');
    fclose(fid_b);
    
    if main_iteration_count == 1
        % CML Air lookup table
%         usergeom_in = usergeom_read(nx,ny);
%         usergeom_peak_in = max(max(usergeom_in));
%         [~,closestIndex] = min(abs(CMLdata_usergeom_peak-usergeom_peak_in));
        closestIndex = 1;
        PressLoad01 = CMLdata_pressload01{1,closestIndex};
        FH01 = CMLdata_FH01{1,closestIndex};
        FH_save1_first= CMLdata_FH_save1{1,closestIndex}(1);
        imf01 = CMLdata_imf01{1,closestIndex};
        
        T_slider=importdata('slider_temp_initial.dat');
        T_disk=importdata('disk_temp_initial.dat');
        [Ce_ANSYS,PressLoad_ANSYS_orig,PressLoad_ANSYS,IMF_ANSYS,FH_ANSYS_mod,FH_ANSYS_orig,QF_ANSYS] = map_CML_ANSYS_avg(PressLoad01,FH01,imf01,U,surf_ele_num,mu);
        ele_conv_input = map_CML_ANSYS_htc(T_slider,T_disk,FH_ANSYS_mod,FH_ANSYS_orig,PressLoad_ANSYS,PressLoad_ANSYS_orig,Ce_ANSYS,k1_phon,k2_phon,k3_phon,b_phon,k_bulk,lambda0_bulk,sigma,gamma,Pr,max_htc);
                
        fid = fopen('Convection_coef_initial_BC.dat', 'w');
        fprintf(fid,'%15.5f\t%15.5f\t%15.5f\n',ele_conv_input');
        fclose(fid);
%         disp([' Usergeom peak for first CML Air solve: ',num2str(usergeom_peak_in),' nm'])
%        disp('First CML Air solve done to get initial htc guess')
        toc;
    else
        ele_conv_input = importdata('Convection_coef_initial_BC.dat');
    end

    %% Slider & Disk
    slider_disk_iteration_flag = 1;
    slider_disk_iteration_count = 1;
    while slider_disk_iteration_flag
        system('set kmp_stacksize=2048K & "C:\Program Files\ANSYS Inc\v140\ansys\bin\winx64\ansys140.exe" -b -i "C:\Users\siddhesh_sakhalkar\Desktop\Siddhesh Codes\Flying touchdown\Model 2b\SolDeformMag.mac" -o "C:\Users\siddhesh_sakhalkar\Desktop\Siddhesh Codes\Flying touchdown\Model 2b\ANSYS_output.txt" -np 7');
        if slider_disk_iteration_count == 1
            [slider_disk_iteration_flag,x_old_disk,f_old_disk,B_old_disk] = disk_temp_compute(slider_disk_iteration_count,0,0,0,U,surf_ele_num); 
        else
            [slider_disk_iteration_flag,x_old_disk,f_old_disk,B_old_disk] = disk_temp_compute(slider_disk_iteration_count,x_old_disk,f_old_disk,B_old_disk,U,surf_ele_num);
        end
        
        slider_disk_iteration_count = slider_disk_iteration_count + 1;
    end
    toc;
    
    %% CML Air
    usergeom_in = usergeom_read(nx,ny);
    usergeom_peak_in = max(max(usergeom_in));
    usergeom_ANSYS_save{i}(main_iteration_count) = usergeom_peak_in;
    disp(['Iteration number: ',num2str(main_iteration_count),' Usergeom peak: ',num2str(usergeom_peak_in),' nm'])
%     [~,closestIndex] = min(abs(CMLdata_usergeom_peak-usergeom_peak_in));   
    closestIndex = 1;
    PressLoad01 = CMLdata_pressload01{1,closestIndex};
    FH01 = CMLdata_FH01{1,closestIndex};
    FH_save1{i}(main_iteration_count) = CMLdata_FH_save1{1,closestIndex}(1);
    
    %% Check for convergence in FH
    
    if main_iteration_count>1
        error_save_FH(main_iteration_count) = abs(FH_save1{i}(main_iteration_count) - FH_save1{i}(main_iteration_count-1));
    else
        error_save_FH(main_iteration_count) = abs(FH_save1{i}(main_iteration_count) - FH_save1_first);
    end
        
    disp(['Iteration number: ',num2str(main_iteration_count),' Abs Error in FH: ',num2str(error_save_FH(main_iteration_count)),' nm'])
    toc;
    
    %% Get "Convection_coef_initial_BC.dat" in ANSYS mesh using "FH01.dat" and "Pressload01.dat" from CML Air mesh and T_slider, T_disk from ANSYS mesh
    T_slider=importdata('slider_temp.dat');
    T_disk=importdata('disk_temp.dat');
    ele_conv_output = map_CML_ANSYS_htc(T_slider,T_disk,FH_ANSYS_mod,FH_ANSYS_orig,PressLoad_ANSYS,PressLoad_ANSYS_orig,Ce_ANSYS,k1_phon,k2_phon,k3_phon,b_phon,k_bulk,lambda0_bulk,sigma,gamma,Pr,max_htc);
    
    %% Broyden to get next ele_conv_input
    
    % Broyden
    x = ele_conv_input(:);
    f = ele_conv_output(:);
    if main_iteration_count == 1
        B = 0.5*speye(surf_ele_num*3);
        x_next = x - B*(x - f);
    else
        s = x - x_old;
        y = (x - f) - (x_old - f_old);
        B = B_old + (s-B_old*y)*y'/(y'*y); 
        x_next = x - B*(x - f);
    end
    
    error = (x-f)./x;
    error_save_conv_ideal(main_iteration_count) = max(abs(error));
    
    ttt = 8; % ele_x = 831 and ele_y = 360 
    error1 = abs((x(ttt)-f(ttt))./x(ttt));
    [~,imax]=max(x);
    error2 = abs((x(imax)-f(imax))./x(imax));
    error_save_conv(main_iteration_count) = max(abs(error));
    conv_TDS_input{i}(main_iteration_count)=x(ttt);
    conv_TDS_output{i}(main_iteration_count)=f(ttt);
    
    disp(['Iteration number: ',num2str(main_iteration_count),' Rel Error in conv: ',num2str(error_save_conv(main_iteration_count))])
    x_old = x;
    f_old = f;
    B_old = B;
    
    %% Check for convergence in ele_conv
    
    if error_save_conv(main_iteration_count) < tol_conv && error_save_FH(main_iteration_count) < tol_FH
        main_iteration_flag = 0; 
    end   
    
    main_iteration_count = main_iteration_count + 1;
    if main_iteration_count > max_iteration
       main_iteration_flag = 0; 
    end 
    
    if main_iteration_flag == 1
        ele_conv_next = reshape(x_next,[surf_ele_num,3]);
        fid_1 = fopen('Convection_coef_initial_BC.dat', 'w');
        fprintf(fid_1,'%15.5f\t%15.5f\t%15.5f\n',ele_conv_next');
        fclose(fid_1);
    end
end

%% Save Output

fid2=fopen('TDStemp.dat','r');
TDStemp{i}=fscanf(fid2,'%f');
fclose(fid2);
fid2b=fopen('TDStemp.dat','w');
fprintf(fid2b,'%15.11f\n',zeros(1,120));
fclose(fid2b);

fid3=fopen('TFCprotrusion.dat','r');
TFCprotrusion{i}=fscanf(fid3,'%f');
fclose(fid3);
fid3b=fopen('TFCprotrusion.dat','w');
fprintf(fid3b,'%15.11f\n',zeros(1,120));
fclose(fid3b);

a=importdata('Usergeom01.dat');
b = a';
c = b(:);
d = rmmissing(c(9:end));
e = reshape(d,nx,ny);
Usergeom01{1,i} = e';

FH01_save{1,i}=FH01;
pressload01{1,i}=PressLoad01;
imf01_save{1,i}=CMLdata_imf01{1,closestIndex};

slider_temp_output{1,i}=importdata('slider_temp.dat');
disk_temp_output{1,i}=importdata('disk_temp.dat');
disk_temp_eff_output{1,i}=importdata('disk_temp_eff.dat');

Convection_coef_input{1,i}=importdata('Convection_coef_initial_BC.dat');
Convection_coef_output{1,i}=importdata('Convection_coef_slider_out.dat');

error_array_conv{i}=error_save_conv;
error_array_FH{i}=error_save_FH;

save('save_data.mat')

end
disp(['Usergeom peak is ',num2str(usergeom_peak_in),' nm'])
ECS_final_temp = slider_temp_output{1,1}(72,1)
disk_final_temp = disk_temp_output{1,i}(72,1)
min_FH(i) = min(min(CMLdata_FH01{1,i}(2:end,2:end)))-ys
toc