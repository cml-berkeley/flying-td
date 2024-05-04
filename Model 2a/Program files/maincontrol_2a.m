tic 
clc
clear
%% Input Paramters
usergeom_peak = [11.7:0.3:14.1];

% Slider FEM model parameters
xmin = 0.6;
xmax = 1;
ymin = 0.25;
ymax = 0.75;
length_slider = 0.8435; % in mm
width_slider = 0.7; % in mm
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
load('usergeom01_first_load.mat')
usergeom_peak_first = max(max(Usergeom01_first));
usergeom_next = Usergeom01_first*usergeom_peak(1)/usergeom_peak_first;
usergeom_write(usergeom_next,xmin,xmax,ymin,ymax,length_slider,width_slider,nx,ny);

for i=1:length(usergeom_peak)
    
disp(['TFC Usergeom peak: ',num2str(usergeom_peak(i))])
        
%% CML Air
[x, y] = system('quick.exe');
    
%% FH solution

fid00=fopen('result.dat','r');
for i_count = 1:35
tline = fgetl(fid00);
end
if tline ==-1
    disp('Weird CML Air error!!')
    weird_error = 1;
    fclose(fid00);
else
    weird_error = 0;
    
    read0 = str2num(tline);
    fclose(fid00);
    FH_save1{i}=read0(1,1);
    FH_save2{i}=read0(1,2);
    FH_save3{i}=read0(1,3);
    FH_save4{i}=read0(1,4);

    disp(['FH at TDS: ',num2str(FH_save1{i}(1)-ys*1e9)])

    if abs(FH_save1{i}(1))>0
        fid0=fopen('result.dat','r');
        for i_count = 1:32
        tline = fgetl(fid0);
        end
        read = str2num(tline);
        fclose(fid0);
        run_save{i}(1,:) = read;

        if read(1,1)==-1
            disp('Slider Crashed!!')
        end
    end
end
%% Save Output

a=importdata('Usergeom01.dat');
b = a';
c = b(:);
d = rmmissing(c(9:end));
e = reshape(d,nx,ny);
Usergeom01{1,i} = e';

FH01{1,i}=importdata('FH01.dat');
FilmCoefficient01{1,i}=importdata('FilmCoefficient01.dat');
pressload01{1,i}=importdata('PressLoad01.dat');
imf01{1,i}=importdata('imf01.dat');
cprss01{1,i}=importdata('cprss01.dat');

%% Initialize
if i < length(usergeom_peak)
    load('usergeom01_first_load.mat')
    usergeom_peak_first = max(max(Usergeom01_first));
    usergeom_next = Usergeom01_first*usergeom_peak(i+1)/usergeom_peak_first;
    usergeom_write(usergeom_next,xmin,xmax,ymin,ymax,length_slider,width_slider,nx,ny);
end
if weird_error == 0
    if abs(FH_save1{i}(1))>0
        fid_s=fopen('result.dat','r');
        for i_count = 1:51
        tline = fgetl(fid_s);
        end
        stiff_read = zeros(3,3);
        newStr = extractAfter(tline,"LOAD(G)");
        stiff_read(1,:) = str2num(newStr);
        tline = fgetl(fid_s);
        newStr = extractAfter(tline,"P-TORQUE(uN-M)");
        stiff_read(2,:) = str2num(newStr);
        tline = fgetl(fid_s);
        newStr = extractAfter(tline,"R-TORQUE(uN-M)");
        stiff_read(3,:) = str2num(newStr);
        fclose(fid_s);
        stiff_save{1,i} = stiff_read;
    end
end

save('save_data.mat')

end
save('save_data_final.mat')

toc

%% Plot
for i=1:length(usergeom_peak)
    min_FH(i) = min(min(FH01{1,i}(2:end,2:end)))-ys;
    FH_ECS(i) = FH_save1{i}(1)-ys*1e9;
    max_pres(i) = max(max(pressload01{1,i}(2:end,2:end))); 
    max_cprss(i) = max(max(cprss01{1,i}))*9.81/1e6;
    max_vdw(i) = min(min(imf01{1,i}))*9.81/1e6;
    max_vdw_pos(i) = max(max(imf01{1,i}))*9.81/1e6;
end

figure
plot(usergeom_peak,min_FH*1e9,'-o','LineWidth',2)
xlabel('TFC Protrusion (nm)')
ylabel('Minimum FH (nm)')
grid on

figure
plot(usergeom_peak,FH_ECS,'-o','LineWidth',2)
xlabel('TFC Protrusion (nm)')
ylabel('FH at ECS (nm)')
grid on