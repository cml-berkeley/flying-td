function [slider_disk_iteration_flag,x_old,f_old,B_old] = disk_temp_compute(slider_disk_iteration_count,x_old,f_old,B_old,U,surf_ele_num)
T_eff = zeros(surf_ele_num*3,1);

rho = 2700;
Cp = 960;
k = 117;
cc = 1.6e7;

elex = importdata('ele_x.dat');
eley = importdata('ele_y.dat');
xx = elex*1e-6;
yy = eley*1e-6;

xx_ele = repmat(sum(elex,2)/3*1e-6,1,3);
yy_ele = repmat(sum(eley,2)/3*1e-6,1,3);

x = xx(:);
y = yy(:);
z = 0;

Area = abs((xx(:,1).*(yy(:,2)-yy(:,3)) + xx(:,2).*(yy(:,3)-yy(:,1)) ...
    + xx(:,3).*(yy(:,1)-yy(:,2)))/2);

htc = importdata('Convection_coef_initial_BC.dat');
Ts = importdata('slider_temp_initial.dat');
slider_temp_final = importdata('slider_temp.dat');
Td = importdata('disk_temp_initial.dat');
Q = htc.*(Ts-Td);

for e = 1:surf_ele_num*3
%% Numerical integration
    r = sqrt((x(e)-xx_ele).^2+(y(e)-yy_ele).^2+z.^2);
    I = Q./(2*pi*k)./r.*exp((U*rho*Cp/2/k).*(x(e)-xx_ele-r));
    
    II = sum(I,2).*Area/3;
    T_eff(e) = 25 + sum(II);

end

T = T_eff+Q(:)/cc;

%% Override to test

% disk_temp_final_eff = 25*ones(surf_ele_num,3);
% disk_temp_final = 25*ones(surf_ele_num,3);

%% Do not override

disk_temp_final_eff = reshape(T_eff,[surf_ele_num,3]);
fid = fopen('disk_temp_eff.dat', 'w');
fprintf(fid,'%15.11f\t%15.11f\t%15.11f\n',disk_temp_final_eff');
fclose(fid);

disk_temp_final = reshape(T,[surf_ele_num,3]);
fid2 = fopen('disk_temp.dat', 'w');
fprintf(fid2,'%15.11f\t%15.11f\t%15.11f\n',disk_temp_final');
fclose(fid2);

%% Check for convergence
x = [Ts(:); Td(:)];
f = [slider_temp_final(:); disk_temp_final(:)];
error = max(max(abs((x - f)./x)));

if error < 1e-3
    slider_disk_iteration_flag = 0;
    disp(['Slider-disk temp computation Iteration number: ',num2str(slider_disk_iteration_count),' Relative Error: ',num2str(error)])
elseif slider_disk_iteration_count > 50
    slider_disk_iteration_flag = 0;
    disp(['Slider-disk temp computation Max Iteration number reached: ',num2str(slider_disk_iteration_count),' Relative Error: ',num2str(error)])
else
    slider_disk_iteration_flag = 1;
    
    %% Broyden
    if slider_disk_iteration_count == 1
        B = 0.1*speye(2*surf_ele_num*3);
        x_next = x - B*(x - f);
    else
        s = x - x_old;
        y = (x - f) - (x_old - f_old);
        B = B_old + (s-B_old*y)*y'/(y'*y);
        x_next = x - B*(x - f);
    end
    x_old = x;
    f_old = f;
    B_old = B;
    
    slider_temp_next_in = reshape(x_next(1:surf_ele_num*3),[surf_ele_num,3]);
    fid3 = fopen('slider_temp_initial.dat', 'w');
    fprintf(fid3,'%15.11f\t%15.11f\t%15.11f\n',slider_temp_next_in');
    fclose(fid3); 
    
    disk_temp_next_in = reshape(x_next(surf_ele_num*3+1:2*surf_ele_num*3),[surf_ele_num,3]);
    fid3 = fopen('disk_temp_initial.dat', 'w');
    fprintf(fid3,'%15.11f\t%15.11f\t%15.11f\n',disk_temp_next_in');
    fclose(fid3); 
    
end
end