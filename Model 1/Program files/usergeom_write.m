function [] = usergeom_write(usergeom_mat,xmin,xmax,ymin,ymax,length_slider,width_slider,nx,ny)
fid = fopen('Usergeom01.dat','w');
% L1 = [0.50610    0.84350    0.17500    0.52500];
% L2 = [401   401];
L1 = [xmin*length_slider	xmax*length_slider     ymin*width_slider	ymax*width_slider];
L2 = [nx		ny];

fprintf(fid,'%10.5f',L1);
fprintf(fid,'\n');

fprintf(fid,'     %d      %d',L2);
fprintf(fid,'\n');

for i = 1:nx
fprintf(fid,'%10.5f',usergeom_mat(i,:));
fprintf(fid,'\n');
end
fclose(fid);
end