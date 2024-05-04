tf = 3; % Chose maximum TFC Power index for plotting

figure
plot(TFCPower(1:tf),ECS_final_temp(1:tf)-25,'-o','LineWidth',2)
xlabel('TFC Power (mW)')
ylabel('ECS temperature increase (^oC)')
grid on

figure
plot(TFCPower(1:tf),disk_final_temp(1:tf)-25,'-o','LineWidth',2)
xlabel('TFC Power (mW)')
ylabel('Disk temperature increase (^oC)')
grid on

figure
plot(TFCPower(1:tf),min_FH(1:tf)*1e9,'-o','LineWidth',2)
xlabel('TFCpower (mW)')
ylabel('FH (nm)')
grid on

%% Disk and Slider durface temperature profile plotter
tt = 3; % Chose TFC Power index for plotting

elex = importdata('ele_x.dat');
eley = importdata('ele_y.dat');

nodes_e(:,1) = elex(:);
nodes_e(:,2) = eley(:);

[nodes,ia,~] = unique(nodes_e,'rows'); % nodes = nodes_e(ia)

elems = zeros(1505,3);
for i = 1:1505

elems(i,1) = find(nodes(:,1)==elex(i,1) & nodes(:,2)==eley(i,1));
elems(i,2) = find(nodes(:,1)==elex(i,2) & nodes(:,2)==eley(i,2));
elems(i,3) = find(nodes(:,1)==elex(i,3) & nodes(:,2)==eley(i,3));

end

% Plot slider temperature profile on ABS
slider_temp_plot = slider_temp_output{1,tt};
slider_temp_plot_e = slider_temp_plot(:);

slider_temp_plot_nodes = slider_temp_plot_e(ia);

figure
patch('Vertices',[nodes slider_temp_plot_nodes],'Faces',elems,'FaceVertexCdata',slider_temp_plot_nodes,'FaceColor','interp','linestyle','none');
grid on
xlim([0,843.5])
xlabel('x (\mum)')
ylabel('y (\mum)')
zlabel('T_s')
colorbar
colormap jet

% Plot disk temperature profile
disk_temp_plot = disk_temp_output{1,tt};
disk_temp_plot_e = disk_temp_plot(:);

disk_temp_plot_nodes = disk_temp_plot_e(ia);

figure
patch('Vertices',[nodes disk_temp_plot_nodes],'Faces',elems,'FaceVertexCdata',disk_temp_plot_nodes,'FaceColor','interp','linestyle','none');
grid on
xlim([0,843.5])
xlabel('x (\mum)')
ylabel('y (\mum)')
zlabel('T_d')
colorbar
colormap jet