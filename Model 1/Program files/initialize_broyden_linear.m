function [] = initialize_broyden_linear(x,x_1,x_2,i1,i2,Usergeom01,xmin,xmax,ymin,ymax,length_slider,width_slider,nx,ny)

usergeom_2 = Usergeom01{1,i2};
usergeom_1 = Usergeom01{1,i1};

usergeom_next = (x-x_1)/(x_2-x_1)*(usergeom_2-usergeom_1) + usergeom_1;
usergeom_write(usergeom_next,xmin,xmax,ymin,ymax,length_slider,width_slider,nx,ny);

end