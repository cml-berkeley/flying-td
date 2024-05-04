function usergeom_mat = usergeom_read(nx,ny)
a=importdata('Usergeom01.dat');
b = a';
c = b(:);
d = rmmissing(c(9:end));
e = reshape(d,nx,ny);
usergeom_mat = e';
end