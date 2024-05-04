fid0=fopen('result.dat','r');
for i_count = 1:32
tline = fgetl(fid0);
end
read = str2num(tline);
fclose(fid0);

fid=fopen('Run.dat','rt');
X = fread(fid) ;
fclose(fid) ;
X = char(X.') ;

Y = X;
Y(278:278+40) = num2str([read(2)*1e-9,read(3)*1e-6,read(4)*1e-6],'%+ .4e    ');

fid2 = fopen('Run.dat','wt') ;
fwrite(fid2,Y) ;
fclose(fid2) ;