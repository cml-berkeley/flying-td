function temp_QF = getQ(imf,u,mu)
% imf is in MPa
NNCF = imf*1e6;
NNCF(NNCF<0)=0;
temp_QF = 0.5*mu*NNCF*u;

end