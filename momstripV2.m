%constants
netta = 377;
kk = 2*pi;
gamma = 1.781072418;

%define cells
NumNodes =150; %number of nodes in 1 direction
stripWidth  = 5;%lambda
delx = stripWidth/NumNodes;
% figure;
phi_inc = pi/2;

%solve
Einc = zeros(NumNodes,1);
zz = zeros(NumNodes);

for mm = 1:NumNodes
    xm = delx*mm;
    xmprev = delx*(mm-1/2);
    xmnext = delx*(mm+1/2);

    
    Einc(mm) = -1j*tan(phi_inc)*(netta/kk) * ( exp(-1j*kk*xmnext*cos(phi_inc))-exp(-1j*kk*xmprev*cos(phi_inc)) );
    for nn = 1:NumNodes
        xn = delx*nn;
        xnprev = delx*nn - delx;
        xnnext = delx*nn + delx;
        
        Rmn = abs(xm-xn);

        if mm==nn
            int1 = delx*(1-2j/pi *(log(gamma*kk*delx/8)-1));
            
            zz(mm,nn) = delx*(kk*netta/4)*(int1);
        else
            int1 = delx/2*besselh(0,2,kk*abs(xm-xn));
            int2 = delx/2*besselh(0,2,kk*abs(xm-(xn+delx/2)));
            int3 = besselh(0,2,kk*abs(xmnext - xn));
            int4 = -besselh(0,2,kk*abs(xmnext - xnnext));
            int5 = -besselh(0,2,kk*abs(xmprev-xn));
            int6 = besselh(0,2,kk*abs(xmprev-xnnext));
            
            zz(mm,nn) = delx*(kk*netta/4)*(int1+int2) ...
                      + netta/(4*kk)*(int3+int4+int5+int6);
        end
        
                   
    end
end


jj =zz\Einc;
plot(abs(jj))
jj= jj.';
% jj(1)=0;
%ok, so now I have the coefficients, jn, but because of my choice of basis
%functions, these are not directly JJ because triangles are fun I guess.
numPointsinCell = 200;
JJ = getTriangleBasis( jj, numPointsinCell,delx);
tt = linspace(0,stripWidth,length(JJ));
figure;plot(tt,abs(JJ))
