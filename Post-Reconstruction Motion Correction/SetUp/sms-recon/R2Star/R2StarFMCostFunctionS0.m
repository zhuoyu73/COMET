function [f,g] = R2StarFMCostFunctionS0(x, data, sen,rInfo,Rr2,Rwe,RS0,maskTE,mask,nc_svd,tt_ext)

FOV= rInfo.FOV/10;
J= rInfo.J;
L= rInfo.L;
N = rInfo.N;
num_coils = rInfo.nCoils;
TEtp = (rInfo.TE - rInfo.TE(1)).*1e-6;
TE = TEtp(maskTE);
numechoes = sum(maskTE);
ttp = rInfo.timingVec(:,:,1,1,1,1,maskTE);
N = rInfo.N;
han_flag = 1;

npm=N*N;
nx=2;
ny=2;
mskpenalwe = logical(ones(N,N));
mskpenalr2 = logical(mask);

nSize = length(x);

r2now = x(1:nSize/4);
wenow = x(nSize/4+1:(2*nSize/4));
S0real =  x((2*nSize/4+1):(3*nSize/4));
S0imag =  x((3*nSize/4+1):end);
S0 = S0real + 1j*S0imag;

for kk = 1 : rInfo.nEchoes
    A{kk}= fast_mr_wt_quad(col(rInfo.kx(:,:,1,1,1,1,kk)),col(rInfo.ky(:,:,1,1,1,1,kk)),FOV,N,2*N,J,col(ttp(:,:,1,1,1,1,kk)),han_flag,wenow,r2now,mask,1,L,1,[],0,nx,ny,1,0);
end

A_full = multi_echo_mri(A,numechoes);
A_fulls=sense_svd(A_full,sen,nc_svd);
datas = prepData(A_fulls,col(data));

resid = datas-A_fulls*S0;
gg = (A_fulls'*((tt_ext).*resid));
ggS0 = (A_fulls'*(resid));

roughr2 = Rr2.penal(Rr2,r2now(mskpenalr2));
roughwe = Rwe.penal(Rwe,wenow);

f = 1/2*(real(resid'*resid))+roughr2+roughwe;

if 1
    subplot(3,1,1), im(reshape(r2now,N,N),[0,40])
    subplot(3,1,2), im(reshape(wenow,N,N),[-500,500])
    subplot(3,1,3), im(reshape(S0,N,N))
    drawnow
end

if nargout > 1
    gradr2 = real(conj(S0).*gg(mskpenalr2))+double(Rr2.cgrad(Rr2,r2now(mskpenalr2)));
    gradwe = (imag(conj(S0).*gg)+double(Rwe.cgrad(Rwe,wenow)));
    gradS0real = -1*real(ggS0) + double(RS0.cgrad(RS0,real(S0)));
    gradS0imag = -1*imag(ggS0) + double(RS0.cgrad(RS0,imag(S0)));
    g = vertcat(gradr2,gradwe,gradS0real,gradS0imag);
    
end
end