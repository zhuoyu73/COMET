function [f,g] = R2StarFMCostFunction(x, data, sen,rInfo,Rr2,Rwe,maskTE,mask,nc_svd,tt_ext,S0)

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

r2now = x(1:nSize/2);
wenow = x(nSize/2+1:end);

for kk = 1 : rInfo.nEchoes
    A{kk}= fast_mr_wt_quad(col(rInfo.kx(:,:,1,1,1,1,kk)),col(rInfo.ky(:,:,1,1,1,1,kk)),FOV,N,2*N,J,col(ttp(:,:,1,1,1,1,kk)),han_flag,wenow,r2now,mask,1,L,1,[],0,nx,ny,1,0);
end

A_full = multi_echo_mri(A,numechoes);
A_fulls=sense_svd(A_full,sen);
datas = prepData(A_fulls,col(data));

resid = datas-A_fulls*S0;
gg = conj(S0).*(A_fulls'*((tt_ext).*resid));

roughr2 = Rr2.penal(Rr2,r2now(mskpenalr2));
roughwe = Rwe.penal(Rwe,wenow);

f = 1/2*(real(resid'*resid))+roughr2+roughwe;

if 1
    subplot(2,1,1), im(reshape(r2now,N,N),[0,40])
    subplot(2,1,2), im(reshape(wenow,N,N),[-500,500])
    drawnow
end

if nargout > 1
    gradr2 = real(gg(mskpenalr2))+double(Rr2.cgrad(Rr2,r2now(mskpenalr2)));
    gradwe = (imag(gg)+double(Rwe.cgrad(Rwe,wenow)));
    g = vertcat(gradr2,gradwe);
    
end
end