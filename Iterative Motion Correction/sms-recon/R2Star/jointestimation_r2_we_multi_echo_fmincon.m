function [r2final, wefinal,s0final] = jointestimation_r2_we_multi_echo_fmincon(data,maskTE,rInfo, FM, sen, r2map,han_flag, numloop,numloop_ext,mask)
%jointestimation_r2_we_multi_echo
%   data: shoudl be reaodut_length X shot X slice X coil
%   maskTE = binary list of echo time to keep
%   slice = slice to process
%   rInfo from recoInfo
%   FM: inital field map
%   sen: SENSE map
%   r2map: initial estimate of R2* map
%   han_flag: use hanning for time segmentation
%   numloop: number of iterations on we and R2*
%   numloop_ext: number of external loop on I0 image
%   mask

%%% k-space 1st echo
ReturnComplete = false;
FOV= rInfo.FOV/10;
J= rInfo.J;
L= rInfo.L;
N = rInfo.N;
num_coils = rInfo.nCoils;
TEtp = (rInfo.TE - rInfo.TE(1)).*1e-6;
TE = TEtp(maskTE)
numechoes = sum(maskTE)
ttp = rInfo.timingVec(:,:,1,1,1,1,maskTE);

% Preparation for saving intermediate results
npm=N*N;
nx=2;
ny=2;
costresult = zeros(1,numloop*numloop_ext);
r2result = zeros(N*N,numloop,numloop_ext);
weresult = zeros(N*N,numloop,numloop_ext);
gradr2_tot=zeros(N*N,numloop,numloop_ext);
gradwe_tot=zeros(N*N,numloop,numloop_ext);
gg_tot=zeros(N*N,numloop,numloop_ext);
imgresult=zeros(N*N,numloop_ext);


% Initialization
C = C2sparse('leak', ones(N,N), 8, 0);
C = sqrt(2^(4)) * C;
xinit=zeros(N*N,1);

sen = reshape(squeeze(sen),N*N,num_coils);
r2temp = zeros(npm,numloop+1,numloop_ext+1);
wetemp = zeros(npm,numloop+1,numloop_ext+1);
r2temp(:,1,1) = col(r2map);
wemap = FM;
wetemp(:,1,1) = col(wemap);

mskpenalwe = logical(ones(N,N));
mskpenalr2 = logical(mask);
powbetar2 = -15; 
%powbetar2 = -4; % AMC
betr2 = 2^powbetar2;
powbetawe = -30; 
%powbetawe = -6; % AMC
betwe = 2^powbetawe;
betS0 = 2^-4;
Rr2 = Robject(mskpenalr2,'edge_type','tight','order',2,'beta',betr2,'potential','quad','delta',0.1);
Rwe = Robject(mskpenalwe,'edge_type','tight','order',2,'beta',betwe,'potential','quad','delta',0.1);
RS0 = Robject(mskpenalwe,'edge_type','tight','order',2,'beta',betS0,'potential','quad','delta',0.1);

%Rr2 = Robject(mskpenalr2,'edge_type','tight','order',2,'beta',betr2,'potential','huber','delta',0.1);
%Rwe = Robject(mskpenalwe,'edge_type','tight','order',2,'beta',betwe,'potential','huber','delta',0.1);
A = cell(numechoes,1);

data = squeeze(data);
data_t0 = squeeze(data(:,:,1,:));
data = col(data);


% ttp=[];
% for kk = 1:numechoes
%     ttp = [ttp;tt+TE(kk)];
% end


for ii=1:numloop_ext
    
    A_t0   = fast_mr_wt_quad(col(rInfo.kx(:,:,1,1,1,1,1)),col(rInfo.ky(:,:,1,1,1,1,1)),FOV,N,2*N,J,col(ttp(:,:,1,1,1,1,1)),han_flag,col(wemap),col(r2map),mask,1,L,1,[],0,nx,ny,1,0);
    A_t0s  = sense_svd(A_t0,sen);
    data0s = prepData(A_t0s,col(data_t0));
    
    if ii ==1
        im_rec = qpwls_pcg(xinit(:),A_t0s,1,col(data0s),0,C,1,45,[],0);
    else
        im_rec = qpwls_pcg(xinit(:),A_t0s,1,col(data0s),0,C,1,15,[],0);
    end
    
    xinit=im_rec(:,end);
    B=im_rec(:,end).*mask(:);
    
    imgresult(:,ii)=B;
    
    nc_svd = A_t0s.rank_svd;
    tt_ext=col(repmat(col(ttp),nc_svd,1));
    
    
    fobj = @(x) (R2StarFMCostFunction(x,data,sen,rInfo,Rr2,Rwe,maskTE,mask,nc_svd,tt_ext,B));
    options = optimoptions('fmincon','Algorithm','interior-point','Display','iter-detailed','SpecifyObjectiveGradient',true,'HessianApproximation','lbfgs','MaxIter',10,'MaxFunEvals',300,'TolFun',1E-3,'TolX',1E-3);
    %options = optimoptions('fmincon','Algorithm','interior-point','Display','iter-detailed','SpecifyObjectiveGradient',true,'Hessian','fin-diff-grads','SubproblemAlgorithm','cg','MaxIter',10,'MaxFunEvals',300,'TolFun',1E-3,'TolX',1E-3);
    xinitFmincon = vertcat(col(r2temp(:,1,ii)),col(wetemp(:,1,ii)));
    lBounds = vertcat(col(zeros(size(r2temp(:,1,ii)))),col(-10000*ones(size(wetemp(:,1,ii)))));
    uBounds = vertcat(col(500*ones(size(r2temp(:,1,ii)))),col(10000*ones(size(wetemp(:,1,ii)))));
    
    [x] = fmincon(fobj,xinitFmincon,[],[],[],[],lBounds,uBounds,[],options);
    xlen = length(x);
    r2temp(:,1,ii+1)=x(1:xlen/2);
    wetemp(:,1,ii+1)=x(xlen/2+1:end);
    
    r2result(:,1,ii+1)=x(1:xlen/2);
    weresult(:,1,ii+1)=x(xlen/2+1:end);
    s0result(:,1,ii+1)=B;
    
    
    r2map = reshape(r2result(:,1,ii+1),N,N);
    wemap = reshape(weresult(:,1,ii+1),N,N);
    s0map = reshape(s0result(:,1,ii+1),N,N);
    
%     for jj = 1:numloop
%         
%         powalp = 16;
%         
%         r2now = r2temp(:,jj,ii);
%         wenow = wetemp(:,jj,ii);
%         
%         for kk = 1 : numechoes
%             A{kk}= fast_mr_wt_quad(col(rInfo.kx(:,:,1,1,1,1,kk)),col(rInfo.ky(:,:,1,1,1,1,kk)),FOV,N,2*N,J,col(ttp(:,:,1,1,1,1,kk)),han_flag,wenow,r2now,mask,1,L,1,[],0,nx,ny,1,0);
%         end
%         
%         A_full = multi_echo_mri(A,numechoes);
%         if (jj == 1)
%             A_fulls=sense_svd(A_full,sen);
%             datas = prepData(A_fulls,col(data));
%         else
%             A_fulls.A = A_full;
%         end
%         
%         if (jj == 1) % first step: steepest descent
%             resid = datas-A_fulls*B;
%             gg = conj(B).*(A_fulls'*((tt_ext).*resid));
%             
%             gradr2 = real(gg(mskpenalr2))+Rr2.cgrad(Rr2,r2temp(mskpenalr2,jj,ii));
%             gradwe = (imag(gg)+Rwe.cgrad(Rwe,wetemp(:,jj,ii)));
%             roughr2 = Rr2.penal(Rr2,r2temp(mskpenalr2,jj,ii));
%             roughwe = Rwe.penal(Rwe,wetemp(:,jj,ii));
%             
%             cost_old = 1/2*(real(resid'*resid))+roughr2+roughwe
%             cost_new = cost_old;
%             
%             newdirr2=-gradr2;
%             newgradr2=gradr2;
%             newdirwe=-gradwe;
%             newgradwe=gradwe;
%         else % after: conjugate gradient
%             cost_old = cost_new;
%             gg = conj(B).*(A_fulls'*((tt_ext).*resid));
%             gradr2 = real(gg(mskpenalr2))+Rr2.cgrad(Rr2,r2temp(mskpenalr2,jj,ii));
%             gradwe= (imag(gg)+Rwe.cgrad(Rwe,wetemp(:,jj,ii)));
%             
%             newgradr2=gradr2;
%             newdirr2 = (-newgradr2 +(newgradr2'*newgradr2)/(oldgradr2'*oldgradr2)*olddirr2);
%             newgradwe=gradwe;
%             newdirwe = (-newgradwe +(newgradwe'*newgradwe)/(oldgradwe'*oldgradwe)*olddirwe);
%         end
%         
%         %        residtot(:,jj,ii)=resid;
%         grad=col(newgradr2);
%         dir=col(newdirr2);
%         
%         if dir' * grad >0
%             warning 'reset loop steepest descent'
%             newdirr2=-gradr2;
%             newgradr2=gradr2;
%             newdirwe=-gradwe;
%             newgradwe=gradwe;
%         end
%         
%         % Find the step size
%         
%         while ((cost_new>=cost_old)|(isnan(cost_new)))
%             if (powalp <= -10)
%                 sprintf('cannot decrease cost function')
%                 break % keyboard
%             end
%             
%             powalp = powalp -1;
%             alp = 2^powalp;
%             r2temp(mskpenalr2,jj+1,ii) = r2temp(mskpenalr2,jj,ii)+alp*(newdirr2);
%             wetemp(:,jj+1,ii) = wetemp(:,jj,ii)+alp*(newdirwe);
%             
%             if 1 % update the model with the new R2 map
%                 r2_tmp=r2temp(:,jj+1,ii);
%                 we_tmp=wetemp(:,jj+1,ii);
%                 A = cell(numechoes,1);
%                 for kk = 1 : numechoes
%                     A{kk}= fast_mr_wt_quad(col(rInfo.kx(:,:,1,1,1,1,kk)),col(rInfo.ky(:,:,1,1,1,1,kk)),FOV,N,2*N,J,col(ttp(:,:,1,1,1,1,kk)),han_flag,we_tmp,r2_tmp,mask,1,L,1,[],0,nx,ny,1,0);
%                 end
%                 A_full = multi_echo_mri(A,numechoes);
%                 %A_fulls=sense_svd(A_full,sen,'CoilRank',3);
%                 A_fulls.A = A_full;
%             end
%             resid = datas-A_fulls*B;
%             roughr2 = Rr2.penal(Rr2,r2temp(mskpenalr2,jj+1,ii))
%             roughwe = Rwe.penal(Rwe,wetemp(:,jj+1,ii))
%             cost_new = 1/2*(real(resid'*resid))+roughr2+roughwe
%         end % end while
%         im(reshape(wenow,128,128)); drawnow
%         powalp
%         if (powalp <= -10)
%             sprintf('cannot decrease cost function')
%             break % keyboard
%         end
%         
%         oldgradr2=newgradr2;
%         olddirr2=newdirr2;
%         oldgradwe=newgradwe;
%         olddirwe=newdirwe;
%         
%         r2map = r2temp(:,jj+1,ii);
%         wemap = wetemp(:,jj+1,ii);
%         gg_tot(:,jj,ii)=gg;
%         gradr2_tot(mskpenalr2,jj,ii)=alp*newdirr2;
%         gradwe_tot(:,jj,ii)=alp*newdirwe;
%         r2result(:,jj,ii) = r2map;
%         weresult(:,jj,ii) = wemap;
%         costresult(:,jj+(ii-1)*numloop) = cost_new;
%         
%         % save intermediate results
%         save intermediate_res r2map wemap B
%     end % end step size
%   r2temp(:,1,ii+1)=r2temp(:,jj+1,ii);
%   wetemp(:,1,ii+1)=wetemp(:,jj+1,ii);
end

% Prepare results for returning to calling function
r2=reshape(r2result,N,N,numloop,numloop_ext);
we=reshape(weresult,N,N,numloop,numloop_ext);
if (ReturnComplete == true)
    r2final = r2;
    wefinal = we;
    s0final = B;
else
    r2final = r2(:,:,1,end);
    wefinal = we(:,:,1,end);
    s0final = B;
end
% save results
if 0
    name_folder = sprintf('R2_FM_Slice_%02d',slice);
    mkdir(name_folder)
    cd(name_folder)
    save costresult costresult
    save gradr2_tot gradr2_tot
    save gradwe_tot gradwe_tot
    save imgresult imgresult
    %save residtot residtot
    r2=reshape(r2result,N,N,numloop,numloop_ext);
    we=reshape(weresult,N,N,numloop,numloop_ext);
    save r2 r2
    save we we
    save powbetawe powbetawe
    save powbetar2 powbetar2
    save TE TE
    
end
end
