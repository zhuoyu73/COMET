function FMestim = FM_jointestimation(rawData,FMInit,N,kx,ky,L,tt,TE,numechoes,numloop_ext,numloop)
% RawData : k-space data
% FMinit: initial estimate field map
% kx,ky,tt : k-space trajectories and timing vector (without TE)
% TE : echo time
% numechoes: number of echoes
% numloop_ext: number of loop outside
% numloop: number of loop inside

%Saving intermediate results
npm=N*N;
costresult = zeros(1,numloop*numloop_ext);
imgresult=zeros(N*N,numloop_ext);

% Initialize image
C = C2sparse('leak', ones(N,N), 8, 0);
C = sqrt(2^(8)) * C;
xinit=zeros(N*N,1);
wetemp = zeros(npm,numloop+1,numloop_ext+1);
wetemp(:,1,1) = col(FMInit);

mskpenalwe = logical(ones(N,N));
powbetawe = -13;
betwe = 2^powbetawe;
Rwe = Robject(mskpenalwe,'edge_type','tight','order',2,'beta',betwe,'potential','quad','delta',0.1);
R = Robj(true(N,N),'edge_type','tight','order',2,'beta',0,'type_denom','matlab','potential','quad');

for ii=1:numloop_ext
    G0 = NUFFT(kx,ky,zeros(size(kx)),N,N,1);
    A0 = TimeSegmentation(G0,squeeze(tt(:,1)+TE(1)),wetemp(:,1,ii),L); 
    im_rec= solve_pwls_pcg(col(xinit),A0, 1, col(rawData(:,1)), R, 'niter',10);
    xinit=im_rec(:,end);
    
    imgresult(:,ii)=im_rec(:,end);
    
    
    for jj = 1:numloop
        
       
        powalp = 0;
        
        A = cell(numechoes,1);
        for kk = 1 : numechoes
            A{kk} = TimeSegmentation(G0,squeeze(tt(:,1)+TE(kk)),wetemp(:,jj,ii),L); 
        end
        
        A_full = multi_echo_mri(A,numechoes);
        
        
        % Find the direction
        
        if (jj == 1) % first step: steepest descent
            resid = rawData(:)-A_full*im_rec(:,end);
            gg = conj(im_rec(:,end)).*(A_full'*((tt(:)).*resid));
            
            gradwe = (imag(gg)+Rwe.cgrad(Rwe,wetemp(:,jj,ii)));
            roughwe = Rwe.penal(Rwe,wetemp(:,jj,ii));
            
            cost_old = 1/2*(real(resid'*resid)) + roughwe
            cost_new = cost_old;
            
            newdirwe=-gradwe;
            newgradwe=gradwe;
            
        else % after: conjugate gradient
            cost_old = cost_new;
            gg = conj(im_rec(:,end)).*(A_full'*((tt(:)).*resid));
            gradwe= (imag(gg)+Rwe.cgrad(Rwe,wetemp(:,jj,ii)));
            
            newgradwe=gradwe;
            newdirwe = (-newgradwe +(newgradwe'*newgradwe)/(oldgradwe'*oldgradwe)*olddirwe);
            
        end
        
        
        grad=col(newgradwe);
        dir=col(newdirwe);
        
        if dir' * grad >0
            warning 'reset loop steepest descent'
            newdirwe=-gradwe;
            newgradwe=gradwe;
        end
        
        % Find the step size
        
        while ((cost_new>=cost_old)|(isnan(cost_new)))
            if (powalp <= -10)
                sprintf('cannot decrease cost function')
                break % keyboard
            end
            
            powalp = powalp -1;
            alp = 2^powalp;
            
            wetemp(:,jj+1,ii) = wetemp(:,jj,ii)+alp*(newdirwe);
            
            
            if 1 % update the model with the new R2 map

                A = cell(numechoes,1);
                for kk = 1 : numechoes
                    A{kk} = TimeSegmentation(G0,squeeze(tt(:,1)+TE(kk)),wetemp(:,jj+1,ii),L); 

                end
                A_full = multi_echo_mri(A,numechoes);
            end
            
            resid = rawData(:)-A_full*im_rec(:,end);
            roughwe= Rwe.penal(Rwe,wetemp(:,jj+1,ii));
            
            cost_new = 1/2*(real(resid'*resid))+roughwe
            
            
            
        end % end while
        
        
        if (powalp <= -10)
            sprintf('cannot decrease cost function')
            break % keyboard
        end
        
        oldgradwe=newgradwe;
        olddirwe=newdirwe;
        costresult(:,jj+(ii-1)*numloop) = cost_new;
        
        if 1
%             subplot(2,1,1), im(reshape(r2map,N,N),[0 40])
         im(reshape(wetemp(:,jj+1,ii),N,N))
%             %             subplot(3,1,3), im(reshape(B,N,N))
             drawnow
        end
    end % end step size
    
    wetemp(:,1,ii+1)=wetemp(:,jj+1,ii);
    
    
end


% name_folder=sprintf('Simu_%s',name);
% 
% 
% mkdir(name_folder)
% cd(name_folder)
% 
% % end
% 
% save costresult costresult
% save imgresult imgresult
% we=reshape(weresult,N,N,numloop,numloop_ext);
% save we we
% save powbetawe powbetawe
% save powbetar2 powbetar2
% 
% % save SNR SNR
% save A_full A_full
% save resid resid
% save B B
% save TE TE
FMestim = wetemp;

end

