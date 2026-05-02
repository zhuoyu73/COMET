function status = spiral_recon_run(recoInfo)
%function status = spiralraw(recoInfo)
%   recoInfo from setup_recon.m
  

sense_reco_list = [7;8];
FM_reco_list = [3;4;5;6;8];

flag_sense = ismember(recoInfo.rec_type,sense_reco_list);
flag_FM = ismember(recoInfo.rec_type,FM_reco_list);

if flag_sense
   if exist('sen.mat','file')
       load sen
   elseif exist('../sen.mat','file')
       load ../sen
   else
       sprintf('Need sen.mat in top level directory for sense reconstruction \n')
       return
   end

  if ~(size(sen,1) == recoInfo.N)
      sen_new = zeros(recoInfo.N,recoInfo.N,recoInfo.nsl,recoInfo.num_coils);
     %RESAMPLE SEN for matrix size of current reconstruction
      for coilIndex = 1:recoInfo.num_coils
          sen_new(:,:,:,coilIndex) = resample_map_resolution(sen(:,:,:,coilIndex),recoInfo.N,recoInfo.nsl);   
      end
      sen=sen_new;
      clear sen_new
      save sen_interp sen
      %keyboard
   end

end

if flag_FM
   if exist('FM.mat','file')
      load FM
   elseif exist('../FM.mat','file')
      load ../FM
   else
      sprintf('Need FM.mat in top level directory for field corrected reconstruction \n')
      return
   end
   if ~(size(FM,1) == recoInfo.N)
       %RESAMPLE FM FOR matrix size of current reconstruction
       FM_new = resample_map_resolution(FM,recoInfo.N,recoInfo.nsl);
       FM = FM_new;
       clear FM_new
       save FM_interp FM
   end
 
end
%keyboard
kx = reshape(recoInfo.kx,[recoInfo.nro recoInfo.nl]);
ky = reshape(recoInfo.ky,[recoInfo.nro recoInfo.nl]);
nro = recoInfo.nro;
nroa = recoInfo.nroa;


skp = 128;


fid = fopen(recoInfo.fname_data);
bytes_to_skip_hdr = fread(fid,1,'uint32');
fseek(fid,bytes_to_skip_hdr,'bof');

bytesPerImage = recoInfo.num_coils*(((skp/2) * 2) + ((2*recoInfo.nroa) * 4)); 

%%%checking file size
%%fsz = (nroa*2*4+128)*nl*nsl*ntp*num_echo;


mkdir('Recon')
cd('Recon')


tmpvol = zeros(recoInfo.N,recoInfo.N,recoInfo.nsl);
dat = zeros(recoInfo.nroa*recoInfo.num_adc,recoInfo.nl,recoInfo.nsl,recoInfo.num_coils); 
dat = zeros(recoInfo.nro,recoInfo.nl,recoInfo.nsl,recoInfo.num_coils); 


for imageIndex = 1:(recoInfo.ntp)
   for shotIndex = 1:recoInfo.nl %% load each interleaf's data into the buffer
        for sliceIndex = 1:recoInfo.nsl
            for adcIndex = 1:recoInfo.num_adc
                for coilIndex = 1:recoInfo.num_coils %% read in each coil data
 
                    %% Read and parse interleaf points
                   sMDH = ice_read_mdh_va21(fid);
                    [raw,rcount] = fread(fid,2*nroa,'float32');

                    if (rcount*recoInfo.num_adc < 2*nro)
                        sprintf('Reached end of file before running out of indices \n')
                         return
                    end

            %% Put it in the correct buffer position
                tmp = raw(1:2:2*nroa)+i*raw(2:2:2*nroa);
                start_ind = (adcIndex-1)*nroa + 1;
                end_ind = min(adcIndex*nroa,nro);
                num_ind = end_ind - start_ind +1;
                tmp2 = col(tmp(1:num_ind)).*exp(-i*2*pi*ky(start_ind:end_ind,shotIndex)*recoInfo.shft_FOVy).*exp(-i*2*pi*kx(start_ind:end_ind,shotIndex)*recoInfo.shft_FOVx);
                dat(start_ind:end_ind,shotIndex,sliceIndex,coilIndex) = tmp2; 
                   
                end % end coilIndex
            end  % end adcIndex
        end  % end sliceIndex
    end  % end shotIndex


% Call correct recon and save output

    switch recoInfo.rec_type
%--------------------------------------------------------------
        case {1,2} %gridding, individual coil images or sos (k2image)
%--------------------------------------------------------------
           for sliceIndex = 1:recoInfo.nsl
              for coilIndex = 1:recoInfo.num_coils
                   imgtmp(:,:,sliceIndex,coilIndex) = k2image(col(kx),col(ky),col(dat(:,:,sliceIndex,coilIndex)),col(recoInfo.ww),recoInfo.N,2.5); %uncorrected for inhomogeneity
          	end
            end
           
            if (recoInfo.rec_type == 1)
                for coilIndex = 1:recoInfo.num_coils
                  img = squeeze(imgtmp(:,:,:,coilIndex));
                  save(sprintf('c%02d_%05d',coilIndex,imageIndex),'img')
                end
            else
                    img = sos_image(imgtmp);
                    save(sprintf('sos_%05d',imageIndex),'img')
            end

 
           
%--------------------------------------------------------------
        case {3,4} %gridding with field map (FM), individual coils or sos	(k2image_we)
%--------------------------------------------------------------
           for sliceIndex = 1:recoInfo.nsl
              for coilIndex = 1:recoInfo.num_coils
                   imgtmp(:,:,sliceIndex,coilIndex) = k2image_we(col(kx),col(ky),col(dat(:,:,sliceIndex,coilIndex)),col(recoInfo.ww),recoInfo.N,col(FM(:,:,sliceIndex)),recoInfo.L, col(recoInfo.tt),2.5); %corrected for inhomogeneity
          	end
            end
           
           if (recoInfo.rec_type == 3)
               for coilIndex = 1:recoInfo.num_coils
                  img = squeeze(imgtmp(:,:,:,coilIndex));
                  save(sprintf('c%02d_%05d',coilIndex,imageIndex),'img')
               end
           else
               img = sos_image(imgtmp);
               save(sprintf('sos_%05d',imageIndex),'img')
           end

%--------------------------------------------------------------
        case {5,6} %iterative with field map (FM), individual coils  or sos  	(fast_mr)
%--------------------------------------------------------------

             
              xinit = col(zeros(recoInfo.N));
              mask = ones(recoInfo.N);
              for sliceIndex = 1:recoInfo.nsl
                  A = fast_mr(col(kx), col(ky), recoInfo.FOV/10, recoInfo.N, 2*recoInfo.N, recoInfo.J, col(recoInfo.tt), col(FM(:,:,sliceIndex)), 0, recoInfo.L, 1);
             
                 for coilIndex = 1:num_coils
                     recon = qpwls_pcg(xinit, A, 1, col(dat(:,:,sliceIndex,coilIndex)), 0, recoInfo.C_regularize, 1, recoInfo.num_iters, mask, 0);
                     imgtmp(:,:,sliceIndex,coilIndex) = reshape(recon(:,end), recoInfo.N, recoInfo.N);
                  end 
              end   

          if (recoInfo.rec_type == 5)
               for coilIndex = 1:recoInfo.num_coils
                  img = squeeze(imgtmp(:,:,:,coilIndex));
                  save(sprintf('c%02d_%05d',coilIndex,imageIndex),'img')
               end
           else
               img = sos_image(imgtmp);
               save(sprintf('sos_%05d',imageIndex),'img')
           end

%--------------------------------------------------------------
        case {7,8} %iterative with SENSE, and SENSE and FM		(fast_mr)
%--------------------------------------------------------------
               if (recoInfo.rec_type == 7)
                     FM = zeros(recoInfo.N,recoInfo.N,recoInfo.nsl);
               end
               tt = reshape(recoInfo.tt,[length(recoInfo.tt)/recoInfo.nl recoInfo.nl]);               
             
               shot_samp_fact = recoInfo.nl/recoInfo.shots_to_use;

              sen = reshape(sen, [recoInfo.N*recoInfo.N recoInfo.nsl recoInfo.num_coils]);

              xinit = col(zeros(recoInfo.N));
              mask = ones(recoInfo.N);
              for sliceIndex = 1:recoInfo.nsl
                   for shotIndex = 1:shot_samp_fact
                        start_shot = (shotIndex-1)*recoInfo.shots_to_use+1;
                        end_shot = shotIndex*recoInfo.shots_to_use;
                        A = fast_mr(col(kx(:,start_shot:end_shot)), col(ky(:,start_shot:end_shot)), recoInfo.FOV/10, recoInfo.N, 2*recoInfo.N, recoInfo.J, col(tt(:,start_shot:end_shot)), col(FM(:,:,sliceIndex)), 0, recoInfo.L, 1);

                        As = sense(A,squeeze(sen(:,sliceIndex,:)));
                        recon = qpwls_pcg(xinit, As, 1, col(squeeze(dat(:,start_shot:end_shot,sliceIndex,:))), 0, recoInfo.C_regularize, 1, recoInfo.num_iters, mask, 0);
                        imgtmp(:,:,sliceIndex,shotIndex) = reshape(recon(:,end), recoInfo.N, recoInfo.N);
                  end 
              end   

               for shotIndex = 1:shot_samp_fact
                   img = squeeze(imgtmp(:,:,:,shotIndex));
                   save(sprintf('sen_%05d',shot_samp_fact*(imageIndex-1)+shotIndex),'img')
               end
         



%--------------------------------------------------------------
        otherwise
           sprintf('Recon Type %d not implemented yet \n',recoInfo.rec_type)
           return
    end
 



end % end imageIndex

cd ..

fclose('all')



