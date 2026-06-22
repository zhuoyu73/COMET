function status = spiral_recon_run(recoInfo)
%function status = spiralraw(recoInfo)
%   recoInfo from setup_recon.m
  


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
        otherwise
           sprintf('Recon Type %d not implemented yet \n',recoInfo.rec_type)
           return
    end
 



end % end imageIndex

cd ..

fclose('all')



