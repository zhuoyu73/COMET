function stat = parsasc ( filename )
% stat = parsasc ( filename )
% A MATLAB replacement for the parsasc script
% Reads a meas.asc file and produces seqinfo.m;  a pre-requisite for spiralraw_* scripts

if findstr(filename,'.dat')
    input_fd = fopen(filename,'rt');
    size_header = fread(input_fd,1,'uint32')
    flag_vb = 1;
elseif findstr(filename,'.asc')
    %input_fd = fopen([filename '.asc'],'rt');
    input_fd = fopen(filename,'rt');
    fseek(input_fd, 0, 'eof');
    size_header = ftell(input_fd);
    fseek(input_fd, 0, 'bof');
    flag_vb = 0;
else
    fname = dir(sprintf('%s*',filename));
    ii = 1;
    while (fname(ii).isdir)
        ii = ii+1;
    end    
    stat = parsasc(fname(ii).name)
    %stat = parsasc([filename '.asc']);
    return
end


output_fd = fopen('seqinfo.m','wt');

keys = {'sSliceArray.asSlice\[0\].dReadoutFOV',
        'alTR\[0\]',
        'alTE\[0\]',
        'sKSpace.lBaseResolution',
        'sKSpace.lPhaseEncodingLines',
        'm_ushSlicePerConcat',
        'sSliceArray.lConc',
        'lRepetitions',
        'm_flReadoutOSFactor',
        'adFlipAngleDegree\[0\]',
        'm_lEchoSpacing',
        'sSliceArray.asSlice\[0\].dThickness',
        'm_tSequenceVariant'};
names = {'fov',
         'TR',
         'TE',
         'nx',
         'ny',
         'nslpc',
         'nconc',
         'ntp',
         'ro_os_fact',
         'FAdeg',
         'Echospus',
         'slthick',
         'sequenceidentity'};
numFound = 0;

while (numFound < length(keys)) %% Stop looking if you've found all the keys

	aLine = fgetl(input_fd);  
    if (ftell(input_fd) > size_header)
        		sprintf('Hit end of file looking for a keyword')
                numFound = numFound+1;
                fseek(input_fd, 0, 'bof'); % rewind
     
    elseif (aLine == -1) 
		sprintf('Hit end of file looking for a keyword')
                numFound = numFound+1;
                fseek(input_fd, 0, 'bof'); % rewind
    else

       indexA = regexp(aLine, [keys{numFound+1} '\s*= ']);
	
      if (indexA)
		
        indexB = regexp(aLine, '= ');
		
        if isempty(str2num(aLine(indexB+2:end))) %% If you can't parse the value into a number
            
            fprintf(output_fd, '%s = ''%s''\n', names{numFound+1}, aLine(indexB+2:end));
            
        else
            
            fprintf(output_fd, '%s = %s\n', names{numFound+1}, aLine(indexB+2:end));
            
        end
		
        %eval([names{numFound+1} ' = ' aLine(indexB+2:end) ';']);
		
        numFound = numFound + 1;
		fseek(input_fd, 0, 'bof'); % rewind
      end

   end

end


%% OTHER PARAMETERS CHANGED IN VB
if flag_vb
keys = {'ushSlicePerConcat',
        'Concatenations',
        'flReadoutOSFactor'    };
names = {'nslpc',
         'nconc',
         'ro_os_fact' };
   fseek(input_fd,0,'bof');
   
numFound = 0;

while (numFound < length(keys)) %% Stop looking if you've found all the keys
	aLine = fgetl(input_fd);  
    if (ftell(input_fd) > size_header)
                sprintf('Did not find extra VB parameter')
                numFound = numFound+1;
                fseek(input_fd, 0, 'bof'); % rewind
     elseif (aLine == -1) 
                sprintf('Did not find extra VB parameter')
                 numFound = numFound+1;
                fseek(input_fd, 0, 'bof'); % rewind
    else
        
       indexA = regexp(aLine, keys{numFound+1});
       %%%indexA = regexp(aLine, 'ReadoutOSFactor');
       if ~isempty(indexA)
            %keyboard
           indexB = regexp(aLine,'\d');
           if isempty(indexB)
               %keyboard
                aLine = fgetl(input_fd);
               aLine = fgetl(input_fd);
               indexC = regexp(aLine,'Precision');
               if ~isempty(indexC)
                   aLine = fgetl(input_fd);
               end     
               fprintf(output_fd, '%s = %d \n', names{numFound+1}, str2num(aLine));
               numFound = numFound + 1;
		       fseek(input_fd, 0, 'bof'); % rewind
             
       end
    end
  end
	
end

% NOW TO SEARCH FOR NUMBER OF COILS
fseek(input_fd,0,'bof')
while (ftell(input_fd) < size_header)
    aLine = fgetl(input_fd);
    indexA = regexp(aLine,'lRxChannelConnected');
    if ~isempty(indexA)
        indexB = regexp(aLine, '= ');
        if ~isempty(str2num(aLine(indexB+2:end))) %% If you can't parse the value into a number
            num_coils = str2num(aLine(indexB+2:end)); 
        end
    end
    
end
fprintf(output_fd, 'num_coils = %d \n', num_coils);



end
%% Export extra lines
fprintf(output_fd, 'if (exist(''nslpc'',''var'') & exist(''nconc'',''var'')) \n');
fprintf(output_fd, 'nsl = nslpc*nconc\n');
fprintf(output_fd,'end \n');

  fprintf(output_fd, 'if exist(''ntp'',''var'') \n');
  fprintf(output_fd, '      ntp = ntp+1 \n');
  fprintf(output_fd, 'end \n');


% NOTE FOLLOWING OVERRIDES SEARCHING FOR RFRES
fprintf(output_fd,'rf_res = 5000; \n');


  %% Close up
fclose(input_fd);
fclose(output_fd);

display('Done reading metadata from meas.asc: Wrote seqinfo.m');

stat = 1;
