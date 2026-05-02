function ascconv = parsasc_full(filename)
%
%  function ascout = parsasc_full(fname)
%  filename is file name of meas.dat
%
%  written by Scott Hoge, 2011_04_25


fid = fopen(filename,'r');
bytes_to_skip_hdr = fread(fid,1,'uint32');
asc = fread(fid,bytes_to_skip_hdr,'char');
asc = char(asc');

ascconv = evalascconv( asc );
save ascconv ascconv


fclose(fid);

function ASCCONV = evalascconv(asc)
%
% this function converts variables expressed in C-code style in to a form
% that Matlab can interpreted and store them as variables.
%

t1 = strfind(asc,'### ASCCONV BEGIN ###');
t2 = strfind(asc,'### ASCCONV END ###');

% extract the ASCCONV section, append each line w/ a semi-colon, and 
% add in a top variable name (here, 'ASCONV') that will be returned.
tmp = regexprep( asc(t1(1):t2(1)-2), '\n', ';\nASCCONV.' ); 
tmp = [ tmp ';' ];
tmp = regexprep( tmp,' # '' ''',''); % Fixes the error in ucWaterSat
tmp = regexprep( tmp,'\[(\d+)\]', '[$1+1]' ); % add '1' to each specified array index,
tmp = regexprep( tmp,'[','(');                % convert C-style array braces
tmp = regexprep( tmp,']',')');                % ... to Matlab-style
tmp = regexprep( tmp,'#','%');                % convert comment lines
tmp = regexprep( tmp,'"','''');               % convert C-strings to Matlab-strings
tmp = regexprep( tmp,'0x([a-f\d]+);','hex2dec(''$1'');'); % convert all 0xDD numbers to hex2dec

tmp = regexprep( tmp,'0x20 % '' ', '3');    %added by ryan and mahshid, 3/4/2014 to read in data with water saturation
tmp = regexprep( tmp,'0x40 % ''@''', '4');    %added by Joe to read in data with no RF pulse

eval(tmp);                              % eval the string as a list of variable declarations

return;

