function [gradients] = readDVSFile(filename)
% readDVSFile Reads Siemens diffusion vector set (.dvs) files for VD/VE 
% Scanners

fid = fopen(filename,'r');
gradients = [];
% Skip commented characters at the beginning
while(1)
    line = fgetl(fid);
    if line(1) ~= '#'
        break;
    end
end

line = fgetl(fid);
while ischar(line)
    if length(line) < 6
        
    elseif strcmpi(line(1:6),'vector')
        
        if line(1) == 'V'
            res = textscan(line,'Vector[%u] = ( %f, %f, %f )');
        else
             res = textscan(line,'vector[%u] = ( %f, %f, %f )');
        end
        gradients = vertcat(gradients,[res{2}, res{3}, res{4}])
    end
    
    line = fgetl(fid);

end

end

