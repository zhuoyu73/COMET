function writeDVSFile(filename,u_x,u_y,u_z)
fid = fopen(filename,'w');

fprintf(fid,'[directions=%i]\n',length(u_x));
fprintf(fid,'CoordinateSystem = xyz\n');
fprintf(fid,'Normalisation = none\n');

for ii = 1:length(u_x)
    fprintf(fid,'Vector[%i] = ( %1.6f, %1.6f, %1.6f )\n',ii-1,u_x(ii),u_y(ii),u_z(ii));
end
fclose(fid);
