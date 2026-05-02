function [read_shift phase_shift slice_shift] = image_shift( ascconv,bAllSlices)
%code that gets the the shifts to be used in the reconstruction of spiral
%images
% Original code from siguhani2
% Edited by JHoltrop December 4, 2012

% inputs
% ascconv: sequncce data created from parsasc
% bAllSlice: Set this value to return a seperate value for each slab(3D)/slice(2D)
%
% outputs:
% read_shift: shift in fraction of the FOV in the read direction
% phase_shift: shift in fraction of the FOV in the phase direction
% slice_shift: shift in fraction of the FOV in the slice direction
%2014_10_03: JLH - Give a different slice shift for each slice
%2014_10_03: JLH - merged the 2d and 3d versions and added the bAllSlices
%flag

warning('MRFIL: This function (%s) is deprecated. Do not use this function for new development.',mfilename);

if nargin == 1
    bAllSlices = 0;
end

if bAllSlices
    nsl = ascconv.sGroupArray.asGroup.nSize;
else
    nsl = 1;
end

read_shift = zeros(nsl,1);
phase_shift = zeros(nsl,1);
slice_shift = zeros(nsl,1);

for nn = 1:nsl
    if isfield(ascconv.sSliceArray.asSlice(nn), 'sPosition')
        if isfield(ascconv.sSliceArray.asSlice(nn).sPosition, 'dSag')
            offset_X = ascconv.sSliceArray.asSlice(nn).sPosition.dSag;
        else
            offset_X = 0;
        end

        if isfield(ascconv.sSliceArray.asSlice(nn).sPosition, 'dCor')
            offset_Y = ascconv.sSliceArray.asSlice(nn).sPosition.dCor;
        else
            offset_Y = 0;
        end

        if isfield(ascconv.sSliceArray.asSlice(nn).sPosition, 'dTra')
            offset_Z = ascconv.sSliceArray.asSlice(nn).sPosition.dTra;
        else
            offset_Z = 0;
        end
    else
        offset_X = 0;
        offset_Y = 0;
        offset_Z = 0;
    end


    if isfield(ascconv.sSliceArray.asSlice(nn).sNormal, 'dSag')
        dc_X = ascconv.sSliceArray.asSlice(nn).sNormal.dSag;
    else
        dc_X = 0;
    end

    if isfield(ascconv.sSliceArray.asSlice(nn).sNormal, 'dCor')
        dc_Y = ascconv.sSliceArray.asSlice(nn).sNormal.dCor;
    else
        dc_Y = 0;
    end

    if isfield(ascconv.sSliceArray.asSlice(nn).sNormal, 'dTra')
        dc_Z = ascconv.sSliceArray.asSlice(nn).sNormal.dTra;
    else
        dc_Z = 0;
    end

    if isfield(ascconv.sSliceArray.asSlice(nn), 'dInPlaneRot')
        InPlaneRot = ascconv.sSliceArray.asSlice(nn).dInPlaneRot;
    else
        InPlaneRot = 0;
    end

    offset_vec = [offset_X;offset_Y;offset_Z];
    norm_vec = [dc_X;dc_Y;dc_Z];
    rot = InPlaneRot;

    maxnrm = max(abs(norm_vec(:)));
    orient_ind = find(abs(norm_vec) == maxnrm);
    % if (orient_ind == 3)
    %    orient_str = 'axial';
    % elseif (orient_ind == 2)
    %     orient_str = 'coronal';
    % elseif (orient_ind == 1)
    %     orient_str = 'sagittal';
    % else
    %     sprintf('Did not find max')
    %     keyboard
    % end

    % shift_vec = offset_vec - (offset_vec'*norm_vec)*norm_vec
    x_vec = [1;0;0];
    y_vec = [0;1;0];
    z_vec = [0;0;1];


    % Phase and Read directions for each orientation determined by cross products

    if (orient_ind == 3) % axial

        if (dc_Z > 0)
            phase_ornt = cross(norm_vec,x_vec);
        elseif (dc_Z < 0)
            phase_ornt = cross(x_vec,norm_vec);
        end

        read_ornt = cross(-y_vec,norm_vec);

    elseif (orient_ind == 2)  % coronal

        if (dc_Y > 0)
            read_ornt = cross(norm_vec,x_vec);
        elseif (dc_Y < 0)
            read_ornt = cross(norm_vec,-x_vec);
        end

        phase_ornt = cross(norm_vec,z_vec);

    elseif (orient_ind == 1)  % sagittal

        if (dc_X > 0)
            read_ornt = cross(-y_vec,norm_vec);
        elseif (dc_X < 0)
            read_ornt = cross(y_vec,norm_vec);
        end

        phase_ornt = cross(z_vec,norm_vec);

    end

    read_vec = read_ornt/norm(read_ornt)    ;    % Make read and phase vectors unit vectors
    phase_vec = phase_ornt/norm(phase_ornt);

    read_in = dot(offset_vec, read_vec);
    phase_in = dot(offset_vec, phase_vec);

    read_shift_tmp = read_in*cos(rot) + phase_in*sin(rot);
    phase_shift_tmp = -read_in*sin(rot) + phase_in*cos(rot);
    slcdir_shift_tmp = dot(offset_vec, norm_vec);

    % lFOV = ascconv.sSliceArray.asSlice(1).dReadoutFOV;

    %Divide by the FOV
    phase_shift(nn) = phase_shift_tmp/ascconv.sSliceArray.asSlice(nn).dPhaseFOV;
    read_shift(nn) = read_shift_tmp/ascconv.sSliceArray.asSlice(nn).dReadoutFOV;
    slice_shift(nn) = slcdir_shift_tmp/(ascconv.sSliceArray.asSlice(nn).dThickness);

end

end
