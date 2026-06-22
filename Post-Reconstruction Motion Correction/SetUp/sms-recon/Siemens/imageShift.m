function [read_shift, phase_shift, slice_shift, rotMatrix] = imageShift(rInfo,bAllSlices)
%code that gets the the shifts to be used in the reconstruction of spiral
%images
% Original code from siguhani2
% Edited by JHoltrop December 4, 2012 
% Updated by Alex Cerjanic to support mapVBVD for VB/VD/VE support

% inputs
% rInfo: reconInfo object initialized with mapVBVD object
% bAllSlice: Set this value to return a separate value for each slab(3D)/slice(2D)
%
% outputs:
% read_shift: shif t in fraction of the FOV in the read direction
% phase_shift: shift in fraction of the FOV in the phase direction
% slice_shift: shift in fraction of the FOV in the slice direction
% rotMatrix: rotation matrix from gradient system to logical system
%2014_10_03: JLH - Give a different slice shift for each slice
%2014_10_03: JLH - merged the 2d and 3d versions and added the bAllSlices
%flag

if nargin == 1
    bAllSlices = 0;
end

if bAllSlices
    nsl = rInfo.DataObject.hdr.MeasYaps.sGroupArray.asGroup{1}.nSize;
else
    nsl = 1;
end

read_shift = zeros(nsl,1);
phase_shift = zeros(nsl,1);
slice_shift = zeros(nsl,1);

for nn = 1:nsl
    
    if isfield(rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}, 'sPosition')
        if isfield(rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.sPosition, 'dSag')
            offset_Sag = rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.sPosition.dSag;
        else
            offset_Sag = 0;
        end

        if isfield(rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.sPosition, 'dCor')
            offset_Cor = rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.sPosition.dCor;
        else
            offset_Cor = 0;
        end

        if isfield(rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.sPosition, 'dTra')
            offset_Tra = rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.sPosition.dTra;
        else
            offset_Tra = 0;
        end
    else
        offset_Sag = 0;
        offset_Cor = 0;
        offset_Tra = 0;
    end


    if isfield(rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.sNormal, 'dSag')
        dc_Sag = rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.sNormal.dSag;
    else
        dc_Sag = 0;
    end

    if isfield(rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.sNormal, 'dCor')
        dc_Cor = rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.sNormal.dCor;
    else
        dc_Cor = 0;
    end

    if isfield(rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.sNormal, 'dTra')
        dc_Tra = rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.sNormal.dTra;
    else
        dc_Tra = 0;
    end

    if isfield(rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}, 'dInPlaneRot')
        InPlaneRot = rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.dInPlaneRot;
    else
        InPlaneRot = 0;
    end

    offset_vec = [offset_Sag;offset_Cor;offset_Tra];
    norm_vec = [dc_Sag;dc_Cor;dc_Tra];
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

        if (dc_Tra > 0)
            phase_ornt = cross(norm_vec,x_vec);
        elseif (dc_Tra < 0)
            phase_ornt = cross(x_vec,norm_vec);
        end

        read_ornt = cross(-y_vec,norm_vec);

    elseif (orient_ind == 2)  % coronal

        if (dc_Cor > 0)
            read_ornt = cross(norm_vec,x_vec);
        elseif (dc_Cor < 0)
            read_ornt = cross(norm_vec,-x_vec);
        end

        phase_ornt = cross(norm_vec,z_vec);

    elseif (orient_ind == 1)  % sagittal

        if (dc_Sag > 0)
            read_ornt = cross(-y_vec,norm_vec);
        elseif (dc_Sag < 0)
            read_ornt = cross(y_vec,norm_vec);
        end

        phase_ornt = cross(z_vec,norm_vec);

    end

    read_vec = read_ornt/norm(read_ornt)    ;    % Make read and phase vectors unit vectors
    phase_vec = phase_ornt/norm(phase_ornt);

    read_in = dot(offset_vec, read_vec);
    phase_in = dot(offset_vec, phase_vec);

    quaternion = reshape(rInfo.DataObject.image.slicePos(4:7, 1),[1,4]);
    SCTRotMatrix = quaternion2rotm(quaternion);
    
%     read_shift_tmp = read_in*cos(rot) + phase_in*sin(rot);
%     phase_shift_tmp = -read_in*sin(rot) + phase_in*cos(rot);
%     slcdir_shift_tmp = dot(offset_vec, norm_vec);
    PRSShiftmm = SCTRotMatrix'*offset_vec;
    
    % lFOV = rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice(1).dReadoutFOV;

    %rInfo.DataObject.image.slicePos(1, nn)

    %Divide by the FOV
    %phase_shift(nn) = rInfo.DataObject.image.slicePos(1, nn)/rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.dPhaseFOV;
    %read_shift(nn) = rInfo.DataObject.image.slicePos(2, nn)/rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.dReadoutFOV;
    %slice_shift(nn) = rInfo.DataObject.image.slicePos(3, nn)/(rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.dThickness);
    phase_shift(nn) = PRSShiftmm(1)/rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.dPhaseFOV;
    read_shift(nn) = PRSShiftmm(2)/rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.dReadoutFOV;
    slice_shift(nn) = PRSShiftmm(3)/(rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{nn}.dThickness);

end
    % Identifying rotation matrix correctly! - AMC

    DicomOrientation = rInfo.DataObject.hdr.Dicom.DICOM.tPatientPosition;
    XYZtoSCTMatrix = zeros(3,3); %Matrix to hold translation from [Phase,Read,Slice] to [Sag,Cor,Tra]
    
    if strcmp('HFS',DicomOrientation)
        XYZtoSCTMatrix = [1,0,0;0,-1,0;0,0,-1];
    elseif strcmp('HFP',DicomOrientation)
        XYZtoSCTMatrix = [-1,0,0;0,1,0;0,0,-1];
    elseif strcmp('HFDR',DicomOrientation)
        XYZtoSCTMatrix = [1,0,0;0,1,0;0,0,-1];
    elseif strcmp('HFDL',DicomOrientation)
        XYZtoSCTMatrix = [0,-1,0;-1,0,0;0,0,-1];
    elseif strcmp('FFS',DicomOrientation)
        XYZtoSCTMatrix = [-1,0,0;0,-1,0;0,0,-1];
    elseif strcmp('FFP', DicomOrientation)
        XYZtoSCTMatrix = [1,0,0;0,1,0;0,0,1];        
    elseif strcmp('FFDR',DicomOrientation)
        XYZtoSCTMatrix = [0,-1,0;1,0,0;0,0,1];      
    elseif strcmp('FFDL',DicomOrientation)
        XYZtoSCTMatrix = [0,1,0;-1,0,0;0,0,1];
    else
        error('Unrecognized DICOM Patient Positioning String Included');
    end
        
    % Complete rotation matrix requires the post multiplication of the
    % PRStoSCTMatrix by the rotation matrix from the quaternion
    quaternion = reshape(rInfo.DataObject.image.slicePos(4:7, 1),[1,4]);
    SCTRotMatrix = quaternion2rotm(quaternion);
    sliceInfo = rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal; 
    rotVector = zeros(3,1);
    if ~isfield(sliceInfo,'dSag')
        rotVector(1) = 0;
    else
        rotVector(1) = sliceInfo.dSag;
    end
    
    if ~isfield(sliceInfo,'dCor')
        rotVector(2) = 0;
    else
        rotVector(2) = sliceInfo.dCor;
    end
    
    if ~isfield(sliceInfo,'dTra')
        rotVector(3) = 0;
    else
        rotVector(3) = sliceInfo.dTra;
    end
    
    sliceInfo = rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{1};
    
%     if ~isfield(sliceInfo,'dInPlaneRot')
%         rotAngle = - pi/2;
%     else
%         rotAngle = sliceInfo.dInPlaneRot - pi/2;
%     end
    
    %rotAngle = rInfo.DataObject.hdr.MeasYaps.sSliceArray.asSlice{1}.dInPlaneRot;
%     SCTRotMatrix = rotationVectorToMatrix(rotAngle*rotVector);
    rotMatrix = XYZtoSCTMatrix*SCTRotMatrix;
    %rotMatrix = PRSRotMatrix;

end

