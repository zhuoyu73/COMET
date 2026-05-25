function [b] = parsePGSEManifestXML(xmlFile)
%PARSESEQUENCEMANIFESTXML Summary of this function goes here
%   Detailed explanation goes her

xmlStruct = xml2struct(xmlFile);

TE = str2num(xmlStruct.SequenceManifest.SequenceParameters.Attributes.TE)/1E6;
TR = str2num(xmlStruct.SequenceManifest.SequenceParameters.Attributes.TR)/1E6;
GradientReversal = @(t,TE) Pulse(TE/2 - t) - Pulse(t - TE/2);

% Parse Slice Select Details
GSliSel = [str2num(xmlStruct.SequenceManifest.GradientPulse{1,1}.Attributes.Read); ...
    str2num(xmlStruct.SequenceManifest.GradientPulse{1,1}.Attributes.Phase); ...
    str2num(xmlStruct.SequenceManifest.GradientPulse{1,1}.Attributes.Slice);]./1E3;
Duration_SliSel = str2num(xmlStruct.SequenceManifest.GradientPulse{1,1}.Attributes.Duration)./1E6;
RiseTime_SliSel = str2num(xmlStruct.SequenceManifest.GradientPulse{1,1}.Attributes.RiseTime)./1E6;

tSliSelStart =  -(Duration_SliSel + RiseTime_SliSel)/2;

% Construct function for this gradient pulse
SliSel = @(t) Trap(t, GSliSel, tSliSelStart, Duration_SliSel, RiseTime_SliSel);

% Parse Slice Select Rewinder Details
GSliSelRew = [str2num(xmlStruct.SequenceManifest.GradientPulse{1,2}.Attributes.Read); ...
    str2num(xmlStruct.SequenceManifest.GradientPulse{1,2}.Attributes.Phase); ...
    str2num(xmlStruct.SequenceManifest.GradientPulse{1,2}.Attributes.Slice);]./1E3;
Duration_SliSelRew = str2num(xmlStruct.SequenceManifest.GradientPulse{1,2}.Attributes.Duration)./1E6;
RiseTime_SliSelRew = str2num(xmlStruct.SequenceManifest.GradientPulse{1,2}.Attributes.RiseTime)./1E6;

tSliSelRewStart = tSliSelStart + Duration_SliSel + RiseTime_SliSel;

% Construct function for this gradient pulse
SliSelRew = @(t) Trap(t, GSliSelRew, tSliSelRewStart, Duration_SliSelRew, RiseTime_SliSelRew);

TEFill1 = str2num(xmlStruct.SequenceManifest.Delay{1,1}.Attributes.DelayTime)./1E6;
TEFill2 = str2num(xmlStruct.SequenceManifest.Delay{1,2}.Attributes.DelayTime)./1E6;
% Process Diffusion Module
LengthOfDiffusionArray = str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.Attributes.Length);


DiffModuleStart = tSliSelRewStart + Duration_SliSelRew + RiseTime_SliSelRew + TEFill1;

for ii = 1:LengthOfDiffusionArray
    
    PreDiffFill(ii) = str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.Delay{1,1}.Attributes.Duration)./1E6;
    DiffFill1(ii) = str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.Delay{1,2}.Attributes.Duration)./1E6;
    DiffFill2(ii) = str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.Delay{1,3}.Attributes.Duration)./1E6;
    PostDiffFill(ii) = str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.Delay{1,4}.Attributes.Duration)./1E6;
    
    GDiff1(:,ii) = [str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.GradientPulse{1,1}.Attributes.Read); ...
        str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.GradientPulse{1,1}.Attributes.Phase); ...
        str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.GradientPulse{1,1}.Attributes.Slice);]./1E3;
    Duration_Diff1(ii) = str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,1}.GradientPulse{1,1}.Attributes.Duration)./1E6;
    RiseTime_Diff1(ii) = str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,1}.GradientPulse{1,1}.Attributes.RiseTime)./1E6;
    
    GSliSel2(:,ii) = [str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,1}.GradientPulse{1,2}.Attributes.Read); ...
        str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.GradientPulse{1,2}.Attributes.Phase); ...
        str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.GradientPulse{1,2}.Attributes.Slice);]./1E3;
    Duration_SliSel2(ii) = str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.GradientPulse{1,2}.Attributes.Duration)./1E6;
    RiseTime_SliSel2(ii) = str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.GradientPulse{1,2}.Attributes.RiseTime)./1E6;
    
    GDiff2(:,ii) = [str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.GradientPulse{1,3}.Attributes.Read); ...
        str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.GradientPulse{1,3}.Attributes.Phase); ...
        str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.GradientPulse{1,3}.Attributes.Slice);]./1E3;
    Duration_Diff2(ii) = str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.GradientPulse{1,3}.Attributes.Duration)./1E6;
    RiseTime_Diff2(ii) = str2num(xmlStruct.SequenceManifest.DiffusionGradientPulseArray.DiffusionArrayEntry{1,ii}.GradientPulse{1,3}.Attributes.RiseTime)./1E6;
    
end
Diff1 = @(t,ii) Trap(t, GDiff1(:,ii), DiffModuleStart + PreDiffFill(ii), Duration_Diff1(ii), RiseTime_Diff1(ii));
SliSel2 = @(t,ii) Trap(t, GSliSel2(:,ii), DiffModuleStart + PreDiffFill(ii) + Duration_Diff1(ii) + RiseTime_Diff1(ii) + DiffFill1(ii), Duration_SliSel2(ii), RiseTime_SliSel2(ii));
Diff2 = @(t,ii) Trap(t, GDiff2(:,ii), DiffModuleStart + PreDiffFill(ii) + Duration_Diff1(ii) + RiseTime_Diff1(ii) + DiffFill1(ii) + Duration_SliSel2(ii) + RiseTime_SliSel2(ii) + DiffFill2(ii), Duration_Diff2(ii), RiseTime_Diff2(ii));


t3DTabStart = DiffModuleStart + PreDiffFill(1) + DiffFill1(1) + Duration_Diff1(1) + RiseTime_Diff1(1) + Duration_SliSel2(1) + RiseTime_SliSel2(1) + DiffFill2(1) + Duration_Diff2(1) + RiseTime_Diff2(1) + PostDiffFill(1) + TEFill2;

% Parse 3D Tab Details
G3DTab = [str2num(xmlStruct.SequenceManifest.GradientPulse{1,3}.Attributes.Read); ...
    str2num(xmlStruct.SequenceManifest.GradientPulse{1,3}.Attributes.Phase); ...
    str2num(xmlStruct.SequenceManifest.GradientPulse{1,3}.Attributes.Slice);]./1E3;

Duration_3DTab = str2num(xmlStruct.SequenceManifest.GradientPulse{1,3}.Attributes.Duration)./1E6;
RiseTime_3DTab = str2num(xmlStruct.SequenceManifest.GradientPulse{1,3}.Attributes.RiseTime)./1E6;

% Construct function for this gradient pulse
Tab3D = @(t) Trap(t, G3DTab, t3DTabStart, Duration_3DTab, RiseTime_3DTab);
f = @(t,ii) (SliSel(t) + SliSelRew(t) + Diff1(t,ii) + SliSel2(t,ii) + Diff2(t,ii) + Tab3D(t));
deltaT = 1E-7;
gamma = getGyromagneticRatio('1H');
for ii = 1:LengthOfDiffusionArray
    ff = @(t) f(t,ii).*GradientReversal(t,TE);
    FG = @(t) deltaT.*trapz(ff(0:deltaT:t).').';
    
    FG2 = @(t) (FG(t))*(FG(t)');
    b(:,:,ii) = gamma.^2.*quadv(FG2,0,TE,'ArrayValued',true)/1E6;
end

end

