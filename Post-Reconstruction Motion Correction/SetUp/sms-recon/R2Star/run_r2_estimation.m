load rInfo
load FM
load sen
load mask
TEmask = logical([1 0 0 1 1 1]);

[r2,we] = createR2starFM(rInfo, FM, sen,mask, TEmask);

