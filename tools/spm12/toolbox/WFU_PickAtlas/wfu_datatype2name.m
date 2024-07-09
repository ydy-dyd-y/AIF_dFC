function  [dataname,bitdepth] = wfu_datatype2name(datatype)
dataname=spm_type(datatype);
bitdepth=spm_type(datatype,'bits');


