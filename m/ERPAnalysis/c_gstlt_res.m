function [fnames,bandpows,bndpwmxs] = ...
    c_gstlt_res(choose,conv,fnames,bandpows,bndpwmxs,outdir)
% C_GSTLT_RES analysis of CENT gestalt protocol data, by conditions which
% are selected via the 'choose' parameter
% INPUT:
%   - 
if nargin < 6,  outdir='';  end
if ~isempty(outdir) && ~isdir(outdir),	mkdir(outdir);    end
% choose some files and do some ctapping over them
files=dir(choose);
numf=numel(files);
[d,f,~]=fileparts(choose);
f=[strrep(f(1),'*','') f(2:end)];
structname=strrep(f,'*','_');
fnames.(structname)={files.name};
% choice of post-processing functions
switch conv
    case 'powers'
        fid=fopen(fullfile(d,'IAPFs.txt'));
        textscan(fid,'%s%s%s',1);
        subjiapf=textscan(fid,'%d%s%f','Delimiter','\t');
        subjs=subjiapf{:,1};
%         conds=subjiapf{:,2};
        iapfs=subjiapf{:,3};
        bandpows.(structname)=zeros(numf,5);
        bndpwmxs.(structname)=zeros(numf,5);
        for i=1:numf
            eeg=pop_loadset('filename',files(i).name,'filepath',d);
            fsubj=str2double(eeg.CTAP.subject(1:4));
            iaf=iapfs(subjs==fsubj);
            if isscalar(iaf)
                iep = [0.5 iaf*0.4 iaf*0.75 iaf*1.2 iaf*2 40];
            else
                disp('** Warning: IAF not read **'); continue;
            end
            nbands = length(iep)-1;
            for b=1:nbands
                [pow,powmx]=genBandPow(eeg,...
                    strcmp({eeg.chanlocs.type},'EEG'),iep(b:b+1));
                bandpows.(structname)(i,b)=pow;
                bndpwmxs.(structname)(i,b)=powmx;
            end
        end
end

end % c_gstlt_res()