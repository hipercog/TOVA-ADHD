function sig=c_gstlt_ersp(ALLEEG,chidx,outdir,conds)
% C_GSTLT_ERSP ersp analysis of CENT gestalt protocol data, by conditions
% INPUT:
%   - 

if nargin<4,    conds={'*','*'};  end;
if nargin<3,    outdir='';  end
if ~isempty(outdir) && ~isdir(outdir),	mkdir(outdir);    end
if nargin<2,    chidx=bionimi2num({'C5' 'D5'}); end;


datadd=find(~cellfun(@isempty,strfind({ALLEEG.setname},conds{1})));
datsub=find(~cellfun(@isempty,strfind({ALLEEG.setname},conds{2})));
nda=numel(datadd);
nsb=numel(datsub);
if nda>nsb,     datadd=datadd(randsample(nda,nsb));
elseif nsb>nda,	datsub=datsub(randsample(nsb,nda));
end
[eegc1,histc1]=pop_grandaverage(ALLEEG,'datasets',datadd);
[eegc2,histc2]=pop_grandaverage(ALLEEG,'datasets',datsub);

lbwh=get(0,'ScreenSize');   lbwh(3)=lbwh(3)*0.9;    lbwh(4)=lbwh(4)*0.9;
imgtitle=[conds{1} ' vs ' conds{2} ' ' bionum2nimi(chidx)];

[ersp,itc,powbase,times,freqs,sig,itcboot,itcphase] = ...
newtimef(...
    {eegc1.data(chidx,:,:) eegc2.data(chidx,:,:)},...
    ALLEEG(1).pnts,...
    [ALLEEG(1).times(1) ALLEEG(1).times(end)],...
    ALLEEG(1).srate,0,...
    'alpha',0.05,...
    'title',conds);
ax=axes('Units','Normal','Position',[.075 .09 .85 .85],'Visible','Off');
set(get(ax,'Title'),'Visible','On');
title(['ERSP: ' imgtitle],...
    'color','m','FontName','Courier','FontWeight','bold','FontSize',15);
set(gcf,'Position',lbwh);
save(fullfile(outdir,['ersp_' strrep(imgtitle,' ','')]),gcf,'-dpng');
% print(gcf,'-dpng',fullfile(outdir,['ersp_' strrep(imgtitle,' ','')]));
% pop_preclust
% std_erspplot
close all
end % c_gstlt_ersp()