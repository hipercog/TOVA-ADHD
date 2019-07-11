function sig=c_gstlt_erp(ALLEEG,chidx,outdir,conds,chlocs)
% C_GSTLT_ERP erp analysis of CENT gestalt protocol data, by conditions
% selected via the 'conds' param
% INPUT:
%   - 

if nargin<5
    [thisdir,~,~]=fileparts(mfilename('fullpath'));
    chlocs=fullfile(thisdir,'chanlocs128_eeg.elp');
end
if nargin<4,    conds={'*','*'};  end;
if nargin<3,    outdir='';  end
if ~isempty(outdir) && ~isdir(outdir),	mkdir(outdir);    end
if nargin<2,    chidx=bionimi2num('A1'); end;

lbwh=get(0,'ScreenSize');
imgtitle=[conds{1} ' vs ' conds{2} ' ' bionum2nimi(chidx)];
locs=readlocs(chlocs);
% baseshift=0.5;  % count in secs of time to shift baseline
alpha=0.05;
datadd=find(~cellfun(@isempty,strfind({ALLEEG.filename},conds{1})));
datsub=find(~cellfun(@isempty,strfind({ALLEEG.filename},conds{2})));
nda=numel(datadd);
nsb=numel(datsub);
if nda>nsb,     datadd=datadd(randsample(nda,nsb));
elseif nsb>nda,	datsub=datsub(randsample(nsb,nda));
end
[erp1,erp2,erpsub,~,sig]=pop_comperp(ALLEEG,1,datadd,datsub,...
    'alpha',alpha,...
    'chans',chidx,...
    'geom','scalp',...
    'title',imgtitle,...
    'std','on',...
    'addavg','on',...
    'subavg','on',...
    'diffavg','on',...
    'tplotopt',{'chanlocs',locs,'legend',conds,'showleg','on','vert',500%,...
%     'frames',[ALLEEG(1).xmin-baseshift ALLEEG(1).xmax-baseshift].*1000
    });
% indicate corrected significance, resize and print to disk
% set(gcf,'Position',lbwh);
% print(gcf,'-dpng',fullfile(outdir,['erp_' strrep(imgtitle,' ','') '_']));
close all
% Bonferroni-Holm significance correction
[corp,h]=bonf_holm(sig,alpha);
sig=corp.*h;
% new image based on erp outputs, better visuals
if sum(h)>0 && sum(h)~=numel(h)
    figure;
    % baseline and smooth
    wndw=ALLEEG(1).pnts/100*10;
    x_eq_z=find(ALLEEG(1).times==0);
    erp1=smooth(erp1(chidx,:)-mean(erp1(chidx,1:x_eq_z)),wndw);
    erp2=smooth(erp2(chidx,:)-mean(erp2(chidx,1:x_eq_z)),wndw);
    erpsub=smooth(erpsub(chidx,:)-mean(erpsub(chidx,1:x_eq_z)),wndw);
    % plot
    bar(erpsub.*h','FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
    hold on
    plot(erp1,'color','k','LineWidth',2);
    plot(erp2,'color','r','LineWidth',2);
    plot(erpsub,'color','g','LineWidth',2);
    set(gca,'YDir','Reverse');
    ylim=get(gca,'ylim');
    line([x_eq_z x_eq_z],ylim,'color','m','LineWidth',2,'LineStyle','--');
    % label
    ylabel('Potential (uV)');
    xlabel('Time (ms)');
    xTicks = get(gca, 'xtick');
%     newlbl=round((((0:100:ALLEEG(1).pnts+100).*1000)./ALLEEG(1).srate)-...
%         baseshift*1000);
    newlbl=round(ALLEEG(1).times([1 xTicks(2:end-1) end]));
    set(gca,'XTickLabel',cellfun(@num2str,num2cell(newlbl),'UniformOutput',0));
%% WIP - plot one p text for each sig group line
    ptx={'***','**','*'};
    pvl=[0.001 0.01 0.05];
    up=find(diff(h)==1);
    down=find(diff(h)==-1);
    for p=1:numel(up)
        if p<=numel(down)
            [pmin,~]=min(corp(up(p):down(p)));
            pidx=find(pmin==corp);
            offset=(erpsub(pidx)/10)*erpsub(pidx)/erpsub(pidx);
            text(pidx,double(erpsub(pidx)+offset),...
                ptx(find(pmin<pvl,1,'first')),'FontSize',15);
        end
    end
    legend(gca,'* p<.05, ** p<.01, *** p<.001',conds{1},conds{2},'difference',...
        'Location','NorthWest');
    title(gca,imgtitle);
    % resize and print to disk
    set(gcf,'Position',lbwh);
    save(fullfile(outdir,['erp_' strrep(imgtitle,' ','')]),gcf,'-dpng');
    print(gcf,'-dpng',fullfile(outdir,['erp_' strrep(imgtitle,' ','')]));
end
close all
fclose all %#ok<PRTCAL>
end % c_gstlt_erp()