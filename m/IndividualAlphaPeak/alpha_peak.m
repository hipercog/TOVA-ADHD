% Calculates IAF and plots stuff
% iafcalcutron(fileEC,fileEO,ch,outputfolder)
% fileEC - filename of the eyes closed recording
% fileEO - filename of the eyes open recording
% ch - vector containing channel indexes ([] for all)
% outputfolder - folder to dump the filter config files

function [iaf,p,f]=alpha_peak(fileEC,fileEO,chIAF)	
	iaf=generateFilters(fileEC,fileEO,chIAF);
end


%% GENERATE FILTERS
function iaf = generateFilters(ecfile,eofile,iafch)

    % Using this for the IAF calculation:
    % M. Lansbergen et. al. 2010
    % The increase in theta/beta ratio on resting-state EEG in boys with
    % attention-deficit/hyperactivity disorder is mediated by slow alpha
    % peak frequency
    XLIM = [0 25];
    iaflimits = [8 12];
    nfft = 2048*2;
	eeg = pop_loadset(ecfile);
	if isempty(iafch), iafch=1:size(eeg.data,1); end;

	P_avg_EC = 0;

	figure();
	subplot(3,1,1);
		hold on;
		for k=1:length(iafch)
			tmp = psd(spectrum.welch('Hamming',1024,50), ...
                eeg.data(iafch(k),:),'NFFT',nfft,'Fs',eeg.srate);    
			plot(tmp.Frequencies,20*log10(tmp.Data),'k','linewidth',1);
			P_avg_EC = P_avg_EC + tmp.Data;
			
		end
		P_avg_EC = P_avg_EC./length(iafch);
		[val,idx]=max(P_avg_EC(tmp.Frequencies>iaflimits(1) & ...
                               tmp.Frequencies<iaflimits(2)));
		wx = tmp.Frequencies(tmp.Frequencies>iaflimits(1) & ...
                             tmp.Frequencies<iaflimits(2));
		iaf = wx(idx);
       
		plot(tmp.Frequencies,20*log10(P_avg_EC),'g','linewidth',2);
		hold off;
		set(gca,'xlim',XLIM,'ylim',[-3 40],'box','on');
		line([iaf iaf],[-3 40],'color','c');
		line([0 15],[20*log10(val) 20*log10(val)],'color','c');
		line([iaflimits(1) iaflimits(1)],[-3 40],'color','r');
        line([iaflimits(2) iaflimits(2)],[-3 40],'color','r');
        
		text(wx(idx)+0.1,double(20*log10(val))+1, ...
		sprintf('Peak = %0.2fHz',wx(idx)), ...
		'FontSize',6,'color','k');
		xlabel('Hz'); ylabel('dB');
		title('Eyes closed');


	% EYES OPEN
	eeg = pop_loadset(eofile);
	if isempty(iafch), iafch=1:size(eeg.data,1); end;

	P_avg_EO = 0;

	subplot(3,1,2);
		hold on;
		for k=1:length(iafch)
			tmp = psd(spectrum.welch('Hamming',1024,50), ...
                eeg.data(iafch(k),:),'NFFT',nfft,'Fs',eeg.srate);    
			plot(tmp.Frequencies,20*log10(tmp.Data),'k','linewidth',1);
			P_avg_EO = P_avg_EO + tmp.Data;
			
		end

		P_avg_EO = P_avg_EO./length(iafch);
		[val,idx]=max(P_avg_EO(tmp.Frequencies>iaflimits(1) & ...
                               tmp.Frequencies<iaflimits(2)));
		wx = tmp.Frequencies(tmp.Frequencies>iaflimits(1) & ...
                             tmp.Frequencies<iaflimits(2));
		iaf = wx(idx);


		plot(tmp.Frequencies,20*log10(P_avg_EO),'g','linewidth',2);
		hold off;
		set(gca,'xlim',XLIM,'ylim',[-3 40],'box','on');
		line([iaf iaf],[-3 40],'color','c');
		line([0 15],[20*log10(val) 20*log10(val)],'color','c');
		line([iaflimits(1) iaflimits(1)],[-3 40],'color','r');
        line([iaflimits(2) iaflimits(2)],[-3 40],'color','r');
		text(wx(idx)+0.1,double(20*log10(val))+1, ...
		sprintf('Peak = %0.2fHz',wx(idx)), ...
		'FontSize',6,'color','k');
		xlabel('Hz'); ylabel('dB');
		title('Eyes open');

	% DIFFERENCE
	P_diff = P_avg_EC-P_avg_EO;
	[val,idx]=max(P_diff(tmp.Frequencies>iaflimits(1) & ...
                         tmp.Frequencies<iaflimits(2)));
	wx = tmp.Frequencies(tmp.Frequencies>iaflimits(1) & ...
                         tmp.Frequencies<iaflimits(2));
	iaf = wx(idx);

	subplot(3,1,3);
		hold on;
		plot(tmp.Frequencies,20*log10(P_diff),'g','linewidth',2);
		plot(tmp.Frequencies,20*log10(P_avg_EO),'k');
		plot(tmp.Frequencies,20*log10(P_avg_EC),'k');
		hold off;
		set(gca,'xlim',XLIM,'ylim',[-3 40],'box','on');
		line([iaf iaf],[-3 40],'color','c');
		line([0 15],[20*log10(val) 20*log10(val)],'color','c');
        line([iaflimits(1) iaflimits(1)],[-3 40],'color','r');
        line([iaflimits(2) iaflimits(2)],[-3 40],'color','r');
		text(wx(idx)+0.1,double(20*log10(val))+1, ...
		sprintf('IAF = %0.2fHz',wx(idx)), ...
		'FontSize',6,'color','k');
		xlabel('Hz'); ylabel('dB');
		title('Difference');

end

