% Calculates IAF and plots stuff
% iafcalcutron(fileEC,fileEO,ch,outputfolder)
% fileEC - filename of the eyes closed recording
% fileEO - filename of the eyes open recording
% ch - vector containing channel indexes ([] for all)
% outputfolder - folder to dump the filter config files




%% GENERATE FILTERS
function [p,f] = spectrumAverage(eegfile,ch)

    % Using this for the IAF calculation:
    % M. Lansbergen et. al. 2010
    % The increase in theta/beta ratio on resting-state EEG in boys with
    % attention-deficit/hyperactivity disorder is mediated by slow alpha
    % peak frequency
%     XLIM = [0 25];
    
    nfft = 2048*2;
	eeg = pop_loadset(eegfile);
	if isempty(ch), ch=1:size(eeg.data,1); end;

	P_avg_EC = 0;


	

		for k=1:length(ch)
			tmp = psd(spectrum.welch('Hamming',1024,50), ...
                eeg.data(ch(k),:),'NFFT',nfft,'Fs',eeg.srate);    
% 			plot(tmp.Frequencies,20*log10(tmp.Data),'k','linewidth',1);
			P_avg_EC = P_avg_EC + tmp.Data;
			
		end
		P_avg_EC = P_avg_EC./length(ch);
		p = P_avg_EC;
		f = tmp.Frequencies;
        figure();
            plot(f,20*log10(p));
            set(gca,'xlim',[0 25],'ylim',[-3 40]);
            
            