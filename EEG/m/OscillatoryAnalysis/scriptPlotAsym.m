%% Band Power asymmetry script

pathPT = 'D:\CENT_data\batch\Lab1_EEG_bdfs_logs\LateralityADHD';
pathHC = 'D:\CENT_data\batch\Healthy_Controls_labEEG\DONE_events_+filter_+badchans_+ICA_+cleaned_+interpd_+ECEO';

% Start with EC/EO data
oldf = cd( pathPT );
if ~isdir( 'bandPowAnalysis' )
    load( 'P3C3F3P4C4F4.mat' );
    preproCENT( 7, pwd, 'outdir','bandPowAnalysis', 'rois',roiLR );
end

cd bandPowAnalysis;

% derive the asymmetry data
asymPT = analysebandpwrs( 0 );

% plot the asymmetry pics
plotBandPowAsym( asymPT );

asymPT = analysebandpwrs( 1 );

% Start with EC/EO data
cd( pathHC );
if ~isdir( 'bandPowAnalysis' )
    load( 'P3C3F3P4C4F4.mat' );
    preproCENT( 7, pwd, 'outdir','bandPowAnalysis', 'rois',roiLR );
end

cd bandPowAnalysis;

% derive the asymmetry data
asymHC = analysebandpwrs( 0 );

% plot the asymmetry pics
plotBandPowAsym( asymHC );

asymHC = analysebandpwrs( 1 );

cd( oldf );

% Test the asyms
grprksp = zeros(3,5,9);
grprksh = zeros(3,5,9);

for rdx = 1:3
    for cdx = 1:5
        pt = asymPT.all(rdx,cdx,:); pt = pt(:)';
        hc = asymHC.all(rdx,cdx,:); hc = hc(:)';
        [grprksp(rdx,cdx,1), grprksh(rdx,cdx,1)] = ranksum( pt, hc );
        pt = asymPT.allec(rdx,cdx,:); pt = pt(:)';
        hc = asymHC.allec(rdx,cdx,:); hc = hc(:)';
        [grprksp(rdx,cdx,2), grprksh(rdx,cdx,2)] = ranksum( pt, hc );
        pt = asymPT.alleo(rdx,cdx,:); pt = pt(:)';
        hc = asymHC.alleo(rdx,cdx,:); hc = hc(:)';
        [grprksp(rdx,cdx,3), grprksh(rdx,cdx,3)] = ranksum( pt, hc );
        pt = asymPT.gbyec(rdx,cdx,:); pt = pt(:)';
        hc = asymHC.gbyec(rdx,cdx,:); hc = hc(:)';
        [grprksp(rdx,cdx,4), grprksh(rdx,cdx,4)] = ranksum( pt, hc );
        pt = asymPT.gbyeo(rdx,cdx,:); pt = pt(:)';
        hc = asymHC.gbyeo(rdx,cdx,:); hc = hc(:)';
        [grprksp(rdx,cdx,5), grprksh(rdx,cdx,5)] = ranksum( pt, hc );
        pt = asymPT.helec(rdx,cdx,:); pt = pt(:)';
        hc = asymHC.helec(rdx,cdx,:); hc = hc(:)';
        [grprksp(rdx,cdx,6), grprksh(rdx,cdx,6)] = ranksum( pt, hc );
        pt = asymPT.heleo(rdx,cdx,:); pt = pt(:)';
        hc = asymHC.heleo(rdx,cdx,:); hc = hc(:)';
        [grprksp(rdx,cdx,7), grprksh(rdx,cdx,7)] = ranksum( pt, hc );
        pt = asymPT.vigec(rdx,cdx,:); pt = pt(:)';
        hc = asymHC.vigec(rdx,cdx,:); hc = hc(:)';
        [grprksp(rdx,cdx,8), grprksh(rdx,cdx,8)] = ranksum( pt, hc );
        pt = asymPT.vigeo(rdx,cdx,:); pt = pt(:)';
        hc = asymHC.vigeo(rdx,cdx,:); hc = hc(:)';
        [grprksp(rdx,cdx,9), grprksh(rdx,cdx,9)] = ranksum( pt, hc );
    end
end


save( 'asymtTests', 'grpth', 'grptp' );