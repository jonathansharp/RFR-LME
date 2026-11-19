function paths = rfr_lme_paths()
% Central locations for RFR-LME inputs, outputs, and downloaded datasets.
paths.data = fullfile('/data2', 'RFR-LME', 'Data');
paths.figures = fullfile('/data2', 'RFR-LME', 'Figures');
paths.models = fullfile('/data2', 'RFR-LME', 'Models');
paths.socat = fullfile('/data2', 'SOCATv2026');
paths.dpath = [fullfile('/med2', 'RFR-LME-Data') filesep];
end
