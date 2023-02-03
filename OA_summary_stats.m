% Save OA indicator summary statistics
% 
% Written by J.D. Sharp: 2/1/23
% Last updated by J.D. Sharp: 2/1/23
% 

% this script defines the bounds of the eighteen LMEs
define_regions
% region labels
reg_lab = {'CCS' 'GA' 'AI' 'EBS' 'BS' 'NBCS' 'NE' 'SE' 'GM' 'CS' ...
           'HI' 'AS' 'JI' 'PK' 'HB' 'JA' 'WI' 'GC'};
% variable information
y_span = [1700,2400;250,550;1800,2500;7.9,8.2;1,5;1.5,7.5;6,16;50,300;8,18];
var_type = {'DIC' 'fCO2' 'TA' 'pH' 'OmA' 'OmC' 'H' 'CO3' 'RF'};
var_lab = {'{\itC}_{T}' '{\itf}CO_{2}' '{\itA}_{T}' 'pH_{T}' '\Omega_{A}' ...
    '\Omega_{C}' '[H^{+}]' '[CO_{3}^{2-}]' 'RF'};
units = {'\mumol kg^{-1}' '\muatm' '\mumol kg^{-1}' '' '' '' 'nmol kg^{-1}' '\mumol kg^{-1}' ''};
rounder = [1 1 1 3 2 2 1 1 2];

% preallocate table
stats = nan(27,18);

% loop through each region
for n = 1:length(region)

    for var_num = 1:9

        % load estimated OA grid
        load(['Data/' region{n} '/gridded_pco2'],'SOCAT_grid');
        load(['Data/' region{n} '/ML_fCO2'],'OAI_grid');
    
        % calculate area-weighted time series
        OAI_grid.(region{n}).var_dom_mean = nan(OAI_grid.(region{n}).dim.z,1);
        area_weights = SOCAT_grid.(region{n}).area_km2.*SOCAT_grid.(region{n}).percent_sea;
        for t = 1:OAI_grid.(region{n}).dim.z
            OAI_grid.(region{n}).var_dom_mean(t) = ...
                squeeze(sum(sum(OAI_grid.(region{n}).(var_type{var_num})(:,:,t).*...
                    area_weights,1,'omitnan'),2,'omitnan'))./...
                    squeeze(sum(sum(area_weights,1,'omitnan'),2,'omitnan'));
        end
    
        % scale H to nanomoles
        if strcmp(var_type{var_num},'H')
            OAI_grid.(region{n}).var_dom_mean = (10^9).*OAI_grid.(region{n}).var_dom_mean;
        end
    
        % calculate trend
        [yf,yr,x] = ...
            leastsq2(OAI_grid.(region{n}).month,...
            OAI_grid.(region{n}).var_dom_mean,0,2,[6 12]);
    
        % determine average modelled climatology
        clim = nan(12,1);
        for m = 1:12
            clim(m) = mean(yf(m:12:end));
        end
    
        % calculate amplitude
        amp = max(clim) - min(clim);
    
        % fill table
        stats(n,(var_num-1)*3+1) = mean(OAI_grid.(region{n}).var_dom_mean);
        stats(n,(var_num-1)*3+2) = x(2)*12;
        stats(n,(var_num-1)*3+3) = amp;
    
        % clean up
        clear SOCAT_grid OAI_grid area_weights t time

    end

end

% save table
writematrix(stats,['Data/SummaryStatistics-' date '.xls']);
