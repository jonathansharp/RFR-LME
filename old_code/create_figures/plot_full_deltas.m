% define regions
define_regions_eiwg
% Set GMM options
set_gmm_options
% Set RFR options
set_rfr_options

%% plot gridded delta values across region
cmap_type = 'cmocean';
cmap_name = 'balance';
zero_piv = 1;
cmap_segs = 17;
plot_delta_mean_full(-17.5,17.5,cmap_type,cmap_name,zero_piv,...
    num_groups,'delta_rfr_grid','\Delta{\itf}_{CO2} (\muatm)',region,lme_shape,lme_idx,rfr_test_idx);

%% plot gridded absolute delta values across region
cmap_type = 'cmocean';
cmap_name = 'amp';
zero_piv = 0;
plot_delta_mean_full(0,50,cmap_type,cmap_name,...
    zero_piv,num_groups,'delta_rfr_grid_abs','|\Delta{\itf}CO_{2}| (\muatm)',...
    region,lme_shape,lme_idx,rfr_test_idx);

%% plot gridded absolute delta values across region
cmap_type = 'cmocean';
cmap_name = 'amp';
zero_piv = 0;
cmap_segs = 9;
plot_delta_mean_full(0,50,cmap_type,cmap_name,cmap_segs,...
    zero_piv,num_groups,'delta_rfr_grid_abs','\Delta{\itf}CO_{2}',...
    region,lme_shape,lme_idx,rfr_test_idx);