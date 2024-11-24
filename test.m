load('E:\OneDrive - Indian Institute of Science\VK Project Dhyaan\V_ProjectDhyaanBK1Programs\data\ftData\064PK\G1_ep_v8.mat');
load('E:\softwares\fieldtrip-20240614\template\layout\rnet64.mat');

goodChannels     = setdiff(data.label,data.label(data.badElecs));

%%%%%%%%%%%%%%%%%%%%%%% Time-frequency representation %%%%%%%%%%%%%%%%%%%%%
windowLenS       = 0.25; % Seconds
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = goodChannels;
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = (1/windowLenS):(1/windowLenS):100; 
cfg.t_ftimwin    = ones(length(cfg.foi),1).* windowLenS;   % length of time window
cfg.toi          = -1:0.05:1.25;

TFR = ft_freqanalysis(cfg,data);

cfg = [];
cfg.baseline     = [-1 0];
cfg.baselinetype = 'db';
cfg.showlabels   = 'yes';
cfg.layout       = layout;
cfg.box          = 'yes';
cfg.colormap     = 'jet';

figure; ft_multiplotTFR(cfg,TFR);
title('Time-frequency representation');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Data segmentation into baseline or stimulus epoch %%%%%%%%%%%%%%%

cfg                 = [];
cfg.channel         = goodChannels;
cfg.toilim          = [-1 + 1/1000 0];
data_bl             = ft_redefinetrial(cfg, data);
cfg.toilim          = [0.25 + 1/1000 1.25]; % Stimulus time period 0.25 to 0.75
data_st             = ft_redefinetrial(cfg, data);

%%%%%%%%%%%%%%%%%%%%%%%% Frequency analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg                 = [];
cfg.channel         = goodChannels;
cfg.output          = 'powandcsd';
cfg.method          = 'mtmfft';
cfg.taper           = 'dpss';
cfg.tapsmofrq       = 2;  % Smoothing over W = +- 4 Hz. For T=0.5, W=4, TW=2, 3 tapers are needed 
cfg.foilim          = [0 100];
cfg.keeptrials      = 'yes';

data_freq_bl        = ft_freqanalysis(cfg, data_bl);
data_freq_st        = ft_freqanalysis(cfg, data_st);

% Convert to Log
cfg = [];
cfg.operation = '10*log10(x1)';
cfg.parameter = 'powspctrm'; % this is represented by x, and x1 refers to the first (and only) input data structure
data_freq_bl_db = ft_math(cfg,data_freq_bl);
data_freq_st_db = ft_math(cfg,data_freq_st);

cfg = [];
cfg.layout        = layout;
cfg.showlabels    = 'yes';
cfg.showoutline   = 'yes';
figure; ft_multiplotER(cfg,data_freq_bl_db,data_freq_st_db);
title('PowerSpectrum for the baseline and stimulus periods');

%%%%%%%%%%%%%%%%%%%%%%%% Connectivity Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = 'ppc'; % choose between coh, imc, wpli or ppc

cfgA = []; % Configuration file for analysis
cgfD = []; % Configuration file for display
cfgD.refchannel    = 'PO3';
cfgD.layout        = layout;
cfgD.showlabels    = 'yes';
cfgD.showoutline   = 'yes';

if strcmp(method,'coh')
    cfgA.method = 'coh';
    cfgA.complex = 'abs';
    cfgD.parameter = 'cohspctrm';
    
elseif strcmp(method,'imc')
    cfgA.method = 'coh';
    cfgA.complex = 'imag';
    cfgD.parameter = 'cohspctrm';
    
elseif strcmp(method,'wpli')
    cfgA.method = 'wpli';
    cfgD.parameter = 'wplispctrm';
    
elseif strcmp(method,'ppc')
    cfgA.method = 'ppc';
    cfgD.parameter = 'ppcspctrm';
end

con_result_bl = ft_connectivityanalysis(cfgA,data_freq_bl);
con_result_st = ft_connectivityanalysis(cfgA,data_freq_st);
ft_multiplotER(cfgD,con_result_bl,con_result_st);

title('PPC for the baseline and stimulus periods');