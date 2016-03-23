% *************************************************************************
%
%   Realtime ECG detection from the Fieldtrip buffer
%
%   Description:
%       This script read blocks of data from the fieldtrip buffer,
%       apply signal processing pipeline (detrend, lowpass filtring, etc..)
%       do threshold detection and play a beep sound each time 
%       the threshold is reached.
%
%   Author : Romain Trachel 
%   Source : ft_realtime_signalviewer.m
%   Date : 04/21/2016
%
% *************************************************************************

% define fieldtrip configuration structure
cfg = struct;

% define blocksize in seconds
cfg.blocksize = 3;
% define overlap between blocs in percent
cfg.overlap = 95 * cfg.blocksize / 100;
% define your channel names 
cfg.channel = {'EX1','EX2'};
% low pass filter settings
cfg.lowpass = 'yes'; 
cfg.l_freq = 15; % low frequency band
cfg.filter_len  = 150;
% high pass filter settings
cfg.highpass = 'yes'; 
cfg.h_freq = 1; % high frequency band
cfg.hp_filter_len  = 250;
% detection settings
cfg.derivative = 'yes';
cfg.deriv_order = 2;
% threshold setting
cfg.threshold = 'std';
cfg.thresh_scale = 3.8;
% or it could be a constant value (e.g. -4*10^6);
cfg.vizualisation = 'yes';

% *************************************************************************
% set the buffer configuration options
% *************************************************************************
% default is detected automatically
if ~isfield(cfg, 'dataformat'),     cfg.dataformat = [];      end
% default is detected automatically
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = [];    end
% default is detected automatically
if ~isfield(cfg, 'eventformat'),    cfg.eventformat = [];     end
 % define blocksize in seconds
if ~isfield(cfg, 'blocksize'),      cfg.blocksize = 1;        end
% define overlap between blocs in seconds
if ~isfield(cfg, 'overlap'),        cfg.overlap = 0.90;          end
% define channel names
if ~isfield(cfg, 'channel'),        cfg.channel = 'all';      end
% take {first, last} data packet in the buffer?
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata = 'last';  end
% capture events {yes, no} ?
if ~isfield(cfg, 'readevent'),      cfg.readevent = 'no';     end
% jump to end of file at initialization
if ~isfield(cfg, 'jumptoeof'),      cfg.jumptoeof = 'no';     end
% connect to the buffer
if ~isfield(cfg, 'dataset') && ~isfield(cfg, 'header') && ~isfield(cfg, 'datafile')
  cfg.dataset = 'buffer://localhost:1972';
end
% *************************************************************************
% set the signal processing configuration options
% *************************************************************************
% enable baseline correction
if ~isfield(cfg, 'demean'),         cfg.demean = 'yes';       end
% enable detrend correction
if ~isfield(cfg, 'detrend'),        cfg.detrend = 'yes';       end
% enable low pass filtering
if ~isfield(cfg, 'lowpass'),         cfg.lowpass = 'yes';   end
% low pass frequency (see ft_preproc_lowpassfilter)
if ~isfield(cfg, 'l_freq'),         cfg.l_freq = 45;       end
% filter type (see ft_preproc_lowpassfilter)
if ~isfield(cfg, 'filter_type'),      cfg.filter_type  = 'fir';   end
% filter length (see ft_preproc_lowpassfilter)
if ~isfield(cfg, 'filter_len'),      cfg.filter_len  = 25;   end
% enable derivative temporal filter 
if ~isfield(cfg, 'derivative'),   cfg.derivative = 'no';    end
% derivative order (ft_preproc_derivative)
if ~isfield(cfg, 'deriv_order'),     cfg.deriv_order = 2;   end
% threshold type {'std', 'no'}
if ~isfield(cfg, 'threshold'),  cfg.threshold = 'std';  end
% threshold scaling factor
if ~isfield(cfg, 'thresh_scale'),  cfg.thresh_scale = 1;    end
% enable vizualisation 
if ~isfield(cfg, 'vizualisation'),  cfg.vizualisation = 'yes';       end

if ~isfield(cfg, 'beep_file');
    % load file from windows (but it wont work for mac or linux)
    cfg.beep_file = 'C:\Windows\Media\ding.wav';
end
% read an audiofile for beep feedback
[beep, beep_fs] = audioread(cfg.beep_file);

% translate dataset into datafile+headerfile
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
cfg = ft_checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear read_header
% start by reading the header from the realtime buffer
hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true, 'retry', true);

chanindx    = match_str(hdr.label, cfg.channel);
nchan       = length(chanindx);
if nchan==0
  error('no channels were selected');
end
 
% determine the size of blocks to process
blocksize = round(cfg.blocksize * hdr.Fs);
overlap   = round(cfg.overlap*hdr.Fs);
 
if strcmp(cfg.jumptoeof, 'yes')
  prevSample = hdr.nSamples * hdr.nTrials;
else
  prevSample  = 0;
end
count       = 0;

% play first beep
sound(beep, beep_fs);

% *************************************************************************
% this is the main task loop where realtime incoming data is processed
% *************************************************************************
while true
 
    % determine number of samples available in buffer
    hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true);

    % see whether new samples are available
    newsamples = (hdr.nSamples*hdr.nTrials-prevSample);

    % *************************************************************************
    %   if there's enought available data: 
    %       read buffer + signal processing + ECG detection
    %   else:
    %       execute the task code (display, servo-commands, etc...)
    % *************************************************************************
    if newsamples>=blocksize
        % determine the samples to process
        if strcmp(cfg.bufferdata, 'last')
            begsample  = hdr.nSamples*hdr.nTrials - blocksize + 1;
            endsample  = hdr.nSamples*hdr.nTrials;
        elseif strcmp(cfg.bufferdata, 'first')
            begsample  = prevSample+1;
            endsample  = prevSample+blocksize ;
        else
            error('unsupported value for cfg.bufferdata');
        end

        % this allows overlapping data segments
        if overlap && (begsample>overlap)
            begsample = begsample - overlap;
            endsample = endsample - overlap;
        end

        % remember up to where the data was read
        prevSample  = endsample;
        count       = count + 1;
        fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

        % read data segment from buffer
        dat = ft_read_data(cfg.datafile, 'header', hdr, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

        % if you want to read some events as well...
        if strcmp(cfg.readevent, 'yes')
            evt = ft_read_event(cfg.eventfile, 'header', hdr, 'minsample', begsample, 'maxsample', endsample);
        end

        % put the data in a fieldtrip-like raw structure
        data.trial{1} = dat;
        data.time{1}  = 1:length(dat)/hdr.Fs;
        data.label    = hdr.label(chanindx);
        data.hdr      = hdr;
        data.fsample  = hdr.Fs;

        % *********************************************************************
        % 1st step : Pre-processing (detrend, artifact detection, etc...)
        % *********************************************************************
        if strcmp(cfg.detrend, 'yes')
            % detrend the data to remove dc offset and slow drifts
            data.trial{1} = ft_preproc_detrend(data.trial{1}, 1, length(dat)-1);
        end

        % *********************************************************************
        % 2nd step : Temporal filtering (lowpass, highpass, bandpass, etc...)
        % *********************************************************************
        if strcmp(cfg.lowpass, 'yes')
            % lowpass filtering to remove 50hz and more...
            data.trial{1} = ft_preproc_lowpassfilter(data.trial{1}, hdr.Fs, ... 
                                           cfg.l_freq, cfg.filter_len, ...
                                           cfg.filter_type, 'onepass');
        end
        if strcmp(cfg.highpass, 'yes')
            % high filtering to remove dc...
            data.trial{1} = ft_preproc_highpassfilter(data.trial{1}, hdr.Fs, ... 
                                           cfg.h_freq, cfg.hp_filter_len, ...
                                           cfg.filter_type, 'onepass');
        end
        % *********************************************************************
        % 3rd step : Pre-thresholding (derivative, baseline, wavelet, etc...)
        % *********************************************************************
        if strcmp(cfg.derivative, 'yes')
            % temporal Nth order derivative
            data.trial{1} = ft_preproc_derivative(data.trial{1}, cfg.deriv_order);
        else
            % or baseline correction to detect the peak amplitude
            data.trial{1} = ft_preproc_baselinecorrect(data.trial{1}, ...
                                        cfg.filter_len+1, length(dat)-1);
        end

        % *********************************************************************
        % 4th step : Thresolding (fixed, std, ci95, etc...)
        % *********************************************************************
        if strcmp(cfg.threshold, 'std')
            % compute std over the current block of data
            data.threshold = cfg.thresh_scale * std(data.trial{1}(cfg.filter_len+1:end));
        else % or use a fixed threhsold
            data.threshold = cfg.threshold;
        end % add more thresholding methods here...

        % see if the last packet of data reached the threshold
        if sum(abs(data.trial{1}(:, end-blocksize+overlap-1:end)) > abs(data.threshold))
            % *****************************
            %   Send servo-command here
            % *****************************
            sound(beep, beep_fs);
        end

        if strcmp(cfg.vizualisation, 'yes')
            % plot the data just like a standard FieldTrip raw data strucute
            plot(data.time{1}(cfg.filter_len+1:end), data.trial{1}(:, cfg.filter_len+1:end));
            xlim([data.time{1}(cfg.filter_len) data.time{1}(end)]);
            ylim([-abs(data.threshold), abs(data.threshold)]);
            % force Matlab to update the figure
            drawnow
        end
    else
        % *********************************************************************
        %   Put your own code here (e.g. display, servo-commands, etc...)
        % *********************************************************************
        
    end % if enough new samples
end % while true
