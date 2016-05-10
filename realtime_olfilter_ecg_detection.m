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
%   Date : 05/30/2016
%
% *************************************************************************

clear all
% define fieldtrip configuration structure
cfg = struct;
% define blocksize in seconds
cfg.blocksize = 0.005;
% define overlap between blocs in percent
cfg.overlap = 90 * cfg.blocksize / 100;
% define your channel names 
cfg.channel = {'1'};
% detection settings
cfg.derivative = 'no';
cfg.deriv_order = 2;
% define time window length
% for visualisation and auto thresholding
cfg.timelen = 1;
% threshold setting
cfg.threshold = 'std'; %10e7;
cfg.thresh_scale = 1.5;
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
% continuous online filter
if ~isfield(cfg, 'olfilter'), cfg.olfilter = 'yes';    end
if ~isfield(cfg, 'olfiltord'), cfg.olfiltord =  3;  end
if ~isfield(cfg, 'olfreq'),      cfg.olfreq = [2 45]; end

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
overlap   = round(cfg.overlap * hdr.Fs);

if strcmp(cfg.jumptoeof, 'yes')
  prevSample = hdr.nSamples * hdr.nTrials;
else
  prevSample  = 0;
end
count    = 0;
prevPeak = prevSample;
prevState = [];

blen = cfg.timelen * hdr.Fs / blocksize;
bdata = zeros([blen, 1]);
mdata = 0;
sdata = 1;
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
        % fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

        % read data segment from buffer
        dat = ft_read_data(cfg.datafile, 'header', hdr, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);
        dat = double(dat);
        
        % if you want to read some events as well...
        if strcmp(cfg.readevent, 'yes')
            evt = ft_read_event(cfg.eventfile, 'header', hdr, 'minsample', begsample, 'maxsample', endsample);
        end

        % put the data in a fieldtrip-like raw structure
        data.trial{1} = dat;
        data.time{1}  = ((begsample:endsample)*1.)/hdr.Fs;
        data.label    = hdr.label(chanindx);
        data.hdr      = hdr;
        data.fsample  = hdr.Fs;
        % *********************************************************************
        % *                 Apply online filter
        % *********************************************************************
        if strcmp(cfg.olfilter, 'yes')
            if ~exist('FM', 'var')
              % initialize online filter
              if cfg.olfreq(1)==0
                fprintf('using online low-pass filter\n');
                [B, A] = butter(cfg.olfiltord, 2 * cfg.olfreq(2)/hdr.Fs);
              elseif cfg.olfreq(2) >= hdr.Fs/2
                fprintf('using online high-pass filter\n');
                [B, A] = butter(cfg.olfiltord, 2 * cfg.olfreq(1)/hdr.Fs, 'high');
              else
                fprintf('using online band-pass filter\n');
                [B, A] = butter(cfg.olfiltord, 2 * cfg.olfreq/hdr.Fs);
              end
              % use one sample to initialize
              [dat, z] = filter(B, A, dat);
            end
            % apply online filter
            [dat, z] = filter(B, A, dat, z', 2);
            FM.z = double(z');
        end

        data.trial{1} = dat; % * 1e-8;
        bdata = circshift(bdata, 1);
        bdata(end) = mean(data.trial{1});
        
        % *********************************************************************
        % 3rd step : Pre-thresholding (derivative, baseline, wavelet, etc...)
        % *********************************************************************
        if strcmp(cfg.derivative, 'yes')
            % temporal Nth order derivative
            s = diff(bdata, cfg.deriv_order);
        else
            s = bdata;
        end
        % *********************************************************************
        % *       Thresolding (fixed, std, ci95, etc...)
        % *********************************************************************
        if strcmp(cfg.threshold, 'std')
            % compute std over the current block of data
            if count > blen
                mdata = mean(s);
                sdata = std(s);
            end
            data.threshold = cfg.thresh_scale * sdata;
        else % or use a fixed threhsold
            data.threshold = cfg.threshold;
        end % add more thresholding methods here...
        
        % see if the last packet of data reached the threshold
        if abs(s(1)) > abs(data.threshold)
            % *****************************
            %   Send servo-command here
            % *****************************
            curPeak = begsample + overlap;
            % if previous peak wasn't in the past 300ms
            if (curPeak - prevPeak) > .33 * hdr.Fs
                sound(beep, beep_fs);
            end
            % save ECG peak sample index
            prevPeak = curPeak;
        end

        if strcmp(cfg.vizualisation, 'yes')

            % plot the data just like a standard FieldTrip raw data strucute
            % plot(data.time{1}, data.trial{1});
            tmin = 0;
            tmax = blen/hdr.Fs;
            plot(linspace(tmin, tmax, blen), s);
            hold on
            plot([tmin, tmax], [data.threshold, data.threshold], '--r');
            plot([tmin, tmax], [-data.threshold, -data.threshold], '--r');
            xlim([tmin, tmax]);
            ylim([-abs(data.threshold) * 1.1, abs(data.threshold) * 1.1]);
            hold off
            % force Matlab to update the figure
            drawnow
        end
    else
        % *********************************************************************
        %   Put your own code here (e.g. display, servo-commands, etc...)
        % *********************************************************************
        
    end % if enough new samples
end % while true
