%% bersim.m
%%Sıla Elif Karaaağaç 23050211070

clear; clc; close all;
set(0,'DefaultFigureVisible','on');
set(groot,'DefaultFigureWindowStyle','normal');

%% User controls
useNoise = true;
snrGrid  = 0:2:30;
mcsRange = 8:10;

maxNumBits    = 1e7;   
maxNumPackets = 1e4;   

ber = zeros(numel(snrGrid), numel(mcsRange));

%% Base HT configuration
cfgHT = wlanHTConfig;
cfgHT.ChannelBandwidth    = 'CBW20';
cfgHT.PSDULength          = 1024;
cfgHT.GuardInterval       = 'Long';
cfgHT.ChannelCoding       = 'BCC';
cfgHT.NumTransmitAntennas = 2;
cfgHT.NumSpaceTimeStreams  = 2;

%% Loop
for m = 1:numel(mcsRange)
    cfgHT.MCS = mcsRange(m);

    ofdmInfo = wlanHTOFDMInfo('HT-Data', cfgHT);
    ind      = wlanFieldIndices(cfgHT);

    fprintf('\n===== MCS %d (2x2) =====\n', cfgHT.MCS);

    for i = 1:numel(snrGrid)
        packetSNR = snrGrid(i) - 10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);

        totalBitErrors       = 0;
        totalBitsTransmitted = 0;
        numPackets           = 0;

        while (totalBitsTransmitted < maxNumBits) && (numPackets < maxNumPackets)

            %% Transmitter
            txPSDU = randi([0 1], cfgHT.PSDULength*8, 1, 'int8');
            tx     = wlanWaveformGenerator(txPSDU, cfgHT);

            %% Channel 
            H        = (randn(2,2) + 1j*randn(2,2)) / sqrt(2);
            rx_clean = tx * H.';

            if useNoise
                rx   = awgn(rx_clean, packetSNR, 'measured');
                nVar = mean(abs(rx_clean(:)).^2) / 10^(packetSNR/10);
            else
                rx   = rx_clean;
                nVar = 1e-12;
            end

            %% Receiver
            htltf      = rx(ind.HTLTF(1):ind.HTLTF(2), :);
            htltfDemod = wlanHTLTFDemodulate(htltf, cfgHT);
            chanEst    = wlanHTLTFChannelEstimate(htltfDemod, cfgHT);

            htdata = rx(ind.HTData(1):ind.HTData(2), :);
            rxPSDU = wlanHTDataRecover(htdata, chanEst, nVar, cfgHT);

            %% BER
            L = min(numel(txPSDU), numel(rxPSDU));
            totalBitErrors       = totalBitErrors + sum(txPSDU(1:L) ~= int8(rxPSDU(1:L)));
            totalBitsTransmitted = totalBitsTransmitted + L;
            numPackets           = numPackets + 1;
        end

        ber(i,m) = totalBitErrors / max(totalBitsTransmitted, 1);
        fprintf('SNR=%.1f dB | BER=%e | packets=%d\n', snrGrid(i), ber(i,m), numPackets);
    end
end

%% Plot 
berPlot = ber;
berPlot(ber == 0) = NaN;  % Hide zero points of BER
figure('Color','w');
semilogy(snrGrid, berPlot, '-o', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('802.11n HT BER (2x2 MIMO, AWGN)');
legend(compose('MCS %d', mcsRange), 'Location', 'southwest');