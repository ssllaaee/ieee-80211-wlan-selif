%% 802.11a GRID OFDM + AWGN : BER vs SNR (Sim vs Theory)
% - Nfft=64, Ncp=16
% - Used subcarriers: -26..-1, +1..+26 (52 used)
% - Pilots: -21, -7, +7, +21 (4 pilots)
% - Data: remaining 48 subcarriers
% - AWGN only (no multipath, no equalization)
% - SNR axis is Es/N0 (per subcarrier symbol energy)

clear; clc; close all;
rng(1);

%% -------------------- Parameters --------------------
Nfft   = 64;
Ncp    = 16;

SNRdB  = 0:1:25;           % x-axis: Es/N0 (dB)
Mlist  = [4 16 64];        % QPSK, 16QAM, 64QAM

Nsym    = 800;             % OFDM data symbols per run
nFrames = 40;              % Monte Carlo runs per SNR (increase for smoother)

names = {'QPSK','16-QAM','64-QAM'};
cols  = {'r','b','g'};

%% -------------------- 802.11a Subcarrier Grid --------------------
% Subcarrier indices (k): -26..-1, +1..+26
sc_all   = [-26:-1 1:26];                 % 52 used
sc_pilot = [-21 -7 7 21];                 % 4 pilots
sc_data  = setdiff(sc_all, sc_pilot);     % 48 data

% MATLAB FFT bin mapping: bin = mod(k, Nfft) + 1
bin_all   = mod(sc_all,   Nfft) + 1;
bin_pilot = mod(sc_pilot, Nfft) + 1;
bin_data  = mod(sc_data,  Nfft) + 1;

nDataSC = numel(sc_data);                 % 48

% Pilot values (simplified fixed pilots)
pilotVal = [1 1 1 -1].';                  % 4x1, BPSK pilots

%% -------------------- Allocate results --------------------
berSim    = zeros(numel(Mlist), numel(SNRdB));
berTheory = zeros(numel(Mlist), numel(SNRdB));

%% -------------------- Main loop --------------------
for mi = 1:numel(Mlist)
    M = Mlist(mi);
    k = log2(M);

    % Theory uses Eb/N0 typically; our x-axis is Es/N0.
    % Eb/N0(dB) = Es/N0(dB) - 10log10(k)
    EbN0dB = SNRdB - 10*log10(k);
    berTheory(mi,:) = ber_theory_gray_qam(EbN0dB, M);

    for si = 1:numel(SNRdB)
        EsN0lin = 10^(SNRdB(si)/10);

        % We normalize constellation to Es=1 (per active subcarrier symbol).
        % With unitary FFT/IFFT scaling, noise variance per subcarrier is N0
        N0 = 1/EsN0lin;

        err = 0; tot = 0;

        for f = 1:nFrames
            %% (1) Bits -> QAM symbols for DATA subcarriers only (48 per OFDM symbol)
            bitsTx  = randi([0 1], nDataSC*Nsym*k, 1);
            symData = qam_gray_mod(bitsTx, M);                 % Es=1
            symDataMat = reshape(symData, nDataSC, Nsym);

            %% (2) Build frequency-domain OFDM symbols with 802.11a grid
            X = zeros(Nfft, Nsym);                             % DC/guards = 0
            X(bin_data, :)  = symDataMat;                      % data
            X(bin_pilot, :) = repmat(pilotVal, 1, Nsym);       % pilots

            %% (3) OFDM modulation: IFFT + CP
            x  = ifft(X, Nfft, 1) * sqrt(Nfft);                % unitary-ish
            xcp = [x(end-Ncp+1:end,:); x];
            tx  = xcp(:);

            %% (4) AWGN
            w  = sqrt(N0/2) * (randn(size(tx)) + 1j*randn(size(tx)));
            rx = tx + w;

            %% (5) Receiver: remove CP + FFT
            rMat  = reshape(rx, Nfft+Ncp, Nsym);
            rNoCP = rMat(Ncp+1:end,:);
            Y = fft(rNoCP, Nfft, 1) / sqrt(Nfft);

            %% (6) Extract data subcarriers -> demod -> BER
            Ydata  = Y(bin_data, :);
            bitsRx = qam_gray_demod(Ydata(:), M);

            err = err + sum(bitsRx ~= bitsTx);
            tot = tot + numel(bitsTx);
        end

        % Log-plot için: BER=0 olursa çizilemez -> floor koy
        if err == 0
            berSim(mi,si) = 0.5/tot;
        else
            berSim(mi,si) = err/tot;
        end
    end
end

%% -------------------- Plot (reference style) --------------------
%% -------------------- Plot (reference style) --------------------
figure('Color','w');
ax = axes; hold(ax,'on');

for mi = 1:numel(Mlist)
    c = cols{mi};
    semilogy(ax, SNRdB, berSim(mi,:),    [c 'o-'], 'LineWidth', 1.2, 'MarkerSize', 4);
    semilogy(ax, SNRdB, berTheory(mi,:), [c '--'], 'LineWidth', 1.2);
end

% FORCE LOG AXIS + LIMITS
ax.YScale = 'log';
ax.YLim   = [1e-5 1];
ax.XLim   = [min(SNRdB) max(SNRdB)];
grid(ax,'on'); grid(ax,'minor');

% Make axes white even if MATLAB is in dark mode
ax.Color  = 'w';
ax.XColor = 'k';
ax.YColor = 'k';

xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
title(sprintf('802.11a Grid OFDM+AWGN (Nfft=%d, Ncp=%d): Sim vs Theory', Nfft, Ncp));
legend('QPSK Sim','QPSK Theory','16-QAM Sim','16-QAM Theory','64-QAM Sim','64-QAM Theory', ...
       'Location','southwest');

%% ===================== Local functions =====================

function Pb = ber_theory_gray_qam(EbN0dB, M)
    EbN0 = 10.^(EbN0dB/10);

    if M == 2 || M == 4
        Pb = qfunc_local(sqrt(2*EbN0));   % BPSK/QPSK exact
        return;
    end

    k = log2(M);
    a = sqrt(3*k./(M-1) .* EbN0);
    Q = qfunc_local(a);

    % Gray square M-QAM BER approx (good match)
    Pb = (4/k)*(1-1/sqrt(M)).*Q - (4/k)*(1-1/sqrt(M)).^2 .* (Q.^2);
end

function y = qfunc_local(x)
    y = 0.5*erfc(x./sqrt(2));
end

function sym = qam_gray_mod(bits, M)
    k = log2(M);
    bits = bits(:);

    if M == 2
        sym = 1 - 2*bits;
        sym = sym / sqrt(mean(abs(sym).^2));
        return;
    end

    L = sqrt(M);
    if mod(k,2)~=0 || round(L)~=L
        error('M must be 2 or square QAM (4,16,64,...)');
    end

    Nsym = numel(bits)/k;
    if Nsym ~= round(Nsym)
        error('bit length must be multiple of log2(M)');
    end
    Nsym = round(Nsym);

    b = reshape(bits, k, Nsym);
    bh = k/2;

    bI = b(1:bh,:).';
    bQ = b(bh+1:end,:).';

    gI = bi2de(bI,'left-msb');
    gQ = bi2de(bQ,'left-msb');

    nI = gray2bin_int(gI);
    nQ = gray2bin_int(gQ);

    aI = 2*nI - (L-1);
    aQ = 2*nQ - (L-1);

    sym = (aI + 1j*aQ).';

    % Normalize to Es=1
    scale = sqrt(2*(L^2 - 1)/3);
    sym = sym / scale;
end

function bits = qam_gray_demod(sym, M)
    k = log2(M);
    sym = sym(:);

    if M == 2
        bits = real(sym) < 0;
        bits = bits(:);
        return;
    end

    L = sqrt(M);
    if mod(k,2)~=0 || round(L)~=L
        error('M must be 2 or square QAM (4,16,64,...)');
    end
    bh = k/2;

    % Undo normalization
    scale = sqrt(2*(L^2 - 1)/3);
    r = sym * scale;

    I = real(r); Q = imag(r);

    nI = round((I + (L-1))/2);
    nQ = round((Q + (L-1))/2);
    nI = min(max(nI,0), L-1);
    nQ = min(max(nQ,0), L-1);

    gI = bin2gray_int(nI);
    gQ = bin2gray_int(nQ);

    bI = de2bi(gI, bh, 'left-msb');
    bQ = de2bi(gQ, bh, 'left-msb');

    bitsMat = [bI bQ];
    bits = reshape(bitsMat.', [], 1);
end

function b = gray2bin_int(g)
    g = uint32(g); b = g;
    while any(g ~= 0)
        g = bitshift(g, -1);
        b = bitxor(b, g);
    end
    b = double(b);
end

function g = bin2gray_int(b)
    b = uint32(b);
    g = bitxor(b, bitshift(b, -1));
    g = double(g);
end