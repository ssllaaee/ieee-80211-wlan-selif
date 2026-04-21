%% tx_rx_chain.m
%% 802.11n HT-only SISO BER simulation
%% TX: PSDU -> scramble -> BCC -> puncture -> interleave -> map -> IFFT/CP
%% RX: HT-LTF -> ch.est -> FFT/equalize -> demod -> deinterleave -> depuncture -> Viterbi -> descramble
%% Packet: [HT-SIG | HT-STF | HT-LTF | HT-Data]

clearvars -except sanityOnly sanityMCSRange sanityPSDULength;
clc; close all;

if ~exist('sanityOnly',    'var'), sanityOnly      = false; end
if ~exist('sanityMCSRange','var'), sanityMCSRange  = 0;     end
if ~exist('sanityPSDULength','var'), sanityPSDULength = 64; end

%% ── Controls ──────────────────────────────────────────────────────────────
useNoise      = true;
useFlatFading = true;
snrGrid       = 0:2:30;
mcsRange      = 0:7;
PSDULength    = 1024;
quickMode     = true;
plotFloorBER  = true;

maxNumBits       = 1e6;  maxNumPackets    = 120;
minNumPackets    = 12;   targetBitErrors  = 150;
zeroErrorPktStop = 24;   progressEvery    = 8;

if quickMode
    maxNumBits = 3e5;  maxNumPackets    = 40;
    minNumPackets = 8; targetBitErrors  = 60;
    zeroErrorPktStop = 12;
end

if sanityOnly
    useNoise = false; useFlatFading = false;
    snrGrid = 0; mcsRange = sanityMCSRange; PSDULength = sanityPSDULength;
    maxNumPackets = 1; maxNumBits = PSDULength*8;
end

%% ── Init ──────────────────────────────────────────────────────────────────
cfg = struct('MCS', 0, 'PSDULength', PSDULength, 'ScramblerInit', int8([1 0 1 1 1 0 1]));
phy = ht20_phy();
ber = zeros(numel(snrGrid), numel(mcsRange));

%% ── Main loop ─────────────────────────────────────────────────────────────
for im = 1:numel(mcsRange)
    cfg.MCS = mcsRange(im);
    mcs = mcs_info(cfg.MCS);
    ref = make_reference_fields(cfg, phy, mcs);
    self_test_chain(cfg, phy, mcs, ref);
    fprintf('\n===== MCS %d (HT-only SISO) =====\n', cfg.MCS);

    for is = 1:numel(snrGrid)
        tSNR   = tic;
        snrDat = snrGrid(is) - 10*log10(phy.NFFT / mcs.Nsd);
        bitErr = 0; bitCnt = 0; pktCnt = 0;
        stopReason = 'max limits';

        while bitCnt < maxNumBits && pktCnt < maxNumPackets
            txPSDU = int8(randi([0 1], cfg.PSDULength*8, 1));
            txPkt  = build_packet(txPSDU, cfg, phy, mcs, ref);

            h  = 1;
            if useFlatFading, h = (randn + 1j*randn)/sqrt(2); end
            rx = txPkt * h;

            if useNoise
                nVar = mean(abs(rx).^2) / 10^(snrDat/10);
                rx   = rx + sqrt(nVar/2)*(randn(size(rx)) + 1j*randn(size(rx)));
            else
                nVar = 1e-12;
            end

            rxPSDU = recover_psdu(rx, cfg, phy, mcs, ref, nVar);
            L      = min(numel(txPSDU), numel(rxPSDU));
            bitErr = bitErr + sum(txPSDU(1:L) ~= rxPSDU(1:L));
            bitCnt = bitCnt + L;
            pktCnt = pktCnt + 1;

            if mod(pktCnt, progressEvery) == 0
                fprintf('  SNR=%.1f dB | pkt=%d | err=%d | BER~=%e\n', ...
                    snrGrid(is), pktCnt, bitErr, bitErr/max(bitCnt,1));
            end
            if pktCnt >= minNumPackets    && bitErr >= targetBitErrors
                stopReason = 'target bit errors reached'; break; end
            if pktCnt >= zeroErrorPktStop && bitErr == 0
                stopReason = 'early stop: zero errors';  break; end
        end

        ber(is,im) = bitErr / max(bitCnt,1);
        fprintf('SNR=%.1f dB | BER=%e | packets=%d | %s | %.2f s\n', ...
            snrGrid(is), ber(is,im), pktCnt, stopReason, toc(tSNR));
    end
end

if sanityOnly, fprintf('\nSanity-only run completed.\n'); end

%% ── Plot ──────────────────────────────────────────────────────────────────
berPlot = ber;
if plotFloorBER
    berPlot(ber == 0) = 0.5 / max(maxNumBits,1);
else
    berPlot(ber == 0) = NaN;
end
figure('Color','w');
semilogy(snrGrid, berPlot, '-o', 'LineWidth', 1.5);
grid on; xlabel('SNR (dB)'); ylabel('BER');
title('802.11n HT BER (HT-only SISO)');
legend(compose('MCS %d', mcsRange), 'Location','southwest');
if plotFloorBER
    subtitle('Zero-BER points shown at measurement floor, not true zero');
end

%% ════════════════════════════════════════════════════════════════════════
%%  PHY / MCS
%% ════════════════════════════════════════════════════════════════════════
function phy = ht20_phy()
    phy.NFFT       = 64;
    phy.Ncp        = 16;
    phy.DataIdx    = [-28:-22 -20:-8 -6:-1 1:6 8:20 22:28];
    phy.PilotIdx   = [-21 -7 7 21];
    phy.PilotVal   = [1 1 1 -1];
    phy.PilotPol   = [1 -1 1 1 -1 1 1 1 1 -1 1 1];
    phy.SigDataIdx = [-26:-22 -20:-8 -6:-1 1:6 8:20 22:26];
    phy.HTSTFIdx   = [-24 -20 -16 -12 -8 -4 4 8 12 16 20 24];
    phy.HTSTFVal   = [1 1 -1 -1 1 1 -1 1 -1 1 1 1] * sqrt(13/6);
    phy.HTLTFSub   = [-28:-1 1:28];
    phy.HTLTFSeq   = repmat([1 -1], 1, numel(phy.HTLTFSub)/2);
    phy.DataFFT    = sub2fft(phy.DataIdx,  phy.NFFT);
    phy.PilotFFT   = sub2fft(phy.PilotIdx, phy.NFFT);
    phy.LTFFFT     = sub2fft(phy.HTLTFSub, phy.NFFT);
end

function mcs = mcs_info(id)
    tbl = [52 1 1 2  26  52;   % MCS0  BPSK   1/2
           52 2 1 2  52 104;   % MCS1  QPSK   1/2
           52 2 3 4  78 104;   % MCS2  QPSK   3/4
           52 4 1 2 104 208;   % MCS3  16QAM  1/2
           52 4 3 4 156 208;   % MCS4  16QAM  3/4
           52 6 2 3 208 312;   % MCS5  64QAM  2/3
           52 6 3 4 234 312;   % MCS6  64QAM  3/4
           52 6 5 6 260 312];  % MCS7  64QAM  5/6
    assert(id>=0 && id<=7, 'Only MCS 0..7 supported.');
    r = tbl(id+1,:);
    mcs = struct('Nsd',r(1),'Nbpscs',r(2),'Rnum',r(3),'Rden',r(4),'Ndbps',r(5),'Ncbps',r(6));
end

%% ════════════════════════════════════════════════════════════════════════
%%  REFERENCE FIELDS & SELF-TEST
%% ════════════════════════════════════════════════════════════════════════
function ref = make_reference_fields(cfg, phy, mcs)
    Nservice = 16; Ntail = 6;
    Nsym = ceil((Nservice + cfg.PSDULength*8 + Ntail) / mcs.Ndbps);
    Npad = Nsym*mcs.Ndbps - (Nservice + cfg.PSDULength*8 + Ntail);
    ref.Nsym     = Nsym;   ref.Npad  = Npad;
    ref.Nservice = Nservice; ref.Ntail = Ntail;
    ref.HTSIG    = htsig_waveform(htsig_bits(cfg), phy);
    ref.HTSTF    = htstf_waveform(phy);
    ref.HTLTF    = htltf_waveform(phy);
    ref.Preamble = [ref.HTSIG; ref.HTSTF; ref.HTLTF];
end

function self_test_chain(cfg, phy, mcs, ref)
    x  = int8(randi([0 1], cfg.PSDULength*8, 1));
    % interleaver round-trip
    assert(all(x(1:mcs.Ncbps) == ht_deinterleave(ht_interleave(x(1:mcs.Ncbps), mcs.Nbpscs, mcs.Nsd), mcs.Nbpscs, mcs.Nsd)), ...
        'Interleaver mismatch MCS %d.', cfg.MCS);
    % full chain round-trip
    assert(all(x == recover_psdu(build_packet(x,cfg,phy,mcs,ref), cfg,phy,mcs,ref, 1e-12)), ...
        'Noise-free chain mismatch MCS %d.', cfg.MCS);
end

%% ════════════════════════════════════════════════════════════════════════
%%  TX CHAIN
%% ════════════════════════════════════════════════════════════════════════
function pkt = build_packet(psdu, cfg, phy, mcs, ref)
    % 1. Scramble (prepend service field)
    bits = ht_scramble([zeros(ref.Nservice,1,'int8'); psdu(:)], cfg.ScramblerInit);
    % 2. Append tail + pad, encode, puncture
    bits = puncture(bcc_encode([bits; zeros(ref.Ntail+ref.Npad,1,'int8')]), mcs.Rnum, mcs.Rden);
    bits = bits(1:ref.Nsym*mcs.Ncbps);
    % 3. Per-symbol: interleave -> modulate -> subcarrier map -> IFFT+CP
    txData = zeros(phy.NFFT+phy.Ncp, ref.Nsym);
    for n = 1:ref.Nsym
        seg = bits((n-1)*mcs.Ncbps + (1:mcs.Ncbps));
        d   = ht_modulate(ht_interleave(seg, mcs.Nbpscs, mcs.Nsd), mcs.Nbpscs);
        pol = phy.PilotPol(mod(n-1, numel(phy.PilotPol))+1);
        F   = zeros(phy.NFFT,1);
        F(phy.DataFFT)  = d;
        F(phy.PilotFFT) = (phy.PilotVal * pol).';
        t = ifft(F)*sqrt(phy.NFFT);
        txData(:,n) = [t(end-phy.Ncp+1:end); t];   % prepend CP
    end
    pkt = [ref.Preamble; txData(:)];
end

%% ════════════════════════════════════════════════════════════════════════
%%  RX CHAIN
%% ════════════════════════════════════════════════════════════════════════
function psdu = recover_psdu(rx, cfg, phy, mcs, ref, ~)
    [rxLTF, rxData] = split_ht_fields(rx, ref);
    H    = estimate_channel(rxLTF, phy);
    Y    = fft_equalize(rxData, phy, ref, H);
    bits = demodulate_symbols(Y, mcs, ref);
    bits = deinterleave_bits(bits, mcs, ref);
    bits = decode_bcc(bits, mcs, ref);
    psdu = descramble_bits(bits, cfg, ref);
end

function [rxLTF, rxData] = split_ht_fields(rx, ref)
    nSig = numel(ref.HTSIG); nStf = numel(ref.HTSTF); nLtf = numel(ref.HTLTF);
    rxLTF  = rx(nSig+nStf+1        : nSig+nStf+nLtf);
    rxData = rx(nSig+nStf+nLtf+1   : end);
end

function H = estimate_channel(rxLTF, phy)
    % Average two HT-LTF symbols, divide by known sequence
    raw = rxLTF(2*phy.Ncp+1:end);
    L1  = fft(raw(1:phy.NFFT))            / sqrt(phy.NFFT);
    L2  = fft(raw(phy.NFFT+1:2*phy.NFFT)) / sqrt(phy.NFFT);
    H   = zeros(phy.NFFT,1);
    H(phy.LTFFFT) = ((L1(phy.LTFFFT)+L2(phy.LTFFFT))/2) ./ phy.HTLTFSeq(:);
end

function Y = fft_equalize(rxData, phy, ref, H)
    % Remove CP, FFT, single-tap equalization
    Y = reshape(rxData, phy.NFFT+phy.Ncp, ref.Nsym);
    Y = fft(Y(phy.Ncp+1:end,:), [], 1) / sqrt(phy.NFFT);
    Y = Y(phy.DataFFT,:) ./ H(phy.DataFFT);
end

function bits = demodulate_symbols(Y, mcs, ref)
    bits = zeros(ref.Nsym*mcs.Ncbps, 1, 'int8');
    for n = 1:ref.Nsym
        id = (n-1)*mcs.Ncbps + (1:mcs.Ncbps);
        bits(id) = ht_demodulate(Y(:,n), mcs.Nbpscs);
    end
end

function bits = deinterleave_bits(bits, mcs, ref)
    for n = 1:ref.Nsym
        id = (n-1)*mcs.Ncbps + (1:mcs.Ncbps);
        bits(id) = ht_deinterleave(bits(id), mcs.Nbpscs, mcs.Nsd);
    end
end

function bits = decode_bcc(bits, mcs, ref)
    bits = viterbi_decode(depuncture(bits, mcs.Rnum, mcs.Rden), ref.Ntail);
end

function psdu = descramble_bits(bits, cfg, ref)
    bits = ht_descramble(bits(1:ref.Nservice+cfg.PSDULength*8), cfg.ScramblerInit);
    psdu = bits(ref.Nservice+1 : ref.Nservice+cfg.PSDULength*8);
end

%% ════════════════════════════════════════════════════════════════════════
%%  PREAMBLE WAVEFORMS
%% ════════════════════════════════════════════════════════════════════════
function bits = htsig_bits(cfg)
    sig1 = [u2b(cfg.MCS,7) 0 u2b(cfg.PSDULength,16)];
    sig2 = [0 1 0 0 0 0 0 0 0 0];
    bits = int8([sig1 sig2 crc8_htsig([sig1 sig2]) zeros(1,6)].');
end

function y = htsig_waveform(bits, phy)
    coded = bcc_encode(bits);
    y = zeros(2*(phy.NFFT+phy.Ncp), 1);
    for n = 1:2
        b = signal_interleave(coded((n-1)*48+(1:48)));
        F = zeros(phy.NFFT,1);
        F(sub2fft(phy.SigDataIdx,phy.NFFT)) = (1-2*double(b))/sqrt(2);
        F(phy.PilotFFT) = (phy.PilotVal * phy.PilotPol(n)).';
        t = ifft(F)*sqrt(phy.NFFT);
        y((n-1)*(phy.NFFT+phy.Ncp)+(1:phy.NFFT+phy.Ncp)) = [t(end-phy.Ncp+1:end); t];
    end
end

function y = htstf_waveform(phy)
    F = zeros(phy.NFFT,1);
    F(sub2fft(phy.HTSTFIdx,phy.NFFT)) = phy.HTSTFVal(:);
    t = ifft(F)*sqrt(phy.NFFT);
    y = repmat([t(end-phy.Ncp+1:end); t], 2, 1);
end

function y = htltf_waveform(phy)
    F = zeros(phy.NFFT,1);
    F(phy.LTFFFT) = phy.HTLTFSeq(:);
    t = ifft(F)*sqrt(phy.NFFT);
    y = [t(end-2*phy.Ncp+1:end); t; t];
end

%% ════════════════════════════════════════════════════════════════════════
%%  BCC ENCODER
%% ════════════════════════════════════════════════════════════════════════
function out = bcc_encode(bits)
    g0=[1 0 1 1 0 1 1]; g1=[1 1 1 1 0 0 1]; s=zeros(1,6);
    out = zeros(2*numel(bits), 1, 'int8');
    for k = 1:numel(bits)
        r = [double(bits(k)) s];
        out(2*k-1) = int8(mod(sum(r.*g0),2));
        out(2*k)   = int8(mod(sum(r.*g1),2));
        s = r(1:6);
    end
end

%% ════════════════════════════════════════════════════════════════════════
%%  PUNCTURING
%% ════════════════════════════════════════════════════════════════════════
function pat = puncture_pattern(Rnum, Rden)
    switch sprintf('%d%d',Rnum,Rden)
        case '12', pat = [];
        case '23', pat = [1 1 0 1];
        case '34', pat = [1 1 0 1 1 0];
        case '56', pat = [1 1 0 1 0 1 0 1 1 0];
        otherwise,  error('Unsupported code rate %d/%d',Rnum,Rden);
    end
end

function out = puncture(bits, Rnum, Rden)
    pat = puncture_pattern(Rnum, Rden);
    if isempty(pat), out=bits; return; end
    pat = repmat(pat(:), ceil(numel(bits)/numel(pat)), 1);
    out = bits(logical(pat(1:numel(bits))));
end

function out = depuncture(bits, Rnum, Rden)
    pat = puncture_pattern(Rnum, Rden);
    if isempty(pat), out=double(bits); return; end
    out = -ones(ceil(numel(bits)/sum(pat))*numel(pat), 1);
    ii  = 1;
    for k = 1:numel(out)
        if pat(mod(k-1,numel(pat))+1) && ii<=numel(bits)
            out(k)=double(bits(ii)); ii=ii+1;
        end
    end
end

%% ════════════════════════════════════════════════════════════════════════
%%  INTERLEAVER
%% ════════════════════════════════════════════════════════════════════════
function out = ht_interleave(bits, Nbpscs, Nsd)
    Ncbps=Nsd*Nbpscs; Ncol=13; Nrow=4*Nbpscs; s=max(Nbpscs/2,1);
    k=(0:Ncbps-1).';
    i=Nrow*mod(k,Ncol)+floor(k/Ncol);
    j=s*floor(i/s)+mod(i+Ncbps-floor(Ncol*i/Ncbps),s);
    tmp=zeros(Ncbps,1,'int8'); out=zeros(Ncbps,1,'int8');
    tmp(i+1)=int8(bits(k+1)); out(j+1)=tmp(i+1);
end

function out = ht_deinterleave(bits, Nbpscs, Nsd)
    Ncbps=Nsd*Nbpscs; Ncol=13; Nrow=4*Nbpscs; s=max(Nbpscs/2,1);
    k=(0:Ncbps-1).';
    i=Nrow*mod(k,Ncol)+floor(k/Ncol);
    j=s*floor(i/s)+mod(i+Ncbps-floor(Ncol*i/Ncbps),s);
    out=zeros(Ncbps,1,'int8'); out(k+1)=int8(bits(j+1));
end

function out = signal_interleave(bits)
    k=(0:47).'; i=3*mod(k,16)+floor(k/16);
    out=int8(bits(i+1));
end

%% ════════════════════════════════════════════════════════════════════════
%%  MODULATION
%% ════════════════════════════════════════════════════════════════════════
function s = ht_modulate(bits, Nbpscs)
    b=double(bits(:));
    switch Nbpscs
        case 1, s=1-2*b;
        case 2, s=(1-2*b(1:2:end)+1j*(1-2*b(2:2:end)))/sqrt(2);
        case 4, s=qam_map(b,4,[-3 -1 3 1],sqrt(10));
        case 6, s=qam_map(b,6,[-7 -5 -1 -3 7 5 1 3],sqrt(42));
        otherwise, error('Unsupported Nbpscs.');
    end
    s=s(:);
end

function bits = ht_demodulate(sym, Nbpscs)
    N=numel(sym); bits=zeros(N*Nbpscs,1,'int8');
    switch Nbpscs
        case 1, bits=int8(real(sym)<0);
        case 2
            z=sym*sqrt(2);
            bits(1:2:end)=int8(real(z)<0);
            bits(2:2:end)=int8(imag(z)<0);
        case 4
            z=sym*sqrt(10);
            for k=1:N, bits((k-1)*4+(1:4))=int8([slice16(real(z(k))) slice16(imag(z(k)))]); end
        case 6
            z=sym*sqrt(42);
            for k=1:N, bits((k-1)*6+(1:6))=int8([slice64(real(z(k))) slice64(imag(z(k)))]); end
    end
end

function s = qam_map(bits, bpSym, map, scale)
    N=numel(bits)/bpSym; s=zeros(N,1);
    for k=1:N
        b=bits((k-1)*bpSym+(1:bpSym));
        if bpSym==4
            s(k)=(map(b(1)*2+b(2)+1)+1j*map(b(3)*2+b(4)+1))/scale;
        else
            s(k)=(map(b(1)*4+b(2)*2+b(3)+1)+1j*map(b(4)*4+b(5)*2+b(6)+1))/scale;
        end
    end
end

function b = slice16(v)
    if v<-2, b=[0 0]; elseif v<0, b=[0 1]; elseif v<2, b=[1 1]; else, b=[1 0]; end
end

function b = slice64(v)
    pts=[-7 -5 -3 -1 1 3 5 7];
    gray=[0 0 0;0 0 1;0 1 1;0 1 0;1 1 0;1 1 1;1 0 1;1 0 0];
    [~,idx]=min(abs(v-pts)); b=gray(idx,:);
end

%% ════════════════════════════════════════════════════════════════════════
%%  SCRAMBLER
%% ════════════════════════════════════════════════════════════════════════
function out = ht_scramble(bits, seed)
    s=seed(:).'; out=bits;
    for k=1:numel(bits)
        fb=xor(s(7),s(4));
        out(k)=int8(xor(double(bits(k)),fb));
        s=[fb s(1:6)];
    end
end

function out = ht_descramble(bits, seed)   % XOR-based: same as scramble
    out = ht_scramble(bits, seed);
end

%% ════════════════════════════════════════════════════════════════════════
%%  VITERBI DECODER
%% ════════════════════════════════════════════════════════════════════════
function out = viterbi_decode(bits, Ntail)
    persistent T
    if isempty(T), T=viterbi_trellis(); end
    if mod(numel(bits),2), bits(end+1)=0; end
    N=numel(bits)/2; INF=1e9;
    pm=INF*ones(T.Ns,1); pm(1)=0;
    prevS=zeros(T.Ns,N,'uint8'); prevB=zeros(T.Ns,N,'uint8');

    for t=1:N
        r0=double(bits(2*t-1)); r1=double(bits(2*t));
        npm=INF*ones(T.Ns,1);
        for ns=1:T.Ns
            s0=double(T.PrevState(ns,1))+1; s1=double(T.PrevState(ns,2))+1;
            m0=pm(s0)+bm(T.Out(ns,1,1),r0)+bm(T.Out(ns,1,2),r1);
            m1=pm(s1)+bm(T.Out(ns,2,1),r0)+bm(T.Out(ns,2,2),r1);
            if m0<=m1
                npm(ns)=m0; prevS(ns,t)=T.PrevState(ns,1); prevB(ns,t)=T.PrevBit(ns,1);
            else
                npm(ns)=m1; prevS(ns,t)=T.PrevState(ns,2); prevB(ns,t)=T.PrevBit(ns,2);
            end
        end
        pm=npm;
    end

    st=1; out=zeros(N,1,'int8');
    for t=N:-1:1, out(t)=prevB(st,t); st=double(prevS(st,t))+1; end
    out=out(1:end-Ntail);
end

function m = bm(e,r)   % branch metric (erasure-aware)
    if r<0, m=0; else, m=double(e~=r); end
end

function T = viterbi_trellis()
    g0=[1 0 1 1 0 1 1]; g1=[1 1 1 1 0 0 1]; Ns=64;
    T=struct('Ns',Ns,'PrevState',zeros(Ns,2,'uint8'), ...
             'PrevBit',zeros(Ns,2,'uint8'),'Out',zeros(Ns,2,2,'uint8'));
    fill=zeros(Ns,1,'uint8');
    for s=0:Ns-1
        for b=0:1
            r=double([b bitget(uint8(s),6:-1:1)]);
            ns=uint8(bi2de(r(1:6),'left-msb'));
            sl=double(fill(ns+1))+1;
            T.PrevState(ns+1,sl)=uint8(s);
            T.PrevBit(ns+1,sl)=uint8(b);
            T.Out(ns+1,sl,:)=uint8([mod(sum(r.*g0),2) mod(sum(r.*g1),2)]);
            fill(ns+1)=uint8(sl);
        end
    end
end

%% ════════════════════════════════════════════════════════════════════════
%%  UTILITIES
%% ════════════════════════════════════════════════════════════════════════
function idx = sub2fft(k, NFFT)
    idx = mod(k,NFFT)+1;
end

function b = u2b(v, n)
    b=zeros(1,n);
    for k=1:n, b(k)=bitget(uint32(v),k); end
end

function crc = crc8_htsig(bits)
    poly=[1 0 0 0 0 0 1 1 1];
    reg=[logical(bits(:).') false(1,8)];
    for k=1:numel(bits)
        if reg(k), reg(k:k+8)=xor(reg(k:k+8),logical(poly)); end
    end
    crc=double(reg(end-7:end));
end