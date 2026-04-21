# TX-RX Chain Notes

This note summarizes the manual IEEE 802.11n HT transmitter-receiver chain used in this repository.

## Main Flow

1. The transmitter starts with PSDU bits.
2. The bit stream is scrambled to randomize long repetitive patterns.
3. BCC encoding and puncturing apply channel coding with the selected coding rate.
4. Interleaving distributes coded bits to reduce the impact of burst errors.
5. QAM mapping converts grouped bits into complex modulation symbols.
6. Data and pilot symbols are placed onto OFDM subcarriers.
7. IFFT and cyclic prefix generation transform the signal into the time domain.
8. The packet passes through the channel and is affected by noise, represented with SNR-dependent BER behavior.
9. The receiver detects the packet, estimates the channel using HT-LTF symbols, equalizes the received signal, and recovers the transmitted bits.
10. Demapping, deinterleaving, depuncturing, Viterbi decoding, and descrambling reconstruct the final PSDU.

## Why This File Matters

The MATLAB code shows the implementation details, but this note gives a fast conceptual reference for the order and purpose of the main PHY blocks. It is intended to make the manual TX-RX chain easier to read inside the repository.
