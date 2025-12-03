# EECE 442 – Digital Communication System Simulation  
Full Physical-Layer Communication Chain (MATLAB)

This repository contains the complete implementation of the EECE 442 course project at the American University of Beirut.  
The work covers both **Part I** (Sampling, Quantization, Source Coding) and **Part II** (Channel Coding, Modulation, AWGN/ISI Channels, Detection).

The implementation follows the official project descriptions:  
- Part I Requirements: Sampling, Quantization, Source Coding  
- Part II Requirements: Channel Coding, Modulation, Equalization, BER Analysis  
- Full Final Report: Figures, results, and analysis

---

## 📌 Overview

The project constructs a **complete digital communication simulation chain**, starting from a real signal and ending with reconstructed data after channel transmission.

The chain includes:
1. Message Sampling  
2. Quantization (Uniform, Lloyd–Max, Two-level)  
3. Lossless Source Coding (Fixed-length, Huffman, Block Coding)  
4. Channel Coding (Repetition)  
5. Modulation (BPSK, QPSK)  
6. AWGN Channel & ISI Channel  
7. Equalization and Detection (Hard Decision, MLSE/Viterbi)  
8. End-to-End BER, SNR, Throughput, and Rate–Distortion Evaluation  

All modules are connected through **`full_chain.m`**.

---

## 📁 Repository Structure

```
.
├── Sampling/
├── Quantization/
├── Coding/
├── ChannelCoding/
├── Modulation/
├── Channel/
├── full_chain.m
├── report.pdf
```

---

## ⚙️ Module Descriptions

### **1. Sampling**
- Uniform time-domain sampling  
- Reconstruction from samples  
- Fourier-series sampling  

---

### **2. Quantization**
- Two-level quantizer  
- Uniform M-level quantization  
- Lloyd–Max quantizer  
- MSE & thresholds/levels analysis  

---

### **3. Source Coding**
- Fixed-length coding  
- Huffman coding  
- Block Huffman coding   

---

### **4. Channel Coding**
- Repetition L = 3  
- Majority-vote decoding  
- Coded vs. uncoded BER  

---

### **5. Modulation**
- BPSK  
- QPSK (Gray mapping)  

---

### **6. Channels & Detection**
- AWGN  
- ISI channel  
- MLSE (Viterbi)   

---

### **7. End-to-End Metrics**
- BER curves  
- Bits/symbol  
- Quantization SNR  
- Modulation performance  

---

## ▶ Running the Code

Run full chain:

```matlab
full_chain
```

Run test modules:

```matlab
test_sampling
test_quantization
test_source_coding
test_modulation
test_channel_awgn
test_mlse
```

---


## 👥 Authors
- **Tony Abi Haidar**
- **Kassem Hassoun**

American University of Beirut – EECE 442: Communication Systems
Professor: Jihad Fahs  
Fall 2025–2026
