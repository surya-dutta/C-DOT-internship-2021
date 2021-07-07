EbNo = 1:1:10;
numLLRRuns = 5;

% Modulation properties
M = 4;
k = log2(M);

% Create a rate 1/2, constraint length 7 comm.ConvolutionalEncoder System
% object. This encoder takes one-bit symbols as inputs and generates 2-bit
% symbols as outputs.
codeRate = 1/2;
constlen = 7;
codegenpoly = [171 133];    
trellis = poly2trellis(constlen, codegenpoly);
enc = comm.ConvolutionalEncoder(trellis);
dSpect = distspec(trellis,14);

% Create a comm.QPSKModulator System object and two comm.QPSKDemodulator
% System objects one each for hard decision and LLR demodulation.
qpskMod = comm.QPSKModulator('BitInput',true);
qpskDemodHD = comm.QPSKDemodulator('BitOutput',true,...
    'DecisionMethod', 'Hard decision');
qpskDemodLLR = comm.QPSKDemodulator('BitOutput',true,...
    'DecisionMethod', 'Log-likelihood ratio');

% Create a comm.AWGNChannel System object. 
chan = comm.AWGNChannel('NoiseMethod',...
    'Signal to noise ratio (Eb/No)',...
    'SignalPower', 1,...
    'SamplesPerSymbol', 1);

% Adjust signal-to-noise ratio for coded bits and multi-bit symbols.
adjSNR = EbNo - 10*log10(1/codeRate) + 10*log10(k);

% Compute Noise Variance which is required by the demodulator in  LLR mode
NoiseVariance = 10.^(-adjSNR/10);

% Create comm.ViterbiDecoder System objects to act as the hard-decision,
% soft-decision decoders.

vitDecHD  = comm.ViterbiDecoder(trellis, 'InputFormat', 'Hard',...
    'TracebackDepth', 32);

vitDecSD  = comm.ViterbiDecoder(trellis, 'InputFormat', 'Soft',...
  'SoftInputWordLength',3, 'TracebackDepth', 32); 

% Create comm.ErrorRate System objects to compare the decoded bits to the
% original transmitted bits. The Viterbi decoder creates a delay in the
% output decoded bit stream equal to the traceback length. To account for
% this delay set the 'ReceiveDelay' property of the comm.ErrorRate objects
% to 32.
errorCalcHD  = comm.ErrorRate('ReceiveDelay', 32);
errorCalcSD  = comm.ErrorRate('ReceiveDelay', 32);

% Before using a comm.ViterbiDecoder object in the 'soft decision' mode,
% the output of the comm.QPSKDemodulator object needs to be quantized.
scalQuant = dsp.ScalarQuantizerEncoder;
scalQuant.Partitioning = 'Unbounded';


% Since the AWGN Channel as well as the RANDI function use the default
% random stream, the following commands are executed so that the results
% will be repeatable
s0 = RandStream.getGlobalStream;
s = RandStream.create('mt19937ar', 'seed',12345);
RandStream.setGlobalStream(s);


% Number of bits per iteration
bitsPerIter = 1.2e4;
% Maximum number of iterations
maxNumIters = 100;
% Maximum number of bit errors to collect
maxNumErrs  = 300;

% Initialize variables for storing BER results
ber_HD = zeros(3,length(EbNo));
ber_SD = zeros(3,numLLRRuns); 

% Set up a figure for visualizing BER results
fig = figure;
grid on;
ax = fig.CurrentAxes;
hold(ax,'on');
ax.YScale = 'log';
xlim(ax, [EbNo(1)-1, EbNo(end)+1]); ylim(ax, [1e-6 1]);
xlabel(ax, 'Eb/No (dB)'); ylabel(ax, 'BER');
title(ax, 'Soft vs. Hard Decision Demodulation');
fig.NumberTitle = 'off';
fig.Renderer = 'zbuffer';
set(fig, 'DefaultLegendAutoUpdate', 'off');

% Use 'bercoding' to calculate theoretical upper bounds on BER
theoryBER_HD = bercoding(adjSNR, 'conv', 'hard', codeRate, dSpect);
theoryBER_Soft = bercoding(adjSNR(1:numLLRRuns), 'conv', 'soft', codeRate, dSpect);

semilogy(EbNo, theoryBER_HD, 'mo-', EbNo(1:numLLRRuns), theoryBER_Soft, 'go-');
legend('Hard Decision: Theoretical Upper Bound','Soft Decision: Theoretical Upper Bound', ...
       'Location', 'SouthWest');
 
% This is the common loop for hard decision and LLR with unquantized and
% soft decision decoding
for idx=1:numLLRRuns    

    % Reset objects with state at the start of an EbNo value
    reset(errorCalcHD)
    reset(errorCalcSD)
    reset(enc)
    reset(vitDecHD)
    reset(vitDecSD)
    iter=1;
    
    % set the channel EbNo value for simulation
    chan.EbNo = adjSNR(idx);
    
    % set the NoiseVariance property, required for LLR demodulation
    qpskDemodLLR.Variance = NoiseVariance(idx);
    
    % Fine tune the quantizer range according to the noise SNR
    scalQuant.BoundaryPoints = (-2.1:.7:2.1)/NoiseVariance(idx);
  
    % Exit loop when either the number of bit errors exceeds 'maxNumErrs'
    % or the maximum number of iterations have completed
    while (iter <= maxNumIters)
        
        data = randi([0 1], bitsPerIter, 1);        % Generate message bits        
        encData = enc(data);                        % Convolutionally encode the data            
        modData = qpskMod(encData);                 % Modulate the encoded data        
        channelOutput = chan(modData);              % Pass the modulated signal
                                                    % through an AWGN channel     
        demodDataHD = qpskDemodHD(channelOutput);   % Hard decision Demodulation        
        demodDataLLR = qpskDemodLLR(channelOutput); % 'LLR' Demodulation
                
        % Hard-decision decoding: Pass the demodulated data through the
        % Viterbi decoder; and compute and accumulate errors     
        decDataHD = vitDecHD(demodDataHD);
        ber_HD(:,idx) = errorCalcHD(data, decDataHD);
        
        % Soft-decision decoding: The demodulated data must pass through a
        % quantizer before being fed to the Viterbi decoder. However the
        % output from the demodulator must be sign-reversed before being
        % fed to the quantizer. This is because in the soft-decision
        % decoding mode, the comm.ViterbiDecoder object assumes that
        % positive numbers correspond to 1s and negative numbers to 0s.
        % Thus the decoding operation consists of feeding in the sign
        % reversed data from the comm.PSKDemodulator object to the
        % dsp.scalarQuantizerEncoder object and feeding in the output from
        % this dsp.scalarQuantizerEncoder object to the
        % comm.ViterbiDecoder. Compute and accumulate errors.
        quantizedValue = scalQuant(-demodDataLLR);
        decDataSD = vitDecSD(double(quantizedValue));
        ber_SD(:,idx) = errorCalcSD(data, decDataSD);
        
        iter = iter+1;
    end
    
    % Plot results
    semilogy(ax, EbNo(1:idx), ber_HD(1,1:idx), 'r*', ...
             EbNo(1:numLLRRuns),ber_SD(1,1:numLLRRuns), 'k*');
    legend('Hard Decision: Theoretical Upper Bound', ...
           'Soft Decision: Theoretical Upper Bound', ...
           'Hard Decision: Simulation' , ...
           'Soft Decision: Simulation',...
           'Location', 'SouthWest');
    drawnow;
end

% Perform curve fitting and plot the results
fitBER_HD  = berfit(EbNo, ber_HD(1,:));
fitBER_SD = berfit(EbNo(1:numLLRRuns), ber_SD(1,:));
semilogy(ax, EbNo(1:numLLRRuns), fitBER_HD, 'r*-', EbNo(1:numLLRRuns), fitBER_SD, 'b*-')
hold(ax,'off');

% Restore default stream
RandStream.setGlobalStream(s0);