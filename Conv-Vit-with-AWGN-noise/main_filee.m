%%%MAIN FILE%%%
clc;clear;
%%%Convolution encoder and soft decision viterbi decoder%%%
%Initialize Simulation Parameters
L = 10;               % message length
R = 1/2;              % code rate
dbs = -1:10;          % SNR per bit, in dB
trials = 1e3;         % number of trials to perform
g = [1 1 1;1 0 1];    % generator sequence
n = length(g);        % convolution Code (1/n) parameter
memory_els = 3;       % number of memory elements
%Initialize and compute Shannon Limit/Uncoded Efficiency
errs   = 0*dbs;
EbN0 = 10.^(dbs/10);
sigs = 1./sqrt(2*R*EbN0);
%ber0 = logspace(-6,-2.1,81);
%ber1 = logspace(-6,-0.99,81);
%db0  = 10*log10((2.^(2*R*(1+log2((ber0.^ber0).*(1-ber0).^(1-ber0))))-1)/(2*R));
%db1  = 20*log10(erfinv(1-2*ber1));
for trial = 1:trials
    m = round(rand(1,L));       % message signal
    %ENCODER: Convolution Encoder
    c = convolutional_coding(m,g);
    %NOISE CHANNEL - Noise is directly added to the encoded sequence%
    noise = randn(1,length(c));
        for i=1:length(dbs)
        r = 2*c - 1 + sigs(i)*noise;
        %DECODER: Viterbi Decoder
        mhat = viterbi_decoder(r,g,1);
        errs(i) = errs(i) + sum(mhat~= m); 
        end
    %Plot Simulated Result
    ber = errs/(L*trial);
    [trial, errs];
    semilogy(dbs, ber,'o-');
    hold off;
    xlabel('SNR per bit, E_b / N_0 (dB)');
    ylabel('Bit-Error Rate');
    axis([-1 10 1e-6 1])
    title(['After ',num2str(trial),' trials (',num2str(L*trial),' msg bits)']);
    grid on;
    drawnow;
end
