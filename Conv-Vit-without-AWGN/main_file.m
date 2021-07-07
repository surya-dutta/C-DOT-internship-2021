%%MAIN FILE%%
%Initialize Simulation Parameters
message=[1 0 1 0 1]; % message
R = 1/2;                 % code rate
dbs = -1:10;             % SNR per bit, in dB
trials = 100;           % number of trials to perform
n=2;                     % Convolution Code parameter
k=length(message);       % Length of the message             
mem = 3;                 % Number of memory elements
%Initialize and compute Shannon Limit/Uncoded Efficiency - taken from paper
errs   = 0*dbs;
EbN0 = 10.^(dbs/10);
sigs = 1./sqrt(2*R*EbN0);
ber0 = logspace(-6,-2.1,81);
ber1 = logspace(-6,-0.99,81);
db0  = 10*log10((2.^(2*R*(1+log2((ber0.^ber0).*(1-ber0).^(1-ber0))))-1)/(2*R));
db1  = 20*log10(erfinv(1-2*ber1));
for trial = 1:trials
    %--------- ENCODER: 1/2 Convolution Encoder -------%
    c = convenc_1_2(message);
    %--------- DECODER: 1/2 Viterbi Decoder -----------%
    dec_op = viterbidec_1_2(c); %%%generalised
    %--------------------------------------------------%
    %%%PLOT SIMULATED RESULTS
    ber = errs/(k*trial); %%%Bit error rate is 0, because no errors/noise are introduced
    %[trial, errs]
    %semilogy(dbs, ber,'o-', db0, ber0,':', db1, ber1,':');
    %hold off;
    %xlabel('SNR per bit, E_b / N_0 (dB)');
    %ylabel('Bit-Error Rate');
    %axis([-1 10 1e-6 1])
    %title(['After ',num2str(trial),' trials (',num2str(k*trial),' msg bits)']);
    %grid on;
    %drawnow;
    %%%Graph doesnt show any error
end
