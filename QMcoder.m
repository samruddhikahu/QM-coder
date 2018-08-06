function [B_out] = QMcoder(Seq)
load('Prob_Est_Table1.mat');
% Initializing the encoder:
Adc = 65536;
Cdc = 0;
CT = 11;
S = 0;
ST = 0;
MPS = 0;
Qe = PE_Table(S+1,2);
B_out = [];
Cs = [CT, ST];
for i = 1:length(Seq)
    D = bin2dec(Seq(i));
    if(D==MPS)
        % Code_MPS(S);
        Adc = Adc - Qe;
        if(Adc<32768)
            if(Adc<Qe)
                Cdc = Cdc + Adc;
                Adc = Qe;
            end
            % Estimate Qe after MPS:
            S = PE_Table(S+1,4);
            Qe = PE_Table(S+1,2);
            % Renormalization:
            [Adc, Cdc, Cs, B_out] = Renorm_e(Adc, Cdc, Cs, B_out);
        end
    else
        % Code_LPS(S):
        Adc = Adc - Qe;
        if(Adc>=Qe)
            Cdc = Cdc + Adc;
            Adc = Qe;
        end
        
        % Estimate Qe after LPS;
        if(PE_Table(S+1,5)==1)
            MPS = 1 - MPS;
        end
        S = PE_Table(S+1,3);
        Qe = PE_Table(S+1,2);
        
        % Renormalization:
        [Adc, Cdc, Cs, B_out] = Renorm_e(Adc, Cdc, Cs, B_out);
    end
    %B_dc = [B_dc B_out];
end

%Flush Procedure:
CT = Cs(1,1);
ST = Cs(1,2);
T = Cdc + Adc - 1;
T = bitand(T,hex2dec('FFFF0000'));
if(T<Cdc)
    T = T + 32768;
end
C = T;
clear T;
while(C>0)
    C = bitshift(C,CT);
    % Byte Out Procedure:
    T = bitshift(C,-19);
    if(T>255)
        %%%%
        B_out(end) = B_out(end) + 1;
        %%% Stuff zeros:
        if(B_out(end)==255)
            B_out = [B_out 0];
        end
        %%% Output Stacked Zeros:
        while(ST>0)
            B_out = [B_out 0];
            ST = ST - 1;
        end
        T = bitand(T,255);
        B_out = [B_out T];
    else
        if(T==255)
            ST = ST + 1;
        else
            %%% Output Stacked FFs;
            while(ST>0)
                B_out = [B_out 255 0];
                ST = ST - 1;
            end
            B_out = [B_out T];
        end
    end
    C = bitand(C,524287);
    CT = 8;
end 
end