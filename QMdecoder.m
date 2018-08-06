function Seq = QMdecoder(B_out, len)
load('Prob_Est_Table1.mat');
% Initdec:-
Adc = 65536;
Cdc = 0;
S = 0;
MPS = 0;
Qe = PE_Table(S+1,2);
BP = 0;
%Seq_dcr = [];
%x = length(Seq_dcr);
% Byte_in:
[BP, Cdc] = Byte_in(BP, Cdc, B_out);
Cdc = bitshift(Cdc,8);
% Byte_in:
[BP, Cdc] = Byte_in(BP, Cdc, B_out);
Cdc = bitshift(Cdc,8);
CT = 0;
for i = 1:len
    %%% Decode(S):
    Adc = Adc - Qe;
    Cx = bitshift(Cdc,-16);
    Clow = bitand(Cdc,hex2dec('0000FFFF'));
    if(Cx<Adc)
        if(Adc<32768)
            % Conditional MPS exchange:
            if(Adc<Qe)
                D = 1 - MPS;
                
                % Estimate Qe after LPS;
                if(PE_Table(S+1,5)==1)
                    MPS = 1 - MPS;
                end
                S = PE_Table(S+1,3);
                Qe = PE_Table(S+1,2);
                
            else
                D = MPS;
                
                % Estimate Qe after MPS:
                S = PE_Table(S+1,4);
                Qe = PE_Table(S+1,2);
                
            end
            % Renorm_d:
            while(Adc<32768)
                if(CT==0)
                    % Byte in:
                    [BP, Cdc] = Byte_in(BP, Cdc, B_out);                    
                    CT = 8;
                end
                
                Adc = Adc*2;
                Cdc = Cdc*2;
                CT = CT - 1;
            end
            
        else
            D = MPS;
        end
    else
        % Conditional LPS exchange:
        if(Adc<Qe)
            D = MPS;
            Cx = Cx - Adc;
            Adc = Qe;
            % Estimate Qe after MPS:
            S = PE_Table(S+1,4);
            Qe = PE_Table(S+1,2);
            
        else
            D = 1 - MPS;
            Cx = Cx - Adc;
            Adc = Qe;
            % Estimate Qe after LPS:
            if(PE_Table(S+1,5)==1)
                MPS = 1 - MPS;
            end
            S = PE_Table(S+1,3);
            Qe = PE_Table(S+1,2);
            
        end
        Cdc = bitshift(Cx,16) + Clow;
        % Renorm_d:
        while(Adc<32768)
            if(CT==0)
                % Byte in:                
                [BP, Cdc] = Byte_in(BP, Cdc, B_out);
                CT = 8;
            end
            
            Adc = Adc*2;
            Cdc = Cdc*2;
            CT = CT - 1;
        end
        
    end
    %Seq_dcr = [Seq_dcr num2str(D)];
    Seq(i) = num2str(D);
    %x = length(Seq_dcr);
end
end