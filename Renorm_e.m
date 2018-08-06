function [A, C, Cs, B_out] = Renorm_e(A, C, Cs, B_out)
% if(CT==11)
%     B_out = [];
% end
CT = Cs(1,1);
ST = Cs(1,2);
while(A<32768)
    A = A*2;
    C = C*2;
    CT = CT - 1;
    if(CT==0)
        % Byte out:
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
    else
        T = [];
        B_out = [B_out T];
    end
end
Cs = [CT ST];
end
