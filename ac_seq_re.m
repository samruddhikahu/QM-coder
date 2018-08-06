function Ir = ac_seq_re(Seq, H, W)
nB = (H*W)/64;
p = 0;
for j = 1:nB
    k = 1;
    while(k<=63)
        p = p + 1;
        D = str2num(Seq(p));
        if(D==1)
            Ir(j,k:63) = 0;
            break;
        else
            p = p + 1;
            D = str2num(Seq(p));
            while(D==0)
                Ir(j,k) = 0;
                k = k + 1;
                p = p + 1;
                D = str2num(Seq(p));
            end
            % Decode V:
            % Decode Sign of V:
            p = p + 1;
            D = str2num(Seq(p));
            if(D==1)
                Sgn = 'neg';
            else
                Sgn = 'pos';
            end
            % Decode log2(Sz):
            n = 0;
            p = p + 1;
            D = str2num(Seq(p));
            while(D)
                n = n + 1;
                p = p + 1;
                D = str2num(Seq(p));
            end
            n = n - 1;
            if(n==0)
                Sz = 1;
            elseif(n>0)
                %base = 2^n;
                Szr = bin2dec(Seq(p+1:p+n));
                Sz = Szr + (2^n);
                p = p + n;
            else
                Sz = 0;
            end
            abV = Sz + 1;
            if(Sgn=='pos')
                V = abV;
            else
                V = -abV;
            end
            
            Ir(j,k) = V;
            k = k + 1;
        end
    end
            
end
end