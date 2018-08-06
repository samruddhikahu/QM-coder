function Seq = ac_seq(I)
[s s1] = size(I);
Seq = [];
for j = 1:s
    clear If;
    clear cnt;
    clear Int;
    clear last;
    count = 0;
    z = 1;
    for i = 2:64
        if(I(j,i)==0)
            count = count + 1;
            last(z) = i;
            z = z + 1;
        else
            count = 0;
            z = 1;
        end
    end
        
    If = I(j,:);
    if(count>0)
        %If(last(1)) = 10;
        if(last(1)<=64)
            If(last(1):64) = [];
        end
    end
    Ig = If(:,2:end);
    len = length(Ig);
    k = 1;
    while(k<=63)
        if(k>len)
            % Code 1:
            Seq = [Seq '1'];
            break;
        else
            % Code 0:
            Seq = [Seq '0'];
            
            V = Ig(k);
            while(V==0)
                % Code 0:
                Seq = [Seq '0'];
                k = k + 1;
                  V = Ig(k);
            end
            
            % Code 1:
            Seq = [Seq '1'];
            % Encode V:
            if(V<0)
                Seq = [Seq '1'];
            else
                Seq = [Seq '0'];
            end
            Sz = abs(V) - 1;
            n = floor(log2(Sz));
            for i = 1:n+1
                Seq = [Seq '1'];
            end
            Seq = [Seq '0'];
            if(n>0)
                Szr = Sz - (2^n);
                Seq = [Seq dec2bin(Szr,n)];
            end
        end
        k = k + 1;
    end
end
end