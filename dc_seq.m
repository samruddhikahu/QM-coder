function Seq = dc_seq(I)
dc_prev = 0;
Seq = [];
for j = 1:length(I)
    diff = I(j) - dc_prev;
    dc_prev = I(j);
    if(diff==0)
        Seq = [Seq '0'];
    else
        Seq = [Seq '1'];
        % Encode V;
        if(diff<0)
            Seq = [Seq '1'];
        else
            Seq = [Seq '0'];
        end
        Sz = abs(diff) - 1;
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
end
end