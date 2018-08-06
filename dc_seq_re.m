function dc = dc_seq_re(Seq, H, W)
dc_prev = 0;
p = 1;
nB = (H*W)/(8*8);
dc = zeros(1,nB);
for i = 1:nB
   D = str2num(Seq(p));
   if(D==0)
       dc_diff = 0;
   else
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
         abdiff = Sz + 1;
         if(Sgn=='pos')
             dc_diff = abdiff;
         else
             dc_diff = -abdiff;
         end
             
       % dc_diff = V;
   end
   dc(i) = dc_diff + dc_prev;
   dc_prev = dc(i);
   p = p + 1;
end

end