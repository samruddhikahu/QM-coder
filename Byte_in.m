function [BP, Cdc] = Byte_in(BP, Cdc, B_out)
%global B_out;
BP = BP + 1;
if(BP>length(B_out))
    B = 0;
else
    B = B_out(BP);
end
if(B==255)
    % Unstuff zeros:
    BP = BP + 1;
    if(B_out(BP)==0)
        Cdc = bitor(Cdc,hex2dec('FF00'));
        %BP = BP - 1;
    else
        % End of code segment
    end
else
    Cdc = Cdc + bitshift(B,8);
end

end