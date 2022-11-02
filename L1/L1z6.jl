# convert number fron binary (IEE 754 standard) to decimal
function BitsToFloat64(s)
    (rs, rc, rm) = (s[1], s[12:-1:2], s[13:end])
    s = nothing
    c = 0 
    m = 1.0
    
    if rs == 1 s = -1 else s = 1 end
    
    for (cnt, b) in enumerate(rc)
        c += (Int16(b)-48)*2^(cnt-1)
    end
   
    for (cnt, b) in enumerate(rm)
        m += (Int16(b)-48)*2.0^(-cnt)
    end
    c = c-1023
  
    return s*2.0^c*m  
end 
