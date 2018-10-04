FUNCTION encode_header_slice, hour, minute, second, millisecond, counter, nslices, pcounter, dof, grey=grey, tas=tas
   ;Encode the DMT header slice
   ;Output will be 8 bytes
   ;AB 7/20216
   
   IF n_elements(grey) eq 0 THEN grey=0
   IF n_elements(tas) eq 0 THEN tas=100
   
   IF grey eq 0 THEN BEGIN
      bslice=bytarr(64)
      bslice[0:6]=(dec2bin8(nslices))[1:7]
      bslice[7]=dof
      bslice[8:12]=(dec2bin8(hour))[3:7]
      bslice[13:18]=(dec2bin8(minute))[2:7]
      bslice[19:24]=(dec2bin8(second))[2:7]
      bslice[25:34]=(dec2bin16(millisecond))[6:15]
      bslice[35:47]=(dec2bin16(counter))[3:15]
      bslice[48:63]=(dec2bin16(pcounter))[0:15]

      out=bytarr(8)  ;Convert 64 bits to 8 bytes
      FOR i=0,7 DO out[i]=bin2dec(bslice[i*8:(i+1)*8-1])
      return,reverse(out)  ;For some reason this is reversed in DMT mono data
     
   ENDIF ELSE BEGIN
      bslice=bytarr(128)
      bslice[120:127]=reverse((dec2bin8(nslices-1))[0:7])  ;The grey probe does not count the trailing slice, so subtract one here
      bslice[115:119]=reverse((dec2bin8(hour))[3:7])
      bslice[109:114]=reverse((dec2bin8(minute))[2:7])
      bslice[103:108]=reverse((dec2bin8(second))[2:7])
      bslice[93:102]=reverse((dec2bin16(millisecond))[6:15])
      ;This probe uses microseconds + 125ns counter rather than larger 125ns counter for mono
      ;Compute that here
      nanosecond=counter * 125d
      microsecond=floor(nanosecond/1000)
      newcounter=counter mod 8  ;Only need last 3 remaining bits (max counter = 8)
      bslice[83:92]=reverse((dec2bin16(microsecond))[6:15])
      bslice[80:82]=reverse((dec2bin16(newcounter))[13:15])
      bslice[64:79]=reverse((dec2bin16(pcounter))[0:15])
      bslice[56:63]=reverse((dec2bin8(tas))[0:7])
      ;bslice[0:55] is all zero

      out=bytarr(64)  ;Convert 128 bits to 64 grey mode (0-3) bytes
      FOR i=0,63 DO out[i]=bslice[i*2]*1 + 2*bslice[i*2+1] 
      return,out
   ENDELSE
END
