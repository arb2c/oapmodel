FUNCTION create_dmt_frame, x, time, pcounter, grey=grey, rlehbflag=rlehbflag, tas=tas, dofflag=dofflag
   ;Create DMT formatted particle, including timing data and the image
   ;AB 7/2016
   IF n_elements(grey) eq 0 THEN grey=0
   IF n_elements(tas) eq 0 THEN tas=100
;IF time gt 2.09 then stop

;print,x,format='(64i1)'
   ;Build frame
   ndiodes=fix(64)
   nslices=n_elements(x)/ndiodes
   nslices=fix(nslices)+1   ;Number slices in particle, add one for the blank trailing slice that will be added
   
   ;Time in DMT form
   sfm=long(time)
   hour=sfm/3600l
   minute=(sfm mod 3600)/60
   second=(sfm mod 3600) mod 60
   millisecond=fix((time-sfm)*1000)
   counter=fix(((time*1000) mod 1) * 1e6 / 125d)  ;125ns increments
   header=encode_header_slice(hour, minute, second, millisecond, counter, nslices, pcounter, dofflag, tas=tas, grey=grey)

   cimage=bytarr(10000)   ;Initialize compressed image
   nn=0        ;Position in cimage
   
   IF grey eq 0 THEN BEGIN
      ;----------Mono probes-----------
      byteform=bytarr(n_elements(x)/8) + 8  ;8 extra bytes for a blank slice
      FOR i=0L,n_elements(byteform)-1 DO byteform[i]=not(byte(bin2dec(reverse(x[i*8:i*8+7]))))
      sync=bytarr(8)+170
      ;byteform=[byteform, sync, header]
      byteform=[sync, header, byteform]
      rlehbflag=bytarr(10000)  ;Keeps track of whether a particular compressed byte is a RLEHB
      
      ;Find chunks of 0s, 1s, or random data
      u=[-1,uniq(byteform)]  ;May want to adjust u to ignore single zero or one bytes
      ;u=[-1,2,3,4,55,66,107,171,200,210,450,480,481]   ;Test of u for code below, tests OK.
      ;Check for counts greater than 32 and adjust u   
      count=u[1:*]-u 
      w=where(count gt 32, nw)
      FOR i=nw-1,0,-1 DO BEGIN
         FOR j=(count[w[i]]-1)/32,1,-1  DO BEGIN  ;Have to consider >64, >96, etc.
            u=[u[0:w[i]], u[w[i]]+32*j, u[w[i]+1:*]]  ;Do in reverse order to avoid indexing errors
         ENDFOR
      ENDFOR

      FOR i=1,n_elements(u)-1 DO BEGIN
         chunk=byteform[u[i-1]+1:u[i]]
         count=byte(n_elements(chunk))
         IF count gt 32 THEN stop,'overcount in create_dmt_frame'
         IF chunk[0] eq 0 THEN BEGIN  ;Chunk of zeros
            rlehb=(128b or count-1)
            cimage[nn]=rlehb   ;Write number of zeros
            rlehbflag[nn]=1
            nn=nn+1
         ENDIF
         IF chunk[0] eq 255 THEN BEGIN  ;Chunk of ones
            rlehb=(64b or count-1)
            cimage[nn]=rlehb   ;Write number of ones
            rlehbflag[nn]=1
            nn=nn+1
         ENDIF
         IF (chunk[0] ne 255) and (chunk[0] ne 0) THEN BEGIN  ;Chunk of data
            rlehb=count-1
            cimage[nn]=rlehb   ;Write rlehb
            rlehbflag[nn]=1
            nn=nn+1
            cimage[nn:nn+count-1]=chunk  ;Write data
            nn=nn+count
         ENDIF
      ENDFOR
      rlehbflag=rlehbflag[0:nn-1]
   ENDIF ELSE BEGIN
      ;---------Grey probes-------------
      goodslices=where(total(x,1) gt 0)  ;Only compress the slices with image data, otherwise can get overrunning count values >256
      image=[[header], [3-x[*,goodslices]]]  ;Add header to front of particle, reverse image so dark=0 and no_shadow=3
      ;Change image into a series of values and counts
      ;This works, but need to account for count>127 which is a PITA
      ;u=uniq(image)
      ;value=image[u]          ;Value of a pixel
      ;count=[u[0]+1, u[1:*]-u]  ;Number of times the pixel repeats

      ;Another way of making the count (brute force)
      i=0L
      curval=-1  ;Current value
      value=bytarr(10000)
      count=bytarr(10000)
      c=-1L
      FOR i =0, n_elements(image)-1 DO BEGIN
         IF (image[i] ne curval) or (count[c] eq 127) THEN BEGIN
            ;New value detected
            c++
            curval=image[i]
            value[c]=curval
            count[c]=1
         ENDIF ELSE count[c]++  ;No new value, add to count
      ENDFOR
      count=count[0:c]
      value=value[0:c]
      
      IF max(count) gt 127 then stop,'Need to account for counts > 127 in grey compression' 
      ;Cheap first try, don't bother with 4 and 6 bit options      
      FOR i=0,n_elements(count)-1 DO BEGIN
         valuebyte=4b+byte(value[i])     ;Adding 4 sets the '2-bit' flag
         countbyte=128b+byte(count[i])-1 ;Adding 128 sets the count flag, subtract 1 since the original write counts as first instance
         cimage[nn]=valuebyte
         nn++
         IF count[i] gt 1 THEN BEGIN   ;Write a count if it is greater than 1
            cimage[nn]=countbyte
            nn++
         ENDIF
      ENDFOR
      cimage[nn]=7b     ;Write 11 to indicate no shadow on trailing slices
      nn++
      cimage[nn]=255b   ;Two trailing slices
      nn++
         
      
;print,image,format='(64i1)'
;d=quick_decompress_dmt_grey(cimage)
;s=size(image,/dim)
;print,''
;print,d[*,0:s[1]+2],format='(64i1)'
;stop


;FOR i=0,n_elements(fn) -1 DO convert_oapmodel2,fn[i],grey=1,maxtime=30,format='DMT',/breakframes
;d=decompress_dmt_grey(cimage)
   ENDELSE

;Debugging output, prints original and compressed/decompressed particle
; zz=quick_decompress_dmt(cimage[0:nn-1])
; print,x,format='(64i1)'
; print,'**'
; print,zz,format='(64i1)'
; print,' '
; print,' '
; if total(x-zz[0:n_elements(x)-1]) ne 0 THEN stop

   return, cimage[0:nn-1]
END
