FUNCTION create_2ds_frame, x, slicecounter, pcounter, array=array
   ;Build frame
   ndiodes=fix(128)
   id=12883s              ;Word 1, '2S'
   nh=0s                  ;This will update as particle is compressed, assume no overload, empty fifo, etc.
   nv=0s                  
   nn=0s                  ;A generic counter, will be set to either nh or nv at the end
   pcount=fix(pcounter)   ;Particle counter  (not Number particle events)
   nslices=n_elements(x)/ndiodes
   nslices=fix(nslices)   ;Number slices in particle
   
   cimage=intarr(10000)   ;Initialize cimage
   FOR j=0,nslices-1 DO BEGIN    ;Write each slice, could be one or more words
      slice=x[*,j]
      numclear=0l
      numshaded=0l
      newslice=1 ;Turn on flag for new slice
      FOR idiode=0,ndiodes-1 DO BEGIN   ;Scroll through diodes, counting nshaded/nclear
         
         ;Decide if need to write a new word
         writeword=0
         IF (slice[idiode] eq 0) and (numshaded gt 0) THEN writeword=1   ;New shade/clear pair ready
         IF (slice[idiode] ne 0) and (idiode eq ndiodes-1) THEN writeword=1  ;End of slice with shaded diode
         IF (slice[idiode] eq 0) and (numclear eq ndiodes-1) THEN writeword=1 ;Empty slice
         
         IF (writeword eq 1) THEN BEGIN   
            ;numshaded not updated in this case, do it here
            IF (slice[idiode] ne 0) and (idiode eq ndiodes-1) THEN numshaded=numshaded+1
            
            f=0s   ;This word
            f=(f or numclear)           ;add in numclear
            f=(f or ishft(numshaded,7)) ;add in numshaded at bits 7-13
            f=(f or ishft(newslice,14)) ;flag new slice
            IF (numclear eq ndiodes-1) and (slice[idiode] eq 0) THEN f='7fff'x  ;Blank slice special case
            IF numshaded eq ndiodes THEN f='4000'x ;Full slice special case
            ;print,newslice, numclear, numshaded                       
            cimage[nn]=f
            nn=nn+1    ;keep count of written words
            
            numclear=0l  ;Reset counters
            numshaded=0l
            newslice=0   ;turn off this flag until get to next slice                     
         ENDIF

         ;Increment nclear/nshaded here
         IF slice[idiode] eq 0 THEN numclear=numclear+1 ELSE numshaded=numshaded+1
      ENDFOR  ;diodes
   ENDFOR   ;slices
   ;print,'-------------------------------------------'
   f=ishft(1,15)        ;Timing word, with bit 15 set to flag counter
   ;fullcounter=long((time+inttime)*tas/res)  ;Total slice count.  Use t+inttime for greater (albeit false) accuracy.
   ;Take care of very short interarrivals from blank slices in an image
   ;IF iframe gt 0 THEN fullcounter=fullcounter+slicestart
   f=(ishft(slicecounter,-16) and 'ffff'x)  ;Write timing word 1, bits 16-31
   cimage[nn]=f
   nn=nn+1
   f=(slicecounter and 'ffff'x)  ;Write timing word 1, bits 0-15
   cimage[nn]=f
   nn=nn+1
   IF array eq 'V' THEN nv=nn ELSE nh=nn
   frame=[id,nh,nv,pcount,nslices,cimage[0:nn-1]] 
   return, frame
END
