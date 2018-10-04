FUNCTION oap_image, data, ind, makeplot=makeplot, threshold=threshold, compress=compress, clip=clip, datadir=datadir
   ;FUNCTION to read in a series of modeled images from the raw data file.
   ;Enter data structure from the main model output file, and the indices
   ;of desired particles.
   ;AB 5/2012
   
   IF n_elements(makeplot) eq 0 THEN makeplot=0
   IF n_elements(threshold) eq 0 THEN threshold=0
   IF n_elements(compress) eq 0 THEN compress=0     ;Remove blocks of empty slices (only for mulitple particles!)
   IF n_elements(clip) eq 0 THEN clip=0             ;Remove empty leading slices, make sure 1 blank slice at end
   IF n_elements(datadir) eq 0 THEN datadir=''     

   openr,lun,datadir+data.imagefile,/get_lun
   
   ;Just read in the whole file if no index specified
   IF n_elements(ind) eq 0 THEN BEGIN
      im=bytarr(data.op.ndiodes, max(data.slicestop))
      readu,lun,im
   ENDIF ELSE BEGIN
      
      im=bytarr(data.op.ndiodes, (data.slicestop[ind[0]]-data.slicestart[ind[0]])+1)
      point_lun,lun,data.slicestart[ind[0]]*data.op.ndiodes
      readu,lun,im
      FOR k=1L,n_elements(ind)-1 DO BEGIN
         im2=bytarr(data.op.ndiodes, (data.slicestop[ind[k]]-data.slicestart[ind[k]])+1)
         point_lun,lun,data.slicestart[ind[k]]*data.op.ndiodes
         readu,lun,im2
         im=[[im],[im2]]
      ENDFOR
   ENDELSE
   free_lun,lun
   IF threshold gt 0 THEN im=((im > 1)-1)<1  ;turn all pixels above threshold of 2 to 1
   
   ;Get rid of leading empty slice and make sure last slice is empty... only works for first slice of ENTIRE ENSEMBLE
   IF clip eq 1 THEN BEGIN
      slicetotal=total(im,1)
      imin=min(where(slicetotal gt 0))
      IF imin ne -1 THEN im=im[*,imin:*]    ;Clip first slice
      IF slicetotal[n_elements(slicetotal)-1] ne 0 THEN im=[[im],[bytarr(data.op.ndiodes)]] ;Add last empty slice if needed
   ENDIF

   ;Get rid of multiple (3+) blank slices in a row
   IF compress eq 1 THEN BEGIN
      slicetotal=total(im,1)
      IF n_elements(slicetotal) gt 2 THEN BEGIN
         good=where((slicetotal+slicetotal[1:*]+slicetotal[2:*]) gt 0, ngood)
         IF ngood gt 0 THEN im=im[*,good]
      ENDIF
   ENDIF
   
   IF makeplot eq 1 THEN BEGIN
      tvlct,r,g,b,/get  ;save color table
       
      ;Grey color table
      r1=[255,   200, 150,  75,   0]
      g1=[255,   200, 150,  75,   0]
      b1=[255,   200, 150,  75,   0]
      tvlct,r1,g1,b1        ;load up new ct
      width=800
      height=600
      window,/free,xsize=width,ysize=height
      s=size(im,/dim)
      i=0
      REPEAT BEGIN
         istart=i*height
         istop=((i+1)*height)<(s[1]-1)
         tv,im[*,istart:istop]+1, i*(s[0]+3), 0
         IF threshold gt 0 THEN tv,im[*,istart:istop]*3+1, i*(s[0]+3), 0
         i=i+1
      ENDREP UNTIL ((i-1)*s[0] gt width) or (i*height gt s[1])
      tvlct,r,g,b  ;restore original color table
   ENDIF
   
   IF makeplot eq 2 THEN BEGIN
      ;Same as above, but with blue particles and outlined background
      tvlct,r,g,b,/get  ;save color table
      
      ;Blue/black color table
      r1=[255, 0, 0, 0, 0]
      g1=[255, 0, 0, 0, 0]
      b1=[255, 255, 250, 75, 0]
      IF threshold eq 0 THEN BEGIN
         r1=[255, 150, 0, 0, 0]
         g1=[255, 150, 0, 0, 0]
         b1=[255, 255, 255, 150, 0]
      ENDIF
      tvlct,r1,g1,b1    ;load up new ct
      width=800
      height=600
      window,/free,xsize=width,ysize=height+6
      s=size(im,/dim)
      i=0
      REPEAT BEGIN
         istart=i*height
         istop=((i+1)*height)<(s[1]-1)
         tv,bytarr(s[0]+2, height+3)+4, i*(s[0]+8), 1   ;outline
         IF threshold eq 0 THEN tv,im[*,istart:istop], i*(s[0]+8)+1, 2
         IF threshold gt 0 THEN tv,im[*,istart:istop], i*(s[0]+8)+1, 2
         i=i+1
    
      ENDREP UNTIL ((i-1)*s[0] gt width) or (i*height gt s[1])
      tvlct,r,g,b  ;restore original color table
   ENDIF
   
   return,im
END
      