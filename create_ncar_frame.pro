FUNCTION create_ncar_frame, x, time, clockMhz, dofreject, ndiodes=ndiodes
   ;Create NCAR formatted particle, including timing data and the image
   ;AB 9/2018
   IF n_elements(tas) eq 0 THEN tas=100
   IF n_elements(ndiodes) eq 0 THEN ndiodes=fix(64)

   ;Number slices in particle, no trailing slice in NCAR data (can add one here if needed)
   ;nslices=fix(n_elements(x)/ndiodes)
   slicetotal=total(x,1)
   nslices=max(where(slicetotal gt 0)) + 1   ;Do it this way to ignore trailing empty slices that often come in
   
   ;Time in NCAR form
   counter = ulong64(double(time)*clockMhz*1e6) 
   
   ;Sync pattern depends on dofreject and the version of acquisition card (flagged by clockMhz)
   IF clockMhz eq 12.0 THEN BEGIN
      sync=ulong64('AAAAAA0000000000'x)
      IF dofreject eq 1 THEN sync=ulong64('AAAAAB0000000000'x)
      timeline=(counter and '000000ffffffffff'x) or sync 
   ENDIF
   IF clockMhz eq 33.0 THEN BEGIN  ;new card in 2018, a.k.a. F2DC_v2
      sync=ulong64('AAAA000000000000'x)
      IF dofreject eq 1 THEN sync=ulong64('AAAA800000000000'x)
      timeline=(counter and '00003fffffffffff'x) or sync
   ENDIF

   ;Write the image as a series of ulong64s
   ntrailing=1   ;Empty slices that trail final slice
   image=ulon64arr(nslices+1+ntrailing) + 'ffffffffffffffff'x
   image[0]=swap_endian(timeline)
   FOR i=1,nslices DO BEGIN
      image[i]=swap_endian(bin2dec64(1-x[*,i-1]))
   ENDFOR
   return,image
END
