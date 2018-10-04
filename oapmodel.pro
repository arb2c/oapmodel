PRO oapmodel, op   ;model, probe, file=file, writeimages=writeimages, n=n, pathlength=pathlength, tas=tas, taserror=taserror,$
                   ;lam=lam, tau=tau, probename=probename, senserror=senserror, sensitivity=sensitivity
   ;Make a model of OAP images, applying Fresnel diffraction and
   ;other effects.
   ;AB 5/2012
   
   oap_update_op,op
   
   file=op.root+'.dat' 
   imagefile=op.root+'.images'
   close,1,2
   IF op.writeimages eq 1 THEN BEGIN
      openw,1,imagefile
      openw,2,imagefile+'_float'
   ENDIF
   ;help,op,/st
   print,'Creating file: ',file
  
   ;-----------------Initialize variables-------------------------
   arraywid=op.ndiodes*op.res        ;meters
   diodehist=lonarr(op.ndiodes)
   n=float(op.n)
   
   ;Make compatible processing options used to call soda2_findsize.pro, soda2_samplearea.pro
;   opsoda={res:op.res*1e6, reconstruct:1, smethod:'fastcircle', armwidth:op.armwidth*100, numdiodes:op.ndiodes}
;   IF n_elements(popsoda) THEN ptr_free,popsoda
;   popsoda=ptr_new(opsoda)
;   miscsoda={yres:op.res*1e6}
;   IF n_elements(pmiscsoda) THEN ptr_free,pmiscsoda
;   pmiscsoda=ptr_new(miscsoda)
   
;   opareasize=opsoda
;   opareasize.smethod='areasize'
;   popareasize=ptr_new(opareasize)
   
   ;------------------Initialize particles-------------------------------
   ;Dimensions of the volume to model, in meters.  x and y plane is seen by probe.  z dimension is length of track
   ;xdim will be array width plus 1 centimeter on each sideinthist=make1dhist(t,intendbins)
   
   modelwidth=arraywid+op.pad*2      ;add an extra centimeter to each side for particle radius
   modelheight=op.armwidth+op.pad*2
   seed=5  ;Use a reproducible seed for now
   
   ;Generate a distribution of particle sizes
   IF op.monodisperse[0] ne 0 THEN BEGIN
      ;Create a monodisperse distribution.  If monodisperserange is non-zero, it makes
      ;it over a range of values about the midpoint.  This is meant to simulate an entire bin size.
      ;Monodisperse can also be an array of sizes, each will have same concentration.
      monoindex=fix(randomu(seed,n)*n_elements(op.monodisperse))  ;Index of which size to assign each particle
      d=fltarr(n) + op.monodisperse[monoindex] + randomu(seed,n)*op.monodisperserange - op.monodisperserange/2.0      
      ;Use n0 (in unit of #/m3) to get length, or defaults to 10/L when n0=0.
      IF op.n0 ne 0 THEN numperm3=op.n0 ELSE numperm3=1.0e4
      modellength=n / (numperm3*modelwidth*modelheight) 
   ENDIF ELSE BEGIN
      monoindex=bytarr(n)  ;Set all to zero, since this will not be used for indexing bitmapfile
      IF op.mu eq 0 THEN BEGIN
         ;Use inverse of n=exp(-lam*d) to make random points
         d=(-alog(randomu(seed,n)))/op.lam  ;Diameter, meters
      ENDIF
      
      ;If it is a gamma function, compute it this way
      IF op.mu ne 0 THEN BEGIN
         xvalues=dindgen(1e5)/1.0e6  ;Size array up to 10 cm
         d=randomf(xvalues,op.n0*exp(-op.lam*xvalues)*xvalues^op.mu,n)
      ENDIF
      
      ;Compute path length based on n and n0.
      ;Solving integral of general gamma function gives n0=n*(lam^(1+mu)) / gamma(1+mu) / volume
      modellength=n * (op.lam^(1+op.mu) / gamma(1+op.mu) / (op.n0*modelwidth*modelheight)) 
   
      ;If mu is negative (superexponential) then integral is undefined.
      IF modellength lt (n/1.0e4) THEN stop,'Model length '+string(modellength)+' too short.  Need to adjust n0 downwards'
      IF finite(modellength) eq 0 THEN stop,'Infinite model length'
   ENDELSE
 
   ;Particle positions
   ;X-posistion
   IF op.spinningdisk eq 1 THEN x=randomu(seed,n,/double)*op.res*2 - op.res  + modelwidth/2.0 $ ;Simulate a spinning disk with 2 pixels of 'wobble' 
           ;IF op.spinningdisk eq 1 THEN x=randomu(seed,n)*op.res*20 + modelwidth/2.0 - op.res*10 $ ;Extra wobble
      ELSE x=randomu(seed,n,/double)*modelwidth  ;randomly across array 
   ;Time dimension
   y=randomu(seed,n,/double)*modellength
   ;Z-position
   z=randomu(seed,n,/double)*modelheight  
   IF op.centerz eq 1 THEN z=fltarr(n)+modelheight/2.0
   IF op.fixedz ne 0 THEN z=fltarr(n)+modelheight/2.0 + op.fixedz  ;Fixed z, relative to center of DoF
   
   modelboxsize=[modelwidth, modellength, modelheight]
   print,'Max time: ',max(y)/op.tas
  
   ;--------Create separate shattering mode--------
   
   ;Shattered particles will be created from particles that would have hit a probe tip, defined 
   ;  as 2 1cm square areas on the upper and lower portions of the modelbox. 
   IF op.shatter THEN BEGIN      
      IF op.pad lt 0.01 THEN stop,'Pad length should be at least 1 cm for shattering simulation.'
      
      ;Find particles that would hit a probe tip
      zcenter=modelheight/2.0
      xcenter=modelwidth/2.0
      dcof=abs(z-zcenter)           ;distance from center of focus
      dx=abs(x-xcenter)
      
      ;Use figure 7 in Vidaurre and Hallet (JTECH 2009) to determine the number of resultant shattered particles
      ;Have fit a line to their graph, n=5.0*diam[mm]
      nshatters=floor(d*1000*10.0)    ;Compute for all particles for now
      tiphits=where((dcof gt op.armwidth/2.0) and (dcof lt (op.armwidth/2.0 + 0.01)) and (dx lt 0.005) and (nshatters gt 1),ntiphits)
      
      ;Only use nshatters where conditions above (in tiphits) were met
      nshatters=nshatters[tiphits]
      
      totalshatters=total(nshatters)
      xs=fltarr(totalshatters)
      ys=fltarr(totalshatters)
      zs=fltarr(totalshatters)
      ds=fltarr(totalshatters)
      
      ;For each particle that hits a tip, create sub-particles
      istart=0
      FOR i=0L,ntiphits-1 DO BEGIN
         
         ;Diameter of shattered particles, randomly exponential, mass is conserved
         IF op.shape eq 'column' THEN BEGIN
            columncorrect=op.columnaspect/sqrt(1+float(op.columnaspect)^2)
            dminor=d[tiphits[i]]*columncorrect/op.columnaspect * 1.0e2  ;cm
            dmajor=d[tiphits[i]]*columncorrect * 1.0e2  ;cm
            mass=0.65*(dminor^2)*dmajor *0.9  ;g
            
            ;Use inverse of n=exp(-lam*d) to make random points
            ;Assume lambda of 60, which is what we saw in NAMMA for the shattered mode
            dshatter=(-alog(randomu(seed,nshatters[i])))/(60.0*1e2)  ;Diameter, meters
            dminor=dshatter*columncorrect/op.columnaspect *1.0e2 ;cm
            dmajor=dshatter*columncorrect *1.0e2 ;cm
            mass_shatter=total(0.65*(dminor^2)*dmajor *0.9) ;g
            factor=(mass/total(mass_shatter))^(0.333)  ;Compute this factor to conserve mass
            dshatter=dshatter*factor
            
         ENDIF
         
         IF op.shape eq 'disc' THEN BEGIN
            mass=!pi/6*((d[tiphits[i]]*1.0e2)^3)  ;g
             ;Use inverse of n=exp(-lam*d) to make random points
            ;Assume lambda of 60, which is what we saw in NAMMA for the shattered mode
            dshatter=(-alog(randomu(seed,nshatters[i])))/(60.0*1e2)  ;Diameter, meters
            mass_shatter=total(!pi/6*(dshatter*1.0e2)^3)
            factor=(mass/total(mass_shatter))^(0.333)  ;Compute this factor to conserve mass
            dshatter=dshatter*factor          
         ENDIF
         
         ;Compute position of shattered particles.
         ;Azimuth is random in a sector of 45 degrees after shattering event
         azmain=randomu(seed,1)*2*!pi  ;has to be done in 2 steps for addition of arrays to work properly
         azimuth=azmain[0] + randomu(seed,nshatters[i])*!pi/8
         ;Distance will be approximately proportional to mass/area, i.e. d^1.
         ;distance=0.005+dshatter*80  ;The 0.005 is to move it to the edge of the tip.  Am assuming now all particles hit the center of the tip.
         
         ;Random distance
         distance=0.005+randomu(nshatters[i])*0.02  ;The 0.005 is to move it to the edge of the tip.
         
         
         ;A more complex way that I have not fully worked out....
         ;fd=0.5*0.5*op.tas^2*0.47*!pi/4*dshatter^2   ;Drag force, [kg*m/s2 ], assuming 0.5 kg/m3 air density.  0.47 is Cd of a sphere.
         ;fi=0.5*op.tas^2*(mass_shatter/1000)         ;Inertial force, assume 20% of TAS
         ;distance=0.5*mass*op.tas^2*0.02  ;Assume travel time of 2cm from tip EDGE (not original z) to sample area
         
         istop=istart+nshatters[i]-1
         xs[istart:istop]=xcenter+distance*cos(azimuth)
         ys[istart:istop]=y[tiphits[i]]+randomu(seed, nshatters[i])*0.02  ;Randomly within 2cm of original time
         zs[istart:istop]=zcenter + op.armwidth/2.0 + 0.005 + distance*sin(azimuth)
         ds[istart:istop]=dshatter      
         istart=istart+nshatters[i]
      ENDFOR

      x=[x,xs]
      y=[y,ys]
      z=[z,zs]
      d=[d,ds]
   ENDIF ;Shattering
    
   ;----------Initialize more variables and sort particles in time-----------------
   n=n_elements(x)    ;Includes both real and shattered particles
   print,n
   slicecount=0L
   slicestart=lonarr(n)
   slicestop=lonarr(n)
   
   fpart={diam50:fltarr(n), diam25:fltarr(n), areasize:fltarr(n), allin50:intarr(n), allin25:intarr(n), centerin:intarr(n),$
         ar50:fltarr(n), ar25:fltarr(n), area:fltarr(n), areafilled:fltarr(n), sumgrey:fltarr(n), dedge:fltarr(n), $
         pixelscleared:fltarr(n), pixelsfirstslice:fltarr(n), maxshadow:fltarr(n)}
   dpart=fpart  ;Delayed particles
   
   s=sort(y)
   x=x[s]
   y=y[s]
   z=z[s]
   d=d[s] > (1.0e-6)  ;Make minimum size of 1um to avoid floating errors
   IF op.minsize gt 0 THEN BEGIN   ;This is to set a minimum size, good for low lambda where don't care about small particles
      bad=where(d lt op.minsize, nbad)
      d[bad]=0  ;These will now be skipped in the main loop
   ENDIF
   imagearea=fltarr(n_elements(d))
   t=y/op.tas  ;time
   inttime=[0,t[1:*]-t] ;inttime
   
   ;---------Make images------------------
   diodeposition=(dindgen(op.ndiodes)*op.res + op.res/2.0) + op.pad
   zcenter=modelheight/2.0
   dcof=abs(z-zcenter)           ;distance from center of focus
   modelres=op.res/10.0          ;meters, model resolution
   modelsperpixel=op.res/modelres
   lastpercentcomplete=-1
   
   FOR i=0L,n-1 DO BEGIN
      ;Diode positions relative to center of image
      deltax=diodeposition-x[i]   
      
      ;Make sure particle is within armwidth.
      ;Also test to see if particle is close enough (r_over_r < 5) to the diode array to be imaged
      IF (dcof[i]*2 lt op.armwidth) and (min(abs(deltax)) lt (2.5*d[i]) and (d[i] gt 0)) THEN BEGIN
      
         ;Create hi-res Fresnel image         
         IF op.shape eq 'disc' THEN BEGIN
            fim=fresnel_image_fast3(d[i], dcof[i], center=cx, res=modelres*1e6, wavelength=op.wavelength) 
            imagearea[i]=!pi/4*d[i]^2
         ENDIF ELSE BEGIN
            ;Use Fourier Transform method for other shapes
            fim=fresnel_fft(d[i], dcof[i], aspect=op.columnaspect, phi=randomu(seed)*!pi, psi=0 ,center=cx, $
                res=modelres, /crop, shape=op.shape, transparency=op.transparency, wavelength=op.wavelength,$
                smoothparticle=op.smoothparticle, idealarea=idealarea, bitmapfile=op.bitmapfile[monoindex[i]])
            ;Modelres is variable for the FFT method
            modelsperpixel=op.res/modelres
            imagearea[i]=idealarea  ;m2, computed from the undiffracted image (for use with complicated shapes)
            
            ;Undiffracted column method below
            ;fim=make_column(d[i]*1e6, op.columnaspect, randomu(seed)*!pi ,center=c, res=modelres*1e6)
         END
         
         IF cx gt 0 THEN BEGIN   ;Skip particles with no shadow       
            xpos=deltax/op.res*modelsperpixel + cx   ;Align to center of the hi-res image
          
            ;NOTE: ypos is not really aligned with the actual slices, random factor added.
            cy=(size(fim,/dim))[1]
            numslices=fix((cy/modelsperpixel/op.taserror)>2)+1
            ypos=findgen(numslices)*modelsperpixel*op.taserror - randomu(seed)*modelsperpixel 
             
            ;Perfect Fresnel image, resampled to OAP resolution
            fim_oap=downsample_grey(fim, xpos, ypos, modelsperpixel*op.diodespacing, modelsperpixel, op.sensitivity*0+1, thresholds=op.greythresholds)     
;IF total(fim_oap[*,numslices-1]) gt 0 THEN stop            

            ;Thresholded
            im_50=((fim_oap.grey > 1)-1)<1  ;turn all pixels above threshold of 2 to 1
            ;Now using maxshadow for this.  IF min(fim_oap.levels)*100 gt (100-op.dofthreshold) THEN im_50=im_50*0  ;Apply DoF threshold, erase particle if not met 
            im_25=fim_oap.grey<1
            x50=soda2_findsize(im_50,op.res*1e6,op.res*1e6)
            ;x50area=soda2_findsize(im_50,popareasize,miscsoda.yres)
            x25=soda2_findsize(im_25,op.res*1e6,op.res*1e6)
            fpart.diam50[i]=x50.diam
            fpart.diam25[i]=x25.diam
            fpart.allin50[i]=x50.allin
            fpart.allin25[i]=x25.allin
            fpart.ar50[i]=x50.ar
            fpart.ar25[i]=x25.ar
            fpart.area[i]=total(im_50)
            fpart.areafilled[i]=total(fillholes(im_50))
            fpart.sumgrey[i]=total(fim_oap.grey)
            fpart.dedge[i]=dedge(im_50,res=op.res*1e6)
            fpart.areasize[i]=x50.areasize
            fpart.centerin[i]=x50.centerin
            fpart.maxshadow[i]=1-min(fim_oap.levels)
        
            ;Image with time delay, sensitivity, skipped slices, etc applied   
            dim=baumgardner_delay(fim,modelres,op.tas,tau=op.tau,xpos=xpos,/clean)
           
            dim_oap=downsample_grey(dim, xpos, ypos, modelsperpixel*op.diodespacing, modelsperpixel, op.sensitivity, thresholds=op.greythresholds)      
            IF op.skipslice THEN BEGIN
               temp50=((dim_oap.grey > 1)-1)<1  ;make a 50% image
               w=where(total(temp50,1),nw)
               IF nw gt 0 THEN dim_oap.grey[*,0:w[0]]=0  ;clear all slices through the first slice with a 50% shadow 
               dpart.pixelscleared[i]=total(temp50 - (((dim_oap.grey > 1)-1)<1))   ;Save this value for testing
               IF nw gt 1 THEN dpart.pixelsfirstslice[i]=total(temp50[*,w[1]])  ;Number of pixels on first recorded slice
            ENDIF
            ;Thresholded
            im_50=((dim_oap.grey > 1)-1)<1  ;turn all pixels above threshold of 2 to 1
            ;Now using maxshadow for this.  IF min(dim_oap.levels)*100 gt (100-op.dofthreshold) THEN im_50=im_50*0  ;Apply DoF threshold, erase particle if not met 
            im_25=dim_oap.grey<1
            x50=soda2_findsize(im_50,op.res*1e6,op.res*1e6)
            ;x50area=soda2_findsize(im_50,popareasize,pmiscsoda)
            x25=soda2_findsize(im_25,op.res*1e6,op.res*1e6)
            dpart.diam50[i]=x50.diam
            dpart.diam25[i]=x25.diam
            dpart.allin50[i]=x50.allin
            dpart.allin25[i]=x25.allin
            dpart.ar50[i]=x50.ar
            dpart.ar25[i]=x25.ar
            dpart.area[i]=total(im_50)
            dpart.areafilled[i]=total(fillholes(im_50))    
            dpart.sumgrey[i]=total(dim_oap.grey)
            dpart.dedge[i]=dedge(im_50,res=op.res*1e6)
            dpart.areasize[i]=x50.areasize
            dpart.centerin[i]=x50.centerin
            dpart.maxshadow[i]=1-min(dim_oap.levels)
            
            diodehist=diodehist+total(im_50,2)
            
            IF op.writeimages eq 1 THEN BEGIN
               ;Crop to size
               w=where(total(dim_oap.grey,1),nw)                          
               ;Now using maxshadow for this.  IF (nw gt 0) and min(dim_oap.levels)*100 le (100-op.dofthreshold) THEN BEGIN
               IF (nw gt 0) THEN BEGIN
                  istart=(min(w)-1) >0  ;Give 1 slice for spacing
                  istop=max(w)
                  writeu,1,dim_oap.grey[*,istart:istop]              
                  writeu,2,dim_oap.levels[*,istart:istop]
                  ;Save start position
                  slicestart[i]=slicecount
                  slicecount=slicecount+istop-istart+1
                  slicestop[i]=slicecount-1      
               ENDIF   
            ENDIF

            percentcomplete=floor(float(i)/n*100)
            IF percentcomplete ne lastpercentcomplete THEN print,percentcomplete
            lastpercentcomplete=percentcomplete
         ENDIF
      ENDIF
   ENDFOR
   
   ;Save data
   data={op:op, fpart:fpart, dpart:dpart, n:n, x:x, y:y, z:z, d:d*1e6, imagearea:imagearea*1e12, t:t, inttime:inttime,$
         modelboxsize:modelboxsize, diodehist:diodehist, $
         slicestart:slicestart, slicestop:slicestop, imagefile:imagefile}
   save,data,file=file,/compress
   ;ptr_free,popsoda,pmiscsoda
   close,1,2
END
