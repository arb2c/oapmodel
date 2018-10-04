FUNCTION fresnel_fft,diam,z,aspectratio=aspectratio,phi=phi,psi=psi,wavelength=wavelength,$
         res=res,center=center,makeplot=makeplot,shape=shape,crop=crop,transparency=transparency,$
         post=post, smoothparticle=smoothparticle, transmission=transmission, holetransmission=holetransmission,$
         areapercent=areapercent, idealarea=idealarea, levels=levels, bitmapfile=bitmapfile, mout=m, aout=a
   ;Make a 1-um resolution Fresnel diffraction image of a column.
   ;Diam, Z, res in meters
   ;Transparency is a fraction of the overall diameter to shadow on the border
   ;Transmission is the I/I_0 value to be assigned to the particle itself.  Usually 0.
   ;Idealarea is an output variable, given in square meters.
   ;This version follows Trester, Computing in Science and Engineering, Sept 1999, p.77
   ;AB 05/2012
   
   IF n_elements(wavelength) eq 0 THEN wavelength=0.65e-6 ;meters
   IF n_elements(makeplot) eq 0 THEN makeplot=0
   IF n_elements(shape) eq 0 THEN shape='column'
   IF n_elements(crop) eq 0 THEN crop=0
   IF n_elements(transparency) eq 0 THEN transparency=0
   IF n_elements(transmission) eq 0 THEN transmission=0
   IF n_elements(holetransmission) eq 0 THEN holetransmission=1
   IF n_elements(aspectratio) eq 0 THEN aspectratio=2
   IF n_elements(phi) eq 0 THEN phi=0
   IF n_elements(psi) eq 0 THEN psi=0
   IF n_elements(smoothparticle) eq 0 THEN smoothparticle=0  ;Smoothing window in microns
   IF n_elements(areapercent) eq 0 THEN areapercent=0
  
   ;N=wavelength*z, keep N constant to maintain reasonable compute times
   n=1000.0  
   center=(n-1)/2.0     
   
   ;Zd=z*wavelength/R^2, keep this constant too by adjusting d
   Zd=z*wavelength/((diam/2.0)^2)  
   
   d_adj=2*sqrt(n/zd)        ;There is an implicit image resolution of 1.0um/pixel, so this is in microns
   res=diam/d_adj            ;Resolution of the final image, microns per pixel
   edge_adj=diam*(1-transparency)/res    ;Adjust edge thickness too based on the transparency
   undiffracted=0

   ;Check for oversize images
   IF (d_adj gt (n*0.8)) or (zd lt 0.01) THEN BEGIN   ;Keep 80% of pixels for padding
      ;The FFT method will not work for Zd less than about 0.01, where image
      ;becomes larger than number of pixels allotted. Just return an 
      ;undiffracted image for this case.
      res=2.0*diam/n    ;meters/pixel
     
      CASE shape OF
         'column'  :a=make_column(diam, aspectratio, phi, center=c, res=res, fixedn=n,edge=diam*(1-transparency), transmission=transmission)
         'dendrite':a=make_dendrite(diam, aspectratio, phi, center=c, res=res, fixedn=n, transmission=transmission)
         'ellipse' :a=make_ellipse(diam, aspectratio, phi, center=c, res=res, fixedn=n,edge=diam*(1-transparency), transmission=transmission)
         'plate'   :a=make_plate(diam,psi,phi,center=c,res=res,fixedn=n,edge=diam*(1-transparency), transmission=transmission)
         'discfft' :a=make_disc(diam,center=c,res=res,fixedn=n,edge=diam*(1-transparency), transmission=transmission)
         'bitmap'  :a=make_bitmap(diam,phi,center=c,res=res,fixedn=n,edge=diam*(1-transparency), transmission=transmission, bitmapfile=bitmapfile)
      ENDCASE    
      IF smoothparticle gt 0 THEN a=smooth(a,smoothparticle/res) 
      imfinal=a
      undiffracted=1
      idealarea=(n^2-total(a)) * res^2
   ENDIF ELSE BEGIN   
      ;Make the undiffracted image
      CASE shape OF
         'plate'  : a=make_plate(d_adj,0,phi,center=c,res=1.0,fixedn=n,edge=edge_adj, transmission=transmission) 
         'discfft': a=make_disc(d_adj,center=c,res=1.0,fixedn=n,edge=edge_adj, transmission=transmission)
         'column' : a=make_column(d_adj, aspectratio, phi, center=c, res=1.0, fixedn=n, edge=edge_adj, transmission=transmission, $
                    holetransmission=holetransmission, areapercent=areapercent)
         'dendrite':a=make_dendrite(d_adj, aspectratio, phi, center=c, res=1.0, fixedn=n, transmission=transmission)
         'ellipse': a=make_ellipse(d_adj, aspectratio, phi, center=c, res=1.0, fixedn=n, edge=edge_adj, transmission=transmission)
         'bitmap':  a=make_bitmap(d_adj,phi,center=c,res=1.0,fixedn=n,edge=edge_adj, transmission=transmission, bitmapfile=bitmapfile)
      ENDCASE
      idealarea=(n^2-total(a)) * res^2    
      IF smoothparticle gt 0 THEN a=smooth(a,smoothparticle/res) 
      ;Gaussian smoothing is very slow, commented out for now
      ;IF smoothparticle gt 0 THEN a=gauss_smooth(a,2/res, /edge_wrap) ;Gaussian smooth with 2 micron window
      
      ;Make distance matrix
      dist=fltarr(n,n)    
      FOR i=0,n-1 DO dist[i,*]=i-center
      cornerdist=(dist+center)^2 + transpose(dist+center)^2
      m=a*exp(complex(0,-!pi)/n*cornerdist)
         
      ;Compute diffracted image
      im=abs(fft(m))^2
      imfinal=im*n^2   ;Normalize image intensity to 1   
   ENDELSE
   ;Cropping and plotting
   IF crop eq 1 THEN BEGIN
      imagedx=where(min(imfinal, dim=2) lt 0.90,nx)
      imagedy=where(min(imfinal, dim=1) lt 0.90,ny)
      IF (nx gt 0) and (ny gt 0) THEN BEGIN
         imfinal=imfinal[min(imagedx):max(imagedx), min(imagedy):max(imagedy)]
         a=a[min(imagedx):max(imagedx), min(imagedy):max(imagedy)]
      ENDIF
      center=(max(imagedx)-min(imagedx))/2.0
      IF undiffracted eq 1 THEN BEGIN 
         ;Undiffracted images will be cropped very tight.  This
         ;adds a little white padding for delay computations, etc.
         s=size(imfinal,/dim)
         pad=2
         paddedimage=fltarr(s[0]+pad*2,s[1]+pad*2)+1
         paddedimage[pad:s[0]+pad-1,pad:s[1]+pad-1]=imfinal
         imfinal=paddedimage
      ENDIF
   ENDIF
  
   IF makeplot eq 2 THEN BEGIN
      psize=600.0
      window,/free,xsize=psize*2,ysize=psize
      isize=size(imfinal,/dim)
      mag=psize/max(isize)
      tv,congrid(rotate(imfinal/max(imfinal),2)*250,isize[0]*mag,isize[1]*mag,/interp)
      tv,congrid(a,isize[0]*mag,isize[1]*mag,/interp)*150,psize,0
      stop  
   ENDIF
   IF makeplot eq 1 THEN BEGIN
      psize=500.0
      IF n_elements(post) ne 0 THEN BEGIN
         set_plot,'ps'
         device,file=post
      ENDIF ELSE window,/free,xsize=psize*2,ysize=psize*2,xpos=1200
      isize=size(imfinal,/dim)
      mag=psize/max(isize)
      imfinal=reverse(reverse(imfinal,2),1)  ;For some unknown reason the diffracted images are inverted
      tv,bytarr(psize*2,psize*2)+200   ;uniform background
      tv,congrid(imfinal*200<255,isize[0]*mag,isize[1]*mag,/interp),psize,psize
      tv,congrid(a,isize[0]*mag,isize[1]*mag,/interp)*200<255,0,psize
      dim=baumgardner_delay(imfinal,1.0e-6,100,/clean)
      tv,congrid(dim*200<255,isize[0]*mag,isize[1]*mag,/cubic),0,0

      mpp=(25.0e-6/res)
      n=isize[0]/mpp
      xpos=findgen(n+1)*mpp
      ypos=xpos
      blocky=downsample_grey(dim,xpos,ypos,mpp-1,mpp-1,fltarr(200)+1.0,show=1)  ;Dim is modified in here to overlay pixels
      dim[where(dim ne 0)]=1 ;Get rid of background dim, only show pixels
      levels=blocky.levels
      tv,congrid(dim*200<255,isize[0]*mag,isize[1]*mag,/interp,center=0),psize,0 ;Turn center on/off if ugly
      xyouts,10,psize+10,'a',/device,charsize=6,charthick=1,color=0,font=1
      xyouts,psize+10,psize+10,'b',/device,charsize=6,charthick=1,color=0,font=1
      xyouts,10,10,'c',/device,charsize=6,charthick=1,color=0,font=1
      xyouts,psize+10,10,'d',/device,charsize=6,charthick=1,color=0,font=1
      plots,[0,psize*2],[psize,psize],color=0,thick=2,/device
      plots,[psize,psize],[0,psize*2],color=0,thick=2,/device
      ;tv,congrid(
      ;stop  
      IF n_elements(post) ne 0 THEN BEGIN
         device,/close
         set_plot,'x'
      ENDIF
      window,4,xsize=1000,ysize=1000
      tv,dim*200
   ENDIF

   return,imfinal
END