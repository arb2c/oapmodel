FUNCTION make_agg, monomerdiam, n, res=res, center=center, fixedn=fixedn, truediam=truediam
BROKEN, don't use yet.

;FUNCTION to return an image an aggregate of n monomers with diameter monomerdiam
   ;Diam and res in microns
   ;AB 5/2012
   
   IF n_elements(res) eq 0 THEN res=1.0
   IF n_elements(fixedn) eq 0 THEN fixedn=0
   IF n_elements(shape) eq 0 THEN shape='discfft'

   ;Make aggregates of particles with some random positions
   a=fltarr(fixedn,fixedn)+1
   FOR i=1,multi DO BEGIN
      CASE shape OF
         'plate'  : b=make_plate(d_adj,psi,randomu(xxx)*2*!pi,center=c,res=1.0,fixedn=n,edge=edge_adj) 
         'discfft': b=make_disc(d_adj,center=c,res=1.0,fixedn=n,edge=edge_adj)
         'column' : b=make_column(d_adj, aspectratio, randomu(xxx)*2*!pi, center=c, res=1.0, fixedn=n, edge=edge_adj)
      ENDCASE
      ;Randomly place
      c=shift(b,randomu(xxx)*n/5-n/10,randomu(xxx)*n/5-n/10)

      ;Make touching aggregates
      IF (touch eq 1) and (i gt 1) THEN BEGIN
         REPEAT c=shift(b,randomu(xxx)*n/5-n/10,randomu(xxx)*n/5-n/10) UNTIL (min(a+c) eq 0)
      ENDIF
      ;Touching aggregates in a chain configuration
      IF (touch eq 2) and (i gt 1) THEN BEGIN
         particlearea=n^2-total(b)
         REPEAT BEGIN
            c=shift(b,randomu(xxx)*n/5-n/10,randomu(xxx)*n/5-n/10) 
            dummy=where(a+c eq 0, ncommon)
         ENDREP UNTIL ((ncommon gt 0) and (ncommon lt particlearea/10.0))
      ENDIF
      a=a*c
   ENDFOR

   FOR i=0,roughen-1 DO BEGIN
      radius=1
      k= SHIFT(DIST(2*radius+1), radius, radius) LE radius
      ;k=[[0,1,0],[1,1,1],[0,1,0]]
      ainv=abs(1-a)  ;Works better on the inverse image so edges don't get dilated/eroded
      b=dilate(ainv,k)
      w=where(ainv ne b, nw)
      r=randomu(xxx, nw)
      good=where(r lt 0.5)     
      a[w[good]]=0   ;Modify the original image, not the inverse
   ENDFOR


   raw=read_png(file,r,g,b)
   iblack=where((r eq 0) and (g eq 0) and (b eq 0))
   ;iwhite=where(r eq 255)
   black=where(raw eq iblack[0])
   im=raw*0      ;Starting out with white image on black background.  Will reverse later.
   im[black]=1
   
   s=size(im)   
   IF s[2] gt 2 THEN im2=erode(im,[[0,1,0],[1,1,1],[0,1,0]]) xor im ELSE im2=im   ;only use outline
   z=where(im2 eq 1, area)  ;these 3 lines get area and indices of each occluded element
   y_ind=z/s[1]             ;remember that y is along airflow, x along array
   x_ind=z mod s[1]         ;area here is area of perimeter, not full particle
   temp=min_circle_fast(x_ind,y_ind)
   imagediam=temp.diam+1    ;plus 1 is due to the fact we're looking a pixels with 0.5 shading threshold, not points
   
   c2c=floor(sqrt(s[1]^2+s[2]^2))+1 ;Corner to corner distance of the image array
   im2=bytarr(c2c,c2c)     ;Create larger image for rotation padding
   xr=[(c2c-s[1])/2, (c2c-s[1])/2 + s[1]]  ;Find x and y indices to place image in center
   yr=[(c2c-s[2])/2, (c2c-s[2])/2 + s[2]]
   im2[xr[0]:xr[1]-1,yr[0]:yr[1]-1]=im       ;Place image into padded image
   imr=rot(im2,phi*180/!pi)                  ;Rotate, don't use cubic which actually makes image worse
   finalsize=floor(c2c*diampixels/imagediam) ;Image size required to get correct resolution   
   a=congrid(imr, finalsize, finalsize)      ;Image sized to correct resolution   
   a=abs(1-a)                                ;Reverse black/white
   
   ;Insert particle into frame of fixedn if specified (for FFT diffraction), otherwise just return 'a'
   IF fixedn ne 0 THEN BEGIN; n=round(diamfloat/res*1.5)>2 ELSE n=fixedn
      img=fltarr(fixedn,fixedn)+1
      center=(fixedn-1)/2.0   ;Centerpoint of image
      ;Insert the bitmap in the center of the image
      istart=long(center)-finalsize/2
      img[istart:istart+finalsize-1, istart:istart+finalsize-1]=a
   ENDIF ELSE img=a
   
   return,img
END
   
   