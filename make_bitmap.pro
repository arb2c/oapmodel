FUNCTION make_bitmap, diam, phi, res=res, center=center, fixedn=fixedn, edge=edge, $
                      transmission=transmission, bitmapfile=bitmapfile, imagediam=imagediam
   ;FUNCTION to return an image based on a specified PNG file.
   ;Only fully black pixels in the PNG file will be counted as part of shaded image.
   ;Diam and res in microns
   ;phimajor is angle of major axis in radians
   ;AB 5/2012
   
   IF n_elements(res) eq 0 THEN res=1.0
   IF n_elements(fixedn) eq 0 THEN fixedn=0
   IF n_elements(edge) eq 0 THEN edge=0
   IF n_elements(transmission) eq 0 THEN transmission=0
   IF n_elements(bitmapfile) eq 0 THEN bitmapfile='pumpkin.png'
   IF n_elements(phi) eq 0 THEN phi=0  ;rotation angle, radians
   diampixels=diam/res  ;Diameter of desired image, in pixels
   
   raw=read_png(bitmapfile,r,g,b)
   IF n_elements(r) gt 0 THEN BEGIN  ;PNG contains a color table
      iblack=where((r eq 0) and (g eq 0) and (b eq 0))
      black=where(raw eq iblack[0])
   ENDIF ELSE BEGIN   ;PNG is 4-layer with standard color table
      raw=reform(raw[0,*,*])  ;Just use first layer
      black=where(raw eq 0)
   ENDELSE
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
   finalsize=floor(c2c*diampixels/imagediam)>1 ;Image size required to get correct resolution   
   a=congrid(imr, finalsize, finalsize)      ;Image sized to correct resolution   
   a=abs(1-a)                                ;Reverse black/white
   
   ;Insert particle into frame of fixedn if specified (for FFT diffraction), otherwise just return 'a'
   IF fixedn ne 0 THEN BEGIN; n=round(diamfloat/res*1.5)>2 ELSE n=fixedn
      img=fltarr(fixedn,fixedn)+1
      center=(fixedn-1)/2.0   ;Centerpoint of image
      ;Insert the bitmap in the center of the image
      istart=long(center)-finalsize/2
      IF istart lt 0 THEN return,0
      img[istart:istart+finalsize-1, istart:istart+finalsize-1]=a
   ENDIF ELSE img=a
   
   return,img
END
   
   