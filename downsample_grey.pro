FUNCTION downsample_grey, image, xpos, ypos, xsize, ysize, sensitivity, show=show, thresholds=thresholds
   ;FUNCTION to downsample a hi-res image into a 4-level
   ;grey-scale image to simulate an OAP.
   ;image is hi-res image of I/I_0 as in Korolev 2007.
   ;xpos=left edge positions (in pixels) of the diodes
   ;ypos=bottom edge positions (in pixels) of the diodes
   ;xsize and ysize number of pixels to count in each dimension

   IF n_elements(show) eq 0 THEN show=0  ;Modify the original image to show shadowed pixels
   IF n_elementS(thresholds) eq 0 THEN thresholds=[0.25, 0.5, 0.75]  ;Default 3-levels thresholds for I/I0

   nx=n_elements(xpos)
   ny=n_elements(ypos)
   psize=size(image,/dim)
   grey=bytarr(nx,ny)
   levels=fltarr(nx,ny)+1.0
   xpos=round(xpos)
   ypos=round(ypos)
   
   ;Find diodes that are fully/partly in bounds
   goodx=where((xpos+xsize gt 0) and (xpos lt psize[0]),ngoodx)
   goody=where((ypos+ysize gt 0) and (ypos lt psize[1]),ngoody)
   IF (ngoodx eq 0) or (ngoody eq 0) THEN return,{grey:grey, levels:levels}
   
   left=xpos[goodx]>0
   right=(xpos[goodx]+(xsize>1)-1) < (psize[0]-1)
   bottom=ypos[goody]>0
   top=(ypos[goody]+(ysize>1)-1) < (psize[1]-1)
   ;totalimaged=floor(xsize)*floor(ysize)

   FOR i=0,ngoodx-1 DO BEGIN
      FOR j=0,ngoody-1 DO BEGIN
         ;Make correction for partly in-bounds diodes
         ;nimaged=(right[i]-left[i]+1)*(top[j]-bottom[j]+1)  ;Number counted
         ;tot=total((image)[left[i]:right[i],bottom[j]:top[j]]) 
         ;m=(tot+(totalimaged-nimaged)*1.0) / totalimaged
         
         m=mean((image)[left[i]:right[i],bottom[j]:top[j]])
         levels[goodx[i],goody[j]]=m
         ;IF nimaged ne totalimaged and (m ne 1) then stop       
         IF m lt (1-thresholds[0])*sensitivity[i] THEN grey[goodx[i],goody[j]]=1
         IF m lt (1-thresholds[1])*sensitivity[i] THEN grey[goodx[i],goody[j]]=2
         IF (show eq 1) and  (m lt (1-thresholds[1])*sensitivity[i]) THEN BEGIN 
            ;Print data and highlight the diode stripes
            print,mean((image)[left[i]:right[i],bottom[j]:top[j]]),left[i],right[i],bottom[j],top[j]  
            grey[goodx[i],goody[j]]=2  
            image[left[i]:right[i],bottom[j]:top[j]]=0;image[left[i]:right[i],bottom[j]:top[j]]+0.3  
         ENDIF
         IF m lt (1-thresholds[2])*sensitivity[i] THEN grey[goodx[i],goody[j]]=3
      ENDFOR
   ENDFOR
   return,{grey:grey, levels:levels}
END
   
   

   