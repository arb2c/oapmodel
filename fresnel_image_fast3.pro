FUNCTION fresnel_image_fast3,diam,z,wavelength=wavelength,center=center, res=res, crop=crop
   ;Make a 1-um resolution Fresnel diffraction image.
   ;Diam, Z, in meters
   ;Res in microns
   ;AB 05 2012

   IF n_elements(crop) eq 0 THEN crop=0
   IF n_elements(res) eq 0 THEN res=1.0
   center=0                        ;Need to define first in case of early return
   ntest=round(diam*1e6*3.0/res)>1   ;6X radius... number of microns/pixels to model
   x=findgen(ntest)*res            ;Little_r... compute diffraction at all these 1-micron distances 
   f=korolev_fresnel(diam,z,littler=x/1.0e6,wavelength=wavelength)
   
   ;Find the last place where intensity is lower than 0.9, add 10, and crop there
   n=max(where(smooth(f.intensity,5<(n_elements(f.intensity)-1)) lt 0.95,nw))
   center=n                ;Send this back to calling program to indicate centerpoint
   IF nw eq 0 THEN return,fltarr(1,1)  ;No shadows at all
   IF nw eq 1 THEN BEGIN & im3=fltarr(3,3)+1 & im3[1,1]=min(f.intensity) & return,im3 & ENDIF  ;Single shadow
   x=findgen(n)            ;Little_r... compute diffraction at all these 1-micron distances 
   d=fltarr(n,n)           ;Distance from center
   im1=fltarr(n,n)+1
   FOR i=0,n-1 DO d[i,*]=sqrt(x[i]^2+x^2)
  
   s=sort(d)
   dround=round(d)  ;to match up with x...
   u=uniq(dround[s])
   ;Make a 'lookup table' for all the distances and fill in the particle
   FOR j=1,n_elements(u)-1 DO BEGIN
      du=dround[s[u[j]]]
      w=s[u[j-1]+1:u[j]]   ;where(dround eq du)
      IF du ge n-1 THEN im1[w]=1 ELSE $   ;n-1 keeps a border pixel with value=1.  This is important for delay computations later.
         im1[w]=f.intensity[du]
   ENDFOR

   im2=[reverse(im1),im1[1:n-1,*]]
   imfinal=[[reverse(im2,2)],[im2[*,1:n-1]]]
   
   IF crop eq 1 THEN BEGIN
      imagedx=where(min(imfinal, dim=2) lt 0.90,nx)
      imagedy=where(min(imfinal, dim=1) lt 0.90,ny)
      IF (nx gt 0) and (ny gt 0) THEN BEGIN
         imfinal=imfinal[min(imagedx):max(imagedx), min(imagedy):max(imagedy)]
       ENDIF
      center=(max(imagedx)-min(imagedx))/2.0
   ENDIF
   return,imfinal
END