pro plotfit

; if plot needs to on a file instead of the screen, uncomment following  
  set_plot,'ps'
  device,/encapsulated,/cm,ysize=30.,xsize=20.,yoffset=1.0,filename='fit.ps'
  !p.font=0
  device,/helvetica,font_index=5

  print,'resolution: ',!d.x_px_cm

; teken twee grafieken op 1 blad
  !p.multi=[0,1,4]

; lees plot titel
  ifile = 'plottitel.txt'
  openu,10,ifile
  titel='abcdefghijk'
  readf,10,titel
  instr='abcdefghijk'
  readf,10,instr
  print,'instr=',instr
  decaytime='abcdefghijk'
  readf,10,decaytime
  decaytimee='abcdefghijk'
  readf,10,decaytimee
  starttime='abcdefghijk'
  readf,10,starttime
  source='abcdefghijk'
  readf,10,source
  mt=100.0
  readf,10,mt
  plindex='abcdefghijk'
  readf,10,plindex
  plindexe='abcdefghijk'
  readf,10,plindexe
  chi='abcdefghijk'
  readf,10,chi
  chi2='abcdefghijk'
  readf,10,chi2
  ib2=1
  readf,10,ib2
  close,10
  print,'decaytime,mt=',decaytime,mt

  openu,10,'ana.plot.out'
  pdat=dblarr(9)
  readf,10,pdat
  pdat2=dblarr(6)
  readf,10,pdat2
  pdat3=dblarr(9)
  readf,10,pdat3
  burstid=-1
  readf,10,burstid
  close,10
     
; 
  logg = 0
  print,'instr=',string(instr,format='(a3)')
  if(string(instr,format='(a3)') eq "PCA") then logg=1

; read data
  ifile = 'tt.out'
  openu,10,ifile
  dsize=lonarr(1)
  readf,10,dsize
  dat=fltarr(7,dsize(0))
  readf,10,dat
  close,10

; determine left and right border of plot
  xmin = min(dat(0,*))
  xmax = max(dat(0,*))
  d1=-10.
  d2=+200.
  xmin=max(xmin,d1)
  xmax=min(xmax,d2)
  dx = xmax-xmin
  xmin = xmin - 0.05*dx
  xmax = xmax + 0.05*dx
  xmax=mt
  dx = xmax-xmin

; determine lower and upper border of plot
  if(logg eq 1) then begin
     ymax = 2.*max(dat(1,*)+dat(2,*))
     ymin = ymax / 1000.
     dy = ymax-ymin
  endif else begin
     ymin = min(dat(1,*))
     ymax = max(dat(1,*)+dat(2,*))
     dy = ymax-ymin
     ymin = ymin - 0.05*dy
     ymax = ymax + 0.05*dy
     dy = ymax-ymin
  endelse

  yplot=fltarr(2,2)

  lc1 = '!5Peak flux: '+string(pdat3(2),format="(f7.3)")+'+/-'+string(pdat3(3),format="(f6.3)")+' c/s/cm!e2'
  lc2 = '!5Fluence  : '+string(pdat3(4),format="(f7.2)")+'+/-'+string(pdat3(5),format="(f7.2)")+' c/cm!e2'
  lc3 = '!5Duration : '+string(pdat3(0),format="(f8.2)")+'+/-'+string(pdat3(1),format="(f7.2)")+' s'
  lc4 = '!5Best fit model: '+string(pdat3(6),format="(f2.0)")
  lc5 = '!5Best fluence model: '+string(pdat3(7),format="(f2.0)")
  lc6 = '!5Best duration model: '+string(pdat3(8),format="(f2.0)")
  lc7 = '!5Grey line: exponential fit'
  lc8 = '!5Black line: power law + gauss fit'
  lc9 = '!5Dashed line: power law'
  lb  = 'Power law index: '+string(pdat2(0),format="(f7.3)")+'+/-'+string(pdat2(1),format="(f7.3)")+' (chi!e2!n-red='+string(pdat2(4),format="(f6.2)")+')'
  lbg = 'Power law index: '+string(pdat2(2),format="(f7.3)")+'+/-'+string(pdat2(3),format="(f7.3)")+' (chi!e2!n-red='+string(pdat2(5),format="(f6.2)")+')'
  lb2 = 'Expon. time: '+string(decaytime,format="(f7.2)")+'+/-'+string(decaytimee,format="(f7.2)")+' s (chi!u2!n-red='+string(chi2,format="(f6.2)")+')'
;  lc  = 'Power law chi: '+chi
;  lc2 = 'Expon. chi: '+chi2
  lbs = 'MJD '+starttime
  print,lb
; plot
  !p.position=[0.12,0.73,0.95,0.95]
  for i = 0, dsize(0) - 1 do begin
     deltat=dat(0,dsize(0)-1)-dat(0,dsize(0)-2)
     if( i lt dsize(0)-1 ) then deltat=dat(0,1)-dat(0,0)
     yplot(0,0) = dat(0,i)+0.5*deltat
     yplot(0,1) = dat(0,i)+0.5*deltat
     yplot(1,0) = dat(1,i)-dat(2,i)
     yplot(1,1) = dat(1,i)+dat(2,i)
     if(yplot(1,0) lt ymin) then yplot(1,0)=ymin
;     print,'dat1,dat2=',dat(1,i),dat(2,i)
     if( i eq 0) then begin
        plot,yplot(0,0:1),yplot(1,0:1),$
             xtype=0,$
             ytype=logg,$
             xstyle=5,$
             ystyle=1,$
             title=titel+'  /  '+lbs+' / burstid='+string(burstid,format='(i4)'),$
             xtitle='!5Time [s]',$
             ytitle='!5Flux',$
             charsize=1.5,$
             xrange=[xmin,xmax],yrange=[ymin,ymax],$
             linestyle=0,thick=1
        vert=fltarr(2,2)
        vert(0,0) = 0.
        vert(1,0) = 0.
        vert(0,1) = ymin
        vert(1,1) = ymax
        oplot,vert(0:1,0),vert(0:1,1),linestyle=6,color=150,thick=5
        oplot,yplot(0,0:1),yplot(1,0:1)
        axis,xstyle=1,xax=1,xticklen=0.02,xrange=[xmin,xmax],charsize=0.01,xtitle='',color=0
;        xyouts,xmin+0.35*dx,ymin+0.94*dy,source+' / '+instr+' / '+lbs,charsize=0.8
;        xyouts,xmin+0.38*dx,ymin+0.87*dy,lbs,charsize=0.8
;        xyouts,xmin+0.45*dx,ymin+0.62*dy,lc,charsize=0.8
;        xyouts,xmin+0.45*dx,ymin+0.55*dy,lc2,charsize=0.8
;        xyouts,xmin+0.5*dx,ymin+1.02*dy,titel,charsize=0.99,align=0.5
     endif else begin
        if( dat(2,i) gt 0.) then oplot,yplot(0,0:1),yplot(1,0:1)
     endelse
;     for i = 0, dsize(0)-1 do begin
;        print,'hh',i+1,dat(0,i),dat(3,i)
;     endfor
;     print,'deltat=',deltat
     yplot(0,0) = dat(0,i)
     yplot(0,1) = dat(0,i)+deltat
     yplot(1,0) = dat(1,i)
     yplot(1,1) = dat(1,i)
     oplot,yplot(0,0:1),yplot(1,0:1)
     oplot,dat(0,*),dat(3,*),color=150,thick=5
     oplot,dat(0,*),dat(5,*)
     oplot,dat(0,*),dat(4,*),linestyle=2
;     if(logg eq 1) then oplot,dat(0,*),dat(5,*)
;     if(logg eq 1) then oplot,dat(0,*),dat(4,*),linestyle=2
;     if ( logg eq 0 ) then begin
;        xyouts,xmin+0.5*dx,ymin+0.9*dy,lc1,charsize=0.8
;        xyouts,xmin+0.5*dx,ymin+0.83*dy,lc2,charsize=0.8
;        xyouts,xmin+0.5*dx,ymin+0.76*dy,lc3,charsize=0.8
;        xyouts,xmin+0.5*dx,ymin+0.69*dy,lc4,charsize=0.8
;     endif else begin
;        exp1=alog10(ymin)+(alog10(ymax)-alog10(ymin))*0.90
;        exp2=alog10(ymin)+(alog10(ymax)-alog10(ymin))*0.83
;        exp3=alog10(ymin)+(alog10(ymax)-alog10(ymin))*0.76
;        exp4=alog10(ymin)+(alog10(ymax)-alog10(ymin))*0.69
;        xyouts,xmin+0.5*dx,10.^exp1,lc1,charsize=1.0,font=-1
;        xyouts,xmin+0.5*dx,10.^exp2,lc2,charsize=0.8,font=1
;        xyouts,xmin+0.5*dx,10.^exp3,lc3,charsize=0.8,font=1
;        xyouts,xmin+0.5*dx,10.^exp4,lc4,charsize=0.8,font=0
;     endelse

  endfor

  !p.position=[0.12,0.51,0.95,0.73]
  ymin = -25
  ymax =  25
  dy=50.
  for i = 0, dsize(0) - 1 do begin
     deltat=dat(0,dsize(0)-1)-dat(0,dsize(0)-2)
     if( i lt dsize(0)-1 ) then deltat=dat(0,1)-dat(0,0)
     yplot(0,0) = dat(0,i)+0.5*deltat
     yplot(0,1) = dat(0,i)+0.5*deltat
     yplot(1,0) = (dat(1,i)-dat(2,i)-dat(3,i))/dat(2,i)
     yplot(1,1) = (dat(1,i)+dat(2,i)-dat(3,i))/dat(2,i)
     if( i eq 0) then begin
        plot,yplot(0,0:1),yplot(1,0:1),$
             xtype=0,$
             ytype=0,$
             xstyle=5,$
             ystyle=1,$
             title='',$
             xtitle='!5Time [s]',$
             ytitle='!5Deviation [sigma]',$
             charsize=1.5,$
             xrange=[xmin,xmax],yrange=[ymin,ymax],$
             linestyle=0,thick=1
        xyouts,xmin+0.05*dx,ymin+0.9*dy,'Exponential fit',charsize=0.8
        xyouts,xmin+0.05*dx,ymin+0.83*dy,lb2,charsize=0.8
        xyouts,xmin+0.6*dx,ymin+1.9*dy,lc1,charsize=0.8
        xyouts,xmin+0.6*dx,ymin+1.83*dy,lc2,charsize=0.8
        xyouts,xmin+0.6*dx,ymin+1.76*dy,lc3,charsize=0.8
        xyouts,xmin+0.6*dx,ymin+1.69*dy,lc4,charsize=0.8
        xyouts,xmin+0.6*dx,ymin+1.62*dy,lc5,charsize=0.8
        xyouts,xmin+0.6*dx,ymin+1.55*dy,lc6,charsize=0.8
        xyouts,xmin+0.35*dx,ymin+1.9*dy,lc7,charsize=0.6
        if ( logg eq 1 ) then begin
           xyouts,xmin+0.35*dx,ymin+1.85*dy,lc8,charsize=0.6
           xyouts,xmin+0.35*dx,ymin+1.8*dy,lc9,charsize=0.6
        endif else begin
           xyouts,xmin+0.35*dx,ymin+1.85*dy,lc8,charsize=0.6
           xyouts,xmin+0.35*dx,ymin+1.8*dy,lc9,charsize=0.6
        endelse
     endif else begin
        if(abs(dat(6,i)-1.) lt 0.01) then oplot,yplot(0,0:1),yplot(1,0:1)
     endelse
        axis,xstyle=1,xax=1,xticklen=0.02,xrange=[xmin,xmax],charsize=0.01,xtitle='',color=0
     yplot(0,0) = dat(0,i)
     yplot(0,1) = dat(0,i)+deltat
     yplot(1,0) = (dat(1,i)-dat(3,i))/dat(2,i)
     yplot(1,1) = (dat(1,i)-dat(3,i))/dat(2,i)
     if(abs(dat(6,i)-1.) lt 0.01) then oplot,yplot(0,0:1),yplot(1,0:1)
;     oplot,dat(0,*),dat(4,*)

  endfor

; power law
  !p.position=[0.12,0.29,0.95,0.51]
  ymin = -25
  ymax =  25
  dy=50.
  for i = 0, dsize(0) - 1 do begin
     deltat=dat(0,dsize(0)-1)-dat(0,dsize(0)-2)
     if( i lt dsize(0)-1 ) then deltat=dat(0,1)-dat(0,0)
     yplot(0,0) = dat(0,i)+0.5*deltat
     yplot(0,1) = dat(0,i)+0.5*deltat
     yplot(1,0) = (dat(1,i)-dat(2,i)-dat(4,i))/dat(2,i)
     yplot(1,1) = (dat(1,i)+dat(2,i)-dat(4,i))/dat(2,i)
     if(abs(dat(6,i)) lt 0.01) then begin
        yplot(1,0) = -100.
        yplot(1,1) = -100.
     endif
     if( i eq 0) then begin
        plot,yplot(0,0:1),yplot(1,0:1),$
             xtype=0,$
             ytype=0,$
             xstyle=5,$
             ystyle=1,$
             title='',$
             xtitle='!5Time [s]',$
             ytitle='!5Deviation [sigma]',$
             charsize=1.5,$
             xrange=[xmin,xmax],yrange=[ymin,ymax],$
             linestyle=0,thick=1
        xyouts,xmin+0.05*dx,ymin+0.9*dy,'Power law fit ',charsize=0.8
        xyouts,xmin+0.05*dx,ymin+0.83*dy,lb,charsize=0.8
     endif else begin
        if(abs(dat(6,i)-1.) lt 0.01) then oplot,yplot(0,0:1),yplot(1,0:1)
     endelse
        axis,xstyle=1,xax=1,xticklen=0.02,xrange=[xmin,xmax],charsize=0.01,xtitle='',color=0
     yplot(0,0) = dat(0,i)
     yplot(0,1) = dat(0,i)+deltat
     yplot(1,0) = (dat(1,i)-dat(4,i))/dat(2,i)
     yplot(1,1) = (dat(1,i)-dat(4,i))/dat(2,i)
     if(abs(dat(6,i)-1.) lt 0.01) then oplot,yplot(0,0:1),yplot(1,0:1)
;     oplot,dat(0,*),dat(4,*)

  endfor

  ymin = -25
  ymax =  25
  dy=50.
  !p.position=[0.12,0.07,0.95,0.29]
  for i = 0, dsize(0) - 1 do begin
     deltat=dat(0,dsize(0)-1)-dat(0,dsize(0)-2)
     if( i lt dsize(0)-1 ) then deltat=dat(0,1)-dat(0,0)
     yplot(0,0) = dat(0,i)+0.5*deltat
     yplot(0,1) = dat(0,i)+0.5*deltat
     yplot(1,0) = (dat(1,i)-dat(2,i)-dat(5,i))/dat(2,i)
     yplot(1,1) = (dat(1,i)+dat(2,i)-dat(5,i))/dat(2,i)
     if( i eq 0) then begin
        plot,yplot(0,0:1),yplot(1,0:1),$
             xtype=0,$
             ytype=0,$
             xstyle=1,$
             ystyle=1,$
             title='',$
             xtitle='!5Time [s]',$
             ytitle='!5Deviation [sigma]',$
             charsize=1.5,$
             xrange=[xmin,xmax],yrange=[ymin,ymax],$
             linestyle=0,thick=1
        ibstr = ''
;        ibstr='(plus Gaussian, centroid='+string(pdat(2),format='(f6.2)')+'+/-'+string(pdat(3),format='(f5.2)')+' s, st.dev.='+string(pdat(4),format='(f7.2)')+'+/-'+string(pdat(5),format='(f7.2)')+' s, fracgs='+string(pdat(7),format='(f6.2)')+'+/-'+string(pdat(8),format='(f5.2)')+')'
        ibstr=' plus Gaussian (st.dev.='+string(pdat(4),format='(f6.2)')+'+/-'+string(pdat(5),format='(f5.2)')+' s, fracgs='+string(pdat(7),format='(f7.2)')+'+/-'+string(pdat(8),format='(f6.2)')+')'
        xyouts,xmin+0.05*dx,ymin+0.9*dy,'Power law fit '+ibstr,charsize=0.8
        xyouts,xmin+0.05*dx,ymin+0.83*dy,lbg,charsize=0.8
     endif else begin
        if(abs(dat(6,i)-1.) lt 0.01) then oplot,yplot(0,0:1),yplot(1,0:1)
     endelse
     yplot(0,0) = dat(0,i)
     yplot(0,1) = dat(0,i)+deltat
     yplot(1,0) = (dat(1,i)-dat(5,i))/dat(2,i)
     yplot(1,1) = (dat(1,i)-dat(5,i))/dat(2,i)
     if(abs(dat(6,i)-1.) lt 0.01) then oplot,yplot(0,0:1),yplot(1,0:1)
;     oplot,dat(0,*),dat(4,*)

  endfor

  device,/close

  end
