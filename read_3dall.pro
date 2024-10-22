PRO read_3dall

  SAT_KEY='3dall'

  ;DATA_DIR = '~/GITM_events/gitm_201303_fac3/Backup/data_0317_1/'$
  ;           +SAT_KEY+'/'
  DATA_DIR = './'

  SAVE_DIR ='./'
  SAVE_DIR = SAVE_DIR + SAT_KEY +'/'

  SAVE_DIR_EXIST=FILE_TEST(SAVE_DIR,/DIRECTORY)                               
                                                                               
  IF (SAVE_DIR_EXIST gt 0) THEN BEGIN                                          
     print, 'SAVE DIR EXISTS !'                                                
                                                                               
  ENDIF ELSE BEGIN                                                             
                                                                               
     FILE_MKDIR, SAVE_DIR                                                      
     print, 'SAVE DIR HAS BEEN CREATED !'                                      
                                                                               
  ENDELSE

  FLIST=FINDFILE(DATA_DIR+'3DUSR*.bin')
  NFILES=N_ELEMENTS(FLIST)                                                     
                                                                               
  print, nfiles

  ;nfiles=1

  FOR IFILE=0,nfiles-1 DO BEGIN                                                
                                                                               
                                                                               
     FILENAME = FLIST(IFILE)                                                   
     PRINT, IFILE;,FILENAME      
                                                                               
     data = 0.0                                                                
                                                                               
     read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $           
                             vars, data, rb, cb, bl_cnt, itime1                
                                                                               
                                                                               
     alt = reform(data(2,*,*,*))                      
     lat = reform(data(1,*,*,*))                                 
     lon = reform(data(0,*,*,*))                        
     
     lon1=lon(*,*,2)
     lat1=lat(*,*,2)
     alt1=alt(2,2,*)

     ;help,lon1,lat1,alt1

     ;help,alt,lat,lon

     ;for i=0,nvars-1 do print, tostr(i)+'. '+vars(i) 
  

     Eden = reform(data(3,*,*,*))
     Wi = reform(data(4,*,*,*))
     ExB = reform(data(5,*,*,*))
     Rho = reform(data(6,*,*,*))


     SAVE_STR=strmid(FILENAME,0,20)
     print, save_str
     SAVE_FN=SAVE_DIR+save_str+'.sav'

     SAVE, itime1, lon1, lat1, alt1, Eden, Wi, ExB, Rho, FILENAME=SAVE_FN

  ENDFOR

END
