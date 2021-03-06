;---------------------------------------------------------------------------
;       Classification Extension-based Cloud Detection (CECD v1.0) method
;
;INPUT:
; 1) Resampled 300m Landsat/Sentinel-2 Image
; 2) Cloud-free Proba-V image
; 3) Parameter: Band index of Blue, Red, NIR, and SWIR bands of Landsat/Sentinel-2 image
;   For instance: the list of bands of Landsat-8 OLI is coastal, blue, green, red, nir, swir1, swir2, so the input band index is 1,3,4,5
;
;OUTPUT:
; 1) THOT image
; 2) THMFCM image, THOT image after MFCM correction
; 3) Training roi, Training Samples
;
;AUTHOR:
; Xidong Chen (chenxd@radi.ac.cn), Liangyun Liu (liuly@radi.ac.cn), Yuan Gao.
;
;---------------------------------------------------------------------------
pro CECD_Tool_Main,poss1,poss2,reproj_cloud_file,ref_file,dir
  compile_opt idl2
  envi_open_file,reproj_cloud_file,r_fid = HIF_fid
  envi_open_file,ref_file,r_fid = RIF_fid
  map_info = envi_get_map_info(fid=HIF_fid)
  THOT,dir,HIF_fid,RIF_fid,poss1,poss2,THOT_file=THOT_file,data_cloudy=data_cloudy
  MFCM,dir,THOT_file,data_cloudy,HIF_fid,THMFCM_file=THMFCM_file
  Sample_Collection,dir,THMFCM_file,data_cloudy
  print,'Finished CECD Sampling!'
end

pro THOT,dir,HIF_fid,RIF_fid,poss1,poss2,THOT_file=THOT_file,data_cloudy=data_cloudy
  compile_opt idl2
  ;-------------------------HOT----------------------;
  envi_file_query,HIF_fid,ns=ns,nl=nl,dims=dims
  n_b = n_elements(poss1) ; number of overlap bands
  data_clear=dblarr(ns*nl,n_b);
  data_cloudy=dblarr(ns*nl,n_b)
  for i_band =0,n_b-1 do begin
    data_clear[*,i_band]=envi_get_data(fid=RIF_fid,dims=dims,pos=poss1[i_band])
    data_cloudy[*,i_band]=envi_get_data(fid=HIF_fid,dims=dims,pos=poss2[i_band])
  endfor
  mask_area=where(data_cloudy[*,0] gt 0, masknum)
  index = where((data_cloudy[*,0]-0.5*data_cloudy[*,1]) gt 0.075,cloud_count)
  result = regress(double(data_clear[mask_area,0]), double(data_clear[mask_area,1]), SIGMA=sigma, CONST=const,MEASURE_ERRORS=measure_errors,CORRELATION=correlation)
  cosin = [sqrt((result*result)/(1+result*result)),sqrt(1/(1+result*result))]
  if result le 1.4 then cosin=[1,0.5]
  ;------------------------THOT----------------------;
  data_cloudy_clear_difference=dblarr(ns*nl,n_b)
  HOT_cloudy = dblarr(ns*nl)
  for i_band =0, n_b-1 do data_cloudy_clear_difference[*,i_band]=mean_filter(reform(data_cloudy[*,i_band],ns,nl),3)-mean_filter(reform(data_clear[*,i_band],ns,nl),3)
  HOT_cloudy[mask_area] = data_cloudy[mask_area,0]*cosin[0] - data_cloudy[mask_area,1]*cosin[1];fix((data_cloudy[*,*,0]*sinx[0] - data_cloudy[*,*,1]* cosx[0])*100.0)
  measure_errors = REPLICATE(0.5, N_ELEMENTS(HOT_cloudy))
  mask_area = where(data_cloudy_clear_difference[*,0] gt 0,masknum)
  THOT_result = REGRESS(double(Transpose(data_cloudy_clear_difference[mask_area,*])), double(HOT_cloudy[mask_area]), $
    SIGMA=sigma, CONST=const,MEASURE_ERRORS=measure_errors,CORRELATION=correlation,/double)
  ;THOT
  THOT=dblarr(masknum,1)
  for i_band=0,n_b-1 do THOT=THOT+data_cloudy_clear_difference[mask_area,i_band]*THOT_result[i_band]
  THOT=THOT+const
  THOT_image=fltarr(ns*nl)
  THOT_image[mask_area] = THOT
  ;-------------refine THOT---------------;
  ;Due to the problem of geometric matching, there may be a deviation of several pixels between a pair of pixels, resulting in some abnormal values in THOT, so use some simple rules to refine the THOT
  Bdiff = data_cloudy-data_clear
  envi_enter_data,reform(HOT_cloudy,ns,nl),r_Fid=HOT_fid
  envi_doit,'class_doit',fid=HOT_fid,pos=0,dims=dims,out_bname='Kmeans_HOT',$
    method=7,r_fid=PCPS_fid,change_thresh=0.05,num_classes=2,iterations=1,/in_memory
  HOT_kmeans = envi_get_data(fid=PCPS_fid,dims=dims,pos=0)
  mean_value = mean(THOT_image[where(HOT_kmeans eq 1 and data_cloudy[*,0] gt 0)])
  THOT_image[where(data_cloudy[*,3] lt 0.09 or Bdiff[*,3] lt 0 or Bdiff[*,0] lt -0.05 or (data_cloudy[*,0] lt 0.09 and data_cloudy[*,0] gt 0))] = mean_value
  THOT_file = dir+'THOT.tif'
  envi_file_mng,id = PCPS_fid,/remove
  envi_file_mng,id = HOT_fid,/remove
  openw,lun,THOT_file,/get_lun
  writeu,lun,THOT_image
  free_lun,lun
  envi_setup_head,fname=THOT_file,ns=ns,nl=nl,nb=1,data_type=size(THOT_image,/type),$
    interleave=0,offset=0,map_info=map_info,/write
End

pro MFCM,dir,THOT_file,data_cloudy,HIF_fid,THMFCM_file=THMFCM_file
  compile_opt idl2
  envi_open_file,THOT_file,r_fid=THOT_Fid
  envi_file_query,THOT_fid,ns=l8_ns,nl=l8_nl,dims=THOT_dims
  THOT_data = envi_get_data(fid=THOT_Fid,dims=THOT_dims,pos=0)
  valid_area = where(data_cloudy[*,0] eq 0)
  kernel = [[0.125,0.125,0.125],[0.125,0,0.125],[0.125,0.125,0.125]]
  ;-------------------parameter init----------------------;
  m = 2
  max_ite = 200
  neighbor_inf = 0.3
  obj_the = 500
  sigma = 0.001
  C_0_back = 0.0
  C_1_back = 0.0
  J_back = 0.0
  betaa = fltarr(l8_ns*l8_nl)
  k=0
  J = fltarr(max_ite)
  C_0 = fltarr(max_ite)
  C_1 = fltarr(max_ite)
  ;----------------------MFCM--------------------;
  HOT_empiral = data_cloudy[*,0]-0.5*data_cloudy[*,1]
  mask = bytarr(l8_ns,l8_nl)
  mask[where(HOT_empiral gt 0.04)] = 1
  envi_enter_data,mask,r_fid=mask_fid
  envi_doit,'class_doit',fid=THOT_Fid,pos=0,dims=THOT_dims,out_bname='Kmeans_THOT',m_fid=mask_fid,m_pos=0,$
    method=7,r_fid=PCPS_fid,change_thresh=0.05,num_classes=2,iterations=1,/in_memory
  kmeans_THOT = envi_get_data(fid=PCPS_fid,dims=THOT_dims,pos=0)
  U = fltarr(l8_ns*l8_nl)+0.4
  U[where(kmeans_THOT eq 2)]=1
  envi_file_mng,id = PCPS_fid,/remove
  envi_file_mng,id = mask_fid,/remove
  while(k lt max_ite) do begin
    x = THOT_data - betaa
    x_averNeighbor = reform(convol(reform(THOT_data,l8_ns,l8_nl),kernel,/center),1,l8_ns*l8_nl)
    ;----------update cluster central----------;
    C_0[k] = total((x+neighbor_inf*x_averNeighbor)*(U^m),/NaN) / ((1+neighbor_inf)*total(U^m,/NaN))
    C_1[k] = total((x+neighbor_inf*x_averNeighbor)*((1-U)^m),/NaN) / ((1+neighbor_inf)*total((1-U)^m,/NaN))
    Cmatrix_0 = make_array(l8_ns,l8_nl,value=C_0[k])
    Cmatrix_1 = make_array(l8_ns,l8_nl,value=C_1[k])
    ;----------update betaa------------;
    betaa = THOT_data - ((U^m)*Cmatrix_0 + (1-U)^m*Cmatrix_1)/((U^m) + (1-U)^m)
    ;------------update U----------------;
    U = ((x-Cmatrix_0)^2+neighbor_inf*(x_averNeighbor-Cmatrix_0)^2)^(-1/(m-1))/$
      (((x-Cmatrix_0)^2+neighbor_inf*(x_averNeighbor-Cmatrix_0)^2)^(-1/(m-1))+$
      ((x-Cmatrix_1)^2+neighbor_inf*(x_averNeighbor-Cmatrix_1)^2)^(-1/(m-1)))
    ;------------Calculate  j------------;
    J[k] = total(((x-Cmatrix_0)^2)*(U^m),/NaN)+ total(((x-Cmatrix_1)^2)*((1-U)^m),/NaN) + $
      neighbor_inf*total(((x_averNeighbor-Cmatrix_0)^2)*(U^m),/NaN)+ neighbor_inf*$
      total(((x_averNeighbor-Cmatrix_1)^2)*((1-U)^m),/NaN)
    if((abs(C_0[k]-C_0_back) lt sigma)and (abs(C_1[k]-C_1_back) lt sigma)) then break
    if(abs(J[k]-J_back) lt obj_the) then break
    C_0_back = C_0[k]
    C_1_back = C_1[k]
    J_back = J[k]
    k=k+1
    print,k
  endwhile
  ;---------------------Correct THOT--------------------;
  x_corrected = THOT_data-betaa
  x_corrected[where(THOT_data eq 0)] =0
  THMFCM_file = dir+'THMFCM.tif'   ; the image derived further improving the contrast between clouds and surfaces using MFCM
  openw,lun,THMFCM_file,/get_lun
  writeu,lun,x_corrected
  free_lun,lun
  envi_setup_head,fname=THMFCM_file,ns=l8_ns,nl=l8_nl,nb=1,data_type=size(x_corrected,/type),$
    interleave=0,offset=0,map_info=map_info,/write
end

pro Sample_Collection,dir,THMFCM_file,data_cloudy
  compile_opt idl2
  envi_open_file,THMFCM_file,r_fid=x_fid
  envi_file_query,x_fid,dims=THOT_dims,ns=l8_ns,nl=l8_nl
  ;--------------OTSU Classify cloud------------; For convenience here, we use kmeans instead of otsu
  envi_doit,'class_doit',fid=x_fid,pos=[0],dims=THOT_dims,out_bname='THMFCM_Kmeans',$
    method=7,r_fid=THMFCM_fid,change_thresh=0.05,num_classes=2,iterations=1,/in_memory
  kmeans_HOT = envi_get_data(fid=THMFCM_fid,dims=THOT_dims,pos=0)
  cloud = bytarr(l8_ns,l8_nl)
  cloud[where(kmeans_HOT eq 2 and data_cloudy[*,0] gt 0)] = 1
  confident_zone,cloud,zone=cloud_zone
  bright = bytarr(l8_ns,l8_nl)
  HOT_emperical = data_cloudy[*,0]-0.5*data_cloudy[*,1]
  bright[where(HOT_emperical gt 0.05 and cloud eq 0)] = 1
  confident_zone,bright,zone=bright_zone
  bright = !null
  surface_cover = bytarr(l8_ns,l8_nl)
  surface_cover[where(kmeans_HOT eq 1 and data_cloudy[*,0] gt 0 and cloud eq 0)] = 1
  confident_zone,surface_cover,zone=surface_zone
  surface_cover = !null
  envi_file_mng, id = THMFCM_fid,/remove
  cloud_index = where(cloud_zone eq 1,cloud_count)
  valid_index = where(data_cloudy[*,0] gt 0,valid_count)
  ncloud_count = valid_count-cloud_count
  cloud_number = 300*cloud_count*1.0/valid_count
  noncloud_number = 300-cloud_number
  cloud_number = cloud_number>200
  noncloud_number = noncloud_number <300
  cl_index = []
  nc_index = []
  ;----------cloud-------------;
  cl_start_index = [0.08,0.2,0.3,0.4,0.5,0.6,0.7]
  cl_end_index = [0.2,0.3,0.4,0.5,0.6,0.7,1.8]
  cl_bin_count = n_elements(cl_start_index)
  for i_bin =0, cl_bin_count-1 do begin
    cl_index_ith = (where(cloud_zone eq 1 and data_cloudy[*,0] ge cl_start_index[i_bin] and data_cloudy[*,0] lt cl_end_index[i_bin], cl_count))[long(randomu(random,cl_count*1.0/cloud_count*cloud_number>20,1)*cl_count)]
    if cl_count gt 5 then cl_index = [cl_index,cl_index_ith]
  endfor
  ;----non cloud-------;
  nc_start_index = [0,0.1,0.2,0.3,0.5,0.8]
  nc_end_index = [0.1,0.2,0.3,0.5,0.7,1.8]
  nc_bin_count = n_elements(nc_start_index)
  for i_bin =0, nc_bin_count-1 do begin
    nc_index_ith = (where(cloud_zone ne 1 and data_cloudy[*,0] gt nc_start_index[i_bin] and data_cloudy[*,0] lt nc_end_index[i_bin], nc_count))[long(randomu(random,nc_count*1.0/ncloud_count*noncloud_number>20,1)*nc_count)]
    if nc_count gt 5 then nc_index = [nc_index,nc_index_ith]
  endfor
  nc_index_ith = (where(cloud_zone ne 1 and data_cloudy[*,2] gt 0 and data_cloudy[*,2] lt 0.18 and data_cloudy[*,0] gt 0.25, nc_count))[long(randomu(random,nc_count*1.0/ncloud_count*noncloud_number>40,1)*nc_count)]
  if nc_count gt 5 then nc_index = [nc_index,nc_index_ith]
  ;---------------save roi-----------------;
  cl_xy=array_indices(reform(data_cloudy[*,0],l8_ns,l8_nl),cl_index)
  cl_id= ENVI_CREATE_ROI(color=2,ns=l8_ns,nl=l8_nl,name='cloud')
  ENVI_DEFINE_ROI,cl_id,/point,xpts=cl_xy[0,*],ypts=cl_xy[1,*]
  nc_xy=array_indices(reform(data_cloudy[*,0],l8_ns,l8_nl),nc_index)
  nc_id= ENVI_CREATE_ROI(color=9,ns=l8_ns,nl=l8_nl,name='clear')
  ENVI_DEFINE_ROI,nc_id,/point,xpts=nc_xy[0,*],ypts=nc_xy[1,*]
  roi_ids = [cl_id,nc_id]
  roi_dir = dir+'Training_Roi.roi'
  envi_save_rois,roi_dir,roi_ids
  envi_delete_rois,roi_ids
end

pro confident_zone,data,zone=zone
  compile_opt idl2
  ns = (size(data,/dimension))[0]
  nl = (size(data,/dimension))[1]
  erode_cloud = erode(data,[5,5])
  kernelSize = [5,5]
  kernel_1 = bytarr(kernelSize[0], kernelSize[1])
  kernel_2 = bytarr(kernelSize[0], kernelSize[1])
  kernel_3 = bytarr(kernelSize[0], kernelSize[1])
  kernel_4 = bytarr(kernelSize[0], kernelSize[1])
  kernel_5 = bytarr(kernelSize[0], kernelSize[1])
  kernel_6 = bytarr(kernelSize[0], kernelSize[1])
  kernel_7 = bytarr(kernelSize[0], kernelSize[1])
  kernel_8 = bytarr(kernelSize[0], kernelSize[1])
  kernel_1[0,0] = 1
  kernel_2[4,0] = 1
  kernel_3[0,4] = 1
  kernel_4[4,4] = 1
  kernel_5[2,0] = 1
  kernel_6[2,4] = 1
  kernel_7[0,2] = 1
  kernel_8[4,2] = 1
  LFE1 = erode_cloud-convol(erode_cloud, kernel_1, /CENTER)
  LFE2 = erode_cloud-convol(erode_cloud, kernel_2, /CENTER)
  LFE3 = erode_cloud-convol(erode_cloud, kernel_3, /CENTER)
  LFE4 = erode_cloud-convol(erode_cloud, kernel_4, /CENTER)
  LFE5 = erode_cloud-convol(erode_cloud, kernel_5, /CENTER)
  LFE6 = erode_cloud-convol(erode_cloud, kernel_6, /CENTER)
  LFE7 = erode_cloud-convol(erode_cloud, kernel_7, /CENTER)
  LFE8 = erode_cloud-convol(erode_cloud, kernel_8, /CENTER)
  LFE_max = max([[abs(reform(LFE1,ns*nl))],[abs(reform(LFE2,ns*nl))],[abs(reform(LFE3,ns*nl))],[abs(reform(LFE4,ns*nl))],$
    [abs(reform(LFE5,ns*nl))],[abs(reform(LFE6,ns*nl))],[abs(reform(LFE7,ns*nl))],[abs(reform(LFE8,ns*nl))]],dimension=2,Min_Subscrip)
  index = where(LFE_max eq 0 and erode_cloud eq 1)
  LFE_max = !null
  zone = bytarr(ns*nl)
  zone[index] = 1
end

pro CECD_SampleGeneration_GUIResponse,event
  compile_opt idl2
  uname_id=widget_info(event.id,/uname)
  case fix(uname_id) of
    1: begin    ;Landsat/Sentinel-2
      L8_filepath=dialog_pickfile(title='Landsat/Sentinel-2 image')
      if L8_filepath eq '' then return
      L8_file_id=widget_info(event.top,find_by_uname='L8_filepath')
      widget_control,L8_file_id,set_value=L8_filepath
    end
    2: begin    ;proba-v
      Proba_filepath=dialog_pickfile(title='Proba-V image')
      if Proba_filepath eq '' then return
      Proba_file_id=widget_info(event.top,find_by_uname='Proba_filepath')
      widget_control,Proba_file_id,set_value=Proba_filepath
    end
    3: begin    ;output
      output=dialog_pickfile(title='Output Path',/directory)
      if output eq '' then return
      output_id=widget_info(event.top,find_by_uname='output')
      widget_control,output_id,set_value=output
    end
    4: begin    ;excuate
      Landsat_band_id=widget_info(event.top,find_by_uname='Landsat_band')
      widget_control,Landsat_band_id,get_value=Landsat_band
      if Landsat_band eq '' then return
      Landsat_band = fix(strsplit(Landsat_band[0],',',/extract))
      ;
      Proba_band_id=widget_info(event.top,find_by_uname='Proba_band')
      widget_control,Proba_band_id,get_value=Proba_band
      if Proba_band eq '' then return
      Proba_band = fix(strsplit(Proba_band[0],',',/extract))
      ;
      L8_file_id=widget_info(event.top,find_by_uname='L8_filepath')
      widget_control,L8_file_id,get_value=L8_filepath
      if L8_filepath eq '' then return
      ;
      proba_file_id=widget_info(event.top,find_by_uname='proba_filepath')
      widget_control,proba_file_id,get_value=proba_filepath
      if proba_filepath eq '' then return
      ;
      output_id=widget_info(event.top,find_by_uname='output')
      widget_control,output_id,get_value=output
      if output eq '' then return
      ;
      widget_control,event.Top,/destroy
      CECD_Tool_Main,Landsat_band,Proba_band,L8_filepath,proba_filepath,output
    end
    5: begin    ;end
      widget_control,event.Top,/destroy
    end
  endcase
end

;save,filename='H:\MODIS_Reprojection_Tool.sav',/routines
pro CECD_SampleGeneration_Tool_V1
  compile_opt idl2
  envi,/restore_base_save_Files
  envi_batch_init
  ;GUI
  Landsat_band ='1,3,4,5' ; the band index of blue, red, nir, and swir band of Landsat/Sentinel-2 image
  Proba_band = '0,1,2,3'
  L8_dir=''
  Proba_dir=''
  outdir = ''
  tlb=widget_base(xsize=420,ysize=260,xoffset=600,yoffset=200,title='CECD SampleGeneration V1',kill_notify='')
  ;CECD inputs
  TileControlsLabel = widget_label(tlb,yoffset=110,xoffset=15, value='Select the index of Red, Blue, NIR, and SWIR band')
  CECD_input_base=widget_base(tlb,xsize=400,ysize=70,xoffset=10,yoffset=25,frame=1,map=1,uname='CECD_input_base')
  ;Landsat-8/Sentinel-2 input
  L8_label=widget_label(CECD_input_base,xsize=150,ysize=30,xoffset=10,yoffset=15,value='Landsat/Sentinel-2 Image')
  L8_filepath=widget_text(CECD_input_base,xsize=25,ysize=1,xoffset=160,yoffset=8,uname='L8_filepath',value=L8_dir)
  L8_input_button=widget_button(CECD_input_base,xsize=30,ysize=30,xoffset=340,yoffset=5,value='...',uname=1,event_pro='CECD_SampleGeneration_GUIResponse')
  ;Proba-V input
  proba_label=widget_label(CECD_input_base,xsize=150,ysize=30,xoffset=10,yoffset=45,value='Cloud-free Proba-V Image')
  proba_filepath=widget_text(CECD_input_base,xsize=25,ysize=1,xoffset=160,yoffset=38,uname='proba_filepath',value=proba_dir)
  proba_input_button=widget_button(CECD_input_base,xsize=30,ysize=30,xoffset=340,yoffset=35,value='...',uname=2,event_pro='CECD_SampleGeneration_GUIResponse')
  ;
  TileControlsLabel2 = widget_label(tlb,yoffset=10,xoffset=160, value='CECD input images') 
  CECD_Set_base=widget_base(tlb,xsize=400,ysize=40,xoffset=10,yoffset=125,frame=1,map=1,uname='CECD_Set_base')
  ;Landsat bands
  Landsat_label=widget_label(CECD_Set_base,xsize=110,ysize=30,xoffset=2,yoffset=15,value='Landsat band index')
  Landsat_band=widget_text(CECD_Set_base,xsize=10,ysize=1,xoffset=120,yoffset=8,uname='Landsat_band',value=Landsat_band,/EDITABLE)
  ;Proba-V bands
  Proba_base=widget_label(CECD_Set_base,xsize=110,ysize=30,xoffset=200,yoffset=15,value='Proba-V band index')
  Proba_band=widget_text(CECD_Set_base,xsize=10,ysize=1,xoffset=320,yoffset=8,uname='Proba_band',value=Proba_band,/EDITABLE)
  ;Output
  TileControlsLabel3 = widget_label(tlb,yoffset=170,xoffset=160, value='Output Path Selection') 
  output_base=widget_base(tlb,xsize=400,ysize=30,xoffset=10,yoffset=185,frame=1,map=1,uname='output_base')
  out_label=widget_label(output_base,xsize=70,ysize=30,xoffset=10,yoffset=10,value='Output path')
  output=widget_text(output_base,xsize=25,ysize=1,xoffset=90,yoffset=5,uname='output',value=outdir)
  output_button=widget_button(output_base,xsize=30,ysize=30,xoffset=260,yoffset=1,value='...',uname=3,event_pro='CECD_SampleGeneration_GUIResponse')
  ;Apply
  apply_button=widget_button(tlb,xsize=60,ysize=30,xoffset=140,yoffset=225,value='OK',uname=4,event_pro='CECD_SampleGeneration_GUIResponse')
  cancel_button=widget_button(tlb,xsize=60,ysize=30,xoffset=240,yoffset=225,value='Cancel',uname=5,event_pro='CECD_SampleGeneration_GUIResponse')
  widget_control,tlb,/realize
  Xmanager,'CECD_SampleGeneration_Tool_V1',tlb
end