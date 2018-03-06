;+
; :description:
;    This function simulates the whole chain of the data simulation for a given source
;    
; :returns:
;   a structure with: 
;   eventlist         : the time ordered stx_sim_detector_eventlist
;   filtered_eventlist : the time filtered time ordered stx_sim_detector_eventlist
;   triggers           : the triggers generated by the time filtering
;   
;
; :keywords:
;    start_time:   in, optional, type="double", default="0"
;                  the timepoint where does the flare beginn
;
;    src_struct:   in, optional, type="stx_sim_source_structure", default=" tx_sim_source_structure()"
;                  the constant flare struct for the event to simulate
;
;    photoncount:  in, optional, type="ulong", default="same number as the source simulation will produce"
;                  the constant flare struct for the event to simulate
;    
;    energy_profile_type: in, optional, type=string, default='power'
;                        type of distribution
;                        available types: 
;                        'thermal' - "thermal" distribution: F(E) = K * exp( -c * E ) , 
;                        'power'  - power law distribution: F(E) = K * E^(-gamma) ,
;                         
;    energy_profile_param: in, optional, type=array
;                        parameters of the distribution
;                        1) when the type is 'thermal' the param should be an array [ K, c ], default=[ 1, 0.1 ] 
;                        2) when the type is 'power' the param should be single number (value of gamma), default = 3.5              
;                  
;    time_profile_type:  in, optional, type=string, default='gauss'
;                        type of distribution
;                        available types: 
;                        'uniform' - uniform distribution, 
;                        'linear'  - linear rising of the elements density beginning from 0, 
;                        'gauss'   - gaussian distribution with defined mean and FWHM (those paramiters are given by param value),  
;                        'exp'     - exponent rising of the elements density beginning from 0, the base of exponent can be given by param value
;    time_profile_param: in, optional, type=double or dblarr
;                        parameters of the distribution
;                        used only when type 'gauss' or 'exp' is set:
;                        1) when the type is 'gauss' the param should be an array [mean, FWHM] (in seconds), default=[ 50.d, 20.d ] 
;                        2) when the type is 'exp' the param should be signe number which will be the base of the exponent, default=e ( exp(1) )
;                  
;    duration:     in, optional, type="double", default="50"
;                  the time the flare should last in seconds        
;                  
;    T_L :         in, type="double", optional, default="2 microseconds"
;                  the latency time for the detector coincidence time filtering
;                   
;    T_R :         in, type="double", optional, default="10 microseconds"
;                  the read out time for the detector coincidence time filtering
;                  
;    plotting:     in, optional, type="flag", default="off"
;                  do some plotting
; 
; :examples:
;   
;   res = STX_DS_DEMO(/plotting)
;   
; :History:
;   23-mar-2014, richard.schwartz@nasa.gov, changed spelling of filterd to filtered, changed ambiguous keywords
;     by adding _type to time_profile and energy_profile
;   10-apr-2014, richard.schwartz@nasa.gov, used stx_contstruct_sim_detector_eventlist on triggers_out to make them
;     more compatible with accumulator routines which look for a start_time
;   19-may-2014 - Laszlo I. Etesi (FHNW), replaced all "event_list" with "eventlist"
;   28-jul-2014 - Laszlo I. Etesi, using new __define for named structures
;   
;-
function stx_ds_demo, $
    
    start_time = start_time, $
    src_struct = src_struct, $
    photoncount = photoncount, $
    time_profile_type = time_profile, $
    duration = duration, $
    time_profile_param = time_profile_param, $
    energy_profile_param = energy_profile_param, $
    energy_profile_type = energy_profile, $
    T_L = T_L, $
    T_R = T_R, $
    plotting = plotting
    
  default, start_time, 0  
  default, time_profile, "gauss" 
  default, duration, 50d
  default, time_profile_param, [25d,10d]
  default, energy_profile, "power"
  default, energy_profile_param, 3.5 
  default, src_struct, stx_sim_source_structure()
  default, T_L, 2d-6 ;2 microseconds
  default, T_R, 10d-6 ;10 microseconds
  
  
  ;create the photon candidates   
  ph_list = stx_sim_flare(src_struct = src_struct)
  
  
  default, photoncount, n_elements(ph_list)
  
  
  
  ;expand the photon list
  if  photoncount gt n_elements(ph_list) then begin
    ph_list = reproduce(ph_list,ceil(photoncount/float(n_elements(ph_list))))
    ph_list = reform(ph_list,n_elements(ph_list))
  end
  
  ;trim the photon list to the number of photoncount
  ph_list = ph_list[0:photoncount-1]
  
  ;TODO: n.h. change time in seconds to discrete timepoints of the hardware frequency
  ;create and assign time and energy points for all photons
  duration_in = duration
  ph_list.time = stx_sim_time_distribution(data_granulation=data_granulation,nofelem=photoncount,type=time_profile,length=duration_in,param=time_profile_param) * data_granulation
  ph_list.energy = stx_sim_energy_distribution(nofelem=photoncount,type=energy_profile)
  
  ;do the hit test for each photon candidate
  ph_hits = stx_sim_energyfilter(ph_list)
    
  
  ;transform the passed photons into detector events
  events = replicate({stx_sim_detector_event},n_elements(ph_hits))  
  events.relative_time = ph_hits.time
  events.detector_index = ph_hits.subc_d_n
  events.pixel_index = ph_hits.pixel_n
  events.energy_ad_channel = stx_sim_energy2ad_channel(ph_hits.energy)
  
  ;create a event list
  eventlist = stx_construct_sim_detector_eventlist(start_time=start_time, detector_events=events, sources=src_struct)
  
  ;timeorder the eventlist
  eventlist.detector_events = stx_sim_timeorder_eventlist(eventlist.detector_events)
  
  ;applay time filtering to the eventlist
  ;t1 = systime(1)
  filtered_events = stx_sim_timefilter_eventlist(eventlist.detector_events, triggers_out=triggers_out, T_L=T_L, T_R=T_R)
  ;print,'Time for stx_sim_timefilter_eventlist ', systime(1) - t1
   ;create the filtered event list
  filtered_eventlist = stx_construct_sim_detector_eventlist(start_time=start_time, detector_events = filtered_events, sources = src_struct)
  trigger_eventlist = stx_construct_sim_detector_eventlist( start_time = start_time, detector_events = triggers_out, sources = src_struct )
  ;prepare the output struct  
;  out = { $
;    type                : "stx_sim_datasimulation_output", $
;    filtered_eventlist  : filtered_eventlist, $
;    eventlist          : eventlist, $
;    triggers            : triggers_out $ 
;  }
   out = { $
    type                : "stx_sim_datasimulation_output", $
    filtered_eventlist  : filtered_eventlist, $
    eventlist          : eventlist, $
    trigger_eventlist  : trigger_eventlist $ 
  }
 
  
  ;do some plotting if desired
  if keyword_set(plotting) then begin
    !P.multi = [0,1,2]
    ;stx_sim_timefilter_eventlist_plot,eventlist.detector_events, filtered_eventlist.detector_events, triggers_out, t_l=t_l, t_r=t_r, timerange=[0,duration]
    
    binsize=0.01
    
    hist_ev = histogram(out.eventlist.detector_events.relative_time,binsize=binsize, locations=loc_ev)
    plot, loc_ev, hist_ev, title="detector coincidence time filtering", xtitle="time in seconds", ytitle="# / "+trim(binsize)+"s", /ylog
    
    hist_fiev = histogram(out.filtered_eventlist.detector_events.relative_time,binsize=binsize, locations=loc_fiev)
    oplot, loc_fiev, hist_fiev, linestyle=1
    
    hist_tr = histogram(triggers_out.relative_time,binsize=binsize, locations=loc_t)
    oplot, loc_t, hist_tr, linestyle=2
   
    al_legend, ["events","passed events","triggers","latency time = "+trim(t_l)+"s","read out time = "+trim(t_r)+"s"], linestyle=[0,1,2,-1,-1]
    
    
    hist = histogram(out.filtered_eventlist.detector_events.energy_ad_channel,binsize=1, locations=loc)
    plot, loc, hist, /ylog, yrange=[10,max(hist)]
    
       
    !P.multi = -1    
  endif
  
  return, out
end