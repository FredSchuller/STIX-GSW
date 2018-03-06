;+
; :description:
;   This function converts random events that were generated by Gordon's
;   'stixfsw_randomtest.pro' routine to stix_sim_detector_eventlists
;   
;   ** Structure <1165f780>, 4 tags, length=16, data length=14, refs=1:
;      DET             INT              2
;      PIXEL           INT             10
;      ADC2            INT            725
;      RELTIME         DOUBLE           1.000000
;
; :categories:
;    flight software simulator, converter
;
; :returns:
;    a stx_sim_detector_eventlist structure containing randomly
;    generated detector events
;
; :examples:
;    eventlist = stx_fsw_convert_rnd_events2eventlist(rnd_events=rnd_events)
;
; :history:
;     03-jul-2015, Laszlo I. Etesi (FHNW), initial release
;
;-
function stx_fsw_convert_rnd_events2eventlist, rnd_events=rnd_events
  detector_events = stx_construct_sim_detector_event($
    relative_time = rnd_events.reltime, $
    detector_index = rnd_events.det, $
    pixel_index = rnd_events.pixel, $
    energy_ad_channel = rnd_events.adc2)
  
  return, stx_construct_sim_detector_eventlist(detector_events=detector_events, start_time=stx_time())
end