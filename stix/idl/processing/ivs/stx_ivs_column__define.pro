;---------------------------------------------------------------------------
; Document name: stx_ivs_column__define.pro
; Created by:    Nicky Hochmuth, 2012/03/05
;---------------------------------------------------------------------------
;+
; PROJECT:
;       STIX
;
; PURPOSE:
;       Performes a recursive interval selection on a spetrogram column
;       each row ist splited on the barycenter
;
; CATEGORY:
;       Stix on Bord Algorithm
;
; CALLING SEQUENCE:
;       obj = stx_ivs_column(start_time, end_time, rows, spectrogram_ptr,split_level, thermalboundary, min_time, min_count)
;       intervals = obj->get_intervals()
;
; HISTORY:
;       2012/03/05, Nicky.Hochmuth@fhnw.ch, initial release
;
;-

;+
; :DESCRIPTION:
;    Initializes the stx_ivs_column object
;    
; :PARAMS:
;   start_time: the index on the time axis of the spectrogram the column starts
;   end_time: the index on the time axis of the spectrogram the column ends
;   rows: an array of all indexes on the energy axis this column is created for
;   spectrogram: a pointer to the spectrogram
;
; :AUTHOR: nicky.hochmuth
;-
function stx_ivs_column::init, start_time, end_time, rows, spectrogram, split_level, thermalboundary, min_time, min_count
   
   ;do some parameter testing
   if ~ppl_typeof(spectrogram, compareto='pointer') then begin
    message, 'Spectrogram has to by a pointer', /continue
    return, 0
   end
   
   self.thermalboundary    = thermalboundary
   self.min_time          = min_time
   self.min_count         = min_count
   
   self.start_time        = start_time
   self.end_time          = end_time
   self.split_level       = split_level
   
   self.spectrogram       = spectrogram
   self.rows              = ptr_new(rows)
   self.ni_values         = ptr_new(bytarr(n_elements((*spectrogram).e_axis.mean)))
   
   ;fill the ni_values for each cell
   self->setNiValues
  return, 1
end

;+
; :description:
;    Cleanup of this class
;-
pro stx_ivs_column::cleanup
  free_pointer, self.rows
  ;free_pointer, self.spectrogram
  free_pointer, self.ni_values
end


;+
; :DESCRIPTION:
;    determines the Ni value for each cell and stores it in the internal 
;
; :AUTHOR: nicky.hochmuth
;-
pro stx_ivs_column::setNiValues
   compile_opt idl2, hidden
    for row_idx=0, n_elements(*self.rows)-1 do begin
      e=(*self.rows)[row_idx]
      
      totalcount = total((*self.spectrogram).data[e, self.start_time:self.end_time],/integer)
      time_span =  stx_time_diff((*self.spectrogram).t_axis.time_start[self.end_time],(*self.spectrogram).t_axis.time_end[self.start_time])
      energy = (*self.spectrogram).e_axis.low[e]
      
      (*self.ni_values)[e]=self->count_lookup(totalcount,energy,time_span)
     endfor
end

;+
; :DESCRIPTION:
;    tests if the split on the barycenter of all ni2 cell will produce ni0 cell on the left or right side if the split
;    if so: change from ni2 state to ni1 state 
;
; :AUTHOR: nicky.hochmuth
;-
pro stx_ivs_column::resetNiValues,idx, split_position
   compile_opt idl2, hidden
   for row_idx=0, n_elements(idx)-1 do begin
     e = idx[row_idx]
     energy = (*self.spectrogram).e_axis.low[e]
     
     start_time_left = self.start_time
     end_time_left = split_position
     totalcount_left = total((*self.spectrogram).data[e, start_time_left:end_time_left])
     time_span_left = stx_time_diff((*self.spectrogram).t_axis.time_end[end_time_left] , (*self.spectrogram).t_axis.time_start[start_time_left])
     ni_value_left = self->count_lookup(totalcount_left,energy,time_span_left)
     
     start_time_right = split_position+1
     end_time_right = self.end_time
     totalcount_right = total((*self.spectrogram).data[e, start_time_right:end_time_right])
     time_span_right = stx_time_diff((*self.spectrogram).t_axis.time_end[end_time_right] , (*self.spectrogram).t_axis.time_start[start_time_right])
     ni_value_right = self->count_lookup(totalcount_right,energy,time_span_right)
     
     old_ni_value = (*self.ni_values)[idx[row_idx]] 
     
     
     switch (old_ni_value) of
       2: begin
         (*self.ni_values)[idx[row_idx]] = (ni_value_left ge 1 && ni_value_right ge 1) $
                                            ;|| (self.split_level eq 0 && (ni_value_left ge 1 || ni_value_right ge 1)) $
                                            ? 2 : 1 
         break
       end
       
       else: begin
       end
     endswitch
     
   end
end




;+
; :DESCRIPTION:
;    in case that there are to less counts for any splitting return the full intervall for merging 
;
; :AUTHOR: nicky.hochmuth
;-
function stx_ivs_column::get_dafault_intervals
  
  low_high = minmax(*self.rows)
  count = total((*self.spectrogram).data[*self.rows,self.start_time:self.end_time], /integer)
  
  return, stx_construct_ivs_interval( (*self.spectrogram).t_axis.time_start[self.start_time], (*self.spectrogram).t_axis.time_end[self.end_time], $
                                      start_time_idx=self.start_time, end_time_idx=self.end_time, $
                                      (*self.spectrogram).e_axis.low[low_high[0]], (*self.spectrogram).e_axis.high[low_high[1]], $
                                      start_energy_idx=(*self.spectrogram).e_axis.HIGH_FSW_IDX[low_high[0]], end_energy_idx=(*self.spectrogram).e_axis.HIGH_FSW_IDX[low_high[1]], $ 
                                      count,trim=0,spectroscopy=0)  
    
  return, intervals
end

;+
; :DESCRIPTION:
;    performes the rekursive interval selection by spliting the time intervals into halves 
;
; :KEYWORDS:
;    split_times: keeps trak of the time split positions for later trimming
;
; :AUTHOR: nicky.hochmuth
;-
function stx_ivs_column::get_intervals, split_times = split_times
  
  ;transform the start and end index to time values 
  start_time = (*self.spectrogram).t_axis.time_start[self.start_time] 
  end_time = (*self.spectrogram).t_axis.time_end[self.end_time]
  
  ;find all cells whitch could split up into halves images 
  idx_2er_img_old = where((*self.ni_values) ge 2 and (*self.ni_values) lt 6, count_2er_img)
  ;idx_2er_img_old = where(*self.ni_values eq 4, count_2er_img)
  
  split_time = -1
  
  if count_2er_img gt 0 then begin
    split_time = self->getbarycenter(idx_2er_img_old)
    ;TODO: n.h check bounds propertly old: split_time+1
    if split_time lt self.end_time then self->resetnivalues, idx_2er_img_old, split_time
    idx_2er_img = where((*self.ni_values) ge 2 and (*self.ni_values) lt 6, count_2er_img)
    ;idx_2er_img = where(*self.ni_values eq 4, count_2er_img)
  endif
  
  ;flag for the case that an cell has a ni_value of 2 but the time range is only a single time_bin 
  ;therefore the cell could not split into halves
  add_2er_cells = 0
   
  intervals = []
   
  ;do ni4 cells exist
  if count_2er_img gt 0 then begin
      ;The remaining time  /energy bins have n1 >1. (orange bins in Figure 5) The counts in all the bins for which Ni>1 are summed over energy.   
      ;Then a common time-division is identified such that there is an equal (or as close as possible to equal) number of energy-summed 
      ;(over the Ni>1 bins) counts in each. (Figure 6) The subalgorithm is then applied to each one of these two bins at each energy. (Figure 7).
      ;->find the barycentre
      
    
    
    ;does a right side exist?
    ;TODO: n.h check bounds propertly old: split_time+1
    if split_time lt self.end_time then begin
      
      ;log the split times 
      split_times = [[split_times], [self.split_level,split_time]]
      
      ;split the column into two halves and get the intervals of each side recursive  
      left = stx_ivs_column(self.start_time,split_time,idx_2er_img,self.spectrogram,level = self.split_level+1,thermalboundary=self.thermalboundary,min_time=self.min_time,min_count=self.min_count)
      right = stx_ivs_column(split_time+1,self.end_time,idx_2er_img,self.spectrogram,level = self.split_level+1,thermalboundary=self.thermalboundary,min_time=self.min_time,min_count=self.min_count)
      
      intervals_left = left->get_intervals(split_times = split_times)
      intervals_right = right->get_intervals(split_times = split_times)
      
      ;concat the left and right intervals to the interval list
      intervals = [intervals,intervals_left,intervals_right]  
       
      ;print, interval_duplicates(merged_intervals) 
      
    endif else begin ;no right side
        add_2er_cells = 1
    end   
  endif


  row_idx = 0
  
  n_rows = n_elements(*self.rows)
  
  ;go thru each cell and decide what doto
  while row_idx lt n_rows do begin

    e = (*self.rows)[row_idx]
    
    ni_value = (*self.ni_values)[e]
    
    ;determine if it is an cadidate for further triming and in what direction
    ;only intervals from the very beginning or to the very end are condidates for trimming
    trim = 0 ; no trimming
    if stx_time_eq(start_time ,(*self.spectrogram).t_axis.time_start[0]) then trim = 1 ;left trimming
    if stx_time_eq(end_time ,  (*self.spectrogram).t_axis.time_end[-1]) then trim = 2 ;right trimming
     
    ;ni=2
    ;all ni=2 cell should have been splited into halves allready but if the column only have a single time bin this it not posible
    if add_2er_cells && ni_value eq 2  then begin
      ;transform the cell to an interval/image 
      intervals = [intervals, stx_construct_ivs_interval(start_time, end_time, $
                                                         start_time_idx=self.start_time , end_time_idx=self.end_time, $
                                                         (*self.spectrogram).e_axis.low[e] ,(*self.spectrogram).e_axis.high[e], $
                                                         start_energy_idx=(*self.spectrogram).e_axis.LOW_FSW_IDX[e], end_energy_idx=(*self.spectrogram).e_axis.HIGH_FSW_IDX[e], $
                                                         total((*self.spectrogram).data[e, self.start_time:self.end_time],/integer),spectroscopy=0,trim=trim)]
                                                         
      ;marc the cell as 0 for later spectroskopy splitting
      (*self.ni_values)[e]=0
      row_idx++
      continue
    endif 
    
    ;ni=1
    if ni_value eq 1  then begin
      intervals = [intervals, stx_construct_ivs_interval(start_time, end_time, $
                                                         start_time_idx=self.start_time , end_time_idx=self.end_time, $
                                                         (*self.spectrogram).e_axis.low[e] ,(*self.spectrogram).e_axis.high[e], $
                                                         start_energy_idx=(*self.spectrogram).e_axis.LOW_FSW_IDX[e], end_energy_idx=(*self.spectrogram).e_axis.HIGH_FSW_IDX[e], $
                                                         total((*self.spectrogram).data[e, self.start_time:self.end_time],/integer),spectroscopy=0,trim=trim)]
      row_idx++
      continue
    endif

    row_idx++
  end
  
;  if interval_duplicates(intervals) then begin
;    print, 1
;  end 
  
  ;animation of finding the intervals
  ;stx_interval_plot, (*self.spectrogram), intervals=intervals,/overplot 
  
  
  return,   intervals
  
end

function stx_ivs_column::get_NiValues
  return, (*self.ni_values)
end

function stx_ivs_column::get_starttime_idx
  return, self.start_time
end

function stx_ivs_column::get_endtime_idx
  return, self.end_time
end

function stx_ivs_column::get_spectrogram
  return, (*self.spectrogram)
end


function stx_ivs_column::getbarycenter, idx
      
      ;the barycentre of the time_axis 
      ;split_t = self.start_time + (self.end_time - self.start_time)/2
    
      ;finding the barycentre: slow
      ;split_time = self.start_time
      ;while split_time lt self.end_time && total((*self.spectrogram).data[self.start_time:split_time,idx]) lt total((*self.spectrogram).data[split_time+1:self.end_time,idx]) do split_time++  
    
      
      ;finding the barycentre: quicker
      split_time = self.start_time
      collumn_mass = n_elements(idx) gt 1 ? total((*self.spectrogram).data[idx,self.start_time:self.end_time],1) : (*self.spectrogram).data[idx,self.start_time:self.end_time] 
      total_mass = total(collumn_mass)
      left_mass = 0
      step=-1
      
      while split_time le self.end_time && left_mass lt total_mass/2 do begin
        step++
        split_time++
        left_mass += collumn_mass[step]
      end
      
      split_time--
      return, split_time
end


;---------------------------------------------------------------------------
; Document name: stx_ivs_count_lookup_img.pro
; Created by:    Nicky Hochmuth, 2012/03/05
;---------------------------------------------------------------------------
;+
; PROJECT:
;       STIX
;
; PURPOSE:
;       defines a count for a good image
;
; CATEGORY:
;       Stix on Bord Algorithm
;
; CALLING SEQUENCE:
;       class = stx_ivs_count_lookup_img(counts.energy,duration)
;
; HISTORY:
;       2012/03/05, Nicky.Hochmuth@fhnw.ch, initial release
;
;-
;+
; :description:
;     a count amount classifier for the quality of producing significat images
;     
;
; :params:
;    count  amount
;    
; :returns:
;   0: if not enough counts for a good image
;   1: if enough counts for a good image
;   2: if enough counts for two or more good images
;-
function stx_ivs_column::count_lookup, counts, energy, time_span, get_thermalboundary=get_thermalboundary, get_min_count_t=get_min_count_t, get_min_count_nt=get_min_count_nt 
  
  default, morethanone, 2
  if keyword_set(get_thermalboundary) then return, self.thermalboundary
  if keyword_set(get_min_count_t) then return, self.min_count[0,0] 
  if keyword_set(get_min_count_nt) then return, self.min_count[0,1]
  
  thermal = energy ge self.thermalboundary
  
  if (counts gt self.min_count[thermal,1]) && (time_span ge (self.min_time[thermal]*morethanone)) then return, 2
  if counts gt self.min_count[thermal,0] then return , 1
  
  return, 0

end

;+
; :description:
;    Constructor
;
; :hidden:
;-
pro stx_ivs_column__define
   compile_opt idl2, hidden
   void = { stx_ivs_column, $
            start_time        : 0,$
            end_time          : 0,$
            split_level       : 0,$
            thermalboundary    : 25.0, $
            min_time          : [10.0,0], $
            min_count         : [[1000l,800],[2000,1600]], $
            rows              : ptr_new(),$
            ni_values         : ptr_new(),$
            spectrogram       : ptr_new()}
end