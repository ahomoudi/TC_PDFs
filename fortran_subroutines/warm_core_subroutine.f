subroutine warm_core(lon,lat,time,temp300,temp500,&
		 temp700,temp850,filter)
	implicit none
	integer, intent(in) :: lon,lat,time
	integer :: x,y,t,mone,none
	double precision, intent(in) :: temp300(lon,lat,time),temp500(lon,lat,time)
	double precision, intent(in) :: temp700(lon,lat,time),temp850(lon,lat,time)
	double precision, intent(out) :: filter(lon,lat,time)

	mone=lon - 3
	none=lat - 3
	do 20,t=1,time 
	  do 10, x=4,mone
	   do 5, y=4,none
		   if(temp300(x,y,t) .GT. 0.00) then
                      if(temp500(x,y,t) .GT. 0.00) then
                        if(temp700(x,y,t) .GT. 0.00) then
                          if(temp850(x,y,t) .GT. 0.00) then
                            if(temp300(x,y,t) .GT. temp850(x,y,t)) then
                              filter(x,y,t) = 1.00
				else
				filter(x,y,t) = 0.00
			     end if 
			   end if
                         end if
                       end if 
		     end if 


5          end do 
10         end do 	
20	end do  
	end subroutine warm_core
