subroutine integrated_temperature(input_array,lon,lat,time,output_array)
	implicit none
	integer, intent(in) :: lon,lat,time
	integer :: x,y,t,mone,none
	double precision, intent(in) :: input_array(lon,lat,time)
	double precision, intent(out) :: output_array(lon,lat,time)
        double precision ::a
	a =0.02040816
	mone=lon - 3
	none=lat - 3
	do 20,t=1,time 
	  do 10, x=4,mone
	   do 5, y=4,none
	   output_array(x,y,t) = ((input_array(x-3,y-3,t)+input_array(x-3,y-2,t)+input_array(x-3,y-1,t) &
		                + input_array(x-3,y,t)  +input_array(x-3,y+1,t)+input_array(x-3,y+2,t) &
				+ input_array(x-3,y+3,t) & 
                                + input_array(x-2,y-3,t)+input_array(x-2,y-2,t)+input_array(x-2,y-1,t) &
		                + input_array(x-2,y,t)  +input_array(x-2,y+1,t)+input_array(x-2,y+2,t) &
				+ input_array(x-2,y+3,t) & 
                                + input_array(x-1,y-3,t)+input_array(x-1,y-2,t)+input_array(x-1,y-1,t) &
		                + input_array(x-1,y,t)  +input_array(x-1,y+1,t)+input_array(x-1,y+2,t) &
				+ input_array(x-1,y+3,t) & 
			        + input_array(x,y-3,t)+input_array(x,y-2,t)+input_array(x,y-1,t) &
		                + input_array(x,y,t)  +input_array(x,y+1,t)+input_array(x,y+2,t) &
				+ input_array(x,y+3,t) & 
                                + input_array(x+1,y-3,t)+input_array(x+1,y-2,t)+input_array(x+1,y-1,t) &
		                + input_array(x+1,y,t)  +input_array(x+1,y+1,t)+input_array(x+1,y+2,t) &
				+ input_array(x+1,y+3,t) & 
                                + input_array(x+2,y-3,t)+input_array(x+2,y-2,t)+input_array(x+2,y-1,t) &
		                + input_array(x+2,y,t)  +input_array(x+2,y+1,t)+input_array(x+2,y+2,t) &
				+ input_array(x+2,y+3,t) & 
                                + input_array(x+3,y-3,t)+input_array(x+3,y-2,t)+input_array(x+3,y-1,t) &
		                + input_array(x+3,y,t)  +input_array(x+3,y+1,t)+input_array(x+3,y+2,t) &
				+ input_array(x+3,y+3,t))/49) 

5          end do 
10         end do 	
20	end do  
	end subroutine integrated_temperature
