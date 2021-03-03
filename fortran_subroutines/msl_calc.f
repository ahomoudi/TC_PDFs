subroutine p850_to_msl(lon,lat,time,zg_array,ta_array,output_array)
	implicit none
	integer, intent(in) :: lon,lat,time
	integer :: x,y,t
	double precision, intent(in) :: zg_array(lon,lat,time),ta_array(lon,lat,time)
	double precision, intent(out) :: output_array(lon,lat,time)
	do 20,t=1,time 
	  do 10, x=1,lon
	   do 5, y=1,lat
	   output_array(x,y,t) = 85000*(1 - ((0.0065* zg_array(x,y,t))&
				/(ta_array(x,y,t)+(0.0065* zg_array(x,y,t)))))** (-5.257)

5          end do 
10         end do 	
20	end do  
	end subroutine p850_to_msl
