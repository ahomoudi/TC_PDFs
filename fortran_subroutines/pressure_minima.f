subroutine minimum_press(lon,lat,time,a,b,ps,box,filter)
	implicit none
	integer, intent(in) :: lon,lat,time
	integer :: x,y,t,mone,none,a,b
	double precision, intent(inout) :: ps(lon,lat,time),box(a,b)
	double precision, intent(out) :: filter(lon,lat,time)

	mone=lon - 3
	none=lat - 3
	do 20,t=1,time 
	  do 10, x=4,mone
	   do 5, y=4,none
		box = reshape((/ps(x-3,y-3,t)  ,ps(x-3,y-2,t),ps(x-3,y-1,t) &
		                , ps(x-3,y,t)  ,ps(x-3,y+1,t),ps(x-3,y+2,t) &
				, ps(x-3,y+3,t) & 
                                , ps(x-2,y-3,t),ps(x-2,y-2,t),ps(x-2,y-1,t) &
		                , ps(x-2,y,t)  ,ps(x-2,y+1,t),ps(x-2,y+2,t) &
				, ps(x-2,y+3,t) & 
                                , ps(x-1,y-3,t),ps(x-1,y-2,t),ps(x-1,y-1,t) &
		                , ps(x-1,y,t)  ,ps(x-1,y+1,t),ps(x-1,y+2,t) &
				, ps(x-1,y+3,t) & 
			        , ps(x,y-3,t),ps(x,y-2,t),ps(x,y-1,t) &
		                , ps(x,y,t)  ,ps(x,y+1,t),ps(x,y+2,t) &
				, ps(x,y+3,t) & 
                                , ps(x+1,y-3,t),ps(x+1,y-2,t),ps(x+1,y-1,t) &
		                , ps(x+1,y,t)  ,ps(x+1,y+1,t),ps(x+1,y+2,t) &
				, ps(x+1,y+3,t) & 
                                , ps(x+2,y-3,t),ps(x+2,y-2,t),ps(x+2,y-1,t) &
		                , ps(x+2,y,t)  ,ps(x+2,y+1,t),ps(x+2,y+2,t) &
				, ps(x+2,y+3,t) & 
                                , ps(x+3,y-3,t),ps(x+3,y-2,t),ps(x+3,y-1,t) &
		                , ps(x+3,y,t)  ,ps(x+3,y+1,t),ps(x+3,y+2,t) &
				, ps(x+3,y+3,t)/), (/a,b/))


		   if(minval(box) .EQ. box(4,4)) then                       
                             filter(x,y,t) = 1.00
		      else
				filter(x,y,t) = 0.00
		     end if 


5          end do 
10         end do 	
20	end do  
	end subroutine minimum_press
