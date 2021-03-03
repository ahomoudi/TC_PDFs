subroutine relative_vorticity(lon,lat,time,u_array,v_array,deltax,deltay,output_array)
implicit none
integer :: lon,lat,time,mone,none
integer :: x,y,t
double precision ::deltax(lat),deltay(lon)
double precision :: u_array(lon,lat,time)
double precision :: v_array(lon,lat,time)
double precision :: output_array(lon,lat,time)
double precision :: a,b,d
a= 0.5
b= 0.15
d= 0.75
mone=lon-1
none=lat-1
do t=1,time 
 do x=2,mone
  do y=2,none
   output_array(x,y,t) = -(a * v_array(x-1,y,t)/deltax(y)) &
    + (a * v_array(x-1,y,t)/deltax(y)) &
    - (a * u_array(x,y-1,t)/deltay(x)) &
    + (a * u_array(x,y+1,t)/deltay(x)) 	   
  end do 
 end do 	
end do  
end subroutine relative_vorticity

