subroutine relative_vorticity(m,n,o,u,v,deltax,deltay,output_array)
implicit none
integer :: m,n,o,mone,none
integer :: x,y,t
double precision ::deltax(n),deltay
double precision :: u(m,n,o)
double precision :: v(m,n,o)
double precision :: output_array(m,n,o)
double precision :: a
a= 0.5
mone=m-2
none=n-2
do t=1,o 
 do x=2,mone
  do y=2,none
   output_array(x,y,t) = (-1.00 *a * v(x,y-1,t)/deltax(y)) &
    + (a * v(x,y+1,t)/deltax(y)) &
    - (-1.00 * a * u(x-1,y,t)/deltay) &
    + (a * u(x+1,y,t)/deltay) 	   
  end do 
 end do 	
end do  
end subroutine relative_vorticity

