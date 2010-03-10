       Program Final
*--------------------------------------------------------------------*
*
*  Final Project
*
*  This program calculates and graphs the orbit of a hypothetical
*  asteroid in the asteroid belt, affected by the gravitational pull
*  of the sun and Jupiter, for two different representative orbital
*  periods of the asteroid (one a resonance of Jupiter's period
*  and one not).
*
*  Scott Bale, College of Arts and Science, Miami University
*
*                                 started  April 21, 1995  scb
*
*                            last revised    May  3, 1995  scb
*
*--------------------------------------------------------------------*
      Double Precision t0, y0(4)
      Double Precision G, m_sun, J_omega
      Double Precision cuberoot, ratio, pi
      integer i

      Parameter (pi = 3.14159 26535 89793)
      Parameter (m_sun = 1.99d30)
      Parameter (J_omega = 2.d0*pi/11.86d0)  !radians/yr
      Parameter (G = 6.6334d-5)     !(kg km/yr^2)*(km^2/kg^2)
      Parameter (cuberoot = 1.d0/3.d0)
*
* Get "ratio" from user
*
*      Write (*,*) 'Please enter the ratio of the orbital'
*      Write (*,*) 'period of the asteroid to that of Jupter.'
*      Write (*,*) '(Keep this ratio between 0.2 and 0.7).'
*      Read (*,*) ratio

      call gINIT
      call VIEWPORT (0.2d0,0.1d0,0.8d0,0.9d0)
      call WINDOW   (-8.d8,-8.d8,8.d8,8.d8)

      DO i=1,2
      call COLOR    (16)
      call LINE     (-8.d8,0.d0,8.d8,0.d0)
      call LINE     (0.d0,-8.d8,0.d0,8.d8)
*
* Initialize:  y0(1) is the initial radius out along the x-axis (km).
* y0(4) is the initial velocity directed in the positive
* y-direction (km/year).  Initially, the y-component of position is
* zero, and the x-component of velocity is zero.
*
      IF (i .eq. 1) ratio = 0.5d0
      IF (i .eq. 2) ratio = 11.d0/30.d0

      y0(4) = ( J_omega * G * m_sun / ratio ) **cuberoot
      y0(1) = G*m_sun/( y0(4)**2.d0 )
      y0(2) = 0.d0
      y0(3) = 0.d0
      t0    = 0.d0

      call R_K_F (t0, y0, ratio)

      call CLEAR
      END DO

      call gEND
      END


*     ________________
      Subroutine R_K_F ( t0, y0, ratio)
*--------------------------------------------------------------------
*
* Time, t, should be in years.
*

      Double Precision epsilon, maxtime, pi, ratio
      Double Precision h, t, y(4), t0, y0(4)
      Double Precision yold(4), yhat(4)
      Double Precision Err, BigErr, MaxErr, hnew, hmin, hmax

      Double Precision f0(4), f1(4), f2(4), f3(4), f4(4), f5(4)

      Double Precision a1, a2, a3, a4, a5
      Double Precision b10, b20,b21, b30,b31,b32, b40,b41,b42,b43,
     +                 b50,b51,b52,b53,b54
      Double Precision c0,c2,c3,c4,c5
      Double Precision d0,d2,d3,d4,d5
      Double Precision J_r, J_omega, J_omega0
      Double Precision J_x, J_xold, J_xnew, J_y, J_yold, J_ynew
      integer i

      Parameter (epsilon = 5.d-5, maxtime = 100.d0)
      Parameter (pi = 3.14159 26535 89793)
      Parameter (a1=0.25d0,a2=0.375d0,a3=1.2d1/1.3d1,a4=1.d0,
     +           a5=0.5d0)

      Parameter (b10=0.25d0, b20=3.d0/3.2d1,b21=9.d0/3.2d1)
      Parameter (b30=1.932d3/2.197d3,b31=-7.2d3/2.197d3,
     +           b32=7.296d3/2.197d3)
      Parameter (b40=4.39d2/2.16d2,b41=-8.d0,b42=3.68d3/5.13d2,
     +           b43=-8.45d2/4.104d3)
      Parameter (b50=-8.d0/2.7d1,b51=2.d0,b52=-3.544d3/2.565d3,
     +           b53=1.859d3/4.104d3,b54=-1.1d1/4.d1)

      Parameter (c0=1.6d1/1.35d2, c2=6.656d3/1.2825d4)
      Parameter (c3=2.8561d4/5.6430d4,c4=-9.d0/5.d1,c5=2.d0/5.5d1)

      Parameter (d0=1.d0/3.6d2,d2=-1.28d2/4.275d3)
      Parameter (d3=-2.197d3/7.524d4,d4=1.d0/5.d1,d5=2.d0/5.5d1)

      Parameter (J_omega = 2.d0*pi/11.86) !radians/year
      Parameter (J_r = 7.783d8)       !km
      Parameter (J_omega0 = 0.d0)     !radians
*
* Initialize "old" variables, which will be used in graphing
*

      J_xold = J_r*cos(J_omega*t0+J_omega0)
      J_yold = J_r*sin(J_omega*t0+J_omega0)
      J_x = J_xold
      J_y = J_yold
      DO i=1,4
         yold(i) = y0(i)
      END DO
*
* Initialize "h", the time increment (in years)
*
      h = 0.5d0

100   call DERS (t0, y0, f0)

200   t = t0 + a1*h
      DO i=1,4
          y(i) = y0(i) + b10*h*f0(i)
      END DO
      call DERS (t,y,f1)

      t = t0 + a2*h
      DO i=1,4
          y(i) = y0(i) + b20*h*f0(i) + b21*h*f1(i)
      END DO
      call DERS (t,y,f2)

      t = t0 + a3*h
      DO i=1,4
          y(i) = y0(i) + b30*h*f0(i) + b31*h*f1(i) + b32*h*f2(i)
      END DO
      call DERS (t,y,f3)

      t = t0 + a4*h
      DO i=1,4
          y(i) = y0(i) + b40*h*f0(i) + b41*h*f1(i) + b42*h*f2(i)
     +                 + b43*h*f3(i)
      END DO
      call DERS (t,y,f4)

      t = t0 + a5*h
      DO i=1,4
          y(i) = y0(i) + b50*h*f0(i) + b51*h*f1(i) + b52*h*f2(i)
     +                 + b53*h*f3(i) + b54*h*f4(i)
      END DO
      call DERS (t,y,f5)
*
* Calculate 5th order approximation (yhat) etc.
*
      BigErr = 0.d0
      DO i=1,4
         Yhat(i) = y0(i)+h*(c0*f0(i) + c2*f2(i) + c3*f3(i)
     +                     +c4*f4(i) + c5*f5(i) )
         ERR     = h * DABS(d0*f0(i)+ d2*f2(i) + d3*f3(i)
     +                     +d4*f4(i)+ d5*f5(i) )
         IF(Err .gt. Big Err)BigErr = Err
      END DO
      MaxErr = h*epsilon
      hnew = 0.90d0 * h * sqrt(sqrt(MaxErr/BigErr) )
*
* Check for increase/decrease in h, as well as h-limits
*
      hmax = 1.d0    !max yrs iterated in one step
      hmin = 0.00001 !min yrs iterated in one step
      IF(hnew .gt. 4.d0*h)hnew = 4.d0*h
      IF(hnew .lt. 0.1d0*h)hnew = 0.1d0*h
      IF(hnew .gt. hmax)hnew = hmax
      IF(hnew .lt. hmin)THEN
         write (*,*) 'h is TOO SMALL'
         write (*,*) t0,y0,yhat,hnew
         stop 'possible problem in rkf integrator'
      ENDIF
*
* If the error is too large, repeat the propogation step using
* hnew
*
      IF(BigErr .gt. MaxErr) THEN
         h=hnew
         goto 200
      ENDIF
*
* Graph Jupiter and the asteroid on the cartesian plane.
* Asteroid:
*
      call COLOR (9)      !erase old line
      call LINE (yold(1), yold(3), y0(1), y0(3) )
      call COLOR (11)     !draw in new position
      call LINE (y0(1), y0(3), yhat(1), yhat(3) )
*
* Jupiter:
*
      !J_x    = J_r*cos(J_omega*t0 + J_omega0)
      !J_y    = J_r*sin(J_omega*t0 + J_omega0)
      J_xnew = J_r*cos(J_omega * (t0+h) + J_omega0)
      J_ynew = J_r*sin(J_omega * (t0+h) + J_omega0)
      call COLOR (4)             !erase old line
      call Line  (J_xold, J_yold, J_x, J_y)
      call COLOR (10)            !draw new position
      call LINE  (J_x, J_y, J_xnew, J_ynew)
*
* Shuffle the graphing variables for the next loop
*
      J_xold = J_x
      J_yold = J_y
      J_x = J_xnew
      J_y = J_ynew
*
* Get the computer's attention, if we need to stop the program
*
      call CURSOR (26,0)
      write (*,500) t0, ratio
500   format (1x,'time=',1x,f10.2,1x,'years',/,1x,'ratio=',1x,f5.3)
*
* The error is small enough, so this step is acceptable.
* Redefine x0 to move the solution along; let the more
* accurate approximation become the initial value for
* the next step.
*
      t0 = t0 + h
      DO i=1,4
          yold(i) = y0(i)
          y0(i) = yhat(i)
      END DO
      h = hnew
*
* Delay
*
*      do i = 1,30
*      end do

*
* Have we gone far enough?  We want to project the solution some
* number of years into the future.  When t0 + h has reached some
* predetermined number of years, then the array yhat(i) contains
* all the information about the asteroid's position and velocity
* that we desire at that time, and the subroutine can end.
*
      IF (t0 + h .lt. maxtime) goto 100



      END


*     _______________
      Subroutine DERS  (t, y, f)
*--------------------------------------------------------------------
      double precision t, y(4), f(4), J_x, J_y
      double precision m_sun, m_Jup, G
      double precision J_omega, J_omega0, J_r
      double precision pi
      Parameter (pi = 3.14159 26535 89793)
      Parameter (m_sun = 1.99d30, m_Jup = 1.90d27)
      Parameter (G = 6.6334d-5)        !(kg km/yr^2)*(km^2/kg^2)
      Parameter (J_omega0 = 0.d0)      !radians
      Parameter (J_omega = 2.d0*pi/11.86d0)  !radians/yr
      Parameter (J_r = 7.783d8)        !km

      J_x = J_r*cos(J_omega*t+J_omega0)
      J_y = J_r*sin(J_omega*t+J_omega0)

      f(1) = y(2)
      f(2) = -(G*m_sun/((y(1)**2.d0+y(3)**2.d0)**1.5d0))*y(1)
     +       +(G*m_Jup/(((J_x-y(1))**2.d0+(J_y-y(3))**2.d0)**1.5d0))
     +       *(J_x-y(1))
      f(3) = y(4)
      f(4) = -(G*m_sun/((y(1)**2.d0+y(3)**2.d0)**1.5d0))*y(3)
     +       +(G*m_Jup/(((J_x-y(1))**2.d0+(J_y-y(3))**2.d0)**1.5d0))
     +       *(J_y-y(3))

      END




