  module onebody_hadron_tensor_sf
    implicit none
    real*8, save :: e,ef,q2,p2,pf2,sina2,wf
  end module onebody_hadron_tensor_sf

! Computes the hadronic response tensor given form factors
! initial and final hadron momentum, energy transfer
! Output is complex 4x4 resp array
  subroutine compute_hadron_tensor(w, wt, xk_x, xk_y, xk_z, q, ff1v, ff2v, ffa, ffp, resp)
    use onebody_hadron_tensor_sf
    use onebody_currents_sf
    implicit none
    integer*4 :: i 
    real*8 :: w, wt, xk, xp, q
    real*8 :: ek,epf
    real*8 :: xk_x,xk_y,xk_z
    real*8 :: p_4(4),pp_4(4),qt_4(4)
    real*8 :: ff1v, ff2v, ffa, ffp
    complex*16 :: resp(4,4)

      ! Define initial and final nucleon 4 momentum
      ! Remember that q points along z
      xk = sqrt(xk_x**2 + xk_y**2 + xk_z**2)
      xp = sqrt(xk_x**2 + xk_y**2 + (xk_z + q)**2)

      ek = sqrt(xmn**2 + xk**2)
      epf = sqrt(xmn**2 + xp**2)

      p_4(1)=ek
      p_4(2)=xk_x
      p_4(3)=xk_y
      p_4(4)=xk_z

      pp_4(1) = epf
      pp_4(2) = xk_x
      pp_4(3) = xk_y
      pp_4(4) = xk_z + q

      qt_4 = pp_4 - p_4
      qt_4(1) = wt
  
      call current_init(p_4,pp_4,qt_4,w)
      call define_spinors()
      call sigccc(resp,ff1v,ff2v,ffa,ffp)
      return
    end subroutine compute_hadron_tensor
    
    subroutine sigccc(resp,ff1v,ff2v,ffa,ffp)
      use onebody_hadron_tensor_sf
      use onebody_currents_sf
      implicit none
      real*8 :: ff1v,ff2v,ffa,ffp
      complex*16 :: resp(4,4)
      complex*16 :: inter(4,4)
      integer*4 :: i,j 
      
      call det_Ja(ff1v,ff2v,ffa,ffp)
      call det_res1b(resp)

      !...Spin average gives a factor of 1/2 to the nuclear response tensor 
      resp = resp*0.5d0
      
      call shift(resp)

      inter = TRANSPOSE(resp)
      resp = inter

      return
    end subroutine sigccc

    !This module shifts the response tensor so that it is in (x,y,z,t), used by
    !TLorentzVector and GENIE convention instead of (t,x,y,z)
    subroutine shift(resp)
      use onebody_hadron_tensor_sf
      use onebody_currents_sf
      implicit none
      complex*16 :: resp(4,4)

      resp = CSHIFT(resp, 1, DIM=1)
      resp = CSHIFT(resp, 1, DIM=2)

      return
    end subroutine shift




  

