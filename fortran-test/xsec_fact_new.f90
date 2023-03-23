  module xsec_fact
    implicit none
    real*8, parameter :: xmp=0.938272046d0,xmn=0.93956563d0,xm=0.5d0*(xmn+xmp)
    real*8, parameter :: hbarc=0.197327053d0
    real*8, parameter :: pi= 4.0d0*atan(1.0d0),c0=1.0d0/16.0d0/pi/pi
    real*8, save :: e,ef,q2,p2,pf2,sina2,wf
  end module xsec_fact


! GENIE edit below, added in nucleon phi (nuphi) variable
  subroutine compute_hadron_tensor(wt, xk_x, xk_y, xk_z, xp_x, xp_y, xp_z, ff1v, ff2v, ffa, ffp, resp)
    use xsec_fact
    use dirac_matrices
    implicit none
    integer*4 :: i 
    real*8 :: wt, xk, xp
    real*8 :: ek,epf
    real*8 :: xk_x,xk_y,xk_z
    real*8 :: xp_x,xp_y,xp_z
    real*8 :: p_4(4),pp_4(4),q_4(4)
    real*8 :: ff1v, ff2v, ffa, ffp
    complex*16 :: resp(4,4)
!     w        energy transfer in MeV
!     wt       corrected energy trasnfer (accounts for binding energy) in MeV
!     xk, xp   initital and final nucleon three-momenta in inverse fm
      
      xk = sqrt(xk_x**2 + xk_y**2 + xk_z**2)
      xp = sqrt(xp_x**2 + xp_y**2 + xp_z**2)

      ek = sqrt(xk**2+xm**2)
      epf = sqrt(xm**2 + xp**2)
      
      !...define the initial and final nucleon 4-momentum
      p_4(1)=ek
      p_4(2)=xk_x
      p_4(3)=xk_y
      p_4(4)=xk_z
      pp_4(1)=epf
      pp_4(2)=xp_x
      pp_4(3)=xp_y
      pp_4(4)=xp_z
 
      q_4 = pp_4 - p_4
      q_4(1) = wt
      
      !write(6,*)"p = (", ek, ", ", p_4(2), ", ", p_4(3), ", ", p_4(4)
      !write(6,*)"pf = (", epf, ", ", pp_4(2), ", ", pp_4(3), ",  ", pp_4(4)
      !write(6,*)"q = (", q_4(1), ", ", q_4(2), ", ", q_4(3), ", ", q_4(4) 

      call current_init(p_4,pp_4,q_4)
      call define_spinors()
      call sigccc(resp,ff1v,ff2v,ffa,ffp)
      return
    end subroutine compute_hadron_tensor
    
    subroutine sigccc(resp,ff1v,ff2v,ffa,ffp)
      use xsec_fact
      use dirac_matrices
      implicit none
      real*8 :: ff1v,ff2v,ffa,ffp
      complex*16 :: resp(4,4)
      complex*16 :: inter(4,4)
      integer*4 :: i,j 
      
      call det_Ja(ff1v,ff2v,ffa,ffp)
      call det_res1b(resp)

      !...Spin average gives a factor of 1/2 to the nuclear response tensor 
      resp = resp*0.5d0
      !..Multiply by 1/mN^2 so that we are in same conventions 
      !resp = resp/xm/xm
      call shift(resp)

      inter = TRANSPOSE(resp)
      resp = inter

      return
    end subroutine sigccc

    !This module shifts the response tensor so that it is in (x,y,z,t), used by
    !TLorentzVector and GENIE convention instead of (t,x,y,z)
    subroutine shift(resp)
      use xsec_fact
      use dirac_matrices
      implicit none
      complex*16 :: resp(4,4)

      resp = CSHIFT(resp, 1, DIM=1)
      resp = CSHIFT(resp, 1, DIM=2)

      return
    end subroutine shift




  

