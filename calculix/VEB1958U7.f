      subroutine umat(amat,iel,iint,kode,elconloc,emec,emec0,
     &        beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,pnewdt,ipkon)

      implicit none
!
      character*80 amat
!
      integer ithermal,icmd,kode,ielas,iel,iint,nstate_,mi(*),iorien,
     &  ipkon(*)
!
      real*8 elconloc(21),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &  time,ttime,pnewdt
!
      real*8 xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*)
!
      real*8 zero, one, two, three, four
      
      integer i, j, k, l, k1
      
      parameter(zero=0.D0, one=1.D0, two=2.D0, three=3.D0, four=4.D0)
!
      
      real*8  cgreen(6),powercg(6),invcg(6),invc(3,3), 
     1 delta(3,3),c(3,3,3,3),cvol(3,3,3,3),ciso(3,3,3,3),
     2 S0isoold(6),H1old(6),H2old(6),H3old(6), 
     3 S0isonew(6),H1new(6),H2new(6),H3new(6),
     4 S0volnew(6),S0new(6),
     5 invc2(3,3),DEVinvc2(3,3),DEVI(3,3),cg(3,3)    
      
      real*8 I1,trC2,I2,dudj,dduddj,dwdi1,ddwddi1,dwdi2,ddwddi2, 
     1 ddwdi1di2,term,term1,term2,term3,term4,term5,term6,term7,term8,
     2 G1,G2,G3,tau1,tau2,tau3,ve,suma,TRinvc2,
     3 barI1,barI2,gamma0,gamma1,gammaneg1,
     4 C10,C20,C30,C01,kappa           
! ------------------------------------------------------------------           
!               CALCULIX UMAT FOR NONLINEAR VISCOELASTIC MODEL
!                              BIDERMAN (1958)
!
!                     WARSAW UNIVERSITY OF TECHNOLOGY
!                       CYPRIAN SUCHOCKI, JULY 2025         
! ------------------------------------------------------------------
!    elconloc(1) - C10
!    elconloc(2) - C20
!    elconloc(3) - C30
!    elconloc(4) - C01
!    elconcoc(5) - G1
!    elconcoc(6) - G2
!    elconcoc(7) - G3
!    elconcoc(8) - tau1
!    elconcoc(9) - tau2
!    elconcoc(10) - tau3
!    elconcoc(11) - kappa
! ------------------------------------------------------------------    
!                                       
!    Material  properties
!
      C10=elconloc(1)
      C20=elconloc(2)
      C30=elconloc(3)
      C01=elconloc(4)          
      G1=elconloc(5) 
      G2=elconloc(6) 
      G3=elconloc(7) 
      tau1=elconloc(8) 
      tau2=elconloc(9) 
      tau3=elconloc(10) 
      kappa=elconloc(11)
!
!    Extract variables from previous increment
!    
      do k1=1,6
         S0isoold(k1)=xstateini(k1,iint,iel)
         H1old(k1)=xstateini(k1+6,iint,iel)
         H2old(k1)=xstateini(k1+12,iint,iel)
         H3old(k1)=xstateini(k1+18,iint,iel)            
      end do                
!   
!    Calculate the right C-G deformation tensor C=2*E+1
!  
      cgreen(1)=two*emec(1)+one
      cgreen(2)=two*emec(2)+one
      cgreen(3)=two*emec(3)+one
      cgreen(4)=two*emec(4)
      cgreen(5)=two*emec(5)
      cgreen(6)=two*emec(6)
!      
! Warning! The following conditional clause can be very helpful
! in allowing some problems to converge to a solution. 
! However, in some cases non-physical solutions can be obtained
! (material overlapping etc., i.e. solutions that normally would not
! be obtained due to lack of convergence). Thus, one should either 
! leave or delete this clause based on his/her reasoning/experience.  !      
      if(vj.le.1.D-30) then 
!       Deformation is reset to zero in order to continue the 
!       calculation           
        cgreen(1)=one
        cgreen(2)=one
        cgreen(3)=one
        cgreen(4)=zero
        cgreen(5)=zero
        cgreen(6)=zero
        vj=one
      end if        
!   
!    Calculate the power of right C-G deformation tensor C*C
! 
      powercg(1)=cgreen(1)**two+cgreen(4)**two+cgreen(5)**two 
      powercg(2)=cgreen(4)**two+cgreen(2)**two+cgreen(6)**two 
      powercg(3)=cgreen(5)**two+cgreen(6)**two+cgreen(3)**two
      powercg(4)=cgreen(1)*cgreen(4)+cgreen(4)*cgreen(2)
     1          +cgreen(5)*cgreen(6)  
      powercg(5)=cgreen(1)*cgreen(5)+cgreen(4)*cgreen(6)
     1          +cgreen(5)*cgreen(3) 
      powercg(6)=cgreen(4)*cgreen(5)+cgreen(2)*cgreen(6)
     1          +cgreen(3)*cgreen(6)  
!    
!    Calculate the invariants of the right C-G tensor
!
      I1=cgreen(1)+cgreen(2)+cgreen(3)
      trC2=powercg(1)+powercg(2)+powercg(3)
      I2=(I1**two-trC2)/two  
!
!    Calculate the inverse of the right C-G tensor 
!    C-1=J-2(C2-I1C+I21)
!
      term=vj**(-two)

      invcg(1)=term*(powercg(1)-I1*cgreen(1)+I2)
      invcg(2)=term*(powercg(2)-I1*cgreen(2)+I2)
      invcg(3)=term*(powercg(3)-I1*cgreen(3)+I2) 
      invcg(4)=term*(powercg(4)-I1*cgreen(4))
      invcg(5)=term*(powercg(5)-I1*cgreen(5))
      invcg(6)=term*(powercg(6)-I1*cgreen(6))                        
!
!    Calculate the 2nd P-K stress
!
! ------------------------------------------------------------------ 
      dudj=kappa/(two*vj)*(vj**two-one)
      dduddj=kappa/two*(one+one/(vj**two))
      
      barI1=vj**(-two/three)*I1
      barI2=vj**(-four/three)*I2
      dwdi1=C10+C20*(barI1-three)+C30*(barI1-three)**two
      ddwddi1=C20+two*C30*(barI1-three)
      dwdi2=C01
      ddwddi2=zero
      ddwdi1di2=zero
      
      gamma0=two*vj**(-two/three)*(dwdi1+barI1*dwdi2)
      gamma1=-two*vj**(-four/three)*dwdi2
      gammaneg1=-two/three*(barI1*dwdi1+two*barI2*dwdi2)
      
      do k1=1,6
          S0volnew(k1)=vj*dudj*invcg(k1)
      end do        
      do k1=1,3
          S0isonew(k1)=gamma0+gamma1*cgreen(k1)+gammaneg1*invcg(k1)
      end do  
      do k1=4,6
          S0isonew(k1)=gamma1*cgreen(k1)+gammaneg1*invcg(k1)
      end do 
      do k1=1,6
         S0new(k1)=S0volnew(k1)+S0isonew(k1)
      end do
!
!    Update overstresses
!
      do k1=1,6
         H1new(k1)=exp(-dtime/tau1)*H1old(k1)
     1     +G1*(1-exp(-dtime/tau1))/(dtime/tau1)*(S0isonew(k1)
     2     -S0isoold(k1))
      end do
      do k1=1,6
         H2new(k1)=exp(-dtime/tau2)*H2old(k1)
     1     +G2*(1-exp(-dtime/tau2))/(dtime/tau2)*(S0isonew(k1)
     2     -S0isoold(k1))
      end do   
      do k1=1,6
         H3new(k1)=exp(-dtime/tau3)*H3old(k1)
     1     +G3*(1-exp(-dtime/tau3))/(dtime/tau3)*(S0isonew(k1)
     2     -S0isoold(k1))
      end do         
!
!    Update total 2nd PK stress
!      
      do k1=1,6 
         stre(k1)=S0new(k1)+H1new(k1)+H2new(k1)+H3new(k1)
      end do                                
!
!    Calculate the stiffness
!

!    Calculate Kronecker delta

      do i=1,3
         do j=1,3
           if(i.eq.j) then
               delta(i,j)=one   
           else 
               delta(i,j)=zero    
           endif
         end do
      end do 

!    Calculate the inverse of the right C-G tensor MATRIX

      invc(1,1)=invcg(1)
      invc(2,2)=invcg(2)
      invc(3,3)=invcg(3)
      invc(1,2)=invcg(4)
      invc(1,3)=invcg(5)
      invc(2,3)=invcg(6)
      invc(2,1)=invc(1,2)
      invc(3,1)=invc(1,3)
      invc(3,2)=invc(2,3)
      
!    Calculate the power of the inverse of the right C-G tensor MATRIX

      do i=1,3
         do j=i,3
            suma=zero
            do k=1,3
               suma=suma+invc(i,k)*invc(k,j)
            end do
            invc2(i,j)=suma
         end do
      end do
      invc2(2,1)=invc2(1,2)
      invc2(3,1)=invc2(1,3)
      invc2(3,2)=invc2(2,3)
      
!    Calculate the right C-G tensor MATRIX

      cg(1,1)=cgreen(1)
      cg(2,2)=cgreen(2)
      cg(3,3)=cgreen(3)
      cg(1,2)=cgreen(4)
      cg(1,3)=cgreen(5)
      cg(2,3)=cgreen(6)
      cg(2,1)=cg(1,2)
      cg(3,1)=cg(1,3)
      cg(3,2)=cg(2,3)      
      
!    Calculate the deviator of invc2 in the reference configuration
      
      TRinvc2=zero             
      do i=1,3
         do j=1,3
            TRinvc2=TRinvc2+cg(i,j)*invc2(i,j)
         end do
      end do
      term8=TRinvc2/three
      do i=1,3
         do j=i,3
            DEVinvc2(i,j)=invc2(i,j)-term8*invc(i,j)
         end do
      end do
      DEVinvc2(2,1)=DEVinvc2(1,2)
      DEVinvc2(3,1)=DEVinvc2(1,3)
      DEVinvc2(3,2)=DEVinvc2(2,3)
      
!    Calculate the deviator of the Kronecker delta in the reference 
!    configuration    
      term1=I1/three
      do i=1,3
         do j=i,3
            DEVI(i,j)=delta(i,j)-term1*invc(i,j)
         end do
      end do
      DEVI(2,1)=DEVI(1,2)
      DEVI(3,1)=DEVI(1,3)
      DEVI(3,2)=DEVI(2,3)

! Calculate viscoelastic factor

      ve=1+G1*(1-exp(-dtime/tau1))/(dtime/tau1)
     1   +G2*(1-exp(-dtime/tau2))/(dtime/tau2)
     2   +G3*(1-exp(-dtime/tau3))/(dtime/tau3)
      if(icmd.ne.3) then
      
      term2=dudj+vj*dduddj
      term3=four/three*vj**(-two/three)*dwdi1
      term4=four*vj**(-four/three)*dwdi2
      term5=four*vj**(-four/three)*ddwddi1
      term6=four*vj**(four/three)*ddwddi2
      term7=four*ddwdi1di2
      
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  cvol(i,j,k,l)=vj*(term2*invc(i,j)*invc(k,l)
     1             -dudj*(invc(i,k)*invc(j,l)
     2             +invc(i,l)*invc(j,k)))
! ------------------------------------------------------------------                    
                  ciso(i,j,k,l)=term3*(I1*((invc(i,k)*invc(j,l)
     -             +invc(i,l)*invc(j,k))/two
     -             -one/three*invc(i,j)*invc(k,l))
     -             -invc(i,j)*DEVI(k,l)-DEVI(i,j)*invc(k,l))
     -             +term4*(delta(i,j)*delta(k,l)
     -             -(delta(i,k)*delta(j,l)
     -             +delta(i,l)*delta(j,k))/two
     -             +two/three*I2*((invc(i,k)*invc(j,l)
     -             +invc(i,l)*invc(j,k))/two
     -             -two/three*invc(i,j)*invc(k,l))
     -             +two/three*vj**two*(invc(i,j)*DEVinvc2(k,l)
     -             +DEVinvc2(i,j)*invc(k,l)))
     -             +term5*DEVI(i,j)*DEVI(k,l)
     -             +term6*DEVinvc2(i,j)*DEVinvc2(k,l)
     -             -term7*(DEVinvc2(i,j)*DEVI(k,l)
     -             +DEVI(i,j)*DEVinvc2(k,l))
                  
                  c(i,j,k,l)=cvol(i,j,k,l)+ve*ciso(i,j,k,l)
               end do
            end do
         end do
      end do 
            
!     
!     stiff stores:
!     1111, 1122, 2222, 1133
!     2233, 3333, 1112, 2212
!
!     1111
      stiff(1)=c(1,1,1,1)
!     1122 = 2211
      stiff(2)=c(1,1,2,2)
!     2222
      stiff(3)=c(2,2,2,2)
!     1133 = 3311
      stiff(4)=c(1,1,3,3)
!     2233 = 3322
      stiff(5)=c(2,2,3,3)
!     3333
      stiff(6)=c(3,3,3,3)
!     1112 = 1121 = 1211 = 2111
      stiff(7)=c(1,1,1,2)
!     2212 = 1222 = 2122 = 2221
      stiff(8)=c(2,2,1,2)
!
!     stiff stores:
!     3312, 1212, 1113, 2213
!     3313, 1213, 1313, 1123
!     
!     3312 = 1233 = 2133 = 3321
      stiff(9)=c(3,3,1,2)
!     1212 = 1221 = 2112 = 2121
      stiff(10)=c(1,2,1,2)
!     1113 = 1131 = 1311 = 3111
      stiff(11)=c(1,1,1,3)
!     2213 = 1322 = 2231 = 3122
      stiff(12)=c(2,2,1,3)
!     
!     3313 = 1333 = 3133 = 3331
      stiff(13)=c(3,3,1,3)  
!     1213 = 1231 = 1312 = 1321 = 2113 = 2131 = 3112 = 3121:
      stiff(14)=c(1,2,1,3)
!     1313 = 1331 = 3113 = 3131
      stiff(15)=c(1,3,1,3)
!     1123 = 1132 = 2311 = 3211
      stiff(16)=c(1,1,2,3)
!
!     stiff stores:
!     2223, 3323, 1223, 1323
!     2323
!
!     2223 = 2232 = 2322 = 3222
      stiff(17)=c(2,2,2,3)
!     3323 = 2333 = 3233 = 3332
      stiff(18)=c(3,3,2,3)
!     1223 = 1232 = 2123 = 2132 = 2312 = 2321 = 3212 = 3221
      stiff(19)=c(1,2,2,3)
!     1323 = 1332 = 2313 = 2331 = 3123 = 3132 = 3213 = 3231
      stiff(20)=c(1,3,2,3)
!     2323 = 2332
      stiff(21)=c(2,3,2,3)
      
      endif 
!
!    Store stresses
!    
      do k1=1,6
         xstate(k1,iint,iel)=S0isonew(k1)
         xstate(k1+6,iint,iel)=H1new(k1)
         xstate(k1+12,iint,iel)=H2new(k1)
         xstate(k1+18,iint,iel)=H3new(k1) 
      end do        
!
      return
      end