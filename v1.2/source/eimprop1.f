c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eimprop1  --  impr. dihedral energy & gradient  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eimprop1" calculates improper dihedral energy and its
c     first derivatives with respect to Cartesian coordinates
c
c
      subroutine eimprop1
      use atmlst
      use atmtyp
      use atoms
      use bound
      use deriv
      use domdec
      use energi
      use group
      use improp
      use inform
      use iounit
      use math
      use torpot
      use usage
      use virial
      implicit none
      integer i,ia,ib,ic,id
      integer ialoc,ibloc,icloc,idloc
      integer iimprop
      real*8 e,dedphi
      real*8 dt
      real*8 ideal,force
      real*8 cosine,sine
      real*8 rcb,angle
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 rt2,ru2,rtru
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 dedxt,dedyt,dedzt
      real*8 dedxu,dedyu,dedzu
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 fgrp
      logical proceed
c
      if (deb_Path) write(iout,*), 'eimprop1 '
c
c
c     zero out energy and first derivative components
c
      eid = 0.0d0
c
c     calculate the improper dihedral angle energy term
c
      do iimprop = 1, niproploc
         i = impropglob(iimprop)
         ia = iiprop(1,i)
         ib = iiprop(2,i)
         ic = iiprop(3,i)
         id = iiprop(4,i)
         ialoc = loc(ia)
         ibloc = loc(ib)
         icloc = loc(ic)
         idloc = loc(id)
c
c     decide whether to compute the current interaction
c
         if (use_group)  call groups (fgrp,ia,ib,ic,id,0,0)
         proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the improper dihedral angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            if (use_polymer) then
               call image (xba,yba,zba)
               call image (xcb,ycb,zcb)
               call image (xdc,ydc,zdc)
            end if
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
c
c     set the improper dihedral parameters for this angle
c
               ideal = vprop(i)
               force = kprop(i)
               if (abs(angle+ideal) .lt. abs(angle-ideal))
     &            ideal = -ideal
               dt = angle - ideal
               do while (dt .gt. 180.0d0)
                  dt = dt - 360.0d0
               end do
               do while (dt .lt. -180.0d0)
                  dt = dt + 360.0d0
               end do
c
c     calculate improper energy and master chain rule term
c
               e = idihunit * force * dt**2
               dedphi = 2.0d0 * radian * idihunit * force * dt
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  dedphi = dedphi * fgrp
               end if
c
c     chain rule terms for first derivative components
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               if (use_polymer) then
                  call image (xca,yca,zca)
                  call image (xdb,ydb,zdb)
               end if
               dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
               dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
               dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
               dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
               dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
               dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c     compute first derivative components for this angle
c
               dedxia = zcb*dedyt - ycb*dedzt
               dedyia = xcb*dedzt - zcb*dedxt
               dedzia = ycb*dedxt - xcb*dedyt
               dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
               dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
               dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
               dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
               dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
               dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
               dedxid = zcb*dedyu - ycb*dedzu
               dedyid = xcb*dedzu - zcb*dedxu
               dedzid = ycb*dedxu - xcb*dedyu
c
c     calculate improper dihedral energy and derivatives
c
               eid = eid + e
               deid(1,ialoc) = deid(1,ialoc) + dedxia
               deid(2,ialoc) = deid(2,ialoc) + dedyia
               deid(3,ialoc) = deid(3,ialoc) + dedzia
c
               deid(1,ibloc) = deid(1,ibloc) + dedxib
               deid(2,ibloc) = deid(2,ibloc) + dedyib
               deid(3,ibloc) = deid(3,ibloc) + dedzib
c
               deid(1,icloc) = deid(1,icloc) + dedxic
               deid(2,icloc) = deid(2,icloc) + dedyic
               deid(3,icloc) = deid(3,icloc) + dedzic
c
               deid(1,idloc) = deid(1,idloc) + dedxid
               deid(2,idloc) = deid(2,idloc) + dedyid
               deid(3,idloc) = deid(3,idloc) + dedzid
c
c     increment the internal virial tensor components
c
               vxx = xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid
               vyx = ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid
               vzx = zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid
               vyy = ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid
               vzy = zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid
               vzz = zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
            end if
         end if
      end do
      return
      end
