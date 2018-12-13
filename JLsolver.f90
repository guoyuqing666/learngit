 program main

    integer i,j,k
    integer status
    integer n,iter,iter0

    character(len=5) :: idtraj

    double precision, dimension(:), pointer:: u,v,v0,a,b,c,d,cp,dp
    double precision, dimension(:), pointer:: T,rho,T0
    double precision, dimension(:,:), pointer:: species,species0,source
    double precision gamma,gamma0,mu,mfluxl,mfluxr,gamma_t
    double precision dz,error,erroru,errory,vlinlet,Tinlet,factor,factor0,dt,hcp
    integer imax,kmax,counter

    common/small/ n
    common/large0/ source,species
    common/large1/ T,rho

    factor0=1.d-3
    dz=0.0001
    n=20000
    gamma=-0.1
    gamma0=gamma
    allocate(u(n),STAT=status)
    allocate(v(n),STAT=status)
    allocate(v0(n),STAT=status)
    allocate(a(n),STAT=status)
    allocate(b(n),STAT=status)
    allocate(c(n),STAT=status)
    allocate(d(n),STAT=status)
    allocate(cp(n),STAT=status)
    allocate(dp(n),STAT=status)
    allocate(T(n),STAT=status)
    allocate(T0(n),STAT=status)
    allocate(rho(n),STAT=status)

    allocate(source(n,7),STAT=status)
    allocate(species(n,6),STAT=status)
    allocate(species0(n,6),STAT=status)

    Tinlet=1000.d0
    hcp=1.2
    !hcp is the heat capacity
    do i=1,n
        rho(i)=300/Tinlet
        T(i)=Tinlet
        T0(i)=T(i)
    enddo
    do i=1,n
        !      species(i,1)=1-(i-1.d0)/(n-1)
        if(i.le.n/2) then
            species(i,1)=1.d0
        else
            species(i,1)=0.d0
        endif
        species(i,2)=1-species(i,1)
        species(i,3)=1.d-4
        species(i,4)=1.d-4
        species(i,5)=0.d0
        species(i,6)=0.d0
    enddo


    !    open(60,file='./u.txt',action='read')
    !    do i=1,n
    !    read(60,'(17(1x,es14.4))') temp,temp,temp,T(i),species(i,1),species(i,2),species(i,3),species(i,4),species(i,5),species(i,6),temp,temp,temp,temp,temp,temp,temp
    !    enddo
    !     close(60)


    do i=1,n
        T0(i)=T(i)
        do j=1,6
            species0(i,j)=species(i,j)
        enddo
    enddo

    mu=5.d-4
    mfluxl=0.1d0
    mfluxr=0.12d0

    do i=1,n
        u(i)=-2.d0*(i-n/2.d0)/n/rho(i)
        !      v(i)=-(-mfluxr-mfluxl)/(2*rho(i)*dz*n)*(1+i*1.d0/n)/2.d0
        v(i)=-(-mfluxr-mfluxl)/(2*rho(i)*dz*n)
        v0(i)=v(i)
    enddo

    !      temp=0.d0
    !      temp=temp-2*rho*v(1)*(dz/2)
    !      do i=2,n-1
    !      temp=temp-2*rho*v(i)*(dz)
    !      enddo
    !      temp=temp-2*rho*v(n)*(dz/2)
    !      do i=1,n
    !      v(i)=v(i)*(-mfluxr-mfluxl)/temp
    !      v0(i)=v(i)
    !      enddo

    u(1)=mfluxl/rho(1)
    u(n)=-mfluxr/rho(n)

    iter=0

    error=1.d0
    do while(error.gt.-1.d-20)
        error=0.d0
        iter=iter+1

        erroru=1.d0
        iter0=0
        do while(erroru.gt.1.d-5)
            iter0=iter0+1
            erroru=0.d0
            vinlet=v0(2)
            a(1)=0.d0
            b(1)=1.d0
            c(1)=0.d0
            d(1)=vinlet
            !      a(1)=0.d0
            !      b(1)=1.d0
            !      c(1)=-1.d0
            !      d(1)=0.d0
            do i=2,n-1
                if(u(i).ge.0.d0) then
                    a(i)=-rho(i)*u(i)/dz-mu/dz/dz
                    b(i)=rho(i)*u(i)/dz+rho(i)*v0(i)+2*mu/dz/dz
                    c(i)=-mu/dz/dz
                    d(i)=-gamma0
                else
                    a(i)=-mu/dz/dz
                    b(i)=-rho(i)*u(i)/dz+rho(i)*v0(i)+2*mu/dz/dz
                    c(i)=rho(i)*u(i)/dz-mu/dz/dz
                    d(i)=-gamma0
                endif
                a(n)=-1.d0
                b(n)=1.d0
                c(n)=0.d0
                d(n)=0.d0
            enddo
            cp(1)=c(1)/b(1)
            dp(1)=d(1)/b(1)
            do i=2,n
                cp(i)=c(i)/(b(i)-a(i)*cp(i-1))
                dp(i)=(d(i)-a(i)*dp(i-1))/(b(i)-a(i)*cp(i-1))
            enddo
            v(n)=dp(n)
            do i=n-1,1,-1
                v(i)=dp(i)-v(i+1)*cp(i)
            enddo

            temp=0.d0
            temp=temp-2*rho(1)*v(1)*(dz/2)
            do i=2,n-1
                temp=temp-2*rho(i)*v(i)*(dz)
            enddo
            temp=temp-2*rho(n)*v(n)*(dz/2)
            do i=1,n
                v(i)=v(i)*(-mfluxr-mfluxl)/temp
            enddo

            do i=2,n-1
                u(i)=(-2*rho(i)*v(i)*dz+rho(i-1)*u(i-1))/rho(i)
            enddo

            gamma_t=0.d0
            do i=1,n
                if(i.eq.1) then
                    gamma_t=gamma_t-rho(i)*u(i)*(v(i+1)-v(i))/dz-rho(i)*v(i)**2+mu*(v(i+2)+v(i)-2*v(i+1))/dz**2
                elseif(i.eq.n) then
                    gamma_t=gamma_t-rho(i)*u(i)*(v(i)-v(i-1))/dz-rho(i)*v(i)**2+mu*(v(i)+v(i-2)-2*v(i-1))/dz**2
                else
                    gamma_t=gamma_t-rho(i)*u(i)*(v(i+1)-v(i-1))/2/dz-rho(i)*v(i)**2+mu*(v(i+1)+v(i-1)-2*v(i))/dz**2
                endif
            enddo
            gamma=gamma_t/n

            erroru=dabs(gamma-gamma0)
            if(iter0.eq.1) error=erroru
            gamma=(gamma+gamma0)/2.d0
            gamma0=gamma

            do i=1,n
                v0(i)=v(i)
            enddo
            !      write(*,*) 'erroru=',erroru
        enddo

        call chem()

        ! this part is for Y and T
        !      do k=1,6
        errory=10.d0
        counter=0
        do while(errory.ge.4.d0)
            dt=1.d-7
            errory=0.d0
            counter=counter+1

            ! this part is for T
            do i=1,n
                dt=min(dt,0.05/dabs(rho(i)*source(i,7)+1.d-20))
                dt=min(dt,1.d-4/dabs(source(i,1)+1.d-20))
                dt=min(dt,1.d-4/dabs(source(i,2)+1.d-20))
                dt=min(dt,1.d-5/dabs(source(i,3)+1.d-20))
                dt=min(dt,1.d-5/dabs(source(i,4)+1.d-20))
                dt=min(dt,1.d-5/dabs(source(i,5)+1.d-20))
                dt=min(dt,1.d-5/dabs(source(i,6)+1.d-20))
            enddo

            a(1)=0.d0
            b(1)=1.d0
            c(1)=0.d0
            d(1)=Tinlet
            do i=2,n-1
                if(u(i).ge.0.d0) then
                    a(i)=(-rho(i)*u(i)/dz-mu/dz/dz)*dt
                    b(i)=rho(i)+(rho(i)*u(i)/dz+2*mu/dz/dz)*dt
                    c(i)=(-mu/dz/dz)*dt
                    !      source(i)=1.d10*species0(i)*(1-species0(i))**4*exp(-1000/T(i))
                    !      d(i)=1.d10*species0(i)*(1-species0(i))**4*exp(-1000/T(i))
                    !      source(i)=7.d8*exp(-((i-n/2.d0)/500.d0)**2)
                    d(i)=rho(i)*T0(i)+source(i,7)/hcp*dt
                else
                    a(i)=(-mu/dz/dz)*dt
                    b(i)=rho(i)+(-rho(i)*u(i)/dz+2*mu/dz/dz)*dt
                    c(i)=(rho(i)*u(i)/dz-mu/dz/dz)*dt
                    !      source(i)=1.d10*species0(i)*(1-species0(i))**4*exp(-1000/T(i))
                    !      d(i)=1.d10*species0(i)*(1-species0(i))**4*exp(-1000/T(i))
                    !      source(i)=7.d8*exp(-((i-n/2.d0)/500.d0)**2)
                    d(i)=rho(i)*T0(i)+source(i,7)/hcp*dt
                endif
                a(n)=0.d0
                b(n)=1.d0
                c(n)=0.d0
                d(n)=Tinlet
            enddo
            cp(1)=c(1)/b(1)
            dp(1)=d(1)/b(1)
            do i=2,n
                cp(i)=c(i)/(b(i)-a(i)*cp(i-1))
                dp(i)=(d(i)-a(i)*dp(i-1))/(b(i)-a(i)*cp(i-1))
            enddo
            temp=0.d0
            T(n)=dp(n)
            do i=n-1,1,-1
                T(i)=dp(i)-T(i+1)*cp(i)
            enddo
            imax=0
            kmax=0
            do i=1,n
                if(errory.lt.dabs(T(i)-T0(i))/Tinlet/dt) then
                    errory=dabs(T(i)-T0(i))/Tinlet/dt
                    imax=i
                    kmax=1
                endif
            enddo

            do k=1,6
                a(1)=0.d0
                b(1)=1.d0
                c(1)=0.d0
                if(k.eq.1) then
                    d(1)=1.d0
                elseif(k.eq.2) then
                    d(1)=0.d0
                else
                    d(1)=0.d0
                endif
                do i=2,n-1
                    if(u(i).ge.0.d0) then
                        a(i)=(-rho(i)*u(i)/dz-mu/dz/dz)*dt
                        b(i)=rho(i)+(rho(i)*u(i)/dz+2*mu/dz/dz)*dt
                        c(i)=(-mu/dz/dz)*dt
                        d(i)=source(i,k)*dt+rho(i)*species0(i,k)
                    else
                        a(i)=(-mu/dz/dz)*dt
                        b(i)=rho(i)+(-rho(i)*u(i)/dz+2*mu/dz/dz)*dt
                        c(i)=(rho(i)*u(i)/dz-mu/dz/dz)*dt
                        d(i)=source(i,k)*dt+rho(i)*species0(i,k)
                    endif
                enddo
                a(n)=0.d0
                b(n)=1.d0
                c(n)=0.d0
                if(k.eq.1) then
                    d(n)=0.d0
                elseif(k.eq.2) then
                    d(n)=1.d0
                else
                    d(n)=0.d0
                endif

                cp(1)=c(1)/b(1)
                dp(1)=d(1)/b(1)
                do i=2,n
                    cp(i)=c(i)/(b(i)-a(i)*cp(i-1))
                    dp(i)=(d(i)-a(i)*dp(i-1))/(b(i)-a(i)*cp(i-1))
                enddo
                species(n,k)=dp(n)
                do i=n-1,1,-1
                    species(i,k)=(dp(i)-species(i+1,k)*cp(i))
                    species(i,k)=max(species(i,k),0.d0)
                    species(i,k)=min(species(i,k),1.d0)
                enddo

                do i=1,n
                    if(errory.lt.dabs(species(i,k)-species0(i,k))/dt) then
                        errory=dabs(species(i,k)-species0(i,k))/dt
                        imax=i
                        kmax=k
                    endif
                enddo

            enddo

            call chem()
            if(mod(counter,1000).eq.0) then
                write(*,'(es14.4,es14.4,es14.4,1x,i6,1x,i6,1x,i8)') errory,dt,species(imax,kmax),imax,kmax,counter
                open(60,file='./u.txt',action='write')
                !      write(60,*) '"coef 4,8,4", "coef 8,16,16", "jpdf1","jpdf2",
                !     . "jpdf3","log jpdf1","log jpdf2","log jpdf3"'
                !      write(60,*) 'ZONE I=',401,', J=',401,', F=POINT'
                do i=1,n
                    write(60,'(17(1x,es14.4))') i*1.d0,u(i),v(i),T(i),species(i,1),species(i,2),species(i,3),species(i,4),species(i,5),species(i,6),source(i,1),source(i,2),source(i,3),source(i,4),source(i,5),source(i,6),source(i,7)
                enddo
                close(60)
            endif
            do i=1,n
                temp=species(i,1)+species(i,2)+species(i,3)+species(i,4)+species(i,5)+species(i,6)
                species(i,1)=species(i,1)/temp
                species(i,2)=species(i,2)/temp
                species(i,3)=species(i,3)/temp
                species(i,4)=species(i,4)/temp
                species(i,5)=species(i,5)/temp
                species(i,6)=species(i,6)/temp
                species0(i,1)=species(i,1)
                species0(i,2)=species(i,2)
                species0(i,3)=species(i,3)
                species0(i,4)=species(i,4)
                species0(i,5)=species(i,5)
                species0(i,6)=species(i,6)
                T0(i)=T(i)
            enddo

        enddo

        if(mod(iter,5).eq.0) then
            open(60,file='./u.txt',action='write')
            !      write(60,*) '"coef 4,8,4", "coef 8,16,16", "jpdf1","jpdf2",
            !     . "jpdf3","log jpdf1","log jpdf2","log jpdf3"'
            !      write(60,*) 'ZONE I=',401,', J=',401,', F=POINT'
            do i=1,n
                write(60,'(17(1x,es14.4))') i*1.d0,u(i),v(i),T(i),species(i,1),species(i,2),species(i,3),species(i,4),species(i,5),species(i,6),source(i,1),source(i,2),source(i,3),source(i,4),source(i,5),source(i,6),source(i,7)
            enddo
            close(60)
        endif

        do i=1,n
            error=max(error,dabs(T(i)-T0(i)))
        enddo

        do k=1,6
            do i=1,n
                rho(i)=300.d0/T(i)
            enddo
        enddo

        write(*,*) error,T(50000),factor
    enddo

    open(60,file='./u.txt',action='write')
    !      write(60,*) '"coef 4,8,4", "coef 8,16,16", "jpdf1","jpdf2",
    !     . "jpdf3","log jpdf1","log jpdf2","log jpdf3"'
    !      write(60,*) 'ZONE I=',401,', J=',401,', F=POINT'
    do i=1,n
        write(60,'(17(1x,es14.4))') i*1.d0,u(i),v(i),T(i),species(i,1),species(i,2),species(i,3),species(i,4),species(i,5),species(i,6),source(i,1),source(i,2),source(i,3),source(i,4),source(i,5),source(i,6),source(i,7)
    enddo
    close(60)

    end

    subroutine chem()
    double precision, dimension(:), pointer:: T,rho
    double precision, dimension(:,:), pointer:: source,species
    double precision k1c,k2c,k3fc,k3bc,k4fc,k1,k2,k3f,k3b,k4f,k4b,keq4
    double precision temp,tt
    double precision ych4,yco,yh2o,yco2,yh2,yo2
    double precision rywch4,rywo2,rywh2o,rywco2,rywco,rywh2
    double precision wch4,wo2,wco2,wh2o,wco,wh2
    double precision cp
    integer i

    common/small/ n
    common/large0/ source,species
    common/large1/ T,rho

    wch4=16
    wo2=32
    wco2=44
    wh2o=18
    wco=28
    wh2=2
    cp=2.d0

    do i=1,n
        if(i.eq.10025) then
            temp=0.d0
        endif
        tt=T(i)
        ych4=species(i,1)
        yo2=species(i,2)
        yco2=species(i,3)
        yh2o=species(i,4)
        yco=species(i,5)
        yh2=species(i,6)

        k1c=4.4d11
        k2c=3.0d8
        k3fc=2.75d9
        k3bc=6.71d10
        k4fc=2.5d16

        k1=k1c*exp(-15095/tt)
        k2=k2c   *exp(-15095/tt)
        k3f=k3fc *exp(-10065/tt)
        k3b=k3bc *exp(-13688/tt)
        k4f=k4fc*tt**(-1)*exp(-20137/tt)
        keq4=10**(-5.9608*(log10(tt))**5+105.04*(log10(tt))**4-749.52*(log10(tt))**3+2716.0*(log10(tt))**2-5019.5*(log10(tt))+3801.3)*(8.3144*10.d0**(-2)*tt)**0.5
        k4b=k4f/keq4

        rywch4=rho(i)*ych4/wch4
        rywo2=rho(i)*yo2/wo2
        rywco2=rho(i)*yco2/wco2
        rywco=rho(i)*yco/wco
        rywh2o=rho(i)*yh2o/wh2o
        rywh2=rho(i)*yh2/wh2

        source(i,1)=-wch4*k1*(rywch4)**0.5*(rywo2)**1.25-wch4*k2*rywch4*rywh2o
        source(i,2)=-wo2*k1*0.5*(rywch4)**0.5*(rywo2)**1.25-wo2*k4f*0.5*(rywh2)**0.5*(rywh2o+1.d-40)**(-1)*(rywo2)**2.25+wo2*k4b*0.5*(rywh2o)
        source(i,3)=wco2*k3f*rywco*rywh2o-wco2*k3b*rywco2*rywh2
        source(i,4)=-wh2o*k2*rywch4*rywh2o-wh2o*k3f*rywco*rywh2o+wh2o*k3b*rywco2*rywh2+wh2o*k4f*(rywh2)**0.5*(rywh2o+1.d-40)**(-1.0)*(rywo2)**2.25-wh2o*k4b*rywh2o
        source(i,5)=wco*k1*(rywch4)**0.5*(rywo2)**1.25+wco*k2*rywch4*rywh2o-wco*k3f*rywco*rywh2o+wco*k3b*rywco2*rywh2
        source(i,6)=wh2*2*k1*(rywch4)**0.5*(rywo2)**1.25+wh2*3*k2*rywch4*rywh2o+wh2*k3f*rywco*rywh2o-wh2*k3b*rywco2*rywh2-wh2*k4f*(rywh2)**0.5*(rywh2o+1.d-40)**(-1)*(rywo2)**2.25+wh2*k4b*rywh2o
        source(i,7)=74870*source(i,1)/wch4+110525*source(i,5)/wco+241830*source(i,4)/wh2o+393590*source(i,3)/wco2
        source(i,7)=source(i,7)/cp
    enddo

    end