!######################################################################!
MODULE tipos
  IMPLICIT NONE
!
  type:: vetor
   real(kind=8):: x,y,z
  end type vetor
!
  type atomoicula
   type(vetor):: pos,vel,forca,ace,posini
   type(vetor):: K1v,K2v,K3v,K4v,K1a,K2a,K3a,K4a
   type(vetor):: pos0,vel0
   type(vetor):: pkn,pkn1,pkn2,pkn3,vkn,vkn1,vkn2,vkn3
   real(kind=8):: massa,raio
   integer(kind=8) :: nv,nv1,tipo,nv_t
   integer(kind=8), dimension(:),allocatable:: vv,vv1,vv_t
   real(kind=8), dimension(:),allocatable:: costhetaijk0
  end type atomoicula
!
END MODULE tipos
!######################################################################!
!
!######################################################################!
MODULE parametros
  IMPLICIT NONE
!
  real(kind=8), parameter :: Pi=4.0D0*DATAN(1.0D0) ! pi
  real(kind=8), parameter :: d2pi=2*Pi,d4pi=4*Pi  ! 2pi e 4Pi
  real(kind=8), parameter :: dt=1.0D-3 , dt2= dt/2 , dt6=dt/6 , dt24=dt/24
!
END MODULE parametros
!######################################################################!

!######################################################################!
MODULE variables
  USE parametros
  USE tipos
  IMPLICIT NONE
!
  type(atomoicula),dimension(:),allocatable:: atomo

  real(kind=8):: tempo,a0, lx, ly, lz, Temp2, Temp1, teta
  integer(kind=8):: tempo_f,ntprint,Np, Np1, Np2, Np3, n1max, n2max, n3max, n1min, n2min, n3min
  type(vetor) :: a1,a2,a3
  real(kind=8):: gamma,massa0
  real(kind=8) :: raiocorte,rc1,rc2
  real(kind=8):: xcm,ycm,zcm
  real(kind=8):: forcaMolaexterna_x,forcaMolaexterna_y,kmolaext,x0mola,y0mola,tpuxar,tempo_max

  real(kind=8),dimension(0:10):: sigij, epiij,Kfene,rfene,drfene,kmolaAng

  
  integer(kind=8) :: NxDep,NyDep,NzDep,iext
  integer(kind=8) :: isig,jepi,lmas
  integer(kind=8) :: Npcell,Ncell,Nnuc,Nfill


  real(kind=8) :: Ecin, Epot
  character :: nomexmake*60, nomeenergia*60, nomeparaview*60
  character :: num1*8,num2*8,num3*8,num4*8
  external :: i3c,writenumber
  character :: i3c*3

!
END MODULE variables
!######################################################################!



PROGRAM simulaCellModelAngPot!14 fev 2020
 USE parametros
 USE tipos
 USE variables

  integer(kind=8) :: i,j,i1,itmp


 call define_GlobalParameters()

!
 do isig=1,1
  do jepi=0,0
    do lmas=1,10,1
    
    call define_PotParameters()

    tempo_f=int(tempo_max/dt)
    call writenumber(isig,num1)
    call writenumber(jepi,num2)
    call writenumber(lmas,num3)

    nomexmake='saidaxmakemol_isig'//TRIM(num1)//'iepi'//TRIM(num2)//'lmas'//TRIM(num3)//'.xyz' 
    nomeenergia= 'energia_tempo_isig'//TRIM(num1)//'iepi'//TRIM(num2)//'lmas'//TRIM(num3)//'.dat'

    tempo=0.0
    itmp=0

    Call Initialize_system()

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i=1,4
        call define_tabelavizinhos()
        call evolucao_RK4
        call saidaxmakemol
        itmp=itmp+1
        CALL saidaParaView(itmp)
        call imprimeresultados()
    end do

    do i=1,5*tempo_f
        call PreditorCorretor()
        if (mod(i,10)==0) call define_tabelavizinhos()
        if (mod(i,ntprint)==0) then
            CALL saidaxmakemol
            itmp=itmp+1
            !CALL saidaParaView(itmp)
            call imprimeresultados()
        end if
!
    end do
!
    deallocate(atomo)
    end do
  end do
 end do
!
END PROGRAM simulaCellModelAngPot !14 fev 2020

!####################################################
! SUB-ROTINA define_posicoesiniciais_cell
!####################################################
SUBROUTINE define_GlobalParameters()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE


call system('rm *.xyz *.dat *.csv')

n1max=100 ; n2max=100 ; n3max=0
n1min=-100 ; n2min=-100 ; n3min=0
a0=1.D0
!
tempo_max=20.0D0
gamma=0.10D0
massa0=1.D0
rc1=a0*0.58
rc2=a0*1.1D0
raiocorte=3*a0

Temp2=0.05d0
ntprint=int(1/dt)/10
lx=30*a0 ; ly=30*a0 ; lz=3*a0 !!

!
END
!######################################################################!

!####################################################
! SUB-ROTINA define_posicoesiniciais_cell
!####################################################
SUBROUTINE define_PotParameters()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE

 real(kind=8) :: depi

kmolaAng(0)=1000.0
kmolaAng(1)=1000.0
kmolaAng(2)=1000.0

depi=(10.0D0)**jepi

sigij(0)=a0*0.5D0 ; epiij(0)=0.01D0
sigij(1)=a0*0.9D0 ; epiij(1)=0.8D0
sigij(2)=a0*0.9D0 ; epiij(2)=0.6D0
sigij(3)=a0*1.0D0 ; epiij(3)=0.03D0*depi
sigij(4)=a0*1.1D0 ; epiij(4)=0.6D0
sigij(6)=a0*1.1D0 ; epiij(6)=0.6D0

sigij(7)=a0*1.1D0 ; epiij(7)=0.6D0
sigij(8)=a0*1.1D0 ; epiij(8)=0.6D0
sigij(10)=a0*1.1D0 ; epiij(10)=0.6D0

Kfene(0)=500.D0 ; Kfene(1)=500.D0 ; Kfene(2)=500.D0
Kfene(3)=100.D0 ; Kfene(4)=100.D0 ; Kfene(6)=100.D0

rfene(0)=a0*0.5    ;rfene(1)=a0/1  ;rfene(2)=a0
rfene(3)=a0 ;rfene(4)=a0 ;rfene(6)=a0

drfene(0)=a0 ;drfene(1)=a0 ;drfene(2)=a0
drfene(3)=a0  ;drfene(4)=a0  ;drfene(6)=a0

!
END !######################################################################!

!####################################################
! SUB-ROTINA define_posicoesiniciais_cell
!####################################################
SUBROUTINE Initialize_system()

 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE

 integer(kind=8) :: itmp


call posicoesiniciais_cell()
print *,'call posicoesiniciais_cell() - ok'
!print*,'Np=',Np

call define_tabelavizinhos1()
print*,'call define_tabelavizinhos1 - ok'

call define_tipo()
print*,'call define_tipo() - ok'

call define_costhetaijk0()
print*,'call define_costhetaijk0() - ok'

call define_tabelavizinhos()
print*,'call define_tabelavizinhos - ok'

call define_vel_inicial
print*,'call define_vel_inicial - ok'

CALL saidaxmakemol()
print*,'call saidaxmakemol  - ok'

itmp=0
CALL saidaParaView(itmp)
print *,'call saidaParaView(0)  - ok'

call define_massa()
print*,'call define_massa() - ok'

call define_ace()
print*,'call define_ace() - ok'

call imprimeresultados()
print*,'call imprimeresultados() - ok'


!
END !######################################################################!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE writenumber(i2,number1)
integer, intent(in)   :: i2
character, intent(out)  :: number1*8

911 FORMAT('_00000'I1'_')
912 FORMAT('_0000'I2'_')
913 FORMAT('_000'I3'_')
914 FORMAT('_00'I4'_')
915 FORMAT('_0'I5'_')
916 FORMAT('_'I6'_')

if (i2.lt.10) then
    WRITE(number1,911) i2
end if
if ((i2.ge.10).and.(i2.lt.100)) then
    WRITE(number1,912) i2
endif
if ((i2.ge.100).and.(i2.lt.1000)) then
    WRITE(number1,913) i2
endif
if ((i2.ge.1000).and.(i2.lt.10000)) then
    WRITE(number1,914) i2
endif
if ((i2.ge.10000).and.(i2.lt.100000)) then
    WRITE(number1,915) i2
endif
if ((i2.ge.100000).and.(i2.lt.1000000)) then
    WRITE(number1,916) i2
endif

return
end subroutine

!####################################################
! SUB-ROTINA define_posicoesiniciais_cell
!####################################################
SUBROUTINE posicoesiniciais_cell()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE
 real(kind=8) :: aleax0,aleay0,aleaz0,phi,theta,delTheta,delphi,raio,raio1,raio2
 real(kind=8):: rx,ry,rz,rx1,ry1,rz1,rx2,ry2,rz2,rx0,ry0,rz0,Rcell,Rnuc,thetacito,raiocito,raioteste
 real(kind=8):: dcelx,dcely
 integer(kind=8) :: j4,j41,j42,k2,k3,k4,k5,j,j1,j5,i,n,n1,Ncx,Ncy,n2,n3,tiporede
 integer,parameter :: seed = 1!8873456

!
j=0
j1=0

k2=0
k3=0

k4=0
k5=0

Rcell=16.D0
Rnuc=Rcell/3

Npcell=int(2*Pi*Rcell*2)
Nnuc=int(2*Pi*Rnuc)
Ncx=5
Ncy=5
Ncell=Ncx*Ncy
Nfill=0
tiporede=4


do i=1,2

 call srand(seed)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !determina o numero de particulas do colageno
 delTheta=2*Pi/Npcell
 do j41=0,Ncx-1
 dcelx=j41*(Rcell+a0)*2.1D0

 do j42=0,Ncy-1
 dcely=j42*(Rcell+a0)*2.1D0
  
 do j4=1,Npcell
  rx=dcelx+(Rcell+a0)*dcos(delTheta*j4) ; ry=dcely+(Rcell+a0)*dsin(delTheta*j4) ; rz=0.D0
  if (i==1) then
   j1=j1+1
  else
   j=j+1
   atomo(j)%pos%x=1*lx/2+rx;  atomo(j)%pos%y=1*ly/2+ry;  atomo(j)%pos%z=1*lz/2+rz
  end if
  rx0=rx ; ry0=ry ; rz0=rz
 end do
 
 end do
 end do
 !
  Np1=j1 !determina o numero de particulas na celula
  k3=j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !determina o numero de particulas do colageno
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !determina o numero de particulas do elastina

SELECT CASE (tiporede)
   CASE (1)
      !   cubica
        a1%x=a0   ; a1%y=0.0D0; a1%z=0.0D0
        a2%x=0.0D0; a2%y=a0   ; a2%z=0.0D0
        a3%x=0.0D0; a3%y=0.0D0; a3%z=0.0D0
        
        !do n3=n3min, n3max
        do n2=n2min, n2max
            do n1=n1min, n1max
            rx=n1*a1%x+n2*a2%x+n3*a3%x
            ry=n1*a1%y+n2*a2%y+n3*a3%y
            rz=n1*a1%z+n2*a2%z+n3*a3%z
            raio=dsqrt(rx**2+ry**2+rz**2)
            if (raio<Rcell) then 
                if (i==1) then
                    k2=k2+1
                else
                    k3=k3+1
                    atomo(k3)%pos%x=lx/2+rx;  atomo(k3)%pos%y=ly/2+ry;  atomo(k3)%pos%z=lz/2+rz
                end if
                rx0=rx ; ry0=ry ; rz0=rz
            end if
            end do
        end do
        !end do

   CASE (2)
      !   triangular
        a1%x=a0   ; a1%y=0.0D0; a1%z=0.0D0
        a2%x=a0/2 ; a2%y=a0*dsqrt(3.0D0)/2   ; a2%z=0.0D0
        a3%x=0.0D0; a3%y=0.0D0; a3%z=0.0D0
        !do n3=n3min, n3max
        do n2=n2min, n2max
            do n1=n1min, n1max
            rx=n1*a1%x+n2*a2%x+n3*a3%x
            ry=n1*a1%y+n2*a2%y+n3*a3%y
            rz=n1*a1%z+n2*a2%z+n3*a3%z
            raio=dsqrt(rx**2+ry**2+rz**2)
            if ((raio<Rcell).and.(raio>Rcell/1.5)) then 
                if (i==1) then
                    k2=k2+1
                else
                    k3=k3+1
                    atomo(k3)%pos%x=lx/2+rx;  atomo(k3)%pos%y=ly/2+ry;  atomo(k3)%pos%z=lz/2+rz
                end if
                rx0=rx ; ry0=ry ; rz0=rz
            end if
            end do
        end do
        !end do

   CASE (3)
      !   Kagome
        a1%x=a0*2   ; a1%y=0.0D0; a1%z=0.0D0
        a2%x=a0*2/2 ; a2%y=a0*2*dsqrt(3.0D0)/2   ; a2%z=0.0D0
        a3%x=0.0D0; a3%y=0.0D0; a3%z=0.0D0
        
        !do n3=n3min, n3max
        do n2=n2min, n2max
            do n1=n1min, n1max
            rx=n1*a1%x+n2*a2%x+n3*a3%x
            ry=n1*a1%y+n2*a2%y+n3*a3%y
            rz=n1*a1%z+n2*a2%z+n3*a3%z

            aleax0=rand()*0.0
            rx1=rx+a1%x/(2+aleax0)
            ry1=ry+a1%y/(2+aleax0)
            rz1=rz+a1%z/(2+aleax0)

            rx2=rx+a2%x/(2+aleax0)
            ry2=ry+a2%y/(2+aleax0)
            rz2=rz+a2%z/(2+aleax0)
    
            raio=dsqrt(rx**2+ry**2+rz**2)
            raio1=dsqrt(rx1**2+ry1**2+rz1**2)
            raio2=dsqrt(rx2**2+ry2**2+rz2**2)
    
            if ( (raio<Rcell).and.(raio1<Rcell).and.(raio2<Rcell) ) then 
     
                if (i==1) then
                    k2=k2+1
                    k2=k2+1
                    k2=k2+1
                else
                    k3=k3+1
                    atomo(k3)%pos%x=lx/2+rx;  atomo(k3)%pos%y=ly/2+ry;  atomo(k3)%pos%z=lz/2+rz
                    k3=k3+1
                    atomo(k3)%pos%x=lx/2+rx1;  atomo(k3)%pos%y=ly/2+ry1;  atomo(k3)%pos%z=lz/2+rz1
                    k3=k3+1
                    atomo(k3)%pos%x=lx/2+rx2;  atomo(k3)%pos%y=ly/2+ry2;  atomo(k3)%pos%z=lz/2+rz2
                end if
                rx0=rx ; ry0=ry ; rz0=rz
            end if
            end do
        end do
        !end do
   CASE (4)
    Nfill=0
    delTheta=2*Pi/Nnuc
    aleax0=rand()*(Rcell-Rnuc-a0)*0.0
    aleay0=rand()*(Rcell-Rnuc-a0)*0.0
    
    dcelx=0.0
    dcely=0.0
    do j41=0,Ncx-1
    dcelx=j41*(Rcell+a0)*2.1D0
    do j42=0,Ncy-1
    dcely=j42*(Rcell+a0)*2.1D0
    
    do j4=1,Nnuc
     j5=int(Rnuc/a0)
     do 
     rx=j5*a0*dcos(delTheta*j4) ; ry=j5*a0*dsin(delTheta*j4) ; rz=0.D0
     raioteste=(rx-aleax0)**2+(ry-aleay0)**2
     
     if (raioteste.ge.Rcell*Rcell*0.99) then
       exit
     end if
     
     if ( (j5 .eq. (int(Rnuc/a0))).and.(mod(j4,3).ne.0) ) then
       j5=j5+1
     end if

     j5=j5+1
     
     if (i==1) then
       k2=k2+1
     else
      Nfill=Nfill+1
      k3=k3+1
      atomo(k3)%pos%x=dcelx+1*lx/2+rx-aleax0;  atomo(k3)%pos%y=dcely+1*ly/2+ry-aleay0;  atomo(k3)%pos%z=1*lz/2+rz
     end if
     rx0=rx ; ry0=ry ; rz0=rz



     end do
    
    end do
    end do
    
    end do

   CASE DEFAULT
      WRITE(*,*)  "--------------"
END SELECT

  !
  Np2=k2  !determina o numero de particulas lestina
  k5=k3
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !determina o numero de particulas do elastina
  Np=Np1+Np2!determina o numero de particulas totais
 !  print*,Np1,Np2,Np
 ! 
! 
  if (i==1) then
   allocate(atomo(Np))
  end if
 
end do
!
atomo%posini%x=atomo%pos%x
atomo%posini%y=atomo%pos%y
atomo%posini%z=atomo%pos%z
!
 RETURN
END

!####################################################
! SUB-ROTINA define_massa
!####################################################
SUBROUTINE define_massa()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE

 integer(kind=8) :: i
!
do i=1,Np1
atomo(i)%massa=massa0
end do!

do i=Np1+1,Np
atomo(i)%massa=massa0
end do

  open(45,file='Define_massa_atomoiculas.dat')
  do i=1,Np
    write(45,*) i,atomo(i)%massa
  end do
  close(45)

 RETURN
END

!####################################################
! SUB-ROTINA define_massa
!####################################################
SUBROUTINE define_tipo()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE

 integer(kind=8) :: i,j,jv
!
!!!! Definicao do tipo 0
do i=1,Np1
atomo(i)%tipo=0
end do!

!!!! Definicao do tipo 1
do i=Np1+1,Np
atomo(i)%tipo=1
end do

!!!! Definicao do tipo 3
do i=Np1+1,Np
 if (atomo(i)%nv1==0) then
   atomo(i)%tipo=3
 end if
end do

 RETURN
END

!####################################################
! SUB-ROTINA define_tabeladevizinhos
!####################################################
SUBROUTINE define_tabelavizinhos1()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE
 integer(kind=8) :: i,j,jv,nvaux,nvaux1
 real(kind=8) :: xij,yij,zij,rij,xi,yi,zi

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do i=1,Np1
  xi=atomo(i)%pos%x ; yi=atomo(i)%pos%y ; zi=atomo(i)%pos%z
  atomo(i)%nv1=0
  atomo(i)%nv_t=0
  do j=1,Np1
   xij=atomo(j)%pos%x-xi ; yij=atomo(j)%pos%y-yi ; zij=atomo(j)%pos%z-zi
   rij= DSQRT(xij**2 + yij**2 + zij**2)
   if(rij.le.rc1) then
        if (j>i) then
            atomo(i)%nv1=atomo(i)%nv1+1
        end if
        if (rij>0.0D0) then
            atomo(i)%nv_t=atomo(i)%nv_t+1
        end if
   end if
  end do


  allocate(atomo(i)%vv1(atomo(i)%nv1))
  allocate(atomo(i)%vv_t(atomo(i)%nv_t))  

  nvaux=0
  nvaux1=0 
  atomo(i)%vv1=0
  atomo(i)%vv_t=0
  
  do j=1,Np1
   xij=atomo(j)%pos%x-xi ; yij=atomo(j)%pos%y-yi ; zij=atomo(j)%pos%z-zi
   rij= DSQRT(xij**2 + yij**2 + zij**2)
   
   if(rij.le.rc1) then
        if (j>i) then
             nvaux=nvaux+1
             atomo(i)%vv1(nvaux)=j
        end if
        if (rij>0.0D0) then
             nvaux1=nvaux1+1
             atomo(i)%vv_t(nvaux1)=j
        end if
   end if
   
  end do
 end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do i=Np1+1,Np
  xi=atomo(i)%pos%x ; yi=atomo(i)%pos%y ; zi=atomo(i)%pos%z
  atomo(i)%nv1=0
  atomo(i)%nv_t=0
  
  do j=Np1+1,Np
   xij=atomo(j)%pos%x-xi ; yij=atomo(j)%pos%y-yi ; zij=atomo(j)%pos%z-zi
   rij= DSQRT(xij**2 + yij**2 + zij**2)
   
   if(rij.le.rc2) then
        if (j>i) then
            atomo(i)%nv1=atomo(i)%nv1+1
        end if
        if (rij>0.0D0) then
            atomo(i)%nv_t=atomo(i)%nv_t+1
        end if
   end if
 
  end do
  
  allocate(atomo(i)%vv1(atomo(i)%nv1))
  allocate(atomo(i)%vv_t(atomo(i)%nv_t))  
  
  nvaux=0
  nvaux1=0 
  atomo(i)%vv1=0
  atomo(i)%vv_t=0
  
  do j=Np1+1,Np
   xij=atomo(j)%pos%x-xi ; yij=atomo(j)%pos%y-yi ; zij=atomo(j)%pos%z-zi
   rij= DSQRT(xij**2 + yij**2 + zij**2)
   
   if(rij.le.rc2) then
        if (j>i) then
             nvaux=nvaux+1
             atomo(i)%vv1(nvaux)=j
        end if
        if (rij>0.0D0) then
             nvaux1=nvaux1+1
             atomo(i)%vv_t(nvaux1)=j
        end if
   end if
   
  end do
 end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  open(45,file='tabeladevizinhos_spring.dat')
  do i=1,Np
   do jv=1,atomo(i)%nv1
    j=atomo(i)%vv1(jv)
     write(45,*) i,atomo(i)%nv1,j
    end do 
   end do
  close(45)

  open(45,file='tabeladevizinhos_angular.dat')
  do i=1,Np
   do jv=1,atomo(i)%nv_t
     j=atomo(i)%vv_t(jv)
     write(45,*) i,atomo(i)%nv_t,j
    end do 
   end do
  close(45)

 RETURN
END

!####################################################
! SUB-ROTINA define_tabeladevizinhos
!####################################################
SUBROUTINE define_tabelavizinhos()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE
 integer(kind=8) :: i,j,jv,nvaux
 real(kind=8) :: xij,yij,zij,rij,xi,yi,zi
!
 do i=1,Np
  xi=atomo(i)%pos%x ; yi=atomo(i)%pos%y ; zi=atomo(i)%pos%z
  atomo(i)%nv=0
  do j=1,Np
   xij=atomo(j)%pos%x-xi ; yij=atomo(j)%pos%y-yi ; zij=atomo(j)%pos%z-zi
   rij= DSQRT(xij**2 + yij**2 + zij**2)
   if ((j>i).and.(rij<=raiocorte)) then
    atomo(i)%nv=atomo(i)%nv+1
   end if
  end do
  
  if (allocated(atomo(i)%vv)) then
   deallocate(atomo(i)%vv)
   allocate(atomo(i)%vv(atomo(i)%nv))
  else
   allocate(atomo(i)%vv(atomo(i)%nv))
  end if
  nvaux=0
  do j=1,Np
    xij=atomo(j)%pos%x-xi ; yij=atomo(j)%pos%y-yi ; zij=atomo(j)%pos%z-zi
    rij= DSQRT(xij**2 + yij**2 + zij**2)
    if ((j>i).and.(rij<=raiocorte)) then
     nvaux=nvaux+1
     atomo(i)%vv(nvaux)=j
    end if
  end do
 end do
!
  open(45,file='tabeladevizinhos_LJ.dat')
  do i=1,Np
   do jv=1,atomo(i)%nv
    j =atomo(i)%vv(jv)
     write(45,*) i,atomo(i)%nv,j
    end do 
   end do
  close(45)


 RETURN
END

!####################################################
! SUB-ROTINA costhetaijk0()
!####################################################
SUBROUTINE define_costhetaijk0()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE

 integer(kind=8) :: i,j,k,jv,kv,naux
 real(kind=8) :: xij,yij,zij,rij2,rij
 real(kind=8) :: xik,yik,zik,rik2,rik
 
!
  do i=1,Np
     naux=0
     do jv=1,atomo(i)%nv_t
     j=atomo(i)%vv_t(jv)
        do kv=jv+1,atomo(i)%nv_t
            k=atomo(i)%vv_t(kv)
            naux=naux+1
        end do
     end do 
     allocate(atomo(i)%costhetaijk0(naux))
  end do

!
  do i=1,Np
     naux=0
     do jv=1,atomo(i)%nv_t
     j=atomo(i)%vv_t(jv)
     xij=atomo(j)%pos%x-atomo(i)%pos%x
     yij=atomo(j)%pos%y-atomo(i)%pos%y
     zij=atomo(j)%pos%z-atomo(i)%pos%z
     rij2=xij**2 + yij**2 + zij**2
     rij=DSQRT(rij2)
        do kv=jv+1,atomo(i)%nv_t
            k=atomo(i)%vv_t(kv)
            xik=atomo(k)%pos%x-atomo(i)%pos%x
            yik=atomo(k)%pos%y-atomo(i)%pos%y
            zik=atomo(k)%pos%z-atomo(i)%pos%z
            rik2=xik**2 + yik**2 + zik**2
            rik=DSQRT(rik2)
            naux=naux+1
            atomo(i)%costhetaijk0(naux)=(xij*xik+yij*yik+zij*zik)/(rij*rik)
        end do
     end do 
  end do
 
  
 open(45,file='define_costhetaijk0.dat')
!
  do i=1,Np
     naux=0
     do jv=1,atomo(i)%nv_t
     j=atomo(i)%vv_t(jv)
     rij=DSQRT(rij2)
        do kv=jv+1,atomo(i)%nv_t
            k=atomo(i)%vv_t(kv)
            naux=naux+1
            write(45,"(3I5,1f15.5)") i,j,k,atomo(i)%costhetaijk0(naux)
        end do
     end do 
  end do
  close(45)
  
  
  
  
  

 RETURN
END

!####################################################
! SUB-ROTINA define_vel_inicial
!####################################################
SUBROUTINE define_vel_inicial()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE
  integer(kind=8) :: i,i1,indexi
 integer :: seed = 86456
 real (kind=8) :: vx, vy, vz, v,vcm0,Tred,Hvx,ri1,xi,yi,zi
!
 seed=86456*lmas
 call srand(seed)
 
 Tred=1.0D-3
 v=1.0D-1
!
 do i=1,Np
    vx=(rand()*2.D0-1.D0)*v
    vy=(rand()*2.D0-1.D0)*v
    vz=0.0D0
    
45 continue
!    fvx=dsqrt(massa0/(d2pi*Tred))*dexp(-massa*0.5D0*vx*vx/Tred)
    Hvx=dexp(-massa0*0.5D0*vx*vx/Tred)
    if (rand().le.Hvx) then
        atomo(i)%vel%x=vx
        atomo(i)%vel%y=vy
        atomo(i)%vel%z=vz
    else
        go to 45
    end if

end do

  xcm=sum(atomo%pos%x)/Np
  ycm=sum(atomo%pos%y)/Np
  zcm=sum(atomo%pos%z)/Np    

  do i=1,Np
    xi=atomo(i)%pos%x-xcm ; yi=atomo(i)%pos%y-ycm ; zi=atomo(i)%pos%z-zcm
    ri1=1.d0/dsqrt(xi**2+yi**2+zi**2)
    atomo(i)%vel%x=atomo(i)%vel%x-0.1*xi
    atomo(i)%vel%y=atomo(i)%vel%y-0.1*yi
    atomo(i)%vel%z=atomo(i)%vel%z-0.1*zi
  end do



 RETURN
END

!####################################################
! SUB-ROTINA saidaxmakemol
!####################################################
SUBROUTINE saidaxmakemol()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE
 integer(kind=8) :: i
 
!
 25 format('Fe  ',1F16.6,1F16.6,1F16.6,' atom_vector ',1F16.6,1F16.6,1F16.6)
 26 format('Al  ',1F16.6,1F16.6,1F16.6,' atom_vector ',1F16.6,1F16.6,1F16.6)
 22 format('Au  ',1F16.6,1F16.6,1F16.6,' atom_vector ',1F16.6,1F16.6,1F16.6)
 23 format('C  ',1F16.6,1F16.6,1F16.6,' atom_vector ',1F16.6,1F16.6,1F16.6)
 24 format('H  ',1F16.6,1F16.6,1F16.6,' atom_vector ',1F16.6,1F16.6,1F16.6)

open(41,file=TRIM(nomexmake),ACCESS='APPEND')
write(41,*) Np
write(41,*) tempo

do i=1,Np
SELECT CASE (atomo(i)%tipo)
   CASE (0)
    write(41,22) atomo(i)%pos%x,atomo(i)%pos%y,atomo(i)%pos%z,atomo(i)%vel%x,atomo(i)%vel%y,atomo(i)%vel%z
   CASE (1)
    write(41,25) atomo(i)%pos%x,atomo(i)%pos%y,atomo(i)%pos%z,atomo(i)%vel%x,atomo(i)%vel%y,atomo(i)%vel%z
   CASE (3)
    write(41,23) atomo(i)%pos%x,atomo(i)%pos%y,atomo(i)%pos%z,atomo(i)%vel%x,atomo(i)%vel%y,atomo(i)%vel%z
END SELECT
end do

close(41)
!
 RETURN
END

!####################################################
! SUB-ROTINA saidaxmakemol
!####################################################
SUBROUTINE saidaParaView(temp_j)
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE
 integer(kind=8) ::  i,temp_j
 character :: num_t*8
!

call writenumber(temp_j,num_t)

28 format(1F16.6,',',1F16.6,',',1F16.6,',',1F16.6,',',1F16.6,',',1F16.6)

nomeparaview='Cell_sur'//TRIM(num_t)//'.csv'

open(21,file=TRIM(nomeparaview),ACCESS='APPEND')
write(21,*) 'x,y,z,vx,vy,vz'
do i=1,Np1
write(21,28) atomo(i)%pos%x,atomo(i)%pos%y,atomo(i)%pos%z,atomo(i)%vel%x,atomo(i)%vel%y,atomo(i)%vel%z
end do

nomeparaview='Cell_fil'//TRIM(num_t)//'.csv'
open(21,file=TRIM(nomeparaview),ACCESS='APPEND')
write(21,*) 'x,y,z,vx,vy,vz'
do i=Np1+1,Np
write(21,28) atomo(i)%pos%x,atomo(i)%pos%y,atomo(i)%pos%z,atomo(i)%vel%x,atomo(i)%vel%y,atomo(i)%vel%z
end do
close(21)


!
 RETURN
END

!####################################################
! SUB-ROTINA calculaenercin
!####################################################
SUBROUTINE calculaenercin()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE
 integer(kind=8) :: i
!
Ecin=0.0D0
!
 do i=1,Np
  Ecin=Ecin+atomo(i)%massa*(atomo(i)%vel%x**2+atomo(i)%vel%y**2+atomo(i)%vel%z**2)/2
 end do
!
 RETURN
END

!####################################################
! SUB-ROTINA renormavel
!####################################################
SUBROUTINE renormavel()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE
 integer(kind=8) :: i,num
 real(kind=8) :: aux,Ecin1
!
num=0
 do i=1,Np1
   num=num+1
   Ecin1=Ecin1+atomo(i)%massa*(atomo(i)%vel%x**2+atomo(i)%vel%y**2+atomo(i)%vel%z**2)/2
 end do

!Temp1=2*Ecin1/(3*num) ! para renormalizar em 3D
Temp1=Ecin1/(num) ! para renormalizar em 2D
aux=dsqrt(Temp2/Temp1)
!
 do i=1,Np1
  atomo(i)%vel%x=atomo(i)%vel%x*aux
  atomo(i)%vel%y=atomo(i)%vel%y*aux
  atomo(i)%vel%z=atomo(i)%vel%z*aux 
end do
!
 RETURN
END

!####################################################
! SUB-ROTINA imprimeresultados
!####################################################
SUBROUTINE imprimeresultados()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE

!
 32 format(10E20.6)
open(21,file=TRIM(nomeenergia),ACCESS='APPEND')
call calculaenercin()
write(21,32) tempo, Ecin/Np, Epot/Np, (Ecin+Epot)/Np
close(21)

!open(21,file='tempo_z_vez_azNp.dat',ACCESS='APPEND')
!write(21,32) tempo, atomo(Np)%pos%z,atomo(Np)%vel%z,atomo(Np)%ace%z
!close(21)

!
!open(21,file=TRIM(nomeangulo),ACCESS='APPEND')
!write(21,32) tempo, teta, atomo(Np)%vel%x,atomo(1)%pos%x,atomo(1)%vel%x
!close(21)
!
 RETURN
END

!####################################################
! SUB-ROTINA define_ace()
!####################################################
SUBROUTINE define_ace()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE
 integer(kind=8) :: i,i1,j,k,jv,kv,naux,indexi,indexj
 real(kind=8) :: xi,yi,zi
 real(kind=8) :: xij,yij,zij,rij2,rij,rij1,forcamola,forcalj,sigrij2,sigrij6,sigrij8,sig,epi
 real(kind=8) :: xik,yik,zik,rik2,rik,rik1,costhetaijk,rijdotrij,costhetaijkMthetaijk0
 real(kind=8) :: rijr02,forcaFENE,auxKfene,auxRfene,auxDRfene,auxrijfene
 !
   atomo%forca%x=0.0D0
   atomo%forca%y=0.0D0
   atomo%forca%z=0.0D0
!
 Epot=0.0D0
!
 do i=1,Np
  xi=atomo(i)%pos%x ; yi=atomo(i)%pos%y ; zi=atomo(i)%pos%z

!  loop LENNARD JONES 
   do jv=1,atomo(i)%nv
    j=atomo(i)%vv(jv)
    xij=atomo(j)%pos%x-xi ; yij=atomo(j)%pos%y-yi ; zij=atomo(j)%pos%z-zi    
    rij2=xij**2 + yij**2 + zij**2

!   rij=DSQRT(rij2)
    sig=sigij(atomo(i)%tipo+atomo(j)%tipo)
    epi=epiij(atomo(i)%tipo+atomo(j)%tipo)
    
    if (atomo(i)%tipo.ne.atomo(j)%tipo) then
        sigrij2=sig*sig/rij2
        sigrij6=sigrij2*sigrij2*sigrij2
        sigrij8=sigrij6*sigrij2
        Epot=Epot+epi*sigrij6*(sigrij6-2)
        forcalj=12*epi*sigrij8*(1-sigrij6)/(sig*sig)
    else
      if (rij2 <= sig) then
        sigrij2=sig*sig/rij2
        sigrij6=sigrij2*sigrij2*sigrij2
        sigrij8=sigrij6*sigrij2
        Epot=Epot+epi*sigrij6*(sigrij6-2)+epi
        forcalj=12*epi*sigrij8*(1-sigrij6)/(sig*sig)
      else
        Epot=0.D0
        forcalj=0.D0
      end if
    end if
    atomo(i)%forca%x=atomo(i)%forca%x+forcalj*xij
    atomo(i)%forca%y=atomo(i)%forca%y+forcalj*yij
    atomo(i)%forca%z=atomo(i)%forca%z+forcalj*zij

    atomo(j)%forca%x=atomo(j)%forca%x-forcalj*xij
    atomo(j)%forca%y=atomo(j)%forca%y-forcalj*yij
    atomo(j)%forca%z=atomo(j)%forca%z-forcalj*zij

   end do 
 
 !  loop FENE
  do jv=1,atomo(i)%nv1
      j=atomo(i)%vv1(jv)
      xij=atomo(j)%pos%x-xi ; yij=atomo(j)%pos%y-yi ; zij=atomo(j)%pos%z-zi
      rij2=xij**2 + yij**2 + zij**2
      rij=DSQRT(rij2)
!
   auxKfene=Kfene(atomo(i)%tipo+atomo(j)%tipo)
   auxRfene=rfene(atomo(i)%tipo+atomo(j)%tipo)
   auxDRfene=drfene(atomo(i)%tipo+atomo(j)%tipo)

   rijr02=1.D0-((rij-auxRfene)/auxDRfene)**2
   
   Epot=Epot-0.5D0*auxKfene*dlog(rijr02)*auxDRfene**2
   forcaFENE=auxKfene*(1.D0-auxRfene/rij)/rijr02
!   print*,i,j,rij2
!   print*,rijr02,dlog(rijr02),forcaFENE
!   
   atomo(i)%forca%x=atomo(i)%forca%x+forcaFENE*xij
   atomo(i)%forca%y=atomo(i)%forca%y+forcaFENE*yij
   atomo(i)%forca%z=atomo(i)%forca%z+forcaFENE*zij

   atomo(j)%forca%x=atomo(j)%forca%x-forcaFENE*xij
   atomo(j)%forca%y=atomo(j)%forca%y-forcaFENE*yij
   atomo(j)%forca%z=atomo(j)%forca%z-forcaFENE*zij

    end do 
    
    !  loop bending
    naux=0
   do jv=1,atomo(i)%nv_t
    j=atomo(i)%vv_t(jv)
    xij=atomo(j)%pos%x-xi
    yij=atomo(j)%pos%y-yi
    zij=atomo(j)%pos%z-zi
    rij2=xij**2 + yij**2 + zij**2
    rij=DSQRT(rij2)
    rij1=1.D0/rij
     do kv=jv+1,atomo(i)%nv_t
      k=atomo(i)%vv_t(kv)
      xik=atomo(k)%pos%x-atomo(i)%pos%x
      yik=atomo(k)%pos%y-atomo(i)%pos%y
      zik=atomo(k)%pos%z-atomo(i)%pos%z
      rik2=xik**2 + yik**2 + zik**2
      rik=DSQRT(rik2)
      rik1=1.D0/rik
      rijdotrij=xij*xik+yij*yik+zij*zik
      costhetaijk=rijdotrij*rij1*rik1
      naux=naux+1
      
      costhetaijkMthetaijk0=costhetaijk-atomo(i)%costhetaijk0(naux)
      
      Epot=Epot+0.5D0*kmolaAng(atomo(i)%tipo+atomo(j)%tipo)*costhetaijkMthetaijk0**2
      forcamola=kmolaAng(atomo(i)%tipo+atomo(j)%tipo)*costhetaijkMthetaijk0*rij1*rik1

      atomo(i)%forca%x=atomo(i)%forca%x+forcamola*(xij+xik-rijdotrij*(rij1*rij1*xij+rik1*rik1*xik) )
      atomo(i)%forca%y=atomo(i)%forca%y+forcamola*(yij+yik-rijdotrij*(rij1*rij1*yij+rik1*rik1*yik) )
      atomo(i)%forca%z=atomo(i)%forca%z+forcamola*(zij+zik-rijdotrij*(rij1*rij1*zij+rik1*rik1*zik) )
! 
      atomo(j)%forca%x=atomo(j)%forca%x-forcamola*(xik-rijdotrij*rij1*rij1*xij)
      atomo(j)%forca%y=atomo(j)%forca%y-forcamola*(yik-rijdotrij*rij1*rij1*yij)
      atomo(j)%forca%z=atomo(j)%forca%z-forcamola*(zik-rijdotrij*rij1*rij1*zij)

      atomo(k)%forca%x=atomo(k)%forca%x-forcamola*(xij-rijdotrij*rik1*rik1*xik)
      atomo(k)%forca%y=atomo(k)%forca%y-forcamola*(yij-rijdotrij*rik1*rik1*yik)
      atomo(k)%forca%z=atomo(k)%forca%z-forcamola*(zij-rijdotrij*rik1*rik1*zik)
! 
     end do
   end do 
    
 end do
 
  do i=1,Np1
  xi=atomo(i)%pos%x ; yi=atomo(i)%pos%y ; zi=atomo(i)%pos%z
 !  loop LENNARD JONES 
    do jv=1,atomo(i)%nv
     j=atomo(i)%vv(jv)
     if (atomo(i)%tipo==atomo(j)%tipo) then
      indexi=int(i/(Npcell+1))+1
      indexj=int(j/(Npcell+1))+1
      if ( indexj > indexi) then
!      write(*,*) i,j,indexi,indexj
! 
           xij=atomo(j)%pos%x-xi ; yij=atomo(j)%pos%y-yi ; zij=atomo(j)%pos%z-zi
           rij2=xij**2 + yij**2 + zij**2
 
       !   rij=DSQRT(rij2)
           sig=sigij(atomo(i)%tipo+atomo(j)%tipo)
           epi=epiij(atomo(i)%tipo+atomo(j)%tipo)
  !     
           sigrij2=sig*sig/rij2
           sigrij6=sigrij2*sigrij2*sigrij2
           sigrij8=sigrij6*sigrij2
           Epot=Epot+epi*sigrij6*(sigrij6-2)
           forcalj=12*epi*sigrij8*(1-sigrij6)/(sig*sig)
  !   
       atomo(i)%forca%x=atomo(i)%forca%x+forcalj*xij
       atomo(i)%forca%y=atomo(i)%forca%y+forcalj*yij
       atomo(i)%forca%z=atomo(i)%forca%z+forcalj*zij
  ! 
       atomo(j)%forca%x=atomo(j)%forca%x-forcalj*xij
       atomo(j)%forca%y=atomo(j)%forca%y-forcalj*yij
       atomo(j)%forca%z=atomo(j)%forca%z-forcalj*zij
     end if
     end if

   end do 
 end do
! 
do i=1,Np
 atomo(i)%ace%x=atomo(i)%forca%x/atomo(i)%massa
 atomo(i)%ace%y=atomo(i)%forca%y/atomo(i)%massa
 atomo(i)%ace%z=atomo(i)%forca%z/atomo(i)%massa
end do

do i=1,Np
 atomo(i)%ace%x=atomo(i)%ace%x-gamma*atomo(i)%vel%x
 atomo(i)%ace%y=atomo(i)%ace%y-gamma*atomo(i)%vel%y
 atomo(i)%ace%z=atomo(i)%ace%z-gamma*atomo(i)%vel%z
end do
!
 RETURN
END

!####################################################
! SUB-ROTINA define_evoluçao_rk4
!####################################################
SUBROUTINE evolucao_RK4()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE
 integer(kind=8) :: i
 real(kind=8):: tempo0

 tempo0=tempo
 call define_ace
 do i=1,Np
 atomo(i)%pos0%x=atomo(i)%pos%x;  atomo(i)%pos0%y=atomo(i)%pos%y;  atomo(i)%pos0%z=atomo(i)%pos%z
 atomo(i)%vel0%x=atomo(i)%vel%x;  atomo(i)%vel0%y=atomo(i)%vel%y;  atomo(i)%vel0%z=atomo(i)%vel%z
 end do

 do i=1,Np
 atomo(i)%k1v%x=atomo(i)%vel%x; atomo(i)%k1v%y=atomo(i)%vel%y; atomo(i)%k1v%z=atomo(i)%vel%z
 atomo(i)%k1a%x=atomo(i)%ace%x; atomo(i)%k1a%y=atomo(i)%ace%y; atomo(i)%k1a%z=atomo(i)%ace%z
 end do
!
 do i=1,Np
 atomo(i)%pos%x=atomo(i)%pos0%x+atomo(i)%k1v%x*dt2
 atomo(i)%pos%y=atomo(i)%pos0%y+atomo(i)%k1v%y*dt2
 atomo(i)%pos%z=atomo(i)%pos0%z+atomo(i)%k1v%z*dt2
 
 atomo(i)%vel%x=atomo(i)%vel0%x+atomo(i)%k1a%x*dt2
 atomo(i)%vel%y=atomo(i)%vel0%y+atomo(i)%k1a%y*dt2
 atomo(i)%vel%z=atomo(i)%vel0%z+atomo(i)%k1a%z*dt2
 end do

 tempo=tempo0+dt2
 call define_ace

 do i=1,Np
 atomo(i)%k2v%x=atomo(i)%vel%x; atomo(i)%k2v%y=atomo(i)%vel%y; atomo(i)%k2v%z=atomo(i)%vel%z
 atomo(i)%k2a%x=atomo(i)%ace%x; atomo(i)%k2a%y=atomo(i)%ace%y; atomo(i)%k2a%z=atomo(i)%ace%z
 end do
!
 do i=1,Np
 atomo(i)%pos%x=atomo(i)%pos0%x+atomo(i)%k2v%x*dt2
 atomo(i)%pos%y=atomo(i)%pos0%y+atomo(i)%k2v%y*dt2
 atomo(i)%pos%z=atomo(i)%pos0%z+atomo(i)%k2v%z*dt2
 
 atomo(i)%vel%x=atomo(i)%vel0%x+atomo(i)%k2a%x*dt2
 atomo(i)%vel%y=atomo(i)%vel0%y+atomo(i)%k2a%y*dt2
 atomo(i)%vel%z=atomo(i)%vel0%z+atomo(i)%k2a%z*dt2
 end do

 tempo=tempo0+dt2
 call define_ace

 do i=1,Np
 atomo(i)%k3v%x=atomo(i)%vel%x; atomo(i)%k3v%y=atomo(i)%vel%y; atomo(i)%k3v%z=atomo(i)%vel%z
 atomo(i)%k3a%x=atomo(i)%ace%x; atomo(i)%k3a%y=atomo(i)%ace%y; atomo(i)%k3a%z=atomo(i)%ace%z
 end do

!
 do i=1,Np
 atomo(i)%pos%x=atomo(i)%pos0%x+atomo(i)%k3v%x*dt
 atomo(i)%pos%y=atomo(i)%pos0%y+atomo(i)%k3v%y*dt
 atomo(i)%pos%z=atomo(i)%pos0%z+atomo(i)%k3v%z*dt
 
 atomo(i)%vel%x=atomo(i)%vel0%x+atomo(i)%k3a%x*dt
 atomo(i)%vel%y=atomo(i)%vel0%y+atomo(i)%k3a%y*dt
 atomo(i)%vel%z=atomo(i)%vel0%z+atomo(i)%k3a%z*dt
 end do

 tempo=tempo0+dt
 call define_ace

 do i=1,Np
 atomo(i)%k4v%x=atomo(i)%vel%x; atomo(i)%k4v%y=atomo(i)%vel%y; atomo(i)%k4v%z=atomo(i)%vel%z
 atomo(i)%k4a%x=atomo(i)%ace%x; atomo(i)%k4a%y=atomo(i)%ace%y; atomo(i)%k4a%z=atomo(i)%ace%z
 end do

!
 do i=1,Np
 atomo(i)%pos%x=atomo(i)%pos0%x+(atomo(i)%k1v%x+2*atomo(i)%k2v%x+2*atomo(i)%k3v%x+atomo(i)%k4v%x)*dt6
 atomo(i)%pos%y=atomo(i)%pos0%y+(atomo(i)%k1v%y+2*atomo(i)%k2v%y+2*atomo(i)%k3v%y+atomo(i)%k4v%y)*dt6
 atomo(i)%pos%z=atomo(i)%pos0%z+(atomo(i)%k1v%z+2*atomo(i)%k2v%z+2*atomo(i)%k3v%z+atomo(i)%k4v%z)*dt6
 end do

 do i=1,Np
 atomo(i)%vel%x=atomo(i)%vel0%x+(atomo(i)%k1a%x+2*atomo(i)%k2a%x+2*atomo(i)%k3a%x+atomo(i)%k4a%x)*dt6
 atomo(i)%vel%y=atomo(i)%vel0%y+(atomo(i)%k1a%y+2*atomo(i)%k2a%y+2*atomo(i)%k3a%y+atomo(i)%k4a%y)*dt6
 atomo(i)%vel%z=atomo(i)%vel0%z+(atomo(i)%k1a%z+2*atomo(i)%k2a%z+2*atomo(i)%k3a%z+atomo(i)%k4a%z)*dt6
 end do

 tempo=tempo0+dt
!

 call define_ace
 
  atomo%pkn3=atomo%pkn2 ; atomo%pkn2=atomo%pkn1 ; atomo%pkn1=atomo%pkn ; atomo%pkn=atomo%vel
  atomo%vkn3=atomo%vkn2 ; atomo%vkn2=atomo%vkn1 ; atomo%vkn1=atomo%vkn ; atomo%vkn=atomo%ace

 RETURN
END

!####################################################
! SUB-ROTINA CALCULA 
!####################################################
SUBROUTINE  PreditorCorretor()
 USE parametros
 USE tipos
 USE variables
 IMPLICIT NONE
 integer(kind=8) :: i
 real(kind=8):: tempo0
!
  tempo0=tempo
  call define_ace

   atomo%pos0=atomo%pos ; atomo%vel0=atomo%vel
   atomo%vkn=atomo%ace  ; atomo%pkn=atomo%vel 
!   
   atomo%pos%x=atomo%pos0%x+(55*atomo%pkn%x-59*atomo%pkn1%x+37*atomo%pkn2%x-9*atomo%pkn3%x)*dt24
   atomo%pos%y=atomo%pos0%y+(55*atomo%pkn%y-59*atomo%pkn1%y+37*atomo%pkn2%y-9*atomo%pkn3%y)*dt24
   atomo%pos%z=atomo%pos0%z+(55*atomo%pkn%z-59*atomo%pkn1%z+37*atomo%pkn2%z-9*atomo%pkn3%z)*dt24
!
   atomo%vel%x=atomo%vel0%x+(55*atomo%vkn%x-59*atomo%vkn1%x+37*atomo%vkn2%x-9*atomo%vkn3%x)*dt24
   atomo%vel%y=atomo%vel0%y+(55*atomo%vkn%y-59*atomo%vkn1%y+37*atomo%vkn2%y-9*atomo%vkn3%y)*dt24
   atomo%vel%z=atomo%vel0%z+(55*atomo%vkn%z-59*atomo%vkn1%z+37*atomo%vkn2%z-9*atomo%vkn3%z)*dt24
!
  tempo=tempo0+dt
  call define_ace
!
   atomo%pos%x=atomo%pos0%x+(9*atomo%vel%x+19*atomo%pkn%x-5*atomo%pkn1%x+atomo%pkn2%x)*dt24
   atomo%pos%y=atomo%pos0%y+(9*atomo%vel%y+19*atomo%pkn%y-5*atomo%pkn1%y+atomo%pkn2%y)*dt24
   atomo%pos%z=atomo%pos0%z+(9*atomo%vel%z+19*atomo%pkn%z-5*atomo%pkn1%z+atomo%pkn2%z)*dt24
!
   atomo%vel%x=atomo%vel0%x+(9*atomo%ace%x+19*atomo%vkn%x-5*atomo%vkn1%x+atomo%vkn2%x)*dt24
   atomo%vel%y=atomo%vel0%y+(9*atomo%ace%y+19*atomo%vkn%y-5*atomo%vkn1%y+atomo%vkn2%y)*dt24
   atomo%vel%z=atomo%vel0%z+(9*atomo%ace%z+19*atomo%vkn%z-5*atomo%vkn1%z+atomo%vkn2%z)*dt24
!
   atomo%pkn3=atomo%pkn2 ; atomo%pkn2=atomo%pkn1 ; atomo%pkn1=atomo%pkn
   atomo%vkn3=atomo%vkn2 ; atomo%vkn2=atomo%vkn1 ; atomo%vkn1=atomo%vkn
!
  tempo=tempo0+dt
!
 RETURN
end subroutine
