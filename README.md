Our simulation was divided in two main codes. The first named,
”simulaCellModelAngPot14FEV2020.f90” and the second “simulaCellModelAngPot22ABRIL2020.f90”.
In the first code we perform an energy minimization given an initial configuration of a system containing 5x5 cells. We solve the classical Newton of coupled equations of movement including a viscous friction that removes energy leading the system to a local minimum of energy. The equations of movement are integrated using a fourth-order runge-kutta method together with a fourth-order predictor corrector method . The configurations of the system are stored in a file with the extension *.xyz, where each position and velocities are saved in regular time steps. After a long time, the energy reaches a local minimum and the simulation is finished.

After that, the second code is started reading the last configurations and these are initial configurations for a simulation where a AFM tip modeled by a spring connected to a base that can be pushed with a constant velocity.  

Our first code named “simulaCellModelAngPot14FEV2020.f90” is composed for some modules and subroutines. These codes are organized as described above.

We create some comum modules:

MODULE tipos # create a structured variable

MODULE parametros # define fixed parameters 

MODULE variables # define global variables

The principal program that share the variables created in module:

PROGRAM simulaCellModelAngPot # the main program

And many subroutines that are called in principal program:

SUBROUTINE define_GlobalParameters() # define global parameters

SUBROUTINE define_PotParameters() # define potential parameters

SUBROUTINE Initialize_system() # define initial system

SUBROUTINE posicoesiniciais_cell() # define initial particles positions

SUBROUTINE define_massa() # define cell particles mass

SUBROUTINE define_tipo() # define type of cells particles

SUBROUTINE define_tabelavizinhos1() # define fixed table of neighbours

SUBROUTINE define_tabelavizinhos() # define variable table of neighbours

SUBROUTINE define_costhetaijk0() # define three-body potentials angle parameters 

SUBROUTINE define_vel_inicial() # define initial particles velocities

SUBROUTINE saidaxmakemol() # define write a xyz file

SUBROUTINE saidaParaView(temp_j) # define write a csv file

SUBROUTINE calculaenercin() # calculate total kinetic energy

SUBROUTINE renormavel() # calculate velocity renormalization factor

SUBROUTINE imprimeresultados() # define write data files

SUBROUTINE define_ace() # calculate forces and accelerations

SUBROUTINE evolucao_RK4() # calculate a fourth-order runge-kutta step 

SUBROUTINE PreditorCorretor() # calculate a fourth-order Preditor Corretor step

SUBROUTINE writenumber(i2,number1) # write a string integer number with 6 digits


The first main program is as follow:

PROGRAM simulaCellModelAngPot!14 fev 2020

USE parametros

USE tipos

USE variables

integer(kind=8) :: i,j,i1,itmp

call define_GlobalParameters()

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

CALL saidaParaView(itmp)

call imprimeresultados()

end if

end do

deallocate(atomo)

end do

end do

end do

END PROGRAM simulaCellModelAngPot !14 fev 2020



Our second code named “simulaCellModelAngPot22ABRIL2020.f90” is composed for some modules and subroutines. These codes are organized as described above.

We create some comum modules:

MODULE tipos # create a structured variable

MODULE parametros # define fixed parameters

MODULE variables # define global variables

The principal program that share the variables created in module:

PROGRAM simulaCellModelAngPot # the main program

And many subroutines that are called in principal program:

SUBROUTINE define_GlobalParameters() # define global parameters

SUBROUTINE define_PotParameters() # define potential parameters

SUBROUTINE Initialize_system() # define initial system

SUBROUTINE posicoesiniciais_cell() # define initial particles positions

SUBROUTINE define_atomo_puxa()                 # define particle atached to AFM

SUBROUTINE define_massa() # define cell particles mass

SUBROUTINE define_tipo() # define type of cells particles

SUBROUTINE define_tabelavizinhos1() # define fixed table of neighbours

SUBROUTINE define_tabelavizinhos() # define variable table of neighbours

SUBROUTINE define_costhetaijk0() # define three-body potentials angle parameters 

SUBROUTINE define_vel_inicial() # define initial particles velocities

SUBROUTINE saidaxmakemol() # define write a xyz file

SUBROUTINE leitura_saidaxmakemol()            # read the xyz file to obtain initial configuration

SUBROUTINE saidaParaView(temp_j) # define write a csv file

SUBROUTINE calculaenercin() # calculate total kinetic energy

SUBROUTINE renormavel() # calculate velocity renormalization factor

SUBROUTINE imprimeresultados() # define write data files

SUBROUTINE define_ace() # calculate forces and accelerations

SUBROUTINE evolucao_RK4() # calculate a fourth-order runge-kutta step 

SUBROUTINE PreditorCorretor() # calculate a fourth-order Preditor Corretor step

SUBROUTINE writenumber(i2,number1) # write a string integer number with 6 digits

The main program is as follow:

program simulaCellModelAngPot!22 Abril 2020

USE parametros

USE tipos

USE variables

integer(kind=8) :: i,j,i1,itmp

call define_GlobalParameters()

do isig=1,1

do jepi=0,0

do lmas=2,10,1

call define_PotParameters()

path='inputData/'

leitura_nomexmake=TRIM(path)//'puxa_saidaxmakemol_isig99iepi99lmas99.xyz'

write(*,*) TRIM(leitura_nomexmake)

nomexmake='puxa_saidaxmakemol_isig'//TRIM('99')//'iepi'//TRIM('99')//'lmas'//TRIM('99')//'.xyz'  

nomeenergia= 'puxa_energia_tempo_isig'//TRIM('99')//'iepi'//TRIM('99')//'lmas'//TRIM('99')//'.dat'

tempo=0.0D0

call leitura_saidaxmakemol()

tempo_f=int(tempo_max/dt)

call writenumber(isig,num1)

call writenumber(jepi,num2)

call writenumber(lmas,num3)

nomexmake='puxa_saidaxmakemol_isig'//TRIM(num1)//'iepi'//TRIM(num2)//'lmas'//TRIM(num3)//'.xyz'  

nomeenergia= 'puxa_energia_tempo_isig'//TRIM(num1)//'iepi'//TRIM(num2)//'lmas'//TRIM(num3)//'.dat'

itmp=0

CALL saidaxmakemol()

print*,'call saidaxmakemol  - ok'

CALL saidaParaView(itmp)

print *,'call saidaParaView(0)  - ok'

call define_ace()

print*,'call define_ace() - ok'

call imprimeresultados()

print*,'call imprimeresultados() - ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

tempo=dt

print*,"tempo=",tempo

print*,"tpuxar=",tpuxar        

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

CALL saidaParaView(itmp)

call imprimeresultados()

end if

end do

deallocate(atomo)

end do

end do

end do

end program simulaCellModelAngPot!22 Abril 2020

Ao executar o script usando os seguinte comandos:

a) chmod +x CompExec1.sh  # Para transformar o script executável.

b) ./CompExec1.sh    # Execute o script .

Que ira compilar o programa usando configurações padrões e executar o código compilado simulaCellModelAngPot22ABRIL2020.x.

A simulação ira iniciar e criar os seguintes arquivos:

Cell_fil_000000_.csv

Cell_sur_000000_.csv

Cell_fil_000001_.csv

Cell_sur_000001_.csv 

Cell_fil_000002_.csv

Cell_sur_000002_.csv

Cell_fil_000003_.csv

Cell_sur_000003_.csv

Cell_fil_000004_.csv

Cell_sur_000004_.csv

saidaxmakemol_isig_000001_iepi_000000_lmas_000001_.xyz

tabeladevizinhos_spring.dat

define_costhetaijk0.dat  

Define_massa_atomoiculas.dat

energia_tempo_isig_000001_iepi_000000_lmas_000001_.dat 

tabeladevizinhos_angular.dat

tabeladevizinhos_LJ.dat

tipos.mod

variables.mod

parametros.mod





