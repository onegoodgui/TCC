
    
module LLS


    
    implicit none
    
    !*****************************************************************
    !m�dulo para dimensionamento de cantoneiras de a�o de abas iguais. 
    !A cantoneira de refer�ncia possui como eixo x o eixo central paralelo a sua aba horizontal, e como eixo y, o eixo central paralelo a sua aba vertical.
    !A aba horizontal est� abaixo do eixo x e a aba vertical est� a esquerda do eixo y.
    !*****************************************************************
    
    
    !*****************************************************************
    !declara��o da vari�vel tipo
    !*****************************************************************
    type LLS_var
    !*****************************************************************
    !    Entradas
    !*****************************************************************
    character(20) :: name = ""    !nome da se��o 
    !type(LLS_var) :: secao        ! secao escolhida/encontrada atraves de fornecimento previo do nome e suas dimensoes caracteristicas
    real(8) :: h = 0.d0           ! Altura da aba da cantoneira  
    real(8) :: t_o = 0.d0         ! Espessura da aba
    real(8) :: A = 0.d0           ! �rea da se��o transversal da cantoneira
    real(8) :: x_g = 0.d0         ! dist�ncia entre a face da aba vertical e o eixo y
    real(8) :: y_g = 0.d0         ! dist�ncia entre a face da aba horizontal e o eixo x
    real(8) :: r_max = 0.d0       ! Raio de giracao maximo
    real(8) :: r_min = 0.d0       ! Raio de giracao minimo
    real(8) :: Ix = 0.d0          ! Momento de in�rcia em rela��o ao eixo x
    real(8) :: Iy = 0.d0          ! Momento de in�rcia em rela��o ao eixo y
    real(8) :: J = 0.d0           ! Momento de in�rcia a tor��o pura
    real(8) :: Cw = 0.d0          ! Constante de empenamento
    real(8) :: r_o= 0.d0          ! raio de giracao polar em rela��o ao centro de cisalhamento
    
    !******************************************************************
    !  Componentes geom�tricos da se��o
    !******************************************************************
    
    
    real(8) :: Ixy = 0.d0         ! Produto de in�rcia em rela��o aos eixos x e y
    real(8) :: rx = 0.d0          ! Raio de gira��o em rela��o ao eixo x
    real(8) :: ry = 0.d0          ! Raio de gira��o em rela��o em eixo y
    real(8) :: Wx_s = 0.d0        ! M�dulo resistente a flex�o em rela��o ao eixo x superior
    real(8) :: Wx_i = 0.d0        ! M�dulo resistente a flex�o em rela��o ao eixo x inferior
    real(8) :: Imax = 0.d0        ! Momento de in�rcia central m�ximo
    real(8) :: Imin = 0.d0
    real(8) :: y_o = 0.d0         !distancia entre o centroide e o centro de cisalhamento da se��o transversal da cantoneira (y_o)
    real(8) :: lambda = 0.d0       !  indice de esbeltez da aba da cantoneira

    
    
    
    end type LLS_var
    !************************************************************************

    
    
    !************************************************************************
    !Declara��o de par�metros
    !************************************************************************
    
    real(8), parameter, private:: pi = 4.d0*atan(1.d0)    !N�mero, se uma circunfer�ncia tem per�metro p e diametro d, pi = p/d
    real(8), parameter, private :: gama_a1 = 1.10         ! coeficiente de redu�ao de resistencia para combinacao normal de a�oes
    
    !************************************************************************
    
    contains
    
    !******************************************************************************************************
    !******************************************************************************************************
    !******************************************************************************************************
    !******************************************************************************************************
    
    !******************************************************************************************************
    ! subrotina que calcula propriedades geom�tricas �teis que n�o constam na tabela AISC
    !******************************************************************************************************
    subroutine LLS_propriedades_geometricas (secao)
    
              
    type(LLS_var),intent(inout) :: secao
  

    ! subrotina para c�lculo do produto de in�rcia dos eixos x e y (Ixy)
    secao%Ixy = secao%h*secao%t_o*(secao%t_o/2 - secao%x_g)*(secao%h/2 - secao%y_g) + secao%t_o*(secao%h-secao%t_o)*(secao%t_o + secao%h/2 - secao%x_g)*(secao%t_o/2 - secao%y_g) 
    
    ! c�lculo do momento de inercia maximo (Imax)
    secao%Imax = (secao%Ix + secao%Iy)/2.d0 + sqrt(((secao%Ix - secao%Iy)/2.d0)**2 + secao%Ixy**2)
    
    
    ! c�lculo do raio de giracao maximo (r_max)
    secao%r_max = sqrt(secao%Imax/secao%A)
    
    
     ! c�lculo do m�dulo resistente a flex�o superior em rela��o ao eixo x (Wx_s)
    secao%Wx_s = secao%Ix/(secao%h-secao%x_g)
    
    
     ! c�lculo do m�dulo resistente a flex�o inferior em rela��o ao eixo x (Wx_i)
    secao%Wx_i = secao%Ix/secao%x_g
         
    
      !c�lculo da distancia entre o centroide e o centro de cisalhamento da se��o transversal da cantoneira (y_o)
     secao%y_o = sqrt((secao%y_g - secao%t_o/2.d0)**2.d0 + (secao%x_g - secao%t_o/2.d0)**2.d0)
     
     
     ! c�lculo da esbeltez da aba (lamba) 
     secao%lambda = secao%h/secao%t_o
     
     
 
    end subroutine
    
     !******************************************************************************************************
     !******************************************************************************************************
    
    
    
     !******************************************************************************************************
     ! c�lculo do fator de reducao para flambagem local - Qs
     !******************************************************************************************************
     
     
     subroutine LLS_Qs(secao, E, fy, Qs) ! c�lculo do fator de reducao para flambagem local
     
     type(LLS_var), intent(inout) :: secao
     real (8) :: E                                  ! M�dulo de elasticidade do a�o
     real (8) :: fy                                 ! tens�o de escoamento do a�o (conforno tabela A1.4 da NBR8800)
     real (8), intent(out) :: Qs                    ! coeficiente de reducao aplicavel a se�oes em que os componentes tem indice de esbeltez superior aos valores da tabela F.1 da NBR8800
    real(8) :: lambda_r                             ! esbeltez limite (tabela F1 da NBR8800)
    real(8) :: lambda_e                             ! esbeltez limite elastica (tabela F1 da NBR8800)
    
    lambda_r = 0.45d0*sqrt(E/fy)                    ! esbeltez limite 'r' da aba da cantoneira. abaixo do qual n�o ocorre flambagem local. Conforme livro Pfeil - Estruturas de Aco, Dimensionamento Pr�tico; pag 133
    lambda_e = 0.91d0*sqrt(E/fy)                    ! esbeltez limite 'e' da aba da cantoneira. acima do qual ocorre flambagem local inel�stica. Conforme livro Pfeil - Estruturas de Aco, Dimensionamento Pr�tico; pag 133

     if( secao%lambda <= lambda_r ) then 
         

         Qs = 1.0d0 ! flambagem local inexistente. apenas flambagem global
         
         return
         
     end if
     
     if (secao%lambda <= lambda_e) then
     
     Qs = 1.340d0 - 0.76d0*(secao%h/secao%t_o)*sqrt(fy/E)  
     
     else
         
     Qs = (0.53d0*E)/(fy*((secao%h/secao%t_o)**2))
     
     end if
     
     end subroutine
     
     !******************************************************************************************************
     !******************************************************************************************************
     
     
     
     
     !******************************************************************************************************
     ! calculo da for�a de flambagem cr�tica - Ncr
     !******************************************************************************************************
     
     
    subroutine LLS_Ncr(secao, G, E, KxLx, KyLy, KzLz, Ncr) ! calculo da for�a de flambagem cr�tica
    
    type(LLS_var),intent(in) :: secao
    real(8), intent(in) :: G          ! modulo de elasticidade transversal do a�o    
    real(8),intent(in):: E            ! m�dulo de elasticidade do a�o 
    real(8),intent(in):: KxLx         ! comprimento de flambagem em rela��o ao eixo principal central x
    real(8),intent(in):: KyLy         ! comprimento de flambagem em rela��o ao eixo principal central y 
    real(8),intent(in):: KzLz         ! comprimento de flambagem em rela��o ao eixo z
    real(8) :: Nex                    ! for�a axial de flambagem elastica por flexao em rela��o ao eixo central principal x
    real(8) :: Ney                    ! for�a axial de flambagem elastica por flexao em rela��o ao eixo central principal y
    real(8) :: Nez                    ! for�a axial de flambagem elastica por tor��o em rela��o ao eixo z
    real(8) :: Neyz                   ! for��o de flambagem el�stica por flexo-tor��o
    real(8),intent(out):: Ncr         ! for�a de flambagem critica
 
    
    

    
     Nex = ((pi**2)*E*secao%A*secao%r_min**2)/(KxLx**2)
     
     Ney = ((pi**2)*E*secao%A*secao%r_max**2)/(KyLy**2)
    
     Nez = (1.d0/secao%r_o**2)*(G*secao%J + (pi**2)*(E*secao%Cw)/(KzLz**2))
     
     Neyz = ((Ney + Nez)/(2.d0*(1.d0 - (secao%y_o/secao%r_o)**2)))*(1.d0 - sqrt(1.d0 - (4.d0*Ney*Nez*(1.d0 - (secao%y_o/secao%r_o)**2))/(Ney + Nez)**2))
     
     Ncr = min(Nex,Neyz)
     
        
     
    end subroutine
    
    !******************************************************************************************************
    !******************************************************************************************************
    
    
    
    
     !******************************************************************************************************
     ! calculo da for�a de flambagem cr�tica - Ncr - CASOS ESPECIAIS - CONECTADO POR UMA ABA - ver pag 123
     !******************************************************************************************************
     
     subroutine LLS_Ncr_exc(secao, G, E, KxLx, Ncr) ! calculo da for�a de flambagem cr�tica
     
    type(LLS_var),intent(in) :: secao
    real(8),intent(inout):: KxLx         ! comprimento de flambagem em rela��o ao eixo principal central x
    real(8), intent(in) :: G             ! modulo de elasticidade transversal do a�o    
    real(8),intent(in):: E               ! m�dulo de elasticidade do a�o 
    real(8) :: Nex                       ! for�a axial de flambagem elastica por flexao em rela��o ao eixo que passa pelo centro geom�trico e � paralelo � aba conectada
    real(8),intent(out):: Ncr            ! for�a de flambagem critica
    
    
    
    if(KxLx/secao%rx <= 80) then
    
    KxLx = 72.d0*secao%rx + 0.75d0*KxLx
    
    else
        
    KxLx = 32.d0*secao%rx + 1.25d0*KxLx    
    
    end if 
    
    Nex = ((pi**2)*E*secao%A*secao%rx**2)/(KxLx**2)
    
    Ncr = Nex
    
    
    end subroutine
    
    !******************************************************************************************************
    ! calculo da for�a axial de compress�o resistente de c�lculo - Ncrd
    !******************************************************************************************************
    
    subroutine LLS_NcRd (secao, G, E, fy, KxLx, KyLy, KzLz, NcRd)
    
    real (8) :: fy                         ! tens�o de escoamento do a�o (conforno tabela A1.4 da NBR8800)
    real (8) :: Qs                         ! coeficiente de reducao aplicavel a se�oes em que os componentes tem indice de esbeltez superior aos valores da tabela F.1 da NBR8800
    real(8) :: lambda_r                    ! esbeltez limite (tabela F1 da NBR8800)
    real(8) :: lambda_e                    ! esbeltez limite elastica (tabela F1 da NBR8800)
    real(8) :: G                           ! modulo de elasticidade transversal do a�o    
    real(8):: E                            ! M�dulo de elasticidade do a�o
    real(8):: KxLx                         ! comprimento de flambagem em rela��o ao eixo principal central x
    real(8):: KyLy                         ! comprimento de flambagem em rela��o ao eixo principal central y
    real(8):: KzLz                         ! comprimento de flambagem em rela��o ao eixo z
    real(8):: Nex                          ! for�a axial de flambagem elastica por flexao em rela��o ao eixo central principal x
    real(8):: Ney                          ! for�a axial de flambagem elastica por flexao em rela��o ao eixo central principal y
    real(8):: Nez                          ! for�a axial de flambagem elastica por tor��o em rela��o ao eixo z
    real(8):: Neyz                         ! for��o de flambagem el�stica por flexo-tor��o
    real(8):: Ncr                          ! for�a de flambagem critica
    real(8):: lambda_0                     ! indice de esbeltez reduzido
    real(8):: qui                            ! parametro adimensional que relaciona tens�o cr�tica e tens�o de escoamento do a�o
    real(8), intent(out) :: NcRd           ! for�a axial de compress�o resistente de c�lculo
    

    
    type(LLS_var), intent(inout) :: secao

    

    call LLS_propriedades_geometricas (secao)
    call LLS_Qs(secao, E, fy, Qs)
    call LLS_Ncr(secao, G, E, KxLx, KyLy, KzLz, Ncr)
    
    
    
    lambda_0 = sqrt((Qs*secao%A*fy)/Ncr)
    
    if(lambda_0 <= 1.50d0) then
        
        qui = 0.658d0**(lambda_0**2)
        
    else
        qui = 0.877d0/(lambda_0)**2
         
    end if
    
       NcRd = (Qs*secao%A*qui*fy)/gama_a1
     
    end subroutine
     
    !******************************************************************************************************
    !******************************************************************************************************   
    
    
    
    
    !******************************************************************************************************
    ! calculo da for�a axial de compress�o resistente de c�lculo - Ncrd - CONECTADO POR UMA ABA - NBR 8800/08 - ver pag 123
    !******************************************************************************************************
    
    subroutine LLS_NcRd_exc (secao, G, E, fy, KxLx, NcRd)
    
    real (8) :: fy                         ! tens�o de escoamento do a�o (conforno tabela A1.4 da NBR8800)
    real (8) :: Qs                         ! coeficiente de reducao aplicavel a se�oes em que os componentes tem indice de esbeltez superior aos valores da tabela F.1 da NBR8800
    real(8) :: lambda_r                    ! esbeltez limite (tabela F1 da NBR8800)
    real(8) :: lambda_e                    ! esbeltez limite elastica (tabela F1 da NBR8800)
    real(8) :: G                           ! modulo de elasticidade transversal do a�o    
    real(8):: E                            ! M�dulo de elasticidade do a�o
    real(8):: KxLx                         ! comprimento de flambagem em rela��o ao eixo principal central x
    real(8):: Nex                          ! for�a axial de flambagem elastica por flexao em rela��o ao eixo que passa pelo centro geom�trico e � paralelo � aba conectada
    real(8):: Ncr                          ! for�a de flambagem critica
    real(8):: lambda_0                     ! indice de esbeltez reduzido
    real(8):: X                            ! parametro adimensional que relaciona tens�o cr�tica e tens�o de escoamento do a�o
    real(8), intent(out) :: NcRd           ! for�a axial de compress�o resistente de c�lculo
    

    
    type(LLS_var), intent(inout) :: secao

    

    call LLS_propriedades_geometricas (secao)
    call LLS_Qs(secao, E, fy, Qs)
    call LLS_Ncr_exc(secao, G, E, KxLx, Ncr)
    
    
    
    lambda_0 = sqrt((Qs*secao%A*fy)/Ncr)
    
    if(lambda_0 <= 1.50d0) then
        
        X = 0.658d0**(lambda_0**2)
        
    else
        X = 0.877d0/(lambda_0)**2
         
    end if
    
       NcRd = (Qs*secao%A*X*fy)/gama_a1
     
    end subroutine
     
    !******************************************************************************************************
    !******************************************************************************************************   
    
    
    
    
    
    
    
    !******************************************************************************************************
    ! subrotina para armazenar as propriedades geometricas de todas as secoes cantoneiras da tabela AISC
    !******************************************************************************************************   

subroutine armazenar_propriedades_cantoneiras_abas_iguais(n_cant, cant)


integer,parameter:: uni_entrada=300
integer :: i                            ! indice que contem o numero de se�oes que antecedem a posi�ao da primeira se��o cantoneira na tabela AISC
integer, intent(out) :: n_cant          ! indice que contem o numero de cantoneiras de abas iguais na tabela AISC
character :: nome
character(20), dimension (74) :: dados  ! string que contem os as propriedades geometricas das se��es cantoneiras da tabela AISC
type(LLS_var), intent(out), allocatable :: cant(:)

rewind(uni_entrada)
nome = '0'
dados='0'
i=0
n_cant=0


!******************************************************************************************************   
! trecho do programa que encontra a posicao da primeira cantoneira na lista de perfis AISC
!******************************************************************************************************   
repete: do while (nome /= 'L')
    
  i=i+1
  
   read (uni_entrada,*) nome

 

end do repete
!******************************************************************************************************   
                                                                                                    

backspace (uni_entrada)
read(uni_entrada, *) dados
backspace (uni_entrada)

!***************************************************************************************************************************   
! trecho do programa que encontra o numero de cantoneiras de abas iguais dentre todas as cantoneiras da lista de perfis AISC  
!***************************************************************************************************************************   
do while (dados(1) == 'L')
    
read(uni_entrada, *) dados

if (dados(6) == dados(14)) then

n_cant=n_cant+1

  end if

    end do
!******************************************************************************************************   

i=0

rewind(uni_entrada)
read (uni_entrada,*) nome
read (uni_entrada,*) dados

backspace (uni_entrada)
backspace (uni_entrada)

repete2: do while (nome /= 'L')
      
  i=i+1
  
   
read (uni_entrada,*) nome
 

end do repete2

backspace (uni_entrada)
read (uni_entrada,*) dados
backspace (uni_entrada) 

!**********************************************************************************************************   
! trecho do programa que armazena as propriedades listadas na tabela AISC de cada cantoneira de abas iguais
!**********************************************************************************************************   
allocate(cant(n_cant))

i=0

do while (dados(1) == 'L')
      
read(uni_entrada, *) dados

if (dados(6) == dados(14)) then
    
i= i+1
read(dados(2), *) cant(i)%name
read(dados(5), *) cant(i)%A
read(dados(6), *) cant(i)%h
read(dados(21), *) cant(i)%t_o
read(dados(27), *) cant(i)%x_g
read(dados(28), *) cant(i)%y_g
read(dados(38), *) cant(i)%Ix
!
cant(i)%Ix = cant(i)%Ix*1000000  ! na tabela AISC, o valor referente ao momento de inercia em rela�ao ai eixo x � Ix/10^6 - unidade [mm]
!
read(dados(41), *) cant(i)%rx
read(dados(42), *) cant(i)%Iy
!
cant(i)%Iy = cant(i)%Iy*1000000  ! na tabela AISC, o valor referente ao momento de inercia em rela�ao ai eixo y � Ix/10^6 - unidade [mm]
!
!
read(dados(45), *) cant(i)%ry
read(dados(47), *) cant(i)%r_min
read(dados(49), *) cant(i)%J

cant(i)%J = cant(i)%J*1000

read(dados(50), *) cant(i)%Cw

cant(i)%Cw = cant(i)%Cw*1000000000

read(dados(58), *) cant(i)%r_o


end if

           end do
!******************************************************************************************************   

end subroutine


!******************************************************************************************************
!******************************************************************************************************



    !******************************************************************************************************
    ! subrotina para encontrar uma secao cantoneira especifica a partir de suas dimens�es caracter�sticas
    !******************************************************************************************************

type(LLS_var) function LLS_propriedades_cantoneira_lista(n_cant, cant, nome_cant)  
 
 integer, intent(in) :: n_cant
 type(LLS_var), intent(in) :: cant(n_cant)
 character(20), intent(in) :: nome_cant

 
 integer,parameter:: uni_entrada=300
 integer :: n
n=1
 
do while(cant(n)%NAME /= nome_cant)

n=n+1

end do

LLS_propriedades_cantoneira_lista = cant(n)

end function LLS_propriedades_cantoneira_lista



    

    end module LLS
    
