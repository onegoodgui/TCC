
    
module LLS


    
    implicit none
    
    !*****************************************************************
    !módulo para dimensionamento de cantoneiras de aço de abas iguais. 
    !A cantoneira de referência possui como eixo x o eixo central paralelo a sua aba horizontal, e como eixo y, o eixo central paralelo a sua aba vertical.
    !A aba horizontal está abaixo do eixo x e a aba vertical está a esquerda do eixo y.
    !*****************************************************************
    
    
    !*****************************************************************
    !declaração da variável tipo
    !*****************************************************************
    type LLS_var
    !*****************************************************************
    !    Entradas
    !*****************************************************************
    character(20) :: name = ""    !nome da seção 
    !type(LLS_var) :: secao        ! secao escolhida/encontrada atraves de fornecimento previo do nome e suas dimensoes caracteristicas
    real(8) :: h = 0.d0           ! Altura da aba da cantoneira  
    real(8) :: t_o = 0.d0         ! Espessura da aba
    real(8) :: A = 0.d0           ! Área da seção transversal da cantoneira
    real(8) :: x_g = 0.d0         ! distância entre a face da aba vertical e o eixo y
    real(8) :: y_g = 0.d0         ! distância entre a face da aba horizontal e o eixo x
    real(8) :: r_max = 0.d0       ! Raio de giracao maximo
    real(8) :: r_min = 0.d0       ! Raio de giracao minimo
    real(8) :: Ix = 0.d0          ! Momento de inércia em relação ao eixo x
    real(8) :: Iy = 0.d0          ! Momento de inércia em relação ao eixo y
    real(8) :: J = 0.d0           ! Momento de inércia a torção pura
    real(8) :: Cw = 0.d0          ! Constante de empenamento
    real(8) :: r_o= 0.d0          ! raio de giracao polar em relação ao centro de cisalhamento
    
    !******************************************************************
    !  Componentes geométricos da seção
    !******************************************************************
    
    
    real(8) :: Ixy = 0.d0         ! Produto de inércia em relação aos eixos x e y
    real(8) :: rx = 0.d0          ! Raio de giração em relação ao eixo x
    real(8) :: ry = 0.d0          ! Raio de giração em relação em eixo y
    real(8) :: Wx_s = 0.d0        ! Módulo resistente a flexão em relação ao eixo x superior
    real(8) :: Wx_i = 0.d0        ! Módulo resistente a flexão em relação ao eixo x inferior
    real(8) :: Imax = 0.d0        ! Momento de inércia central máximo
    real(8) :: Imin = 0.d0
    real(8) :: y_o = 0.d0         !distancia entre o centroide e o centro de cisalhamento da seção transversal da cantoneira (y_o)
    real(8) :: lambda = 0.d0       !  indice de esbeltez da aba da cantoneira

    
    
    
    end type LLS_var
    !************************************************************************

    
    
    !************************************************************************
    !Declaração de parâmetros
    !************************************************************************
    
    real(8), parameter, private:: pi = 4.d0*atan(1.d0)    !Número, se uma circunferência tem perímetro p e diametro d, pi = p/d
    real(8), parameter, private :: gama_a1 = 1.10         ! coeficiente de reduçao de resistencia para combinacao normal de açoes
    
    !************************************************************************
    
    contains
    
    !******************************************************************************************************
    !******************************************************************************************************
    !******************************************************************************************************
    !******************************************************************************************************
    
    !******************************************************************************************************
    ! subrotina que calcula propriedades geométricas úteis que não constam na tabela AISC
    !******************************************************************************************************
    subroutine LLS_propriedades_geometricas (secao)
    
              
    type(LLS_var),intent(inout) :: secao
  

    ! subrotina para cálculo do produto de inércia dos eixos x e y (Ixy)
    secao%Ixy = secao%h*secao%t_o*(secao%t_o/2 - secao%x_g)*(secao%h/2 - secao%y_g) + secao%t_o*(secao%h-secao%t_o)*(secao%t_o + secao%h/2 - secao%x_g)*(secao%t_o/2 - secao%y_g) 
    
    ! cálculo do momento de inercia maximo (Imax)
    secao%Imax = (secao%Ix + secao%Iy)/2.d0 + sqrt(((secao%Ix - secao%Iy)/2.d0)**2 + secao%Ixy**2)
    
    
    ! cálculo do raio de giracao maximo (r_max)
    secao%r_max = sqrt(secao%Imax/secao%A)
    
    
     ! cálculo do módulo resistente a flexão superior em relação ao eixo x (Wx_s)
    secao%Wx_s = secao%Ix/(secao%h-secao%x_g)
    
    
     ! cálculo do módulo resistente a flexão inferior em relação ao eixo x (Wx_i)
    secao%Wx_i = secao%Ix/secao%x_g
         
    
      !cálculo da distancia entre o centroide e o centro de cisalhamento da seção transversal da cantoneira (y_o)
     secao%y_o = sqrt((secao%y_g - secao%t_o/2.d0)**2.d0 + (secao%x_g - secao%t_o/2.d0)**2.d0)
     
     
     ! cálculo da esbeltez da aba (lamba) 
     secao%lambda = secao%h/secao%t_o
     
     
 
    end subroutine
    
     !******************************************************************************************************
     !******************************************************************************************************
    
    
    
     !******************************************************************************************************
     ! cálculo do fator de reducao para flambagem local - Qs
     !******************************************************************************************************
     
     
     subroutine LLS_Qs(secao, E, fy, Qs) ! cálculo do fator de reducao para flambagem local
     
     type(LLS_var), intent(inout) :: secao
     real (8) :: E                                  ! Módulo de elasticidade do aço
     real (8) :: fy                                 ! tensão de escoamento do aço (conforno tabela A1.4 da NBR8800)
     real (8), intent(out) :: Qs                    ! coeficiente de reducao aplicavel a seçoes em que os componentes tem indice de esbeltez superior aos valores da tabela F.1 da NBR8800
    real(8) :: lambda_r                             ! esbeltez limite (tabela F1 da NBR8800)
    real(8) :: lambda_e                             ! esbeltez limite elastica (tabela F1 da NBR8800)
    
    lambda_r = 0.45d0*sqrt(E/fy)                    ! esbeltez limite 'r' da aba da cantoneira. abaixo do qual não ocorre flambagem local. Conforme livro Pfeil - Estruturas de Aco, Dimensionamento Prático; pag 133
    lambda_e = 0.91d0*sqrt(E/fy)                    ! esbeltez limite 'e' da aba da cantoneira. acima do qual ocorre flambagem local inelástica. Conforme livro Pfeil - Estruturas de Aco, Dimensionamento Prático; pag 133

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
     ! calculo da força de flambagem crítica - Ncr
     !******************************************************************************************************
     
     
    subroutine LLS_Ncr(secao, G, E, KxLx, KyLy, KzLz, Ncr) ! calculo da força de flambagem crítica
    
    type(LLS_var),intent(in) :: secao
    real(8), intent(in) :: G          ! modulo de elasticidade transversal do aço    
    real(8),intent(in):: E            ! módulo de elasticidade do aço 
    real(8),intent(in):: KxLx         ! comprimento de flambagem em relação ao eixo principal central x
    real(8),intent(in):: KyLy         ! comprimento de flambagem em relação ao eixo principal central y 
    real(8),intent(in):: KzLz         ! comprimento de flambagem em relação ao eixo z
    real(8) :: Nex                    ! força axial de flambagem elastica por flexao em relação ao eixo central principal x
    real(8) :: Ney                    ! força axial de flambagem elastica por flexao em relação ao eixo central principal y
    real(8) :: Nez                    ! força axial de flambagem elastica por torção em relação ao eixo z
    real(8) :: Neyz                   ! forção de flambagem elástica por flexo-torção
    real(8),intent(out):: Ncr         ! força de flambagem critica
 
    
    

    
     Nex = ((pi**2)*E*secao%A*secao%r_min**2)/(KxLx**2)
     
     Ney = ((pi**2)*E*secao%A*secao%r_max**2)/(KyLy**2)
    
     Nez = (1.d0/secao%r_o**2)*(G*secao%J + (pi**2)*(E*secao%Cw)/(KzLz**2))
     
     Neyz = ((Ney + Nez)/(2.d0*(1.d0 - (secao%y_o/secao%r_o)**2)))*(1.d0 - sqrt(1.d0 - (4.d0*Ney*Nez*(1.d0 - (secao%y_o/secao%r_o)**2))/(Ney + Nez)**2))
     
     Ncr = min(Nex,Neyz)
     
        
     
    end subroutine
    
    !******************************************************************************************************
    !******************************************************************************************************
    
    
    
    
     !******************************************************************************************************
     ! calculo da força de flambagem crítica - Ncr - CASOS ESPECIAIS - CONECTADO POR UMA ABA - ver pag 123
     !******************************************************************************************************
     
     subroutine LLS_Ncr_exc(secao, G, E, KxLx, Ncr) ! calculo da força de flambagem crítica
     
    type(LLS_var),intent(in) :: secao
    real(8),intent(inout):: KxLx         ! comprimento de flambagem em relação ao eixo principal central x
    real(8), intent(in) :: G             ! modulo de elasticidade transversal do aço    
    real(8),intent(in):: E               ! módulo de elasticidade do aço 
    real(8) :: Nex                       ! força axial de flambagem elastica por flexao em relação ao eixo que passa pelo centro geométrico e é paralelo à aba conectada
    real(8),intent(out):: Ncr            ! força de flambagem critica
    
    
    
    if(KxLx/secao%rx <= 80) then
    
    KxLx = 72.d0*secao%rx + 0.75d0*KxLx
    
    else
        
    KxLx = 32.d0*secao%rx + 1.25d0*KxLx    
    
    end if 
    
    Nex = ((pi**2)*E*secao%A*secao%rx**2)/(KxLx**2)
    
    Ncr = Nex
    
    
    end subroutine
    
    !******************************************************************************************************
    ! calculo da força axial de compressão resistente de cálculo - Ncrd
    !******************************************************************************************************
    
    subroutine LLS_NcRd (secao, G, E, fy, KxLx, KyLy, KzLz, NcRd)
    
    real (8) :: fy                         ! tensão de escoamento do aço (conforno tabela A1.4 da NBR8800)
    real (8) :: Qs                         ! coeficiente de reducao aplicavel a seçoes em que os componentes tem indice de esbeltez superior aos valores da tabela F.1 da NBR8800
    real(8) :: lambda_r                    ! esbeltez limite (tabela F1 da NBR8800)
    real(8) :: lambda_e                    ! esbeltez limite elastica (tabela F1 da NBR8800)
    real(8) :: G                           ! modulo de elasticidade transversal do aço    
    real(8):: E                            ! Módulo de elasticidade do aço
    real(8):: KxLx                         ! comprimento de flambagem em relação ao eixo principal central x
    real(8):: KyLy                         ! comprimento de flambagem em relação ao eixo principal central y
    real(8):: KzLz                         ! comprimento de flambagem em relação ao eixo z
    real(8):: Nex                          ! força axial de flambagem elastica por flexao em relação ao eixo central principal x
    real(8):: Ney                          ! força axial de flambagem elastica por flexao em relação ao eixo central principal y
    real(8):: Nez                          ! força axial de flambagem elastica por torção em relação ao eixo z
    real(8):: Neyz                         ! forção de flambagem elástica por flexo-torção
    real(8):: Ncr                          ! força de flambagem critica
    real(8):: lambda_0                     ! indice de esbeltez reduzido
    real(8):: qui                            ! parametro adimensional que relaciona tensão crítica e tensão de escoamento do aço
    real(8), intent(out) :: NcRd           ! força axial de compressão resistente de cálculo
    

    
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
    ! calculo da força axial de compressão resistente de cálculo - Ncrd - CONECTADO POR UMA ABA - NBR 8800/08 - ver pag 123
    !******************************************************************************************************
    
    subroutine LLS_NcRd_exc (secao, G, E, fy, KxLx, NcRd)
    
    real (8) :: fy                         ! tensão de escoamento do aço (conforno tabela A1.4 da NBR8800)
    real (8) :: Qs                         ! coeficiente de reducao aplicavel a seçoes em que os componentes tem indice de esbeltez superior aos valores da tabela F.1 da NBR8800
    real(8) :: lambda_r                    ! esbeltez limite (tabela F1 da NBR8800)
    real(8) :: lambda_e                    ! esbeltez limite elastica (tabela F1 da NBR8800)
    real(8) :: G                           ! modulo de elasticidade transversal do aço    
    real(8):: E                            ! Módulo de elasticidade do aço
    real(8):: KxLx                         ! comprimento de flambagem em relação ao eixo principal central x
    real(8):: Nex                          ! força axial de flambagem elastica por flexao em relação ao eixo que passa pelo centro geométrico e é paralelo à aba conectada
    real(8):: Ncr                          ! força de flambagem critica
    real(8):: lambda_0                     ! indice de esbeltez reduzido
    real(8):: X                            ! parametro adimensional que relaciona tensão crítica e tensão de escoamento do aço
    real(8), intent(out) :: NcRd           ! força axial de compressão resistente de cálculo
    

    
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
integer :: i                            ! indice que contem o numero de seçoes que antecedem a posiçao da primeira seção cantoneira na tabela AISC
integer, intent(out) :: n_cant          ! indice que contem o numero de cantoneiras de abas iguais na tabela AISC
character :: nome
character(20), dimension (74) :: dados  ! string que contem os as propriedades geometricas das seções cantoneiras da tabela AISC
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
cant(i)%Ix = cant(i)%Ix*1000000  ! na tabela AISC, o valor referente ao momento de inercia em relaçao ai eixo x é Ix/10^6 - unidade [mm]
!
read(dados(41), *) cant(i)%rx
read(dados(42), *) cant(i)%Iy
!
cant(i)%Iy = cant(i)%Iy*1000000  ! na tabela AISC, o valor referente ao momento de inercia em relaçao ai eixo y é Ix/10^6 - unidade [mm]
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
    ! subrotina para encontrar uma secao cantoneira especifica a partir de suas dimensões características
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
    
