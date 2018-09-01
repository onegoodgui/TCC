module CargaVento
    
    implicit none
    
   !*****************************************************************************************************!
   !******** M�dulo para determina��o da carga de vento conforme as prescri��es da NBR 6123 (2013) ******!
   !*****************************************************************************************************!
    
    ! Variables -------------------------------------------------------------------------------------------------------------------
    
    type classe
        character :: nome
        real(8) :: b
        real(8) :: p
        
    end type
    
    type tabela_S2
        character(3) :: categoria
        type(classe) :: classe(3)
    
    end type
    !------------------------------
    type Ce
        character(2) :: regiao
        real(8) :: valor
        
    end type
    
    type Cpe
        real(8) :: Cpe_medio
    end type
    
    type angulo
        integer :: theta
        type(Ce) :: Ce(4)
        type(Cpe) :: Cpe(4)
    end type
    
    
    type Coef_pressao
        
        real :: altura_relativa_max
        type(angulo), allocatable :: angulo_theta(:)
        
    end type
    
    type ang_vento
        
        real(8) :: Ce_esquerda   ! corte transversal que intercepta perpendicularmente as regi�es G e E
        real(8) :: Ce_direita   ! corte transversal que intercepta perpendicularmente as regi�es H e F
        
    end type
    
    type Combina_coef
        
        type(ang_vento) :: ang_90
        type(ang_vento) :: ang_0
        
    end type
    
        
    
    
    integer, parameter :: uni_dad_S2 = 101  
    integer, parameter :: uni_dad_pressao = 110 
    real(8), parameter, private:: pi = 4.d0*atan(1.d0)              ! Constante matem�tica pi
    real(8) :: S1                                                   ! Fator topogr�fico S1
    real(8) :: S2                                                   ! Fator topogr�fico S2
    ! -----------------------------------------------------------------------------------------------------------------------------
    contains
    
    ! Subrotinas ------------------------------------------------------------------------------------------------------------------
    
    subroutine Fator_Topografico_S1 (caso, theta, z, d, S1)
    
    character(4), intent(in) :: caso               ! a) Terreno plano ou fracamente acidentado; b) Taludes e morros; c) Vales profundos
    real(8), intent(in) :: theta                   ! Inclina��o m�dia do talude ou enconsta de morro
    real(8), intent(in) :: z                       ! Altura medida a partir da superf�cie do terreno no ponto considerado
    real(8), intent(in) :: d                       ! diferen�a de n�vel entre a base e o topo do talude ou morro
    real(8), intent(out) :: S1                     ! Fator topogr�fico S1

  
    ! ---- Obten��o do valor num�rico de S1 de acordo com os dados de entrada ---- !    
    if (caso == 'a') then
        S1 = 1.0
    else if (caso == 'b') then
        if (theta < 3) then
            S1 = 1.0
        else if (theta < 6) then
            S1 = 1 + ((2.5 - z/d)*tan(3.d0*pi/180))/3 * (theta - 3)
        else if (theta < 17) then
            S1 = 1 + (2.5 - z/d)*tan((theta - 3.d0)*pi/180)
        else if (theta < 45) then
            S1 = (1 + (2.5 - z/d)*tan(17*pi/180)) + (((2.5 - z/d)*(0.31 - tan(17*pi/180)))/28)*(theta - 17)
        else
            S1 = 1 + (2.5 - z/d)*0.31d0
        end if
    else if (caso == 'c') then
        S1 = 0.9
    end if
    
    end subroutine
    
    
    subroutine Transf_dados_S2 (tab_S2)
    
    type(tabela_S2), intent(out) :: tab_S2(5)       ! Tabela com os coeficientes S2
    integer :: junk1
    character :: junk2                              
    integer :: i, n                                 ! Contadores auxiliares
    
! ---- Identifica��o das classes do fator topogr�fico S2 ---- !    
    open (unit=uni_dad_S2 , file="Fator_S2.prn" , action='read')
        do i = 1, 5
            read(uni_dad_S2, *) tab_S2(i)%classe(:)%nome
            if (i == 5) then
                exit
            end if
            rewind(uni_dad_S2)
        end do
   ! ----- Transfer�ncia dos dados da tabela do Fator topogr�fico S2 para vari�vel tab_S2 ----- !
        do i = 1, 5
            read(uni_dad_S2, *) tab_S2(i)%categoria, junk1, junk2, tab_S2(i)%classe(:)%b
            read(uni_dad_S2, *) tab_S2(i)%categoria, junk1, junk2, tab_S2(i)%classe(:)%p
        end do
    close (uni_dad_S2)  

    end subroutine
    
    
    subroutine Fator_Topografico_S2 (cat, cls, z, tab_S2, S2)
    
    character, intent(in) :: cat                ! categoria escolhida de acordo com fator S2
    character, intent(in) :: cls                ! classe escolhida de acordo com fator S2
    real(8), intent(in) :: z                    ! Altura medida a partir da superf�cie do terreno no ponto considerado
    type(tabela_S2), intent(in) :: tab_S2(5)    ! tabela com os valores dos par�metros meteorol�gicos do fator S2
    real(8), intent(out) :: S2                  ! valor de S2
    
    real(8) :: b
    real(8) :: p
    real(8) :: Fr
    integer :: i = 1
    integer :: n = 1

! ---- Obten��o dos par�metros b e p ----
        do while (cat /= tab_S2(i)%categoria)
            i = i+1
        end do
    
        do while (cls /= tab_S2(i)%classe(n)%nome)
            n = n+1
        end do
    
        b = tab_S2(i)%classe(n)%b
        p = tab_S2(i)%classe(n)%p
 
! -- Obten��o do par�metro Fr --
    if (n .EQV. 'a') then
        Fr = 1.d0
    else if (n .EQV. 'b') then
        Fr = 0.98d0
    else
        Fr = 0.95d0
    end if
    
S2 = b*Fr*(z/10)**p

    end subroutine
    
    
    subroutine Coeficiente_pressao ( pressao_vento )
    
    type(Coef_pressao) :: pressao_vento(3)     ! Vari�vel de armazenagem dos coeficientes Ce e Cpe
    integer :: n_linhas(3)                     ! N�mero de linhas de cada caso presente na tabela 5 da NBR 6123:2013
    
    character(10) :: junk1, junk2              ! Auxiliar
    integer :: i, n                            ! Contadores
    
    ! atribui��o dos valores de altura relativa m�xima conforme tabela 5 da NBR 6123 : 2013
    pressao_vento(1)%altura_relativa_max = 1.d0/2.d0
    pressao_vento(2)%altura_relativa_max = 3.d0/2.d0
    pressao_vento(3)%altura_relativa_max = 6.d0
    
    ! defini��o do numero de linhas para cada caso
    n_linhas(1) = 8
    n_linhas(2) = 8
    n_linhas(3) = 9
    
    ! defini��o de tamanho do vetor que guarda os coeficientes
    
    allocate(pressao_vento(1)%angulo_theta(8))
    allocate(pressao_vento(2)%angulo_theta(8))
    allocate(pressao_vento(3)%angulo_theta(9))
    
    open (unit=uni_dad_pressao , file="dados_pressao.prn", action='read')
        
        do i = 1,3
            do n = 1, n_linhas(i)
                read(uni_dad_pressao, *) junk1, pressao_vento(i)%angulo_theta(n)%Ce(:)%regiao
                if (i == 3 .AND. n == 9) then           
                    exit
                end if
                rewind(uni_dad_pressao)
            end do
        end do
    
    
        do i = 1,3
            do n = 1, n_linhas(i)
                read(uni_dad_pressao, *) pressao_vento(i)%angulo_theta(n)%theta, pressao_vento(i)%angulo_theta(n)%Ce(:)%valor, pressao_vento(i)%angulo_theta(n)%Cpe(:)%Cpe_medio
                pressao_vento(i)%angulo_theta(n)%Ce(:)%valor = pressao_vento(i)%angulo_theta(n)%Ce(:)%valor/10
                pressao_vento(i)%angulo_theta(n)%Cpe(:)%Cpe_medio = pressao_vento(i)%angulo_theta(n)%Cpe(:)%Cpe_medio/10
            end do
        end do
        
        close (uni_dad_pressao)
        
    end subroutine
    
    
    subroutine Combinacao_carregamento_vento (h, a, b, ang_cobertura, pressao_vento, coef_max_succao, coef_min_succao)
    
    real(8), intent(in) :: h                        ! comprimento da altura da coluna de sustenta��o + montante imediatamente acima / vista frontal
    real(8), intent(in) :: a                        ! maior dimens�o da estrutura em vista superior
    real(8), intent(in) :: b                        ! menor dimens�o da estrutura em vista superior
    integer, intent(in) :: ang_cobertura            ! inclina��o da viga de cobertura treli�ada entre montante externo e intermedi�rio
    type(Coef_pressao) :: pressao_vento(3)          ! Vari�vel de armazenagem dos coeficientes Ce e Cpe
    real(8), intent(out) :: coef_max_succao(2)
    real(8), intent(out) :: coef_min_succao(2)
    
    !Variaveis internas ----------------------------------------------------
    real(8) :: alturarelativa                       !altura relativa  --> h/b
    real(8) :: Ci(2)
    integer :: i=0, n=0, j=0, k1, k2                  ! contadores
    integer :: n_linhas(3)                          ! N�mero de linhas de cada caso presente na tabela 5 da NBR 6123:2013
    real(8) :: min, max
    real(8), allocatable :: combinacao(:,:)        ! 
    
    ! defini��o do numero de linhas para cada caso
    n_linhas(1) = 8
    n_linhas(2) = 8
    n_linhas(3) = 9
    
    ! defini��o da altura relativa do caso estudado
    alturarelativa = h/b
    
    ! verifica��o do caso em que se encontra a altura relativa
    do i = 1, 3
        if(alturarelativa < pressao_vento(i)%altura_relativa_max) then
            exit
        end if
    end do
    
    ! verifica��o do angulo utilizado no estudo
    do n = 1, n_linhas(i)
        if (ang_cobertura == pressao_vento(i)%angulo_theta(n)%theta) then
            exit
        end if
    end do
    
    allocate(combinacao(6,2))
    
    ! Defini��o dos valores de Cpi (adotados conforme subitem 6.2.5, caso a) - NBR 6123/2013 
    Ci(1) = -0.2d0 ! (-) Suc��o
    Ci(2) = 0.3d0 ! (+) Sobrepress�o
    
    !Combinacao A) [Ce = 90� + Ci = 0,2] -------------------------
    combinacao(1,1) = pressao_vento(i)%angulo_theta(n)%Ce(2)%valor + Ci(1)
    combinacao(1,2) = pressao_vento(i)%angulo_theta(n)%Ce(1)%valor + Ci(1)
    
    !Combinacao B) [Ce = 90� + Ci = -0,3] -------------------------
    combinacao(2,1) = pressao_vento(i)%angulo_theta(n)%Ce(2)%valor + Ci(2)
    combinacao(2,2) = pressao_vento(i)%angulo_theta(n)%Ce(1)%valor + Ci(2)
    
    !Combinacao C) [Ce = 0� + Ci = 0,2] para GE -------------------------
    combinacao(3,1) = pressao_vento(i)%angulo_theta(n)%Ce(3)%valor + Ci(1)
    combinacao(3,2) = pressao_vento(i)%angulo_theta(n)%Ce(3)%valor + Ci(1)
    
    !Combinacao D) [Ce = 0� + Ci = -0,3] para GE -------------------------
    combinacao(4,1) = pressao_vento(i)%angulo_theta(n)%Ce(3)%valor + Ci(2)
    combinacao(4,2) = pressao_vento(i)%angulo_theta(n)%Ce(3)%valor + Ci(2)
    
    !Combinacao C) [Ce = 0� + Ci = 0,2] para FH -------------------------
    combinacao(5,1) = pressao_vento(i)%angulo_theta(n)%Ce(4)%valor + Ci(1)
    combinacao(5,2) = pressao_vento(i)%angulo_theta(n)%Ce(4)%valor + Ci(1)
    
    !Combinacao D) [Ce = 0� + Ci = -0,3] para FH -------------------------
    combinacao(6,1) = pressao_vento(i)%angulo_theta(n)%Ce(4)%valor + Ci(2)
    combinacao(6,2) = pressao_vento(i)%angulo_theta(n)%Ce(4)%valor + Ci(2)
    
    coef_min_succao(:) =   combinacao(1,1)
    coef_max_succao(:) =   combinacao(1,1)
    
    do n = 1, 2
        do j = 1, 6
            if(combinacao(j,n) >= coef_min_succao(1) .AND. combinacao(j,n) >= coef_min_succao(2)) then
                coef_min_succao(:) = combinacao(j,:)
            end if
            if (combinacao(j,n) <= coef_max_succao(1) .AND. combinacao(j,n) <= coef_max_succao(2)) then
                coef_max_succao(:) = combinacao(j,:)
            end if
        end do
    end do
    
  
    end subroutine
    
    !subroutine carga_distribuida_vento(
    
    
    
    
end module CargaVento
    