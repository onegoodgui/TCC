 module Analise_trelica_trapezoidal

    use Matriz_Rigidez
    use CargaVento
    use Estrutura_Trelica
    
    implicit none
    
    ! Variaveis ---------------------------------------------------------------------------------------------------------------------
    integer :: iii
    integer, parameter :: uni_dad_estrut = 150, uni_entrada = 300, uni_dad_tercas = 500
    !Variaveis LLS -----------------------------------------------------------------------------------------------------------------
    integer :: n_cant
    type(LLS_var), allocatable :: cant(:)
    !Variaveis barra ---------------------------------------------------------------------------------------------------------------
    integer :: n_barras
    !integer :: num_nos
    real(8) :: rho = 0.000078d0                     ! Massa espec�fica do a�o [KN/cm�]
    real(8) :: fy = 0.25d0                          ! tens�o de escoamento do a�o [KN/mm^2)
    !real(8) :: E = 200.d0                           ! M�dulo de elasticidade do a�o
    real(8) :: G = 77.d0                            ! modulo de elasticidade transversal do a�o [KN/mm^2]
    ! Variables S1-------------------------------------------------------------------------------------------------------------------
    character(4) :: caso                            ! a) Terreno plano ou fracamente acidentado; b) Taludes e morros; c) Vales profundos
    !integer :: theta                               ! Inclina��o m�dia do talude ou enconsta de morro
    real(8) :: z                                    ! Altura medida a partir da superf�cie do terreno no ponto considerado
    real(8) :: d                                    ! diferen�a de n�vel entre a base e o topo do talude ou morro
    namelist /Fator_Top_S1/ caso, theta, z, d       ! namelist com dados de entrada para avaliar o fator  topogr�fico S1
    ! Variables S2-------------------------------------------------------------------------------------------------------------------   
    character :: cat                                ! categoria referente ao fator S2
    character :: cls                                ! classe referente ao fator S2
    namelist /Fator_Top_S2/ cat, cls                ! namelist com dados de entrada para avaliar o fator  topogr�fico S2
    type(tabela_S2) :: tab_S2(5)                    ! vetor de armazenamento dos dados da tabela num�rica do fator S2
     ! Variables press�o din�mica----------------------------------------------------------------------------------------------------
    real(8) :: Vo                                   ! Velocidade caracter�stica do vento da regi�o, em m/s
    namelist /pressao_dinamica/ Vo                  ! namelist com os dados de entrada para avaliar a press�o din�mica
    ! Variables Coef. press�o--------------------------------------------------------------------------------------------------------
    real(8) :: h                                    ! comprimento da altura da coluna de sustenta��o + montante imediatamente acima / vista frontal
    real(8) :: b                                    ! menor dimens�o da estrutura em vista superior
    real(8) :: a                                    ! maior dimens�o da estrutura em vista superior
    integer :: ang_cobertura                        ! inclina��o da viga de cobertura treli�ada entre montante externo e intermedi�rio
    namelist /Cofc_pressao/ h, a, b, ang_cobertura  ! namelist com os valores de entrada para avaliar Ce e Cpe
    !Variables Press�o Vento---------------------------------------------------------------------------------------------------------
    type(Coef_pressao) :: pressao_vento(3)            ! Vari�vel de armazenamento dos coeficientes Ce e Cpe
    real(8) :: coef_min_succao(2), coef_max_succao(2) ! Ce + Ci
    real(8) :: F(2)                                   ! Carga distribuida do vento
    !Variables Dados Geom�tricos ----------------------------------------------------------------------------------------------------
    namelist /Dados_geometricos/ L, h1, dist_trelica, n_div    ! namelist com os valores de entrada referentes �s imposi��es geom�tricas da estrutura
    !Variables Combina��o de a��es --------------------------------------------------------------------------------------------------
    type(comb_acoes) :: comb_acao(5)
    !Variables Resist�ncia
    real(8), allocatable :: Nrd(:)                   ! For�a resistente de calculo das barras de treli�a [KN] 
    !Variables Otimiza��o
    integer :: caso_vento = 0                           !0 = press�o negativa m�xima; 1 = press�o positiva m�xima
    real(8) :: coeficiente_pp = 1.4d0                   ! coeficiente de pondera��o da carga decorrente do peso pr�prio          
    real(8) :: coeficiente_cobertura = 1.4d0            ! coeficiente de pondera��o da carga decorrente da cobertura          
    real(8) :: coeficiente_sobrecarga = 0.d0            ! coeficiente de pondera��o da carga decorrente da sobrecarga       
    real(8) :: coeficiente_vento = 1.0d0                ! coeficiente de pondera��o da carga decorrente do vento    
    real(8) :: peso_total
    real(8) :: dc
    
    integer, parameter :: uni_dad_vento = 100           ! nome dado para n�mero que identifica arquivo externo
    integer :: i                                        ! variavel de recurso c�clico
    integer :: n                                        ! vari�vel de recurso c�clico 
    integer :: junk1
    character :: junk2
    
    
    contains
    
    !****************************************************************************************************************************************
    subroutine analise_resistencia_axial(barra, fy, E, G, forca_axial, Nrd)
    !****************************************************************************************************************************************
    
    type(barra_trelica), intent(inout), allocatable :: barra(:)     ! Se��o das barras (aplicando simetria), conforme a numera��o e especifica��o da tabela da AISC
    real(8), intent(in) :: fy                                       ! Tens�o de escoamento do a�o [KN/mm^2]
    real(8), intent(in) :: E                                        ! Modulo de elasticidade do a�o [KN/mm^2]
    real(8), intent(in) :: G                                        ! Modulo de elasticidade transversal do a�o [KN/mm^2]
    real(8), intent(in), allocatable :: forca_axial(:)              ! Esfor�o axial nas barras de treli�a [KN]
    real(8) :: Ncrd(n_barras)                                       ! For�a de compress�o resistente de calculo das barras de treli�a [KN]
    real(8) :: Ntrd(n_barras)                                       ! For�a de tra��o resistente de calculo das barras de treli�a [KN]
    real(8), intent(out), allocatable :: Nrd(:)                   ! For�a resistente de calculo das barras de treli�a [KN] 

    
    ! Vari�veis internas ------------------
    integer :: i = 0
    integer :: ii = 0
    
    
    allocate(Nrd(n_barras))
    
 ! trecho do programa que calcula os esfor�os normais de compress�o e tra��o para a carga externa descendente
 !***********************************************************************************************************
    
        do n = 1, n_barras
            if(forca_axial(n)>=0) then
                Ntrd(n) =  barra(n)%s%A*fy/1.1d0
                Nrd(n) = Ntrd(n)
            else
                call LLSt_NcRd (barra(n)%s, E, G, fy, barra(n)%comprimento, barra(n)%comprimento, barra(n)%comprimento, NcRd(n))
                Nrd(n) = Ncrd(n)
            end if
        end do

 
    end subroutine
    
    !**********************************************************************************************************************************
    subroutine funcao_objetivo(nd, nc, Vd, Vc, fob)
    !**********************************************************************************************************************************
    
        integer, intent(in) :: nd                          ! numero de variaveis diescretas
        integer, intent(in) :: nc                          ! numero de variaveis continuas
        integer, intent(in) :: Vd(nd)                      ! valor das variaveis discretas
        real(8), intent(in) :: Vc(nc)                      ! valor das variaveis continuas
        real(8) :: fob                                     ! fun��o objetivo
        real(8) :: func_penal                           ! fun��o penalidade
        real(8) :: g1                                   ! inequa��o de restri��o quanto aos esfor�os axiais
        real(8) :: g1_aux(13)
        real(8) :: g2                                   ! inequa��o de restri��o quanto ao indice de esbeltez das barras
        real(8) :: g2_aux(13)
        real(8) :: g3                                   ! inequa��o de restri��o quanto ao deslocamento da estrutura
        real(8) :: g3_aux
        real(8) :: c = 1                                   ! coeficiente da fun��o de penaliza��o
        integer, allocatable :: i_banzo_superior(:)                  ! armazena os numeros que indicam barras que pertencem ao banzo superior
        integer, allocatable :: i_banzo_inferior(:)                  ! armazena os numeros que indicam barras que pertencem ao banzo inferior
        integer, allocatable :: i_diagonal(:)                        ! armazena os numeros que indicam barras que pertencem � diagonal
        integer, allocatable :: i_montante(:)                        ! armazena os numeros que indicam barras que pertencem ao montante
        character(20), allocatable :: nome_cant_trelica(:)
        
        
        !Corpo da subrotina----------------------------------------------------------------------------
        
        !Atribui coordenadas para os n�s da estrutura
            call node_barras(n_div, h1, num_nos, coord)
            
        !Estabelece o n�mero de barras, conectividades e seus comprimentos
            call barras_trelica (n_div, h1, coord, n_barras, barra)
            
        !Atribui uma se��o transversal para cada barra de acordo com a lista de se��es da AISC
            allocate(nome_cant_trelica(n_barras))
           
        if(h1>0) then
            
            allocate(i_banzo_superior(2*n_div))
            allocate(i_montante(2*n_div+1))
            allocate(i_diagonal(2*n_div))
            allocate(i_banzo_inferior(2*n_div))
            
            i_montante(1) = 1
            i_banzo_inferior(1)= 2
            i_diagonal(1) = 3
            i_banzo_superior(1) = 4
            
            do i = 1, 2*n_div-1
                
                i_montante(i+1) = i_montante(i) + 4
                i_banzo_inferior(i+1) = i_banzo_inferior(i) + 4
                i_diagonal(i+1) = i_diagonal(i) + 4
                i_banzo_superior(i+1) = i_banzo_superior(i) + 4
            
            end do
                i_montante(i+1) = i_montante(i) + 4
                
                nome_cant_trelica(1:n_barras) = ' '
                nome_cant_trelica(i_banzo_superior) = cant(20)%name
                nome_cant_trelica(i_banzo_inferior) = cant(20)%name
                nome_cant_trelica(i_diagonal) = cant(20)%name
                nome_cant_trelica(i_montante) = cant(20)%name
                
                do i = 1, n_barras
                    barra(i)%s%secao = LLS_propriedades_cantoneira_lista(n_cant, cant, nome_cant_trelica(i))
                    call LLS_propriedades_geometricas (barra(i)%s%secao)
                    call LLSt_propriedades_geometricas (barra(i)%s)
                end do
        else
                
                
        end if
    
        !Calcula a a magnitude de carga de peso pr�prio em cada n� da estrutura
            call carga_pp_barras (rho, h1, n_div, num_nos, n_barras, barra, cond_cont)

        !Calcula a carga das ter�as e telhas nos n�s superiores
            call carga_tercas_telhas_barras (h1, dist_trelica, n_div, num_nos, barra, terca, cond_cont)
    
        !Armazena os coeficientes Ce e Cpe da NBR 6123:2013
            call coeficiente_pressao (pressao_vento)
    
        !Combina os valores de Ce e Ci e guarda as combina��es com o valor m�x e min de succ�o
            call combinacao_carregamento_vento (h, a, b, ang_cobertura, pressao_vento, coef_max_succao, coef_min_succao)
    
        !Calcula a carga distribu�da na cobertura devido � a��o do vento // 0 = press�o negativa m�xima; 1 = press�o positiva m�xima
            select case(caso_vento)
                case(0)
                    call carga_distribuida_vento(Vo, S1, S2, coef_max_succao, barra, F)
                case(1)
                    call carga_distribuida_vento(Vo, S1, S2, coef_min_succao, barra, F)
            end select
        !Calcula as cargas nodais da a��o do vento na estrutura
            call carga_vento_sobrecarga_nos(theta, barra, F, cond_cont)
    
            allocate(v_glc(nglc))
            call vetor_vglc (cond_cont, v_glc)
        !Matriz de rigidez - vetor de deslocamentos fixos e livres / vetor de cargas
            allocate(Y(num_nos*ngln,1))
            do i=1, num_nos
                call vetor_Y (i, cond_cont(i), coeficiente_pp*cond_cont(i)%carga_pp, coeficiente_cobertura*cond_cont(i)%carga_telha_telhado, coeficiente_sobrecarga*cond_cont(i)%carga_sobrecarga, coeficiente_vento*cond_cont(i)%carga_vento, Y)
            end do
            
            allocate(MatRigid(ngln*num_nos,ngln*num_nos))
            MatRigid = 0.d0
    
            do i =1, n_barras
                call matriz_rigidez_local(barra(i)%s%A, E, barra(i)%comprimento, k_local)
                call matriz_trasformacao(barra(i)%node(1)%x,barra(i)%node(1)%y, barra(i)%node(2)%x,barra(i)%node(2)%y, barra(i)%comprimento, T)
                call matriz_rigidez_global(k_local, T, k)
                call matriz_rigidez_estrutura(barra(i)%conectividades(1), barra(i)%conectividades(2), k, MatRigid)
            end do 
    
            call matriz_rigidez_nos_restritos (nglc ,v_glc, MatRigid)
            call eliminacao_gauss(num_nos, MatRigid, Y)
            call deslocamentos_matriz_rigidez(nglc, v_glc, MatRigid, Y, deslocamentos)
            call forca_axial_barras(barra, n_barras, deslocamentos, forca_axial)
            call analise_resistencia_axial(barra, fy, E, G, forca_axial, Nrd)
            
            !-----------Peso total da estrutura!-------------------!
            peso_total = SUM(barra(:)%peso)
            
            !------------Deslocamento m�ximo da estrutura!-----------!
            dc = maxval(abs(deslocamentos(:,1)))
            
            
            !-------------trecho em que se calcula a penalidade da restri��o com rela��o aos esfor�os - g1-------------------!
 
            g1 = 1.d0
            
            do i=1, n_barras
        
                g1_aux(i) = ABS(forca_axial(i))/ABS(Nrd(i)) - 1.d0     
     
                if(g1_aux(i) > 0) then
                    g1 = g1 + g1_aux(i)
                end if
      
            end do

    
            !trecho em que se calcula a penalidade da restri��o com rela��o a esbeltez das barras - g2
!****************************************************************************************************    
            g2 = 1.d0
    
            do i=1, n_barras
                g2_aux(i) = barra(i)%comprimento/barra(i)%s%r_min - 200.d0
        
                if(g2_aux(i) > 0) then
                    g2 = g2 + g2_aux(i)
                end if
        
                g2_aux(i) = barra(i)%comprimento/barra(i)%s%r_min - 300.d0
                
                if(g2_aux(i) > 0) then
                    g2 = g2 + g2_aux(i)
                end if
        
            end do

    
            !trecho em que se calcula a penalidade da restri��o com rela��o aos deslocamentos da barra - g3
!****************************************************************************************************     
        
            g3 = 1
            g3_aux = ABS(dc)/(L/250.d0) - 1.d0
    
            if(g3_aux > 0) then
                g3 = g3 + g3_aux
            end if
    
    
     func_penal = 0
     
     if(g1 > 1.d0) func_penal = func_penal + c*g1**2
     
     if(g2 > 1.d0) func_penal = func_penal + c*g2**2
     
     if(g3 > 1.d0) func_penal = func_penal + c*g3**2
    
     fob = peso_total + func_penal
    end subroutine
    
    end module