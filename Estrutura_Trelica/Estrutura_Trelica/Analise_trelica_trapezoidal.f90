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
    real(8) :: rho = 0.000078d0                     ! Massa específica do aço [KN/cm³]
    real(8) :: fy = 0.25d0                          ! tensão de escoamento do aço [KN/mm^2)
    !real(8) :: E = 200.d0                           ! Módulo de elasticidade do aço
    real(8) :: G = 77.d0                            ! modulo de elasticidade transversal do aço [KN/mm^2]
    ! Variables S1-------------------------------------------------------------------------------------------------------------------
    character(4) :: caso                            ! a) Terreno plano ou fracamente acidentado; b) Taludes e morros; c) Vales profundos
    !integer :: theta                               ! Inclinação média do talude ou enconsta de morro
    real(8) :: z                                    ! Altura medida a partir da superfície do terreno no ponto considerado
    real(8) :: d                                    ! diferença de nível entre a base e o topo do talude ou morro
    namelist /Fator_Top_S1/ caso, theta, z, d       ! namelist com dados de entrada para avaliar o fator  topográfico S1
    ! Variables S2-------------------------------------------------------------------------------------------------------------------   
    character :: cat                                ! categoria referente ao fator S2
    character :: cls                                ! classe referente ao fator S2
    namelist /Fator_Top_S2/ cat, cls                ! namelist com dados de entrada para avaliar o fator  topográfico S2
    type(tabela_S2) :: tab_S2(5)                    ! vetor de armazenamento dos dados da tabela numérica do fator S2
     ! Variables pressão dinâmica----------------------------------------------------------------------------------------------------
    real(8) :: Vo                                   ! Velocidade característica do vento da região, em m/s
    namelist /pressao_dinamica/ Vo                  ! namelist com os dados de entrada para avaliar a pressão dinâmica
    ! Variables Coef. pressão--------------------------------------------------------------------------------------------------------
    real(8) :: h                                    ! comprimento da altura da coluna de sustentação + montante imediatamente acima / vista frontal
    real(8) :: b                                    ! menor dimensão da estrutura em vista superior
    real(8) :: a                                    ! maior dimensão da estrutura em vista superior
    integer :: ang_cobertura                        ! inclinação da viga de cobertura treliçada entre montante externo e intermediário
    namelist /Cofc_pressao/ h, a, b, ang_cobertura  ! namelist com os valores de entrada para avaliar Ce e Cpe
    !Variables Pressão Vento---------------------------------------------------------------------------------------------------------
    type(Coef_pressao) :: pressao_vento(3)            ! Variável de armazenamento dos coeficientes Ce e Cpe
    real(8) :: coef_min_succao(2), coef_max_succao(2) ! Ce + Ci
    real(8) :: F(2)                                   ! Carga distribuida do vento
    !Variables Dados Geométricos ----------------------------------------------------------------------------------------------------
    namelist /Dados_geometricos/ L, h1, dist_trelica, n_div    ! namelist com os valores de entrada referentes às imposições geométricas da estrutura
    !Variables Combinação de ações --------------------------------------------------------------------------------------------------
    type(comb_acoes) :: comb_acao(5)
    !Variables Otimização
    
    
    integer, parameter :: uni_dad_vento = 100                     ! nome dado para número que identifica arquivo externo
    integer :: i                                                  ! variavel de recurso cíclico
    integer :: n                                                  ! variável de recurso cíclico 
    integer :: junk1
    character :: junk2
    
    contains
    
    !****************************************************************************************************************************************
    subroutine analise_resistencia_axial(barra, fy, E, G, forca_axial, Nrd)
    !****************************************************************************************************************************************
    
    type(barra_trelica), intent(inout), allocatable :: barra(:)     ! Seção das barras (aplicando simetria), conforme a numeração e especificação da tabela da AISC
    real(8), intent(in) :: fy                                       ! Tensão de escoamento do aço [KN/mm^2]
    real(8), intent(in) :: E                                        ! Modulo de elasticidade do aço [KN/mm^2]
    real(8), intent(in) :: G                                        ! Modulo de elasticidade transversal do aço [KN/mm^2]
    real(8), intent(in), allocatable :: forca_axial(:)              ! Esforço axial nas barras de treliça [KN]
    real(8) :: Ncrd(n_barras)                                       ! Força de compressão resistente de calculo das barras de treliça [KN]
    real(8) :: Ntrd(n_barras)                                       ! Força de tração resistente de calculo das barras de treliça [KN]
    real(8), intent(out), allocatable :: Nrd(:,:)                   ! Força resistente de calculo das barras de treliça [KN] 

    
    ! Variáveis internas ------------------
    integer :: i = 0
    integer :: ii = 0
    
    
    allocate(Nrd(6, n_barras))
    
 ! trecho do programa que calcula os esforços normais de compressão e tração para a carga externa descendente
 !***********************************************************************************************************
    do i = 1, 6
        do n = 1, n_barras
            if(forca_axial(i)>=0) then
                Ntrd(n) =  barra(n)%s%A*fy/1.1d0
                Nrd(i,n) = Ntrd(n)
            else
                call LLSt_NcRd (barra(n)%s, E, G, fy, barra(n)%comprimento, barra(n)%comprimento, barra(n)%comprimento, NcRd(n))
                Nrd(i,n) = Ncrd(n)
            end if
        end do
    end do
 
    end subroutine
    
    !**********************************************************************************************************************************
    subroutine funcao_objetivo(nd, nc, Vd, Vc, fob)
    !**********************************************************************************************************************************
        integer, intent(in) :: nd                          ! numero de variaveis diescretas
        integer, intent(in) :: nc                          ! numero de variaveis continuas
        integer, intent(in) :: Vd(nd)                      ! valor das variaveis discretas
        real(8), intent(in) :: Vc(nc)                      ! valor das variaveis continuas
        real(8) :: fob                                     ! função objetivo
        real(8), allocatable :: Nrd(:,:)                   ! Força resistente de calculo das barras de treliça [KN] 
    
        !Atribui coordenadas para os nós da estrutura
            call node_barras(n_div, h1, num_nos, coord)
            
        !Estabelece o número de barras e suas conectividades
            call barras_trelica (n_div, h1, coord, n_barras, barra)
            
        !Atribui uma seção transversal para cada barra de acordo com a lista de seções da AISC
            call secao_barras(n_cant, cant, n_barras, barra)
    
        !Calcula a a magnitude de carga de peso próprio em cada nó da estrutura
            call carga_pp_barras (rho, h1, n_div, num_nos, n_barras, barra, cond_cont)

        !Calcula a carga das terças e telhas nos nós superiores
            call carga_tercas_telhas_barras (h1, dist_trelica, n_div, num_nos, barra, terca, cond_cont)
    
        !Armazena os coeficientes Ce e Cpe da NBR 6123:2013
            call coeficiente_pressao (pressao_vento)
    
        !Combina os valores de Ce e Ci e guarda as combinações com o valor máx e min de succão
            call combinacao_carregamento_vento (h, a, b, ang_cobertura, pressao_vento, coef_max_succao, coef_min_succao)
    
        !Calcula a carga distribuída na cobertura devido à ação do vento
           ! select case(coeficiente_vento)
           !     case(0)
           !         call carga_distribuida_vento(Vo, S1, S2, coef_max_succao, barra, F)
           !     case(1)
           !         call carga_distribuida_vento(Vo, S1, S2, coef_min_succao, barra, F)
           ! end case
        !Calcula as cargas nodais da ação do vento na estrutura
            call carga_vento_sobrecarga_nos(theta, barra, F, cond_cont)
    
        !Calcula as combinações de ações para Combinações Últimas Normais 
            call combina_acoes (cond_cont, comb_acao)
    
        !Matriz de rigidez - vetor de deslocamentos fixos e livres / vetor de cargas
            !if (coeficiente_vento == 0) then
            !call vetor_Y_vglc (comb_, Y, v_glc)
            
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
            
    end subroutine
    
    end module