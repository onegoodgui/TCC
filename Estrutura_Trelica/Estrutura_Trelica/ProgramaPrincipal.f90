    program Programa_Trelica

    
    use LLS
    use LLSt
    use Estrutura_Trelica
    use CargaVento
    use sislin_pbMEF
    use Matriz_Rigidez
    
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
    real(8) :: rho = 0.000078d0          ! Massa específica do aço [KN/cm³]
    ! Variables S1-------------------------------------------------------------------------------------------------------------------
    character(4) :: caso                            ! a) Terreno plano ou fracamente acidentado; b) Taludes e morros; c) Vales profundos
    !integer :: theta                                ! Inclinação média do talude ou enconsta de morro
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
    
    
    integer, parameter :: uni_dad_vento = 100                     ! nome dado para número que identifica arquivo externo
    integer :: i                                                  ! variavel de recurso cíclico
    integer :: n                                                  ! variável de recurso cíclico 
    integer :: junk1
    character :: junk2
    ! -----------------------------------------------------------------------------------------------------------------------------
  
    
    !Corpo do programa principal --------------------------------
    
    !Abertura de arquivo e leitura de dados---------------------------------------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Fator_Top_S1)
    close (uni_dad_vento)
    
    ! Chamada para cálculo do fator topográfico S1 --------------------------------------------------------------------------------
    call Fator_Topografico_S1 (caso, theta, z, d, S1)
   
      
    !Leitura dos dados do problema referentes ao fator S2 -------------------------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Fator_Top_S2)
    close (uni_dad_vento)
    
    ! Chamada para cálculo do fator topográfico S2 --------------------------------------------------------------------------------
    call Transf_dados_S2 (tab_S2)
    call Fator_Topografico_S2(cat, cls, z, tab_S2, S2)
    
    !Leitura dos dados do problema referentes a pressão dinâmica ------------------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = pressao_dinamica)
    close (uni_dad_vento)
    
    !Leitura dos dados do problema referentes ao coeficiente de pressão -----------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Cofc_pressao)
    close (uni_dad_vento)
    !Leitura dos dados do problema referentes aos dados geométricos de entrada -----------------------------------------------------------
    open (unit = uni_dad_vento, file = 'Vento_Entrada.txt', action='read')        
        read (uni_dad_vento, nml = Dados_geometricos)
    close (unit = uni_dad_vento)
    !Armazenamento dos dados do problema referentes aos vários perfis de terças  -----------------------------------------------------------
    open (unit = uni_dad_tercas, file = 'dados_tercas.txt', action='read')        
        call lista_tercas(uni_dad_tercas, terca)
    close (unit = uni_dad_tercas)
    !Armazenamento das propriedades geométricas referentes aos vários perfis 2L da AISC  -----------------------------------------------------------    
    open (unit = uni_entrada, file = 'lista_perfis_aco.txt', action='read')
        call armazenar_propriedades_cantoneiras_abas_iguais(n_cant, cant)
    close (unit = uni_entrada)
    
    do i = 1, n_cant
     call LLS_propriedades_geometricas (cant(i))
    end do
    
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
    call carga_distribuida_vento(Vo, S1, S2, coef_max_succao, barra, F)
    
    !Matriz de rigidez
    call vetor_Y_vglc (cond_cont, Y, v_glc)
    
    do i = 1, n_barras
        call encontrar_MDC_pbMEF(nne,barra(i)%conectividades(:),MDC)
    end do
    
    call iniciar_MDSB_pbMEF(num_nos,ngln,MDC,MDSB_trelica,err)
    
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
    
    
    end program Programa_Trelica