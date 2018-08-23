
    program Teste
    
    use CargaVento

    implicit none

    ! Variables S1-------------------------------------------------------------------------------------------------------------------
    character(4) :: caso                            ! a) Terreno plano ou fracamente acidentado; b) Taludes e morros; c) Vales profundos
    real(8) :: theta                                ! Inclinação média do talude ou enconsta de morro
    real(8) :: z                                    ! Altura medida a partir da superfície do terreno no ponto considerado
    real(8) :: d                                    ! diferença de nível entre a base e o topo do talude ou morro
    namelist /Fator_Top_S1/ caso, theta, z, d       ! namelist com dados de entrada para avaliar o fator  topográfico S1
    ! Variables S2-------------------------------------------------------------------------------------------------------------------   
    character :: cat                                ! categoria referente ao fator S2
    character :: cls                                ! classe referente ao fator S2
    namelist /Fator_Top_S2/ cat, cls                ! namelist com dados de entrada para avaliar o fator  topográfico S2
    type(tabela_S2) :: tab_S2(5)                    ! vetor de armazenamento dos dados da tabela numérica do fator S2
     ! Variables pressão dinâmica----------------------------------------------------------------------------------------------------
    real(8) :: Vk                                   ! Velocidade característica do vento da região, em m/s
    namelist /pressao_dinamica/ Vk                  ! namelist com os dados de entrada para avaliar a pressão dinâmica
    ! Variables Coef. pressão--------------------------------------------------------------------------------------------------------
    real(8) :: h                                    ! comprimento da altura da coluna de sustentação + montante imediatamente acima / vista frontal
    real(8) :: b                                    ! menor dimensão da estrutura em vista superior
    real(8) :: a                                    ! maior dimensão da estrutura em vista superior
    integer :: ang_cobertura                        ! inclinação da viga de cobertura treliçada entre montante externo e intermediário
    namelist /Cofc_pressao/ h, a, b, ang_cobertura  ! namelist com os valores de entrada para avaliar Ce e Cpe
    type(Coef_pressao) :: pressao_vento(3)          ! Variável de armazenamento dos coeficientes Ce e Cpe
    type(Combina_coef) :: coef_comb
    
    
    integer, parameter :: uni_dad_vento = 100                     ! nome dado para número que identifica arquivo externo
    integer :: i                                                  ! variavel de recurso cíclico
    integer :: n                                                  ! variável de recurso cíclico 
    integer :: junk1
    character :: junk2
    ! -----------------------------------------------------------------------------------------------------------------------------
    
    
    
    ! Abertura de arquivo e leitura de dados---------------------------------------------------------------------------------------
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
    
    
    ! Leitura dos dados de pressão Ce e Cpe ---------------------------------------------------------------------------------------
    call Coeficiente_pressao (pressao_vento)
    call Combinacao_carregamento_vento (h, a, b, ang_cobertura, pressao_vento, coef_comb)
      
    end program Teste

