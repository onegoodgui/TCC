 program Programa_Trelica

    
    use LLS
    use LLSt
    use Estrutura_Trelica
    use CargaVento
    use sislin_pbMEF
    use Matriz_Rigidez
    
    implicit none
    

    ! Variables S1-------------------------------------------------------------------------------------------------------------------
    character(4) :: caso                            ! a) Terreno plano ou fracamente acidentado; b) Taludes e morros; c) Vales profundos
    !integer :: theta                                ! Inclina��o m�dia do talude ou enconsta de morro
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
    !namelist /Comb_acoes/ comb_acao(1)%nome, comb_acao(2)%nome, comb_acao(3)%nome, comb_acao(4)%nome, comb_acao(5)%nome
    
    integer, parameter :: uni_dad_vento = 100                     ! nome dado para n�mero que identifica arquivo externo
    integer :: i                                                  ! variavel de recurso c�clico
    integer :: n                                                  ! vari�vel de recurso c�clico 
    integer :: junk1
    character :: junk2
    ! -----------------------------------------------------------------------------------------------------------------------------
  
    
    !Corpo do programa principal --------------------------------
    
    !Abertura de arquivo e leitura de dados---------------------------------------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Fator_Top_S1)
    close (uni_dad_vento)
    
    ! Chamada para c�lculo do fator topogr�fico S1 --------------------------------------------------------------------------------
    call Fator_Topografico_S1 (caso, theta, z, d, S1)
   
      
    !Leitura dos dados do problema referentes ao fator S2 -------------------------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Fator_Top_S2)
    close (uni_dad_vento)
    
    ! Chamada para c�lculo do fator topogr�fico S2 --------------------------------------------------------------------------------
    call Transf_dados_S2 (tab_S2)
    call Fator_Topografico_S2(cat, cls, z, tab_S2, S2)
    
    !Leitura dos dados do problema referentes a press�o din�mica ------------------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = pressao_dinamica)
    close (uni_dad_vento)
    
    !Leitura dos dados do problema referentes ao coeficiente de press�o -----------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Cofc_pressao)
    close (uni_dad_vento)
    
    !Leitura dos dados do problema referentes aos dados geom�tricos de entrada -----------------------------------------------------------
    open (unit = uni_dad_vento, file = 'Vento_Entrada.txt', action='read')        
        read (uni_dad_vento, nml = Dados_geometricos)
    close (unit = uni_dad_vento)
    
    !Armazenamento dos dados do problema referentes aos v�rios perfis de ter�as  -----------------------------------------------------------
    open (unit = uni_dad_tercas, file = 'dados_tercas.txt', action='read')        
        call lista_tercas(uni_dad_tercas, terca)
    close (unit = uni_dad_tercas)
    
    !Armazenamento das propriedades geom�tricas referentes aos v�rios perfis 2L da AISC  -----------------------------------------------------------    
    open (unit = uni_entrada, file = 'lista_perfis_aco.txt', action='read')
        call armazenar_propriedades_cantoneiras_abas_iguais(n_cant, cant)
    close (unit = uni_entrada)
    
    do i = 1, n_cant
     call LLS_propriedades_geometricas (cant(i))
    end do 
    
    end program Programa_Trelica