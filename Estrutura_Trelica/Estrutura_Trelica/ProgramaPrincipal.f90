 program Programa_Trelica

    
    use LLS
    use LLSt
    use Estrutura_Trelica
    use CargaVento
    use sislin_pbMEF
    use Matriz_Rigidez
    use Analise_trelica_trapezoidal
    use HSAg
    use numRAND
    
    implicit none
    
    integer :: nd                                      ! numero de variaveis diescretas
    integer :: nc                                      ! numero de variaveis continuas
    integer, allocatable :: Vd(:)                      ! valor das variaveis discretas
    real(8), allocatable :: Vc(:)                      ! valor das variaveis continuas
    real(8) :: fob                                     ! função objetivo
    
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
    
    call funcao_objetivo(nd, nc, Vd, Vc, fob)
    
    
    
    
    
   ! do i = 1, n_barras
   !     call encontrar_MDC_pbMEF(nne,barra(i)%conectividades(:),MDC)
   ! end do
    
    !call iniciar_MDSB_pbMEF(num_nos,ngln,MDC,MDSB_trelica,err)
    
    
    
    
    end program Programa_Trelica   